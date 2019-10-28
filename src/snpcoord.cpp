#include <ldmap/ldmap.hpp>
#include <ldmap/genetic_map.hpp>
#include <progress.hpp>
#include <RcppParallel.h>
#include <string>
//[[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "Eigen/src/SparseCore/SparseUtil.h"
#include "Rcpp/Nullable.h"
#include "Rcpp/vector/instantiation.h"
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_sort.h"


//' Find out of vector of SNPs is sorted
//' 
//' @param chr vector of chromosomes	of per region
//' @param pos vector of start positions for each region
//' @return returns true if the vector is sorted
//' @export
//[[Rcpp::export]]
bool sorted_snp_df(const Rcpp::IntegerVector chr,  const Rcpp::IntegerVector pos){


  const size_t p=chr.size();
    //std::vector<std::pair<int,int>> r(p);
  std::pair<int,int> old_p{chr[0],pos[0]};
  std::pair<int,int> new_p;
  for(int i=1; i<p;i++){
    new_p={chr[i],pos[i]};
    if(new_p<old_p){
      return(false);
    }
    old_p=new_p;
  }
  return(true);
}


struct SNPpos{
  std::pair<int,int> c_p;
  SNPpos(const int chrom,const int pos):c_p(std::make_pair(chrom,pos)){};
};

struct GRange{
  std::pair<SNPpos,SNPpos> range;
  GRange(const int chrom,const int start,const int end):range(std::make_pair(SNPpos(chrom,start),SNPpos(chrom,end))){
    if(end<start){
      Rcpp::Rcerr<<"For range: "<<chrom<<":"<<start<<"-"<<end<<std::endl;
      Rcpp::stop("Invalid range! end<start");
    }};
};
bool operator<(const SNPpos p, const GRange g){
  return(p.c_p<g.range.first.c_p);
}
bool operator<=(const SNPpos p, const GRange g){
  return(p.c_p<=g.range.first.c_p);
}
bool operator>(const SNPpos p, const GRange g){
  return(p.c_p>=g.range.second.c_p);
}
// bool operator>=(const SNPpos p, const GRange g){
//   return(p.c_p>=g.range.second.c_p);
// }
bool operator==(const SNPpos p, const GRange g){
  return((p.c_p>=g.range.first.c_p) && (p.c_p<g.range.second.c_p));
}

// 
// struct pairhash {
// public:
//   std::size_t operator()(const std::pair<uint32_t, uint32_t> &x) const
//   {
//     return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
//   }
// };

//' Assign snps to regions of the genome, breaking up regions based on number of SNPs
//' 
//' @param ld_chr vector of chromosomes	of per region
//' @param ld_start vector of start positions for each region
//' @param ld_stop vector of end positions for each region
//' @return returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
//' @export
//[[Rcpp::export]]
Rcpp::StringVector set_ld_region(const Rcpp::IntegerVector ld_chr,
				    const Rcpp::IntegerVector ld_start,
				    const Rcpp::IntegerVector ld_stop,
				    const Rcpp::IntegerVector ld_region_id,
				    const Rcpp::IntegerVector chr,
				    const Rcpp::IntegerVector pos,
				    uint32_t max_size=0,
				    int min_size=1,
				    const bool assign_all=true){

  const size_t ld_size=ld_chr.size();

  if(!std::is_sorted(ld_chr.begin(),ld_chr.end())){
    Rcpp::stop("break regions must be sorted by chromosome!");
  }

  const size_t snp_size = chr.size();
  Rcpp::StringVector ret_region(snp_size);
  using pair_i = uint64_t;
  std::vector< pair_i > tret_region(snp_size);
  using ref_p =std::pair<std::vector<pair_i>::iterator,std::vector<pair_i>::iterator >;
  std::unordered_map<uint32_t,
                     int> ref_map(ld_size*2);
  if(!std::is_sorted(chr.begin(),chr.end())){
    Rcpp::stop("SNPs must be sorted by chromosome!");
  }
  size_t idx1=0;
  size_t idx2=0;
  GRange gr(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
  max_size = max_size == 0 ? snp_size + 1 : max_size;
  uint32_t grc=0;
  
  Progress p(snp_size, true);
  
  while(idx1<snp_size){
    if (Progress::check_abort() ){
      Rcpp::stop("aborted");
    }
    const SNPpos s_pos(chr[idx1],pos[idx1]);
    if(s_pos==gr){
      //point in interval
      grc++;
      pair_i rel= (ld_region_id[idx2] << 16) | (grc/max_size);
      tret_region[idx1]=rel;
      ref_map[rel]++;
      idx1++;
      p.increment();
    }else{
      if(s_pos<gr){
        // point before interval
        if(assign_all){
          Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
          Rcpp::Rcerr<<"Can't map to "<<ld_chr[idx2]<<":"<<ld_start[idx2]<<"-"<<ld_stop[idx2]<<" (ahead)"<<std::endl;
          if(idx2>0){
            Rcpp::Rcerr<<"Can't map to "<<ld_chr[idx2-1]<<":"<<ld_start[idx2-1]<<"-"<<ld_stop[idx2-1]<<" (behind)"<<std::endl;
          }
          Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
        }else{
          tret_region[idx1]=(NA_INTEGER << 16) | NA_INTEGER;
        }
        idx1++;
        p.increment();
      }
      if(s_pos>gr){
        //point after interval
        if((idx2+1)<ld_size){
          grc=0;
          idx2++;
          gr=GRange(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
        }else{
          if(assign_all){
            Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
            Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
          }else{
            tret_region[idx1]=(NA_INTEGER << 16) | NA_INTEGER;
          }
          idx1++;
          p.increment();
        }
      }
    }
  }

  std::transform(tret_region.begin(),tret_region.end(),ret_region.begin(),[&ref_map,min_size](pair_i eli){
    int high16 = eli >> 16;
    int low16 = eli & 0xFFFF;
    low16 = (ref_map[eli] < min_size) ? std::max(0,low16-1) : low16;
      return((std::to_string(high16)+"."+std::to_string(low16))); 
  });
  return(ret_region);
}











std::unordered_map<int,int> reverse_index(Rcpp::IntegerVector index_a){
    const size_t p_a=index_a.size();
    std::unordered_map<int,int> ind_a(p_a);
    auto iap = ind_a.begin();
    for(int i=0; i<p_a; i++){
      iap = ind_a.emplace_hint(iap,index_a(i),i);
    }
    return ind_a;
}


//[[Rcpp::export]]
Rcpp::List join_ids(Rcpp::IntegerVector index_a,Rcpp::IntegerVector index_b){
  const size_t p_a=index_a.size();
  const size_t p_b=index_b.size();

  const auto map_index_a=reverse_index(index_a);
  std::vector<int> sub_ia;
  sub_ia.reserve(std::min(p_a,p_b));
  std::vector<int> sub_ib;
  sub_ib.reserve(std::min(p_a,p_b));

  for(int i=0; i<p_b; i++){
    auto ipm = map_index_a.find(index_b(i));
    if(ipm != map_index_a.end()){
      sub_ia.push_back(ipm->second);
      sub_ib.push_back(i);
    }
  }
  using namespace Rcpp;

  const auto ns = sub_ia.size();
  auto m_ret_l = List::create(_["idx_a"]=Rcpp::wrap(sub_ia),
                              _["idx_b"]=Rcpp::wrap(sub_ib));
  m_ret_l.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  m_ret_l.attr("row.names") = seq(1,ns);

  return(m_ret_l);
}






//' Creation of new ldmap_snps
//' 
//' @param chrom an integer vector of chromosomes
//' @param pos a double vector of positions 
//' @param ascii_ref an optional integer vector of reference allele (see `?utf8ToInt`)
//' @param ascii_ref an optional integer vector of alternate allele (see `?utf8ToInt`)
//' 
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector new_ldmap_snp(Rcpp::IntegerVector chrom=Rcpp::IntegerVector::create(),Rcpp::NumericVector pos=Rcpp::NumericVector::create(), Rcpp::IntegerVector ascii_ref=Rcpp::IntegerVector::create(),Rcpp::IntegerVector ascii_alt=Rcpp::IntegerVector::create()){

  const size_t p=chrom.size();
  Rcpp::NumericVector ret(p);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");


  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");

  Snp mp;

  for(int i=0; i<p; i++){
    mp.str={static_cast<unsigned char>(chrom(i)),
           static_cast<uint64_t>(pos(i)),
           ascii2Ref(static_cast<int>(ascii_ref(i))),
            ascii2Ref(static_cast<int>(ascii_alt(i)))};
    // uint64_t tchrom=chrom(i);
    // uint64_t tpos=static_cast<uint64_t>(pos(i));
    // uint64_t tref=asciiascii_ref(i)-64;
    // uint64_t talt=ascii_alt(i)-64;
    // mp  = {tchrom,tpos,tref,talt};
    ret(i)=mp.dat;
  }
  return ret;
}






//' Formatting method for ldmap snps
//'
//' @param x a vector of ldmap_snps
//' @method format ldmap_snp
//' @export
//' @export format.ldmap_snp
//[[Rcpp::export("format.ldmap_snp")]]
Rcpp::StringVector format_ldmap_snp(Rcpp::NumericVector x){
  Snp mp;
  //  double *bs = reinterpret_cast<double*>(&mp);
  const size_t p=x.size();
  Rcpp::StringVector result(p);
  std::transform(x.begin(),x.end(),result.begin(),[](const SNP s){
                                                    Rcpp::String st=s.char_rep();
                                                    return(st);
                                                  });
  // for(int i=0; i<p; i++){
  //   double bs=x(i);
  //   mp.dat=bs;
  //   result(i)=
  // }
  return(result);
}



// Rcpp::NumericVector sort_snps(Rcpp::NumericVector struct_vec){

//   Rcpp::NumericVector ret(struct_vec.begin(),struct_vec.end());


//   std::sort(ret.begin(),ret.end(),[](const SNP a,const SNP b){
//                                     return (a<b);
//                                   });
//   ret.attr("class") = struct_vec.attr("class");
//   ret.attr("sorted")=Rcpp::LogicalVector::create(true);
//   return(ret);
// }



//' sorting method for ldmap snps
//'
//' @param struct_vec the vector of SNPs
//'
//' @method order ldmap_snp
//' @export
//' @export order.ldmap_snp
//[[Rcpp::export(order.ldmap_snp)]]
Rcpp::IntegerVector order_snps(Rcpp::NumericVector struct_vec){
  using namespace Rcpp;
  const size_t p= struct_vec.size();
  IntegerVector ret=seq(1, p);
  if(ret.size()!=p){
    Rcpp::stop("you can't make seqs right");
  }
  RcppParallel::RVector<double> input_a(struct_vec);

  std::sort(ret.begin(),ret.end(),[&](const int i,const int j){
                                    return (SNP(input_a[i-1])<SNP(input_a[j-1]));
                                  });


  return(ret);
}

//' ranking method for ldmap snps
//'
//' @param struct_vec the vector of SNPs
//'
//' @method rank ldmap_snp
//' @export
//' @export rank.ldmap_snp
//[[Rcpp::export(rank.ldmap_snp)]]
Rcpp::IntegerVector rank_snps(Rcpp::NumericVector struct_vec){
  using namespace Rcpp;
  const size_t p= struct_vec.size();
  // std::vector<int> ret=seq(1, p);
  std::vector<size_t> idx(p);
  std::iota(idx.begin(), idx.end(), 0);


  std::sort(idx.begin(),idx.end(),[&](const int i,const int j){
                                    return (SNP(struct_vec[i])<SNP(struct_vec[j]));
                                  });
  int i=1;
  IntegerVector sret=no_init(p);
  std::for_each(idx.begin(),idx.end(),[&](const int j) mutable{
                                        sret(j)=i++;
                                      });
  return(sret);
}







//' get chroms from a ldmap_snp
//'
//' @param struct_vec the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector chromosomes(Rcpp::NumericVector struct_vec){
  using namespace Rcpp;
  IntegerVector ret = no_init(struct_vec.size());


  std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](const SNP a){
                                                                   return(static_cast<int>(a.snp.str.chrom));
                                                                 });
  return ret;
}




//' get positions from a ldmap_snp
//'
//' @param struct_vec the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector positions(Rcpp::NumericVector struct_vec){
  using namespace Rcpp;
  NumericVector ret = no_init(struct_vec.size());


  std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](const SNP a){
                                                                   return(static_cast<double>(a.snp.str.pos));
                                                                 });
  return ret;
}




//' get ref alleles from a ldmap_snp
//'
//' @param struct_vec the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
SEXP ref_alleles(Rcpp::NumericVector struct_vec,const bool as_ascii_int=true){
  
  if(as_ascii_int){
    using namespace Rcpp;
  Rcpp::IntegerVector ret = no_init(struct_vec.size());
  std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](const SNP a){
                                                                   return(Ref2int(a.snp.str.ref));
                                                                 });
  return ret;
  }else{

    Rcpp::StringVector ret = Rcpp::StringVector(struct_vec.size());
    std::array<SEXP,5> all_chr  {{
                                  Rf_mkChar(Ref2ascii(Ref::A)),
                                  Rf_mkChar(Ref2ascii(Ref::C)),
                                  Rf_mkChar(Ref2ascii(Ref::G)),
                                  Rf_mkChar(Ref2ascii(Ref::T)),
                                  Rf_mkChar(Ref2ascii(Ref::N))}};
    std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[&](const SNP a){
                                                                     return all_chr[static_cast<int>(a.snp.str.ref)];
                                                                   });
    return ret;
  }
}




//' get ref alleles from a ldmap_snp
//'
//' @param struct_vec the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
SEXP alt_alleles(Rcpp::NumericVector struct_vec,const bool as_ascii_int=true){

  if(as_ascii_int){
    using namespace Rcpp;
  IntegerVector ret = no_init(struct_vec.size());
  std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](const SNP a){
                                                                   return(Ref2int(a.snp.str.alt));
                                                                 });
  return ret;
  }else{

    Rcpp::StringVector ret = Rcpp::StringVector(struct_vec.size());
    std::array<SEXP,5> all_chr  {{
                                Rf_mkChar(Ref2ascii(Ref::A)),
                                Rf_mkChar(Ref2ascii(Ref::C)),
                                Rf_mkChar(Ref2ascii(Ref::G)),
                                Rf_mkChar(Ref2ascii(Ref::T)),
                                Rf_mkChar(Ref2ascii(Ref::N))}};

    std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[&](const SNP a){
                                                                     return all_chr[static_cast<int>(a.snp.str.alt)];
                                                                   });
    return ret;
  }
}









//' Formatting method for ldmap snps
//'
//' @param struct_vec the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
SEXP ldmap_snp_2_dataframe(Rcpp::NumericVector struct_vec){

  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");
  using namespace Rcpp;
  const size_t p=struct_vec.size();
  IntegerVector chrom(p);
  NumericVector pos(p);
  IntegerVector ascii_ref(p);
  IntegerVector ascii_alt(p);

  Snp mp;
  //  double *bs = reinterpret_cast<double*>(&mp);
  for(int i=0; i<p; i++){
    double bs=struct_vec(i);
    mp.dat=bs;
    chrom(i)=mp.str.chrom;
    pos(i)=static_cast<double>(mp.str.pos);
    ascii_ref(i)=Ref2int(mp.str.ref);
    ascii_alt(i)=Ref2int(mp.str.alt);
  }
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["pos"]=pos,
                                _["ascii_ref"]=ascii_ref,
                                _["ascii_alt"]=ascii_alt);
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") = seq(1, p);
  return(dfl);
}


//' Find SNPs that match (and find out what to do with them)
//'
//' @param query a vector of (sorted) ldmap_snps SNPs
//' @param reference a vector of (sorted) ldmap_snps SNPs
//'
//' @export
//[[Rcpp::export]]
Rcpp::List join_snp(Rcpp::NumericVector query,Rcpp::NumericVector reference,Rcpp::IntegerVector rsid=Rcpp::IntegerVector::create()){



  Rcpp::Nullable<Rcpp::LogicalVector> n_sorted_a = query.attr("sorted");
  Rcpp::Nullable<Rcpp::LogicalVector> n_sorted_b = reference.attr("sorted");

  const bool sorted_a = n_sorted_a.isNotNull() ? Rcpp::as<bool>(n_sorted_a.get()) : false;
  if(!sorted_a){
    Rcpp::Rcerr<<"Assuming query is sorted (attribute 'sorted' not detected in `query`)"<<std::endl;
  }
    const bool sorted_b = n_sorted_b.isNotNull() ? Rcpp::as<bool>(n_sorted_b.get()) : false;

  if(!sorted_b){
    Rcpp::Rcerr<<"Assuming query is sorted (attribute 'sorted' not detected in `reference`)"<<std::endl;
  }
  RcppParallel::RVector<double> input_a(query);
  RcppParallel::RVector<double> input_b(reference);




  const size_t q_size =query.size();
  Rcpp::IntegerVector ret_index(q_size);
  Rcpp::IntegerVector ret_match(q_size,1);
  Rcpp::NumericVector ret_ref(q_size);

  auto o_ref_b = input_b.begin();
  auto ref_b = input_b.begin();
  auto ref_e = input_b.end();


  for(int i=0; i<q_size; i++){
    SNP q(query(i));
    auto [low,hi] = std::equal_range(ref_b,ref_e,q,[](const SNP a,const SNP b){
                                                     return (a<b);
                                                   });

    for(auto &tl = low; tl!=hi; tl++){

      auto qm = static_cast<int>(q.allele_match(SNP(*tl)))+1;
      if(qm>ret_match(i)){
        ret_index(i)=std::distance(o_ref_b,tl)+1;
        ret_match(i)=qm;
        ret_ref(i)=*tl;
      }
    }
    ref_b=hi;
    if(ref_b==ref_e){
      break;
    }
  }
  ret_ref.attr("class")=reference.attr("class");

  using namespace Rcpp;
  std::vector<std::string> levs({"no_match",
                                 "mismatched_alleles" ,      // mismatch
                                 "one_mono_allelic",
                                 "reverse_ref_alt",
                                 "strand_flip_match",
                                 "strand_flip_reverse_match",
                                 "snp_match"});
  if(levs.size()!=7){
    Rcpp::stop("Ugh! it's "+std::to_string(levs.size()));
  }
  ret_match.attr("levels") =wrap(levs);
  ret_match.attr("class") = Rcpp::StringVector::create("factor");

  auto dfl = List::create(_["query"]=query,
                             _["index"]=ret_index,
                             _["match_type"]=ret_match,
                             _["match"]=ret_ref);
  if(rsid.size()==reference.size()){
    Rcpp::IntegerVector ret_rsid = no_init(q_size);
    std::transform(ret_index.begin(),ret_index.end(),ret_rsid.begin(),[&](const int i){
        return rsid[i];
      });
    dfl["rsid"]=ret_rsid;
  }
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") = seq(1, q_size);
  return dfl;

}
