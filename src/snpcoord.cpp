


#include "R_ext/Arith.h"
#include "ldmap/genetic_map.hpp"
#include "Rinternals.h"
#include "alleles.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <progress.hpp>
#include <RcppParallel.h>
#include <string>
#include <unordered_map>
//[[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "Eigen/src/SparseCore/SparseUtil.h"
#include "Rcpp/Nullable.h"
#include "Rcpp/vector/instantiation.h"
#if __has_include("tbb/parallel_for.h")
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_sort.h"
#endif


#if __has_include(<charconv>)
#include <charconv>
#else
#include <cstdio>
#endif




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







//' Creation of new ldmap_ranges
//'
//' @param chrom an integer vector of chromosomes
//' @param start an integer vector of start positions
//' @param stop an integer vector of stop positions
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector new_ldmap_range(Rcpp::IntegerVector chrom=Rcpp::IntegerVector::create(),
                                    Rcpp::IntegerVector start=Rcpp::IntegerVector::create(),
                                    Rcpp::IntegerVector end=Rcpp::IntegerVector::create()){

  const size_t p=chrom.size();
  Rcpp::NumericVector ret=Rcpp::no_init(p);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");


  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");
  static_assert(0b0000000001==1L);
  static_assert(0b00010011==19L);


  for(int i=0; i<p; i++){
    bed_range br={.str={.end=static_cast<uint64_t>(end(i)),
                        .start=static_cast<uint64_t>(start(i)),
                        .chrom=static_cast<unsigned char>(chrom(i))}};
      ret(i)=br.flt;
  }
  return ret;
}




int snp_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end){

  SNP sx{.snp={.flt=x}};
  auto xb=std::lower_bound(begin,end,sx,[](double  range_a,const SNP &snp_b){
                                          return(Region{.br={.flt=range_a}}< snp_b);
                                        });
  if( (xb==end ) or !(Region{.br={.flt=*xb}} ==sx))
    return NA_INTEGER;
  return std::distance(begin,xb)+1;
}


//' formatting of ldmap_ranges
//'
//' @param ldmap_snp vector of ldmap_snps
//' @param ldmap_range vector of ldmap_snps
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector snp_in_range(Rcpp::NumericVector ldmap_snp,Rcpp::NumericVector ldmap_range){

  const size_t p=ldmap_snp.size();
  RcppParallel::RVector<double> input_snp(ldmap_snp);
  RcppParallel::RVector<double> input_range(ldmap_range);
  Rcpp::IntegerVector ret = Rcpp::no_init(p);
  RcppParallel::RVector<int> output_range(ret);
  auto irb = input_range.begin();
  auto ire = input_range.end();
  static_assert(bed_range{.str={.end=100,.start=50,.chrom=1}}.str.end == 100,"bed_range is having issues");
  // static_assert(bed_range{.str={.end=2147483647ull,.start=247344518ull,.chrom=1}}.str.end == 100,"bed_range is having issues");


  #if __has_include("tbb/parallel_for.h")
  tbb::parallel_for(tbb::blocked_range<size_t>(0,p),
                    [&input_snp=std::as_const(input_snp),&output_range,&irb=std::as_const(irb),&ire=std::as_const(ire)](const tbb::blocked_range<size_t> &r){
                      std::transform(
                                     input_snp.begin()+r.begin(),
                                     input_snp.begin()+r.end(),
                                     output_range.begin()+r.begin(),
                                     [=](const double x){
                                       return(snp_in_range(x,irb,ire));
                                     });
                    });
  #else
  std::transform(
                 input_snp.begin(),
                 input_snp.end(),
                 output_range.begin(),
                 [=](const double x){
                   return(snp_in_range(x,irb,ire));
                 });
#endif
  return ret;

}





//' formatting of ldmap_ranges
//'
//' @param x an ldmap_range
//' @export
//[[Rcpp::export]]
Rcpp::StringVector format_ldmap_range(Rcpp::NumericVector x){

  const size_t p=x.size();
  Rcpp::StringVector ret=Rcpp::no_init(p);

  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");
  static_assert(0b0000000001==1L);
  static_assert(0b00010011==19L);


  for(int i=0; i<p; i++){
    bed_range br={.flt=x(i)};
    auto [end,start,chrom] = br.str;
    if((end==0) or (start==0) or (chrom==0)){
      ret(i)=NA_STRING;
    }else{
      ret(i)=Rcpp::String("chr"+std::to_string(static_cast<int>(chrom))+":"+std::to_string(start)+"_"+std::to_string(end));
    }
  }
  return ret;
}


//' convert ldmap_range to dataframe
//'
//' @param ldmap_range an ldmap_range
//' @export
//[[Rcpp::export]]
SEXP ldmap_range_2_data_frame(Rcpp::NumericVector ldmap_range){

  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");

  using namespace Rcpp;
  const size_t p=ldmap_range.size();
  IntegerVector chrom = Rcpp::no_init(p);
  NumericVector start = Rcpp::no_init(p);
  NumericVector end = Rcpp::no_init(p);


  //  double *bs = reinterpret_cast<double*>(&mp);
  for(int i=0; i<p; i++){
    Region sp{.br={.flt=ldmap_range(i)}};
    chrom(i)= sp.br.str.chrom;
    start(i)= sp.br.str.start;
    end(i)=sp.br.str.end;
  }
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["start"]=start,
                                _["end"]=end);
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") = seq(1, p);
  return(dfl);
}








template<typename Vec_t>
Rcpp::NumericVector new_ldmap_snp_s(Rcpp::IntegerVector chrom,
                                    Rcpp::NumericVector pos,
                                    Vec_t ref,
                                    Vec_t alt,
                                    bool NA2N
                                    ){

  const size_t p=chrom.size();
  const size_t allele_size=ref.size();
  if(allele_size!=alt.size())
    Rcpp::stop("ref and alt must be the same size");


  Rcpp::NumericVector ret=Rcpp::no_init(p);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");

  if(NA2N){
    for(int i=0; i<p; i++){
      SNP mp=SNP::make_snp<true>(static_cast<unsigned char>(chrom(i)),
                                 static_cast<uint64_t>(pos(i)),
                                 Nuc{ascii2Nuc(ref(i))},
                                 Nuc{ascii2Nuc(alt(i))});
      ret(i)=mp.to_double();
    }
  }else{
    for(int i=0; i<p; i++){
      SNP mp=SNP::make_snp<false>(static_cast<unsigned char>(chrom(i)),
                                  static_cast<uint64_t>(pos(i)),
                                  Nuc{ascii2Nuc(ref(i))},
                                  Nuc{ascii2Nuc(alt(i))}
                                  );
      ret(i)=mp.to_double();
    }
  }
  return ret;
}



template<typename Vec_t>
Rcpp::NumericVector new_ldmap_snp_s(Rcpp::IntegerVector chrom,
                                    Rcpp::NumericVector pos,
                                    Vec_t ref,
                                    bool NA2N
                                    ){

  const size_t p=chrom.size();
  const size_t allele_size=ref.size();
  if(allele_size!=p)
    Rcpp::stop("ref and chrom/pos must be the same size");


  Rcpp::NumericVector ret=Rcpp::no_init(p);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");

  if(NA2N){
    for(int i=0; i<p; i++){
      SNP mp=SNP::make_snp<true>(static_cast<unsigned char>(chrom(i)),
                                 static_cast<uint64_t>(pos(i)),
                                 Nuc{ascii2Nuc(ref(i))});
      ret(i)=mp.to_double();
    }
  }else{
    for(int i=0; i<p; i++){
      SNP mp=SNP::make_snp<false>(static_cast<unsigned char>(chrom(i)),
                                 static_cast<uint64_t>(pos(i)),
                                 Nuc{ascii2Nuc(ref(i))});
      ret(i)=mp.to_double();
    }
  }
  return ret;
}



Rcpp::NumericVector new_ldmap_snp_s(Rcpp::IntegerVector chrom,
                                    Rcpp::NumericVector pos,const bool NA2N){

  const size_t p=chrom.size();
  Rcpp::NumericVector ret=Rcpp::no_init(p);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");

  if(NA2N){
    for(int i=0; i<p; i++){
      SNP mp=SNP::make_snp<true>(static_cast<unsigned char>(chrom(i)),
                                 static_cast<uint64_t>(pos(i)));
      ret(i)=mp.to_double();
    }
  }else{
    for(int i=0; i<p; i++){
      SNP mp=SNP::make_snp<false>(static_cast<unsigned char>(chrom(i)),
                                 static_cast<uint64_t>(pos(i)));
      ret(i)=mp.to_double();
    }
  }

  return ret;
}



//' Creation of new ldmap_snps
//'
//' @param chrom an integer vector of chromosomes
//' @param pos a double vector of positions
//' @param ref an optional vector of reference allele (see `?new_ldmap_allele`)
//' @param alt an optional vector of alternate allele (see `?new_ldmap_allele`)
//' @param NA2N an optional boolean specifying whether missing/NA alleles should be treated as "N"
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector new_ldmap_snp(Rcpp::IntegerVector chrom=Rcpp::IntegerVector::create(),
                                  Rcpp::NumericVector pos=Rcpp::NumericVector::create(),
                                  Rcpp::RObject ref=Rcpp::IntegerVector::create(),
                                  Rcpp::RObject alt=Rcpp::IntegerVector::create(),
                                  const bool NA2N=false){




  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");





  const size_t allele_size = Rf_length(SEXP(ref));
  if(Rf_length(SEXP(alt))!=allele_size){
    Rcpp::stop("ref and alt must be the same size or can both be empty/NULL or alt can be empty/NULL");
  }
  bool use_ref=allele_size!=0;
  bool use_alt=Rf_length(SEXP(alt))!=0;

  if(use_ref and use_alt and (ref.sexp_type() != alt.sexp_type())){
    Rcpp::stop("ref and alt must both be of the same type if both are to be used");
  }


  if(use_ref){
    switch (ref.sexp_type()) {
    case INTSXP:{
      if(use_alt)
        return new_ldmap_snp_s<Rcpp::IntegerVector>(chrom,pos,SEXP(ref),SEXP(alt),NA2N);
      return new_ldmap_snp_s<Rcpp::IntegerVector>(chrom,pos,SEXP(ref),NA2N);
    }
    case STRSXP:{
      if(use_alt)
        return new_ldmap_snp_s<Rcpp::StringVector>(chrom,pos,SEXP(ref),SEXP(alt),NA2N);
      return new_ldmap_snp_s<Rcpp::StringVector>(chrom,pos,SEXP(ref),NA2N);
    }
    case RAWSXP:{
      if(use_alt)
        return new_ldmap_snp_s<Rcpp::RawVector>(chrom,pos,SEXP(ref),SEXP(alt),NA2N);
      return new_ldmap_snp_s<Rcpp::RawVector>(chrom,pos,SEXP(ref),NA2N);
    }
    default:
      Rcpp::stop("ref and alt must be coercable to ldmap_allele (i.e must be integer, character or raw)");
    }
  }
  return new_ldmap_snp_s(chrom,pos,NA2N);
}



//' Which SNPs are strand ambiguous
//'
//' @param x the ldmap_snp struct
//'
//' @export
//[[Rcpp::export]]
Rcpp::LogicalVector is_strand_ambiguous(Rcpp::NumericVector x){

  Rcpp::LogicalVector retvec = Rcpp::LogicalVector::import_transform(x.begin(),x.end(),[](double x){
                                                                                         return SNP{.snp={.flt=x}}.is_strand_ambiguous();
                                                                                      }

    );
  return retvec;
}


Rcpp::RawVector new_ldmap_allele_i(Rcpp::IntegerVector allele){
  return Rcpp::RawVector::import_transform(allele.begin(),allele.end(),[&](const int i){
                                                                         return static_cast<Rbyte>(ascii2Nuc(static_cast<char>(i)));});
}

Rcpp::RawVector new_ldmap_allele_s(Rcpp::StringVector allele){
  return Rcpp::RawVector::import_transform(allele.begin(),allele.end(),[&](Rcpp::String i){
                                                                         const char* ip = i.get_cstring();
                                                                         return static_cast<Rbyte>(ascii2Nuc(ip));});
}

Rcpp::RawVector new_ldmap_allele_r(Rcpp::RawVector allele){
  return Rcpp::RawVector::import_transform(allele.begin(),allele.end(),[&](const Rbyte i){
                                                                         return static_cast<Rbyte>(ascii2Nuc(i));});
}


//' Convert allele info to ldmap_alleles
//'
//' @param allele vector of alleles, coded either as character, integer or raw
//'
//' @export
//[[Rcpp::export]]
Rcpp::RawVector new_ldmap_allele(Rcpp::RObject allele=Rcpp::IntegerVector::create()){
  if(allele.isNULL()){
    Rcpp::RawVector ret(0);
    ret.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
    return(ret);
  }
  const size_t p = Rf_length(SEXP(allele));
  Rcpp::RawVector ret = Rcpp::no_init(p);
  switch (allele.sexp_type()) {
  case INTSXP:{
    ret=new_ldmap_allele_i(SEXP(allele));
    break;
  }
  case STRSXP:{
    ret=new_ldmap_allele_s(SEXP(allele));
    break;
  }
  case RAWSXP:{
    ret=new_ldmap_allele_r(SEXP(allele));
    break;
  }
  default:
    Rcpp::stop("`allele` must be of  string or integer type");
  }
  ret.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
  return ret;
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
                                    return SNP{.snp={.flt=input_a[i-1]}}<SNP{.snp={.flt=input_a[j-1]}};
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
                                    return SNP{.snp={.flt=struct_vec[i]}}<SNP{.snp={.flt=struct_vec[j]}};//(SNP(struct_vec[i])<SNP(struct_vec[j]));
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

  std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](double x){
                                                                   return(static_cast<int>(Snp{.flt=x}.str.chrom));
                                                                 });
  return ret;
}




//' get positions from a ldmap_snp
//'
//' @param ldmap_snp the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector positions(Rcpp::NumericVector ldmap_snp){
  using namespace Rcpp;
  NumericVector ret = no_init(ldmap_snp.size());


  std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[](double x){
                                                                   return(static_cast<double>(Snp{.flt=x}.str.pos));
                                                                 });
  return ret;
}

class GlobalNucs{
  mutable bool init_vec;
  std::array<Rcpp::String,16> global_nucs;
  std::array<int,16> global_ints;
  mutable Rcpp::StringVector ascii_strings;


public:
  GlobalNucs():init_vec(false){
    global_nucs[0]=NA_STRING;
    global_ints[0]=0;
    for(char i=1; i<global_nucs.size(); i++){
      global_nucs[i]=Rcpp::String(Nuc2string(Nuc{i}));
    }
  }
  void init_ascii_strings() const{
    ascii_strings=Rcpp::StringVector(global_nucs.begin(),global_nucs.end());
    init_vec=true;
  }
  Rcpp::String get_nuc(const char* i) const{
    return global_nucs[ascii2Nuc(i)];
  }
  Rcpp::String get_nuc(const int i) const{
    if(i>=global_nucs.size()|| i<0){
      return global_nucs[0];
    }
    return global_nucs[ascii2Nuc(i)];
  }
  Rcpp::String get_nuc(const Nuc &i) const{
    if(i.let>=global_nucs.size() || i.let<0){
      return global_nucs[0];
    }
    return global_nucs[i.let];
 }
  Rbyte ascii_int(const Nuc &i) const{
    if(i.let>=global_nucs.size() || i.let<0){
      return global_ints[0];
    }
    return global_ints[i.let];
  }
  Rcpp::StringVector ascii_factors() const{
    if(!init_vec){
      init_ascii_strings();
    }
    return ascii_strings;
  }
  Rcpp::String string_rep(const SNP& s) const{
    auto [alt,ref,pos,chrom] = s.snp.str;
    // Rcpp::Rcerr<<"chrom: "<<static_cast<int>(chrom)<<" pos: "<<pos<<" ref: "<<static_cast<int>(ref)<<" alt: "<<static_cast<int>(alt)<<std::endl;
    if(pos==0 && chrom==0){
      return Rcpp::String(NA_STRING);
    }else{
      return Rcpp::String("chr"+
                          std::to_string(static_cast<int>(chrom))+
                          ":" +
                          std::to_string(static_cast<uint64_t>(pos)) +
                          "_"+
                          std::string(get_nuc(Nuc{static_cast<char>(ref)}).get_cstring())+
                          "_"+
                          std::string(get_nuc(Nuc{static_cast<char>(alt)}).get_cstring()));
    }
  }
  friend Nuc;
  friend SNP;
};

GlobalNucs globalnucs;


//' get ref alleles from a ldmap_snp
//'
//' @param ldmap_snp the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
SEXP ref_alleles(Rcpp::NumericVector ldmap_snp,const bool as_ldmap_allele=true){

  if(as_ldmap_allele){
    using namespace Rcpp;
    Rcpp::RawVector ret = no_init(ldmap_snp.size());
    std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[](double x){
                                                                     return(static_cast<Rbyte>(Snp{.flt=x}.str.ref));
                                                                   });
    ret.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
    return ret;
  }else{
    Rcpp::StringVector ret = Rcpp::StringVector(ldmap_snp.size());
    std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[&](double x){

                                                                     return globalnucs.get_nuc(Nuc{Snp{.flt=x}.str.ref});
                                                                   });
    return ret;
  }
}





//' get ref alleles from a ldmap_snp
//'
//' @param x ldmap allele vec
//'
//' @export
//[[Rcpp::export]]
Rcpp::StringVector format_ldmap_allele(Rcpp::RawVector x){


  return   Rcpp::StringVector::import_transform(x.begin(),x.end(),[&](const Rbyte a){
                                                                    return globalnucs.get_nuc(Nuc{static_cast<char>(a)});
                                                                  });
}





//' coerce ldmap allele to integer
//'
//' @param x ldmap_allele vec
//'
//[[Rcpp::export]]
Rcpp::IntegerVector as_integer_ldmap_allele(Rcpp::RawVector x){


  return   Rcpp::IntegerVector::import_transform(x.begin(),x.end(),[&](const Rbyte a){
                                                                     return static_cast<int>(a);
                                                                   });
}



//' coerce ldmap range to integer(64)
//'
//' @param x ldmap_range vec
//'
//[[Rcpp::export]]
Rcpp::NumericVector as_integer_ldmap_range(Rcpp::NumericVector x){
  Rcpp::NumericVector y = Rcpp::clone(x);
  y.attr("class") = Rcpp::StringVector::create("integer64");
  return y;
}


//' coerce ldmap snp to integer(64)
//'
//' @param x ldmap_snp vec
//'
//[[Rcpp::export]]
Rcpp::NumericVector as_integer_ldmap_snp(Rcpp::NumericVector x){
  Rcpp::NumericVector y = Rcpp::clone(x);
  y.attr("class") = Rcpp::StringVector::create("integer64");
  return y;
}






//' get ref alleles from a ldmap_snp
//'
//' @param ldmap_snp the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
SEXP alt_alleles(Rcpp::NumericVector ldmap_snp,const bool as_ldmap_allele=true){
   if(as_ldmap_allele){
    using namespace Rcpp;
    Rcpp::RawVector ret = Rcpp::no_init_vector(ldmap_snp.size());
    std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[&](double x){
                                                                     return(static_cast<Rbyte>(Snp{.flt=x}.str.alt));
                                                                   });
    ret.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
    return ret;
  }else{
    Rcpp::StringVector ret = Rcpp::StringVector(ldmap_snp.size());
    std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[&](double x){
                                                                     return globalnucs.get_nuc(Nuc{Snp{.flt=x}.str.alt});
                                                                   });
    return ret;
   }
}




SEXP ldmap_snp_2_char_vector(Rcpp::NumericVector ldmap_snp){

  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");

  using namespace Rcpp;
  const size_t p=ldmap_snp.size();
  IntegerVector chrom = Rcpp::no_init(p);
  NumericVector pos = Rcpp::no_init(p);
  StringVector ascii_ref(p);
  StringVector ascii_alt(p);


  //  double *bs = reinterpret_cast<double*>(&mp);
  for(int i=0; i<p; i++){
    SNP sp{.snp={.flt=ldmap_snp(i)}};
    chrom(i)=sp.snp.str.chrom;
    pos(i)=static_cast<double>(sp.snp.str.pos);
    ascii_ref(i)=globalnucs.get_nuc(Nuc{sp.snp.str.ref});
    ascii_alt(i)=globalnucs.get_nuc(Nuc{sp.snp.str.alt});
  }
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["pos"]=pos,
                                _["ascii_ref"]=ascii_ref,
                                _["ascii_alt"]=ascii_alt);
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") = seq(1, p);
  return(dfl);
}






//' Convert ldmap snp back to dataframe
//'
//' @param ldmap_snp the vector of SNPs
//'
//' @export
//[[Rcpp::export]]
SEXP ldmap_snp_2_dataframe(Rcpp::NumericVector ldmap_snp,bool alleles_to_character=false){

  if(alleles_to_character){
    return(ldmap_snp_2_char_vector(ldmap_snp));
  }
  static_assert(sizeof(Snp)==sizeof(double),"packed structure is the size of a double");

  using namespace Rcpp;
  const size_t p=ldmap_snp.size();
  IntegerVector chrom(p);
  NumericVector pos(p);
  RawVector ascii_ref(p);
  RawVector ascii_alt(p);


  //  double *bs = reinterpret_cast<double*>(&mp);
  for(int i=0; i<p; i++){
    SNP sp{.snp={.flt=ldmap_snp(i)}};
    chrom(i)=sp.snp.str.chrom;
    pos(i)=static_cast<double>(sp.snp.str.pos);
    ascii_ref(i)=static_cast<Rbyte>(sp.snp.str.ref);
    ascii_alt(i)=static_cast<Rbyte>(sp.snp.str.alt);
  }
  ascii_ref.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
  ascii_alt.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["pos"]=pos,
                                _["ref"]=ascii_ref,
                                _["alt"]=ascii_alt);
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



  // if(!sorted_b){
  //   Rcpp::Rcerr<<"Assuming query is sorted (attribute 'sorted' not detected in `reference`)"<<std::endl;
  // }
  RcppParallel::RVector<double> input_a(query);
  RcppParallel::RVector<double> input_b(reference);

  const size_t q_size =query.size();
  Rcpp::IntegerVector ret_index(q_size,NA_INTEGER);
  Rcpp::IntegerVector ret_match(q_size,NA_INTEGER);
  Rcpp::NumericVector ret_ref(q_size);

  auto o_ref_b = input_b.begin();
  auto ref_b = input_b.begin();
  auto ref_e = input_b.end();

  for(int i=0; i<q_size; i++){
    double qf=query(i);
    SNP q{.snp={.flt=qf}};
    int rmi=std::numeric_limits<int>::max();
    auto lb = std::lower_bound(ref_b,ref_e,qf,[](double  a,double b){
                                               return (SNP{.snp={.flt=a}}<SNP{.snp={.flt=b}});
                                             });
    for(auto ilb = lb; ilb!=ref_e; ilb++){
      if(auto qm = q.allele_match(SNP{.snp={.flt=(*ilb)}})){
        int tqm=static_cast<int>(*qm);
        if(tqm<rmi){
          rmi=tqm;
          ret_index(i)=std::distance(o_ref_b,ilb)+1;
          ret_match(i)=tqm;
          ret_ref(i)=*ilb;
        }
      }else{
        ref_b=lb;
        break;
      }
    }
    if(ref_b==ref_e){
      break;
    }
  }
  ret_ref.attr("class")=reference.attr("class");
  using namespace Rcpp;
  std::vector<std::string> levs({"perfect_match",
                                 "reverse_match",                  // match -> reverse
                                 "complement_match",               // match -> complement
                                 "reverse_complement_match",       // match -> reverse
                                 "ambig_match",                    // match on ref and alt
                                 "reverse_ambig_match",            // match on ref and alt (ambig)
                                 "complement_ambig_match",         // match on ref and alt (ambig)
                                 "reverse_complement_ambig_match"});
  if(levs.size()!=8){
    Rcpp::stop("factor has incorrect number of levels in join_snp,  it's "+std::to_string(levs.size()));
  }
  ret_match.attr("levels") =wrap(levs);
  ret_match.attr("class") = Rcpp::StringVector::create("factor");
  auto dfl = List::create(_["index"]=ret_index,
                          _["match_type"]=ret_match,
                          _["match"]=ret_ref);
  if(rsid.size()==reference.size()){
    Rcpp::IntegerVector ret_rsid = no_init(q_size);
    std::transform(ret_index.begin(),ret_index.end(),ret_rsid.begin(),[&](const int i){
                                                                        if(i==NA_INTEGER)
                                                                          return NA_INTEGER;
                                                                        return rsid[i-1];
                                                                      });
    dfl["rsid"]=ret_rsid;
  }
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") = seq(1, q_size);
  return dfl;
}





//' Extract alternate allele from
//'
//' @param ref reference sequence (as obtained from a reference genome fasta file)
//' @param alleles_as_ambig IUPAC ambiguity codes representing alleles
//'
//' @export
//[[Rcpp::export]]
Rcpp::StringVector extract_alt(Rcpp::StringVector ref, Rcpp::StringVector alleles_as_ambig){

  const size_t p=ref.size();
  if(p!=alleles_as_ambig.size()){
    Rcpp::stop("ref and alleles must be the same size");
  }
  int i=0;
  Rcpp::StringVector ret = Rcpp::StringVector::import_transform(
                                                                alleles_as_ambig.begin(),
                                                                alleles_as_ambig.end(),
                                                                [&](const Rcpp::String el){
                                                                  const char* cst_el=el.get_cstring();
                                                                  const Rcpp::String ref_el=ref(i++);
                                                                  const char* cst_ref=ref_el.get_cstring();
                                                                  Nuc ret_nuc = Nuc{ascii2Nuc(cst_el)}^Nuc{ascii2Nuc(cst_ref)};
                                                                  return(globalnucs.get_nuc(ret_nuc));
                                                                });
  return ret;
}




//[[Rcpp::export]]
Rcpp::IntegerVector fast_str2int(Rcpp::StringVector input,int offset=0,const std::string prefix="",const int na_val=NA_INTEGER){

  SEXP str_sxp(input);
  const size_t p =input.size();
  Rcpp::IntegerVector ret(p);
  auto pp = get_string_ptr(str_sxp);
  offset = std::max(static_cast<size_t>(offset),prefix.size());
  const bool use_prefix = prefix.size()>0;
  int tresult;
  std::transform(pp,pp+p,ret.begin(),[&](SEXP x){
                                       if(x==R_NaString)
                                         return na_val;
                                       const size_t strl=LENGTH(x);
                                       const char* chp = CHAR(x);
                                       const char* beg=chp+offset;
                                       const char* end=chp+LENGTH( x );
                                       if(offset>=strl){
                                         // Rcpp::Rcerr<<"For string"<<CHAR(x)<<std::endl;
                                         // Rcpp::Rcerr<<"LEN(x): "<<LENGTH(x)<<std::endl;
                                         // Rcpp::Rcerr<<"string offset greater than string length"<<std::endl;
                                         return na_val;
                                       }
                                       if(use_prefix){
                                         if(std::string_view(chp,offset)!=prefix)
                                           return na_val;
                                       }


#if __has_include(<charconv>)
                                       if(auto [p, ec] = std::from_chars(beg, end, tresult);
                                          ec == std::errc())
                                         return tresult;
#else
                                       if(sscanf (beg,"%d",&tresult)>0)
                                         return(tresult);
#endif
                                       return (na_val);


                                     });
  return(ret);
}



//' Create ambiguity codes from two alleles
//'
//' @param A1 allele 1
//' @param A2 allele 2
//'
//' @export
//[[Rcpp::export]]
Rcpp::StringVector make_ambig(Rcpp::StringVector A1, Rcpp::StringVector A2){

  const size_t p=A1.size();
  if(p!=A2.size()){
    Rcpp::stop("A1 and A2 must be the same size");
  }

  static_assert(('A'_N | 'A'_N) == 'A'_N,"combination operator doesn't work!");
  static_assert(('T'_N | 'T'_N) == 'T'_N,"combination operator doesn't work!");
  static_assert(('G'_N | 'G'_N) == 'G'_N,"combination operator doesn't work!");
  static_assert(('A'_N | 'C'_N) == 'M'_N,"combination operator doesn't work!");
  static_assert(('A'_N | 'C'_N) == 'M'_N,"combination operator doesn't work!");
  int i=0;
  Rcpp::StringVector ret = Rcpp::StringVector::import_transform(
                                                                A1.begin(),
                                                                A1.end(),
                                                                [&](const Rcpp::String el){
                                                                  const char* cst_el=el.get_cstring();
                                                                  const Rcpp::String ref_el=A2(i++);
                                                                  const char* cst_ref=ref_el.get_cstring();
                                                                  Nuc ret_nuc = Nuc{ascii2Nuc(cst_el)}|Nuc{ascii2Nuc(cst_ref)};
                                                                  //Rcpp::Rcout<<Nuc2string(ret_nuc)<<std::endl;
                                                                  return(globalnucs.get_nuc(ret_nuc));
                                                                });
  return ret;
}





//' Formatting method for ldmap snps
//'
//' @param x a vector of ldmap_snps
//[[Rcpp::export]]
Rcpp::StringVector format_ldmap_snp(Rcpp::NumericVector x){

  return( Rcpp::StringVector::import_transform(x.begin(),x.end(),[&](double x){
                                                                   return globalnucs.string_rep(SNP{.snp={.flt=x}});
                                                                 }));
}
