


#include "ldmap/genetic_map.hpp"
#include "Rinternals.h"
#include <numeric>
#include "alleles.hpp"
#include <random>
#include <algorithm>
#include <iterator>
#include <limits>
#include <progress.hpp>
#include <RcppParallel.h>
#include <string>
#include <unordered_map>
//[[Rcpp::depends(RcppParallel)]]
// //[[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]

#include <boost/icl/split_interval_set.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/spirit/home/x3.hpp>
#include "Rcpp/Nullable.h"
#include "Rcpp/vector/instantiation.h"
#if __has_include("tbb/parallel_for.h")
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/concurrent_unordered_map.h"


#include "tbb/parallel_sort.h"
#endif
#include <type_traits>
#include <variant>


#if __has_include(<charconv>)
#include <charconv>
#else
#include <cstdio>
#endif



using vector_variant =
  std::variant<Rcpp::NumericVector, Rcpp::RawVector, Rcpp::IntegerVector,Rcpp::StringVector,Rcpp::LogicalVector,Rcpp::List>;

vector_variant get_vector_variant(SEXP x){

  Rcpp::RObject ro(x);
  auto attr = ro.attributeNames();
  switch (TYPEOF(x)) {
  case INTSXP:{
    return Rcpp::IntegerVector(x);
  }
  case RAWSXP:{
    return Rcpp::RawVector(x);
  }
  case REALSXP:{
    return Rcpp::NumericVector(x);
  }
  case STRSXP:{
    return Rcpp::StringVector(x);
  }
  case VECSXP:{
    return Rcpp::List(x);
  }
  case LGLSXP:{
    return Rcpp::LogicalVector(x);
  }
  default:
    Rcpp::stop("x must be a vector type in get_vector_variant (it also can't be a string type...) (this is likely not a user error)");
  }
}


template <typename A, typename B>
const B* indirect_index(const A* ab, const B* bb,const A* a_it){
  return bb+std::distance(ab,a_it);
}


class Indirector{
  size_t offset;
  size_t size;
public:
  Indirector(const size_t offset_,const size_t size_):offset(offset_),size(size_){}

  SEXP operator()(const Rcpp::StringVector input)const{
    const Rcpp::StringVector::const_iterator ib = input.begin();
    auto ret = Rcpp::StringVector(ib+offset,ib+offset+size);
    ret.attr("class")=input.attr("class");
    return ret;
  }
  template<int RTYPE>
  SEXP operator()(const Rcpp::Vector<RTYPE> input)const{
    const auto ib = input.begin();
    auto ret = Rcpp::Vector<RTYPE>(ib+offset,ib+offset+size);
    ret.attr("class")=input.attr("class");
    return ret;
  }
};



template <typename A>
SEXP indirect_index_df(const A* ab, const Rcpp::DataFrame df,const A* a_itb, const A* a_ite){


  auto dfc=df.size();
  const Indirector idr(std::distance(ab,a_itb),std::distance(a_itb,a_ite));
  std::vector<vector_variant> vv;
  vv.reserve(dfc);
  std::transform(df.begin(),df.end(),std::back_inserter(vv),[](auto elem){return get_vector_variant(elem);});
  Rcpp::List ret = Rcpp::List::import_transform(vv.begin(),vv.end(),[&](const vector_variant &elem){
                                                                      return std::visit(idr,elem);
                                                                    });
  ret.attr("names") = df.names();
  ret.attr("class") = Rcpp::StringVector::create("tbl_df","tbl","data.frame");
  ret.attr("row.names") = Rcpp::seq(1, std::distance(a_itb,a_ite));
  return ret;
}



template <typename Iterator>
bool parse_Region(Iterator first, Iterator last, Region &c){
  using boost::spirit::x3::int_;
  using boost::spirit::x3::_attr;
  using boost::spirit::x3::lit;
  using boost::spirit::x3::phrase_parse;
  using boost::spirit::x3::ascii::space;

  int chrom = 0;
  int start = 0;
  int end = 0;

  auto fchrom = [&](auto& ctx){ chrom = _attr(ctx); };
  auto fstart = [&](auto& ctx){ start = _attr(ctx); };
  auto fend = [&](auto& ctx){ end = _attr(ctx); };

  bool r = phrase_parse(first, last,
                       //  Begin grammar

                        lit("chr") >> int_[fchrom] >>':'>>int_[fstart]>>'_'>>int_[fend],
                        space);

 if (!r || first != last) // fail if we did not get a full match
   return false;
 c.br.str={.end=static_cast<uint64_t>(end),
           .start=static_cast<uint64_t>(start),
           .chrom=static_cast<unsigned char>(chrom)};
 return r;
}



//[[Rcpp::export]]
Rcpp::NumericVector parse_ldmap_range(Rcpp::StringVector input){

  Region reg;
  Rcpp::NumericVector ret=Rcpp::NumericVector::import_transform(
                                                                input.begin(),
                                                                input.end(),[&reg](SEXP inp){
                                                                              const size_t p= LENGTH(inp);
                                                                              const char* charp=CHAR(inp);
                                                                              std::string_view sv(charp,p);
                                                                              parse_Region(sv.begin(),sv.end(),reg);
                                                                              return(reg.br.flt);
                                                                            });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");
  return ret;

}



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


Rcpp::NumericVector snps_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end){

  const auto sas= Region::make_Region(x).start_SNP().snp.flt;
  const auto saf= Region::make_Region(x).last_SNP().snp.flt;
  auto xbl =std::lower_bound(begin,end,sas,[](const double& snp_a, const double &snp_b){
                                             return(SNP{.snp={.flt=snp_a}}<SNP{.snp={.flt=snp_b}});
                                           });
  auto xbu =std::upper_bound(xbl,end,saf,[](const double& snp_a, const double &snp_b){
                                           return(SNP{.snp={.flt=snp_a}}<SNP{.snp={.flt=snp_b}});
                                           });
  if(xbl==end){
    auto ret=Rcpp::NumericVector::create();
    ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
    return ret;
  }
  Rcpp::NumericVector ret(xbl,xbu);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  return ret;
}



SEXP snps_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end, Rcpp::DataFrame index_with){

  const auto sas= Region::make_Region(x).start_SNP().snp.flt;
  const auto saf= Region::make_Region(x).last_SNP().snp.flt;
  auto xbl =std::lower_bound(begin,end,sas,[](const double& snp_a, const double &snp_b){
                                             return(SNP{.snp={.flt=snp_a}}<SNP{.snp={.flt=snp_b}});
                                           });
  auto xbu =std::upper_bound(xbl,end,saf,[](const double& snp_a, const double &snp_b){
                                           return(SNP{.snp={.flt=snp_a}}<SNP{.snp={.flt=snp_b}});
                                           });
  if(xbl==end){
    auto ret=Rcpp::NumericVector::create();
    ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
    return ret;
  }
  //  Rcpp::List retdf=
  //  Rcpp::NumericVector ret(xbl,xbu);

  // retdf.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  return indirect_index_df<double>(begin,index_with,xbl,xbu);
}



int range_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end,const bool allow_overlap=false){

  // static_assert(Region{.br={.str={.end=1892607,.start=1,.chrom=1}}}<Region{.br={.str={.end=243199374,.start=1,.chrom=2}}},"I can't write relation operators...");
  const auto sx=Region::make_Region(x).start_SNP();
  auto xb=std::lower_bound(begin,end,sx,[](double  range_a,const SNP &reg_b){
                                          return(Region{.br={.flt=range_a}} < reg_b);
                                        });
  if( xb==end or Region{.br={.flt= (*xb)}} > sx)
    return NA_INTEGER;
  if(!allow_overlap){
    if(Region{.br={.flt=*xb}}.last_SNP()<Region::make_Region(x).last_SNP()){
      Rcpp::stop("range_in_range does not support overlap matches");
    }
  }
  return std::distance(begin,xb)+1;
}


int snp_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end){

  SNP sx{.snp={.flt=x}};
  auto xb=std::lower_bound(begin,end,sx,[](double  range_a,const SNP &snp_b){
                                          return(Region{.br={.flt=range_a}} < snp_b);
                                        });
  if( (xb==end ) or !(Region{.br={.flt=*xb}} ==sx))
    return NA_INTEGER;
  return std::distance(begin,xb)+1;
}




// class

// std::vector<int> snps_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end){
//   SNP sx{.snp={.flt=x}};
//   auto [xbb,xbe]=std::equal_range(begin,end,sx,[](double  range_a,const SNP &snp_b){
//                                           return(Region{.br={.flt=range_a}}< snp_b);
//                                         });
//   if( (xbe==end ) or !(Region{.br={.flt=*xbb}} ==sx))
//     return NA_INTEGER;
//   return std::distance(begin,xb)+1;
// }

//' Assign ranges to ranges
//'
//' @param ldmap_range_query vector of ldmap_ranges
//' @param ldmap_range_target vector of *non-overlapping* ldmap_ranges (must be sorted)
//' @param allow_overlap is it alright if a query is only partially inside the target?
//' @return a vector of integers of length `length(ldmap_range_query)` with the index of the `ldmap_range_target`
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector range_in_range(Rcpp::NumericVector ldmap_range_query,Rcpp::NumericVector ldmap_range_target,bool allow_overlap=false){

  const size_t p=ldmap_range_query.size();
  RcppParallel::RVector<double> input_range_query(ldmap_range_query);
  RcppParallel::RVector<double> input_range_target(ldmap_range_target);

  Rcpp::IntegerVector ret = Rcpp::no_init(p);
  RcppParallel::RVector<int> output_range(ret);
  auto irb = input_range_target.begin();
  auto ire = input_range_target.end();
  if(!std::is_sorted(irb,ire,[](const double a,const double b){
                               return Region::make_Region(a)<Region::make_Region(b);
                             }))
    Rcpp::stop("ldmap_range_query must be sorted");


  // static_assert(Region{.br={.str={.end=100ull,.start=1ull,.chrom=1}}} < Region{.br={.str={.end=100ull,.start=1ull,.chrom=1}}},"bed_range is having issues sorting");


  std::transform(
                 input_range_query.begin(),
                 input_range_query.end(),
                 output_range.begin(),
                 [=](const double x){
                   return(range_in_range(x,irb,ire,allow_overlap));
                 });
  return ret;

}




//' Assign SNPs to ranges
//'
//' @param ldmap_snp vector of ldmap_snps (must be sorted)
//' @param ldmap_range vector of non-overlapping ldmap_ranges (must be sorted)
//' @return a vector of integers of length `length(ldmap_snp)` with the index of the `ldmap_range`
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


Rcpp::List vec_from_bitset(tbb::concurrent_vector<std::pair<double,std::bitset<64> >> &bsp,const size_t num_ranges){
  const size_t p=bsp.size();
  //  Rcpp::IntegerVector ret = Rcpp::no_init(p);
  std::vector<int> ch(num_ranges+1);
  std::iota(ch.begin(),ch.end(),0);

  Rcpp::List outr = Rcpp::List::import_transform(ch.begin(),ch.end(),[&](int i){
                                                                       if(i>0){
                                                                         Rcpp::IntegerVector r(p);
                                                                         return(SEXP(r));
                                                                       }
                                                                       return SEXP(Rcpp::NumericVector(p));
                                                                     });
  std::vector<RcppParallel::RVector<int>> output_ranges;
  Rcpp::NumericVector ret_v = outr[0];
  ret_v.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  std::transform(outr.begin()+1,outr.end(),std::back_inserter(output_ranges),[](SEXP el){
                                                                               Rcpp::IntegerVector rel(el);
                                                                               return RcppParallel::RVector<int>(rel);
                                                                             });
  // auto start = tbb::make_zip_iterator(cnt0,bsp.cbegin());
  // auto end = tbb::make_zip_iterator(cntp,bsp.cend());

  // std::for_each(start,end,
  tbb::parallel_for(tbb::blocked_range<size_t>(0,p),
                    [&](const tbb::blocked_range<size_t> &r){
                      for(int i=r.begin(); i<r.end(); i++){
                        const auto &ir=bsp[i];
                        double sp = ir.first;
                        ret_v[i]=sp;
                        std::bitset<64> x=ir.second;
                        for(int j=0; j<num_ranges; j++){
                          if(x[j]){
                            output_ranges[j][i]=1;
                          }
                        }
                      }});
  return outr;
}


std::vector<RcppParallel::RVector<double>> prepare_ranges(Rcpp::ListOf<Rcpp::NumericVector> &ldmap_ranges){
  std::vector<RcppParallel::RVector<double>> input_ranges;
  auto rcmp = [](const double a,const double b){
                return Region::make_Region(a)<Region::make_Region(b);
              };
  std::transform(ldmap_ranges.begin(),ldmap_ranges.end(),
                 std::back_inserter(input_ranges),[&rcmp](SEXP el){
                                                    Rcpp::NumericVector rel(el);
                                                    RcppParallel::RVector<double> ret(rel);
                                                    if(!std::is_sorted(ret.begin(),ret.end(),rcmp)){
                                                      std::sort(ret.begin(),ret.end(),rcmp);
                                                    }
                                                    return ret;
                                                  });
  return input_ranges;
}

//' Assign SNPs to ranges
//'
//' @param ldmap_snp vector of ldmap_snps (must be sorted)
//' @param ldmap_range vector of non-overlapping ldmap_ranges (must be sorted)
//' @return a vector of integers of length `length(ldmap_snp)` with the index of the `ldmap_range`
//' @export
//[[Rcpp::export]]
Rcpp::List snp_in_ranges(Rcpp::NumericVector ldmap_snp,Rcpp::ListOf<Rcpp::NumericVector> ldmap_ranges){

  const size_t p=ldmap_snp.size();
  RcppParallel::RVector<double> input_snp(ldmap_snp);
  auto input_ranges = prepare_ranges(ldmap_ranges);

  const size_t num_ranges = input_ranges.size();

  if(num_ranges>64){
    Rcpp::stop("snp_in_ranges currently ony supports up to 64 ranges at a time");
  }
  tbb::concurrent_vector<std::pair<double,std::bitset<64>>> bsp;
  bsp.reserve(p/10);


  static_assert(bed_range{.str={.end=100,.start=50,.chrom=1}}.str.end == 100,"bed_range is having issues");
  // static_assert(bed_range{.str={.end=2147483647ull,.start=247344518ull,.chrom=1}}.str.end == 100,"bed_range is having issues");

  tbb::parallel_for(tbb::blocked_range<size_t>(0,p),
                    [&input_snp=std::as_const(input_snp),&bsp,&num_ranges,&input_ranges](const tbb::blocked_range<size_t> &r){
                      for(int i=r.begin(); i<r.end(); i++){
                        const double tx=input_snp[i];
                        std::bitset<64> bs;
                        for(int j=0; j<num_ranges; j++){
                          if( !Rcpp::IntegerVector::is_na(snp_in_range(tx,input_ranges[j].begin(),input_ranges[j].end()))){
                            bs[j]=true;
                          }
                        }
                        if(bs.any()){
                          bsp.emplace_back(std::make_pair(tx,bs));
                        }
                      }
                    });
  const size_t retp=bsp.size();
  Rcpp::StringVector new_names=ldmap_ranges.names();
  new_names.push_front("ldmap_snp");

  auto retlist = vec_from_bitset(bsp,num_ranges);
  retlist.names() = new_names;
  retlist.attr("class") = Rcpp::StringVector::create("tbl_df","tbl","data.frame");
  retlist.attr("row.names") = Rcpp::seq(1, retp);
  return retlist;

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



//' Assign SNPs to ranges
//'
//' @param ldmap_snp vector of ldmap_snps (must be sorted)
//' @param ldmap_range vector of potentially overlapping ldmap_ranges (must be sorted)
//' @return a list of integer vectors giving the ranges to which each SNP belongs
//' @export
//[[Rcpp::export]]
Rcpp::List match_ranges_snps(Rcpp::DataFrame df,Rcpp::NumericVector ldmap_range,const std::string snp_col="snp_struct"){


  Rcpp::NumericVector ldmap_snp=df[snp_col];
  const size_t p=ldmap_range.size();
  RcppParallel::RVector<double> input_snp(ldmap_snp);
  RcppParallel::RVector<double> input_range(ldmap_range);


  auto irb = input_snp.begin();
  auto ire = input_snp.end();
  static_assert(bed_range{.str={.end=100,.start=50,.chrom=1}}.str.end == 100,"bed_range is having issues");
  // static_assert(bed_range{.str={.end=2147483647ull,.start=247344518ull,.chrom=1}}.str.end == 100,"bed_range is having issues");


  Rcpp::List ret = Rcpp::List::import_transform(
                 input_range.begin(),
                 input_range.end(),
                 [=](const double x){
                   return(snps_in_range(x,irb,ire,df));
                 });
  ret.names() = format_ldmap_range(ldmap_range);
  return ret;
}



template<typename T>
auto find_window(const T* begin,const T* end,const  T *query,const T diff){

  auto ld = std::lower_bound(begin,query,*query-diff);
  auto ud = std::upper_bound(query,end,*query+diff);
  // if(ud==end){
  //   ud=end-1;
  // }
  return std::make_pair(ld,ud-1);
}






//' Create overlapping regions based on monotonic, point-level annotation
//'
//' @param ldmap_snp a (sorted) ldmap_snp vector (`length(ldmap_snp)` is referred to  as  `p`)
//' @param cm a numeric vector of (length `p`) per-snp annotations (e.g cumulative recombination rate)
//' @param window the window width.
//'
//' @return a vector of `ldmap_range`s of length `p` giving the window for each SNP.  The width of the window
//' is defined for a target snp `ldmap_snp[i]`, as having the chromosome
//' from `ldmap_snp[i]` and including the position of all `ldmap_snp[j]` snps such that `abs(cm[i]-cm[j])<window` for all values of `j`
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector window_ldmap_range(Rcpp::NumericVector ldmap_snp,Rcpp::NumericVector cm,const double window=1.0){

  const size_t p=ldmap_snp.size();
  Rcpp::NumericVector ret(p);
  double *ldmb = ldmap_snp.begin();
  double *cmb = cm.begin();
  double *cme = cm.end();

  for(int i=0; i<p; i++){
    auto [fw_l, fw_u] =
        find_window(cmb, cme, cmb + i, window/2);
    auto lb = indirect_index(cmb,ldmb,fw_l);
    auto ub = indirect_index(cmb,ldmb,fw_u);
    ret(i)=Region::make_Region(SNP::make_snp(*lb),SNP::make_snp(*ub)).br.flt;
  }
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");
  return ret;
}


template<typename T>
std::vector<T> vec_concatenate_sort(T* beg1,T* end1,T* beg2,T* end2){
  const size_t p=std::distance(beg1,end1);
  const size_t q=std::distance(beg2,end2);
  std::vector<T> ret(p+q);
  auto compab= [](const double &a,const double &b){
                                    return Region{.br={.flt=a}}<Region{.br={.flt=b}};
               };
  auto oit=std::partial_sort_copy(beg1,end1,ret.begin(),ret.end(),compab);
  std::partial_sort_copy(beg2,end2,oit,ret.end(),compab);
  std::inplace_merge(ret.begin(),oit,ret.end(),compab);
  return ret;
}


template<typename T>
Rcpp::NumericVector split_ldmap_ranges_set(T* xb, T* xe){
  //  using interval_t = decltype(Region{.br={.dat=0}}.interval());
  std::array<boost::icl::split_interval_set<uint64_t>,23> sets;

  std::for_each(xb,xe,[&sets](const double& tx){
                                    auto reg = Region::make_Region(tx);
                                    auto chr = reg.chrom()-1;
                                    auto &set_it =sets[chr];
                                    set_it.insert(reg.interval());
                                  });

  size_t ret_size = 0;
  for(int i=0; i<23; i++){
    ret_size+=sets[i].iterative_size();
  }

  int i=0;
  Rcpp::NumericVector ret=Rcpp::no_init(ret_size);
  for(unsigned char chrom=1; chrom<24; chrom++){
    const auto& set_el = sets[chrom-1];
    for(auto &el : set_el){
      ret[i++]=Region{.br={.str={.end=el.upper(),
                                 .start=el.lower(),
                                 .chrom=chrom}}}.br.flt;
    }
  }
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");
  return ret;
}


//' Take a vector of (preferably) sorted and (possibly) overlapping ldmap_ranges and create a new range of (sorted) non-overlapping ldmap_ranges
//'
//' @param x a (preferably sorted) ldmap_range vector (`length(x)` is referred to  as  `p`)
//'
//' @export
//' @return a sorted vector ldmap_ranges of length at least `p` and at most `2p`(?) representing the same intervals
//[[Rcpp::export]]
Rcpp::NumericVector split_ldmap_range_overlap(Rcpp::NumericVector x){
  return(split_ldmap_ranges_set(x.begin(),x.end()));
}


//' Merge two ldmap_range vectors
//'
//' @param x a (preferably sorted) ldmap_range vector (`length(x)` is referred to  as  `p`)
//' @param y a (preferably sorted) ldmap_range vector (`length(y)` is referred to  as  `q`)
//'
//' @return a sorted vector ldmap_ranges of length at most `p+q` and at least `max(p,q)` representing the union of the two sets of ranges
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector merge_ldmap_ranges(Rcpp::NumericVector x,Rcpp::NumericVector y){
  const size_t p=x.size();
  const size_t q=y.size();
  double *xb = x.begin();

  auto ret=vec_concatenate_sort(x.begin(),x.end(),y.begin(),y.end());
  auto match_ab= [](const double &a,const double &b){
                   return Region{.br={.flt=a}}.overlap(Region{.br={.flt=b}});
                 };
  auto itb=ret.begin();
  auto match_r = std::adjacent_find(itb,ret.end(),match_ab);
  if(match_r!=ret.end()){
    return(split_ldmap_ranges_set(&(*ret.begin()),&(*ret.end())));
  }
  Rcpp::NumericVector rret(ret.begin(),ret.end());
  rret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");

  return rret;
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
  std::vector<int> ch(23);
  std::iota(ch.begin(),ch.end(),1);
  Rcpp::StringVector levs=Rcpp::StringVector::import_transform(ch.begin(),ch.end(),[](int i){
                                                                                     if(i==23)
                                                                                       return std::string("chrX");
                                                                                     return std::string("chr"+std::to_string(i));
                                                                                   });
  chrom.attr("levels") =levs;
  chrom.attr("class") = Rcpp::StringVector::create("factor");
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



//[[Rcpp::export]]
Rcpp::IntegerVector sample_interval(Rcpp::IntegerVector n,Rcpp::IntegerVector begin,Rcpp::IntegerVector end,const bool replace =false){


  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  const auto ns = n.size();
  const auto bs = begin.size();
  const auto es = end.size();
  const auto p = std::max(std::max(ns,bs),es);
  size_t ret_size=0;

  for(int i=0; i < p; i++){
    auto mn= n(i % ns);
    ret_size +=mn;
  }

  Rcpp::IntegerVector ret = Rcpp::no_init(ret_size);
  auto ret_b = ret.begin();

  if(replace){
    for(int i=0 ; i<p ; i++){
      auto uid = std::uniform_int_distribution<>(begin[i%bs], end(i%es));
      ret_b=std::generate_n(ret_b,n(i%ns),[&](){return uid(gen);});
    }
  }else{
    for(int i=0 ; i<p ; i++){
      auto nret_b = std::sample(boost::counting_iterator<int>(begin[i%bs]),
                                boost::counting_iterator<int>(end[i%es]),
                                ret_b,n(i%ns),gen);
      if(std::distance(ret_b,nret_b) < n(i%ns)){
        Rcpp::stop("n cannot be  larger than region size, when replace=FALSE");
      }
      ret_b=nret_b;
    }
  }

  return(ret);
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
//' @param struct_vec the vector of SNPs (or ldmap_ranges)
//'
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector chromosomes(Rcpp::NumericVector struct_vec){
  using namespace Rcpp;
  IntegerVector ret = no_init(struct_vec.size());
  if(struct_vec.inherits("ldmap_snp"))
    std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](double x){
      return(static_cast<int>(Snp{.flt=x}.str.chrom));
    });
  else
    std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](double x){
      return(static_cast<int>(bed_range{.flt=x}.str.chrom));
    });
  return ret;
}


//' get starting position from a ldmap_range
//'
//' @param ldmap_range the vector of ldmap_ranges
//'
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector starts(Rcpp::NumericVector ldmap_range){
  using namespace Rcpp;
  IntegerVector ret = no_init(ldmap_range.size());
  if(ldmap_range.inherits("ldmap_snp"))
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
      return(static_cast<int>(bed_range{.flt=x}.str.start));
    });
  else
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
      return(static_cast<int>(bed_range{.flt=x}.str.start));
    });
  return ret;
}



//' get end position from a ldmap_range
//'
//' @param ldmap_range the vector of ldmap_ranges
//'
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector ends(Rcpp::NumericVector ldmap_range){
  using namespace Rcpp;
  IntegerVector ret = no_init(ldmap_range.size());
  if(ldmap_range.inherits("ldmap_snp"))
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
      return(static_cast<int>(bed_range{.flt=x}.str.end));
    });
  else
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
      return(static_cast<int>(bed_range{.flt=x}.str.end));
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
