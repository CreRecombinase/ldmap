#include <numeric>
#include "alleles.hpp"
#include "ldmap.hpp"
#include <Rcpp.h>
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
#include <range/v3/core.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/functional/invoke.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/traits.hpp>
#include <range/v3/utility/static_const.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/cycle.hpp>
#include <range/v3/view/take.hpp>


#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/sample.hpp>
#include <range/v3/view/for_each.hpp>
#include <boost/icl/split_interval_set.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/spirit/home/x3.hpp>

#if __has_include("tbb/parallel_for.h")
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/concurrent_unordered_map.h"
#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/transform.hpp>

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



class Indirector{
  size_t offset;
  size_t size;
public:
  Indirector(const size_t offset_,const size_t size_):offset(offset_),size(size_){}

  SEXP operator()(const Rcpp::StringVector input)const{
    const Rcpp::StringVector::const_iterator ib = input.begin();
    if(size>0){
      auto ret = Rcpp::StringVector(ib+offset,ib+offset+size);
      ret.attr("class")=input.attr("class");
      return ret;
    } 
    auto ret = Rcpp::StringVector();
    ret.attr("class")=input.attr("class");
    return ret;
  }
  template<int RTYPE>
  SEXP operator()(const Rcpp::Vector<RTYPE> input)const{
    if(size>0){
      const auto ib = input.begin();
      auto ret = Rcpp::Vector<RTYPE>(ib+offset,ib+offset+size);
      ret.attr("class")=input.attr("class");
      return ret;
    }
    auto ret = Rcpp::Vector<RTYPE>();
    ret.attr("class")=input.attr("class");
    return ret;
  }
};


template <typename A>
SEXP indirect_index_df(const A* ab, const Rcpp::DataFrame df,const A* a_itb, const A* a_ite){

  auto dfc=df.size(); 
  std::vector<vector_variant> vv;
  vv.reserve(dfc);  
  std::transform(df.begin(),df.end(),std::back_inserter(vv),[](auto elem){return get_vector_variant(elem);});
  const Indirector idr(std::distance(ab,a_itb),std::distance(a_itb,a_ite)); 
  Rcpp::List ret = Rcpp::List::import_transform(vv.begin(),vv.end(),[&](const vector_variant &elem){
    return std::visit(idr,elem);
  });
  ret.attr("names") = df.names();
  ret.attr("class") = Rcpp::StringVector::create("tbl_df","tbl","data.frame");
  ret.attr("row.names") = Rcpp::seq_len(std::distance(a_itb,a_ite));
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
 c.br=make_ldmap_region(chrom,start,end);
 return r;
}




template <typename Iterator>
bool parse_SNP(Iterator first, Iterator last, SNP &c){
  using boost::spirit::x3::int_;
  using boost::spirit::x3::char_;
  using boost::spirit::x3::_attr;
  using boost::spirit::x3::lit;
  using boost::spirit::x3::phrase_parse;
  using boost::spirit::x3::repeat;
  using boost::spirit::x3::ascii::space;

  int chrom = 0;
  int pos  = 0;
  char ref = '\0';
  char alt = '\0';

  auto fchrom = [&](auto& ctx){ chrom = _attr(ctx); };
  auto fpos = [&](auto& ctx){ pos = _attr(ctx); };
  auto fref = [&](auto& ctx){ ref = _attr(ctx); };
  auto falt = [&](auto& ctx){ alt = _attr(ctx); };

  bool r = parse(first, last,
                 repeat(0,1)[lit("chr")] >> int_[fchrom] >>char_(":_")>>int_[fpos]>>char_(":_")>>char_("ACGTMRWSYKVHDBN")[fref]>>char_(":_")>>char_("ACGTMRWSYKVHDBN")[falt]);

  if (!r || first != last) // fail if we did not get a full match
    return false;
  c=SNP::make_snp<false>(static_cast<unsigned char>(chrom),
                         static_cast<uint64_t>(pos),
                         Nuc{ascii2Nuc(ref)},
                         Nuc{ascii2Nuc(alt)});
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
                                                                              if(parse_Region(sv.begin(),sv.end(),reg))
                                                                                return(bit_cast<double>(reg.br));
                                                                              return NA_REAL;
                                                                            });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");
  return ret;
}

//[[Rcpp::export]]
Rcpp::NumericVector parse_ldmap_SNP(Rcpp::StringVector input){

  SNP snp;

  Rcpp::NumericVector ret=Rcpp::NumericVector::import_transform(
                                                                input.begin(),
                                                                input.end(),[&snp](SEXP inp){
                                                                              const size_t p= LENGTH(inp);
                                                                              const char* charp=CHAR(inp);
                                                                              std::string_view sv(charp,p);
                                                                              if(parse_SNP(sv.begin(),sv.end(),snp))
                                                                                return(bit_cast<double>(snp.snp));
                                                                              return NA_REAL;
                                                                            });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  return ret;
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

  const size_t csize = chrom.size();
  const size_t ssize = start.size();
  const size_t esize = end.size();
  for(int i=0; i<p; i++){
    ret(i)=bit_cast<double>(make_ldmap_region(chrom(i%csize),start(i%ssize),end(i%esize)));
  }
  return ret;
}


Rcpp::NumericVector snps_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end){

  const auto sas= bit_cast<double>(Region::make_Region(bit_cast<uint64_t>(x)).start_SNP().snp);
  const auto saf= bit_cast<double>(Region::make_Region(bit_cast<uint64_t>(x)).last_SNP().snp);
  auto xbl =std::lower_bound(begin,end,sas,[](const double& snp_a, const double &snp_b){
                                             return(SNP::make_snp(bit_cast<uint64_t>(snp_a))<SNP::make_snp(bit_cast<uint64_t>(snp_b)));
                                           });
  auto xbu =std::upper_bound(xbl,end,saf,[](const double& snp_a, const double &snp_b){
                                           return(SNP::make_snp(bit_cast<uint64_t>(snp_a))<SNP::make_snp(bit_cast<uint64_t>(snp_b)));
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

  const auto sas = bit_cast<double>(Region::make_Region(bit_cast<uint64_t>(x)).start_SNP().snp);
  const auto saf = bit_cast<double>(Region::make_Region(bit_cast<uint64_t>(x)).last_SNP().snp);
  auto cf = [](const double& snp_a, const double &snp_b){
              return(SNP::make_snp(bit_cast<uint64_t>(snp_a))<SNP::make_snp(bit_cast<uint64_t>(snp_b)));
            };
  auto xbl =std::lower_bound(begin,end,sas,cf);
  auto xbu =std::upper_bound(xbl,end,saf,cf);
  return indirect_index_df<double>(begin,index_with,xbl,xbu);
}



int snp_in_range(const double x, RcppParallel::RVector<double>::const_iterator begin,RcppParallel::RVector<double>::const_iterator end){

  auto sx{SNP::make_snp(bit_cast<uint64_t>(x))};
  auto xb=std::lower_bound(begin,end,sx,[](double  range_a,const SNP &snp_b){
                                          return(Region::make_Region(bit_cast<uint64_t>(range_a)) < snp_b);
                                        });
  if( (xb==end ) or !(Region::make_Region(bit_cast<uint64_t>(*xb)).overlap(sx)))
    return NA_INTEGER;
  return std::distance(begin,xb)+1;
}






//' Quickly calculate distance between to ldmap_ranges
//'
//' @param query vector of query ldmap_snps
//' @param target vector of target ldmap_ranges (must be sorted)
//' @return a vector of integers with the length between the two ranges
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector distance_rr(Rcpp::NumericVector query,Rcpp::NumericVector target){

  using namespace ranges;
  const size_t p=query.size();
  const size_t target_p=target.size();
  RcppParallel::RVector<double> input_range_query(query);
  RcppParallel::RVector<double> input_range_target(target);
  Rcpp::IntegerVector retvec = Rcpp::no_init(p);
  RcppParallel::RVector<int> output_range(retvec);
  std::transform(input_range_query.begin(),input_range_query.end(),input_range_target.begin(),output_range.begin(),[](const double mpa,const double mpb){
    return Region::make_Region(mpa).distance(Region::make_Region(mpb));
  });

  return retvec;
}

//' Assign ranges to nearest ranges
//'
//' @param query vector of query ldmap_snps
//' @param target vector of target ldmap_ranges (must be sorted)
//' @param use_begin logical scalar indicating whether to use the start of the target range(TRUE), the end(FALSE), or them min distance of the two (NA_LOGICAL)
//' @return a vector of integers of length `length(ldmap_range_query)` with the index of the `ldmap_range_target` (or NA_INTEGER if there is no overlap in the set of chromosomes)
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector nearest_snp_range(Rcpp::NumericVector query,Rcpp::NumericVector target){


  using namespace ranges;
  const size_t p=query.size();
  const size_t target_p=target.size();
  RcppParallel::RVector<double> input_range_query(query);
  RcppParallel::RVector<double> input_range_target(target);
  Rcpp::IntegerVector ret = Rcpp::no_init(p);
  RcppParallel::RVector<int> output_v(ret);
  int* op = output_v.begin();
  //  auto output_range = view::all(op);
  //  CPP_assert(input_or_output_iterator<decltype(begin(output_range))>);
  double* irtb = target.begin();
  double* irte = target.end();

  double* irqb = query.begin();
  double* irqe = query.end();
  auto msfun = [](double x){
                 return(SNP::make_snp(x));
               };
  auto mrfun = [](double x){
                 return(Region::make_Region(x));
               };


  using namespace Rcpp;

  // auto irtr = subrange(irtb,irte);
  // auto irt= views::transform(irtr,mrfun);
  // auto irq = subrange(irqb,irqe) | view::transform(msfun);


  if(!is_sorted(input_range_target.begin(),input_range_target.end(),[](double &ra, double rb){
                                                                      return Region::make_Region(ra)<Region::make_Region(rb);
                                                                    }))
    Rcpp::stop("target must be sorted");

  Rcpp::IntegerVector result =Rcpp::IntegerVector::import_transform(query.begin(),query.end(),[&](double x){
                                                                                                return(nearest_snp(SNP::make_snp(x),input_range_target.begin(),input_range_target.end()));
                                                                                              });
  return ret;

}






//' Assign ranges to ranges
//'
//' @param ldmap_range_query vector of ldmap_ranges
//' @param ldmap_range_target vector of *non-overlapping* ldmap_ranges (must be sorted)
//' @param allow_overlap is it alright if a query is only partially inside the target?
//' @return a vector of integers of length `length(ldmap_range_query)` with the index of the `ldmap_range_target`
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector range_in_range(Rcpp::NumericVector ldmap_range_query,Rcpp::NumericVector ldmap_range_target,bool allow_overlap=false){


  using namespace ranges;

  const size_t p=ldmap_range_query.size();
  const size_t target_p=ldmap_range_target.size();
  RcppParallel::RVector<double> input_range_query(ldmap_range_query);

  RcppParallel::RVector<double> input_range_target(ldmap_range_target);

  Rcpp::IntegerVector ret = Rcpp::no_init(p);
  RcppParallel::RVector<int> output_v(ret);
  int* op = output_v.begin();
  //  auto output_range = view::all(op);
  //  CPP_assert(input_or_output_iterator<decltype(begin(output_range))>);
  double* irtb = input_range_target.begin();
  double* irte = input_range_target.end();

  double* irqb = input_range_query.begin();
  double* irqe = input_range_query.end();
  auto mrfun = [](double x){
                 return(Region::make_Region(x));
               };


  auto irtr = subrange(irtb,irte);
  auto irt= view::transform(irtr,mrfun);
  auto irq = subrange(irqb,irqe) | view::transform(mrfun);


  if(!is_sorted(irt,[](Region ra, Region rb){
                      return ra<rb;
  }))
    Rcpp::stop("ldmap_range_query must be sorted");

  if(allow_overlap)
    auto result_thing = transform(irq,op,
                   [&](Region x){
                     return(range_overlap_range(x,irt));
                   });
  else
    auto result_thing = transform(irq,op,
                   [&](Region x){
                     return(range_within_range(x,irt));
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

  // #if __has_include("tbb/parallel_for.h")
  // tbb::parallel_for(tbb::blocked_range<size_t>(0,p),
  //                   [&input_snp=std::as_const(input_snp),&output_range,&irb=std::as_const(irb),&ire=std::as_const(ire)](const tbb::blocked_range<size_t> &r){
  //                     std::transform(
  //                                    input_snp.begin()+r.begin(),
  //                                    input_snp.begin()+r.end(),
  //                                    output_range.begin()+r.begin(),
  //                                    [=](const double x){
  //                                      return(snp_in_range(x,irb,ire));
  //                                    });
  //                   });
  // #else
  std::transform(
                 input_snp.begin(),
                 input_snp.end(),
                 output_range.begin(),
                 [=](const double x){
                   return(snp_in_range(x,irb,ire));
                 });
  //#endif
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
                return Region::make_Region(bit_cast<uint64_t>(a))<Region::make_Region(bit_cast<uint64_t>(b));
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
  retlist.attr("row.names") = Rcpp::seq_len(retp);
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
    if((end==0) or (start==0) or (chrom==0) or (chrom==31 and start==528482304 and end==1954)){
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
  if(ud==begin){
    ud=ud+1;
  }
  return std::make_pair(std::distance(begin,ld),std::distance(begin,ud-1));
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
  Rcpp::NumericVector ret=Rcpp::no_init(p);
  double *ldmb = ldmap_snp.begin();
  double *cmb = cm.begin();
  double *cme = cm.end();

  for(int i=0; i<p; i++){
    auto snpit=SNP::make_snp(*(ldmb+i));
    auto cmbi=cmb+i;
    auto [fw_l, fw_u] =
        find_window(cmb, cme, cmbi, window/2);

    auto lb = SNP::make_snp(*(ldmb+fw_l));
    auto ub = SNP::make_snp(*(ldmb+fw_u));
    if( (lb>snpit) or (snpit>ub)){
      Rcpp::Rcerr<<"Error: "<<lb<<"<"<<snpit<<"<"<<ub<<" is violated"<<std::endl;
      Rcpp::Rcerr<<"Window is "<<*(cmb+fw_l)<<"<="<<*cmbi<<"<="<<*(cmb+fw_u)<<std::endl;
      Rcpp::stop("Error in imp of window_ldmap_range");
    }

    auto ret_el = Region::make_Region(lb,ub);
    if(ret_el.start_SNP()!=lb or ret_el.last_SNP()!=ub){
      Rcpp::Rcerr<<"Error: "<<lb<<"<"<<snpit<<"<<"<<ub<<" is violated"<<std::endl;
      Rcpp::Rcerr<<"Window is "<<*(cmb+fw_l)<<"<="<<*cmbi<<"<="<<*(cmb+fw_u)<<std::endl;
      Rcpp::Rcerr<<"Region is "<<ret_el<<"("<<ret_el.start_SNP()<<"-"<<ret_el.last_SNP()<<")"<<std::endl;
      Rcpp::stop("Error in imp of window_ldmap_range");
    }
    ret(i)=bit_cast<double>(ret_el.br);
  }
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr");
  return ret;
}


template<typename T>
std::vector<T> vec_concatenate_sort(T* beg1,T* end1,T* beg2,T* end2){
  const size_t p=std::distance(beg1,end1);
  const size_t q=std::distance(beg2,end2);
  std::vector<T> ret(p+q);
  auto compab = [](const double a,const double b){
                return Region::make_Region(bit_cast<uint64_t>(a))<Region::make_Region(bit_cast<uint64_t>(b));
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
                        auto reg = Region::make_Region(bit_cast<uint64_t>(tx));
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
      ret[i++]=bit_cast<double>(make_ldmap_region(chrom,el.lower(),el.upper()));
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
                   return Region::make_Region(bit_cast<uint64_t>(a)).overlap(Region::make_Region(bit_cast<uint64_t>(b)));
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
  NumericVector end =   Rcpp::no_init(p);

  static_assert(get_end<int>(set_end(0, 0)) == 0);
  static_assert(static_cast<double>(get_end<int>(make_ldmap_region(1,10583,1892607)))==1892607.0);
  static_assert(get_end<int>(288236057858466047)==1892607);
  //  double *bs = reinterpret_cast<double*>(&mp);
  for(int i=0; i<p; i++){
    double ldi=ldmap_range(i);
    uint64_t ldui = bit_cast<uint64_t>(ldi);
    const Region sp=Region::make_Region(ldui);
    chrom(i)= sp.chrom();
    start(i)= static_cast<double>(sp.start());
    end(i)=static_cast<double>(sp.end());
  }
  std::vector<int> ch(24);
  std::iota(ch.begin(),ch.end(),1);
  Rcpp::StringVector levs=Rcpp::StringVector::import_transform(ch.begin(),ch.end(),[](int i){
                                                                                     if(i==23)
                                                                                       return std::string("chrX");
                                                                                     if(i==24)
                                                                                       return std::string("chrY");
                                                                                     return std::string("chr"+std::to_string(i));
                                                                                   });
  chrom.attr("levels") =levs;
  chrom.attr("class") = Rcpp::StringVector::create("factor");
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["start"]=start,
                                _["end"]=end);

  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") =   seq_len( p);
  
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
  if(chrom.size()!=pos.size())
    Rcpp::stop("length(chrom) != length(pos) in new_ldmap_snp");
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
                                    Rcpp::NumericVector pos,
                                    const bool NA2N){

  const size_t p=chrom.size();
  if(chrom.size()!=pos.size())
    Rcpp::stop("length(chrom) != length(pos) in new_ldmap_snp");
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
Rcpp::IntegerVector sample_interval(Rcpp::IntegerVector n,Rcpp::IntegerVector beginv,Rcpp::IntegerVector endv, const bool replace =false){

  using namespace ranges;

  // std::random_device rd;  //Will be used to obtain a seed for the random number engine
  // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  const auto ns = n.size();
  RcppParallel::RVector<int> beginpv(beginv);
  RcppParallel::RVector<int> endpv(endv);
  RcppParallel::RVector<int> npv(n);

  const auto bs = beginv.size();
  const auto es = endv.size();

  auto begin_range = subrange(beginpv.begin(),beginpv.end());
  auto end_range = subrange(endpv.begin(),endpv.end());
  auto n_range = subrange(npv.begin(),npv.end());

  const auto p = size(n_range);

  size_t ret_size=0;

  for(int i=0; i < p; i++){
    auto mn= n(i % ns);
    ret_size +=mn;
  }
  Rcpp::IntegerVector ret=Rcpp::no_init(ret_size);
  auto rd = std::random_device{}();
  int* outip=ret.begin();
  for(int i=0; i<p; i++){
    int tn=n[i];
    auto range_v = views::ints(beginv[i],endv[i]);
    outip=std::sample(begin(range_v),end(range_v),outip,tn,std::mt19937{rd});
  }
  // auto ret_b = ret.begin();

  // if(replace){
  //   for(int i=0 ; i<p ; i++){
  //     if(i % 100 == 0)
  //       Rcpp::checkUserInterrupt();
  //     auto uid = std::uniform_int_distribution<>(beginv[i%bs], endv(i%es));
  //     ret_b=std::generate_n(ret_b,n(i%ns),[&](){return uid(gen);});
  //   }
  // }else{
  //   for(int i=0 ; i<p ; i++){
  //     if(i % 100 == 0)
  //       Rcpp::checkUserInterrupt();
  //     auto views::ints(beginv[i%bs],endv[
  //     auto nret_b = std::sample(boost::counting_iterator<int>(beginv[i%bs]),
  //                               boost::counting_iterator<int>(end[i%es]),
  //                               ret_b,n(i%ns),gen);
  //     if(std::distance(ret_b,nret_b) < n(i%ns)){
  //       Rcpp::stop("n cannot be  larger than region size, when replace=FALSE");
  //     }
  //     ret_b=nret_b;
  //   }
  // }

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
                                                                                         return SNP::make_snp(bit_cast<uint64_t>(x)).is_strand_ambiguous();
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
  IntegerVector ret=  seq_len( p);
  if(ret.size()!=p){
    Rcpp::stop("you can't make seqs right");
  }
  RcppParallel::RVector<double> input_a(struct_vec);

  std::sort(ret.begin(),ret.end(),
            [&](const int i,const int j){
              return(SNP::make_snp(bit_cast<uint64_t>(input_a[i-1])) <
                     SNP::make_snp(bit_cast<uint64_t>(input_a[j-1])));
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
  // std::vector<int> ret=  seq_len( p);
  std::vector<size_t> idx(p);
  std::iota(idx.begin(), idx.end(), 0);


  std::sort(idx.begin(),idx.end(),[&](const int i,const int j){
                                    return SNP::make_snp(bit_cast<uint64_t>(struct_vec[i-1])) <
                                      SNP::make_snp(bit_cast<uint64_t>(struct_vec[j-1]));
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
  if(!(struct_vec.inherits("ldmap_snp") || struct_vec.inherits("ldmap_range"))){
    Rcpp::stop("must inherit from ldmap_snp or ldmap_range");
  }
  if(struct_vec.inherits("ldmap_snp"))
    std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](double x){
                                                                     auto ret = static_cast<int>(SNP::make_snp(x).chrom());
                                                                     return(ret==0?NA_INTEGER:ret);
                                                                   });
  else
    std::transform(struct_vec.begin(),struct_vec.end(),ret.begin(),[](double x){
                                                                     auto ret = static_cast<int>(Region::make_Region(x).chrom());
                                                                     return(ret==0?NA_INTEGER:ret);
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
  if(!(ldmap_range.inherits("ldmap_snp") || ldmap_range.inherits("ldmap_range"))){
    Rcpp::stop("must inherit from ldmap_snp or ldmap_range");
  }
  if(!ldmap_range.inherits("ldmap_snp"))
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
                                                                       auto ret = static_cast<int>(Region::make_Region(x).start());
                                                                       return(ret ==0 ? NA_INTEGER : ret);
                                                                     });
  else
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
                                                                       auto ret = static_cast<int>(SNP::make_snp(x).pos());
                                                                       return(ret ==0 ? NA_INTEGER : ret);
                                                                     });
  return ret;
}



//' Find convex hull of vector of ranges (or SNPs )
//'
//' @param vector of ldmap_range or ldmap_snp of 
//' @return an ldmap range containg all the ranges (or snps)
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector convex_hull(Rcpp::NumericVector x){
  using namespace Rcpp;
  Rcpp::NumericVector ret(1);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_range","vctrs_vctr"); 
  if(x.inherits("ldmap_snp")){  
    auto [min_pt,max_pt] = std::minmax_element(x.begin(),x.end(),
                                               [](const double a,const double b){
                                                 auto sa=SNP::make_snp(bit_cast<uint64_t>(a));
                                                 auto sb=SNP::make_snp(bit_cast<uint64_t>(b));
                                                 if(sa.chrom()!=sb.chrom()){
                                                   Rcpp::stop("all elements must share chromosomes for `convex_hull`");
                                                 }
                                                 return(sa<sb);
                                               });
    auto asnp=SNP::make_snp(bit_cast<uint64_t>(*max_pt));
    auto bsnp=SNP::make_snp(bit_cast<uint64_t>(*min_pt));
    
    double retflt=bit_cast<double>(make_ldmap_region(asnp.chrom(),bsnp.pos(),asnp.pos()+1));
    ret[0]=retflt;
    return(ret);
  }
  if(x.inherits("ldmap_range")){
    auto start = x.begin();
    auto retflt = std::accumulate(x.begin(), x.end(), *start,
                                  [](const double a, const double b){
                                    auto sa=Region::make_Region(bit_cast<uint64_t>(a));
                                    auto sb=Region::make_Region(bit_cast<uint64_t>(b));
                                    if(sa.chrom()!=sb.chrom()){
                                      Rcpp::stop("all elements must share chromosomes for `convex_hull`");
                                    }
                                    double rr = bit_cast<double>(make_ldmap_region(
                                                                                   sa.chrom(),
                                                                                   std::min(sa.start(),sb.start()),
                                                                                   std::max(sa.end(),sb.end())));
                                    return rr;
                                  });
    ret[0]=retflt;
    return(ret);
  }else{
    Rcpp::stop("`convex_hull` is only defined for ldmap_range and ldmap_snp");
  }
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
  if(!(ldmap_range.inherits("ldmap_snp") || ldmap_range.inherits("ldmap_range"))){
    Rcpp::stop("must inherit from ldmap_snp or ldmap_range");
  }
  if(ldmap_range.inherits("ldmap_range"))
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
      return(static_cast<int>(Region::make_Region(x).end()));
    });
  else
    std::transform(ldmap_range.begin(),ldmap_range.end(),ret.begin(),[](double x){
      return(static_cast<int>(SNP::make_snp(x).pos()+1));
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
  if(!(ldmap_snp.inherits("ldmap_snp"))){
    Rcpp::stop("must inherit from 'ldmap_snp'");
  }

  std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[](double x){
    return(static_cast<double>(SNP::make_snp(x).pos()));
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
    auto [chrom,pos,ref,alt] = get_ldmap_snp(bit_cast<uint64_t>(s.snp));
    // Rcpp::Rcerr<<"chrom: "<<static_cast<int>(chrom)<<" pos: "<<pos<<" ref: "<<static_cast<int>(ref)<<" alt: "<<static_cast<int>(alt)<<std::endl;
    if (pos == 0 or chrom == 0 or (chrom == 15 and pos >= 8727373545472ul))
      return Rcpp::String(NA_STRING);

    return Rcpp::String(
                        "chr" + std::to_string(static_cast<int>(chrom)) + ":" +
                        std::to_string(static_cast<uint64_t>(pos)) + "_" +
                        std::string(get_nuc(Nuc{static_cast<char>(ref)}).get_cstring()) +
                        "_" +
                        std::string(get_nuc(Nuc{static_cast<char>(alt)}).get_cstring()));

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
                                                                     return(static_cast<Rbyte>(SNP::make_snp(x).ref()));
                                                                   });
    ret.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
    return ret;
  }else{
    Rcpp::StringVector ret = Rcpp::StringVector(ldmap_snp.size());
    std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[&](double x){

                                                                     return globalnucs.get_nuc(Nuc{static_cast<char>(SNP::make_snp(x).ref())});
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


//' migrate from old representation of ldmap_snp to new representation
//'
//' @param x ldmap_snp vec
//'
//[[Rcpp::export]]
Rcpp::NumericVector migrate_ldmap_snp(Rcpp::NumericVector x){
  Rcpp::NumericVector ret = Rcpp::NumericVector::import_transform(x.begin(),
                                                                  x.end(),
                                                                  [](double d){
                                                                    Snp ts{.flt=d};
                                                                    return(bit_cast<double>(make_ldmap_snp(ts.str.chrom,ts.str.pos,ts.str.ref,ts.str.alt)));
                                                                  });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  return ret;
}




//' migrate from new representation of ldmap_snp to old representation
//'
//' @param x ldmap_snp vec
//'
//[[Rcpp::export]]
Rcpp::NumericVector old_ldmap_snp(Rcpp::NumericVector x){
  Rcpp::NumericVector ret = Rcpp::NumericVector::import_transform(x.begin(),
                                                                  x.end(),
                                                                  [](double d){
                                                                    auto ts =SNP::make_snp(d);
                                                                    Snp rts{.str = {
                                                                        .alt=static_cast<char>(ts.alt()),
                                                                        .ref=static_cast<char>(ts.ref()),
                                                                        .pos=static_cast<uint64_t>(ts.pos()),
                                                                        .chrom=ts.chrom()}};
                                                                    return(rts.flt);
                                                                  });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  return ret;
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
                                                                     return(static_cast<Rbyte>(SNP::make_snp(x).alt()));
                                                                   });
    ret.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
    return ret;
  }else{
    Rcpp::StringVector ret = Rcpp::StringVector(ldmap_snp.size());
    std::transform(ldmap_snp.begin(),ldmap_snp.end(),ret.begin(),[&](double x){
                                                                     return globalnucs.get_nuc(Nuc{static_cast<char>(SNP::make_snp(x).alt())});
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
    SNP sp{SNP::make_snp(bit_cast<uint64_t>(ldmap_snp(i)))};
    chrom(i)=sp.chrom() !=0 ? sp.chrom() : NA_INTEGER;
    pos(i) =
      sp.pos() != 0 ? static_cast<double>(sp.pos()) : R_NaReal;
    ascii_ref(i)=globalnucs.get_nuc(Nuc{static_cast<char>(sp.ref())});
    ascii_alt(i)=globalnucs.get_nuc(Nuc{static_cast<char>(sp.alt())});
  }
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["pos"]=pos,
                                _["ascii_ref"]=ascii_ref,
                                _["ascii_alt"]=ascii_alt);
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") =   seq_len( p);
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
    SNP sp{SNP::make_snp(bit_cast<uint64_t>(ldmap_snp(i)))};
    chrom(i)=sp.chrom() !=0 ? sp.chrom() : NA_INTEGER;
    pos(i)=sp.pos() !=0 ? static_cast<double>(sp.pos()) : NA_REAL;
    ascii_ref(i)=static_cast<Rbyte>(sp.ref());
    ascii_alt(i)=static_cast<Rbyte>(sp.alt());
  }
  ascii_ref.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
  ascii_alt.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");
  auto dfl = Rcpp::List::create(_["chrom"]=chrom,
                                _["pos"]=pos,
                                _["ref"]=ascii_ref,
                                _["alt"]=ascii_alt);
  dfl.attr("class") = StringVector::create("tbl_df","tbl","data.frame");
  dfl.attr("row.names") =   seq_len( p);
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
    SNP q{SNP::make_snp(bit_cast<uint64_t>(qf))};
    int rmi=std::numeric_limits<int>::max();
    auto lb = std::lower_bound(ref_b,ref_e,qf,[](double  snp_a,double snp_b){
                                                return(SNP::make_snp(bit_cast<uint64_t>(snp_a))<SNP::make_snp(bit_cast<uint64_t>(snp_b)));
                                              });
    for(auto ilb = lb; ilb!=ref_e; ilb++){
      if(auto qm = q.allele_match(SNP::make_snp(bit_cast<uint64_t>(*ilb)))){
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
  dfl.attr("row.names") =   seq_len( q_size);
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
                                                                   return globalnucs.string_rep(SNP::make_snp(bit_cast<uint64_t>(x)));
                                                                 }));
}
