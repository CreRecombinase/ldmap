#include "alleles.hpp"
#include "range/v3/detail/range_access.hpp"
#include "range/v3/range_fwd.hpp"
#include <limits>
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
#include <range/v3/view/iota.hpp>
#include <range/v3/view/cycle.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/take_exactly.hpp>
#include <range/v3/view/zip_with.hpp>
#include <range/v3/view/common.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/algorithm/upper_bound.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/algorithm/for_each.hpp> // specific includes
#include <range/v3/view/transform.hpp>
#include <range/v3/view/adaptor.hpp>
#include <range/v3/view/facade.hpp>
#include <range/v3/utility/semiregular_box.hpp>
#include <Rcpp.h>
#include <variant>
#include <R_ext/Altrep.h>




class GT_view{
  Rcpp::RawVector x;
  size_t p;
  size_t N;
  Rbyte* begin_x;
  mutable std::optional<double> c_af;
public:
  GT_view(Rcpp::RawVector x_);
  GT_view(SEXP x_);
  //  GT_view(const int N_):x(N_*4),p(N_/4),N(N_),begin_x(x.begin()){}
  //  double get_gt(int i) const;
  //  Rcpp::String get_fmt(int i) const;
  Rbyte operator()(int i) const noexcept;
  Rcpp::RawVector operator()(Rcpp::IntegerVector i) const;
  Rcpp::RawVector operator()(Rcpp::CharacterVector i) const;
  Rcpp::RawVector operator()(Rcpp::LogicalVector i) const;
  Rcpp::RawVector operator()(Rcpp::NumericVector i) const;
  Rcpp::RawVector dosage(const bool na_rm) const;
  double cov(const GT_view &other)const;

  int sum(bool na_rm) const;
  double af(bool na_rm) const;
  //  double allele_frequency() const;
  //  double sum() const;
};



class rcpp_doubles : public ranges::view_facade<rcpp_doubles>
{
  friend ranges::range_access;
  Rcpp::Vector<REALSXP> vec;
  struct cursor
  {
  private:
    Rcpp::Vector<REALSXP>::const_iterator iter;
    public:
        cursor() = default;
    cursor(Rcpp::Vector<REALSXP>::const_iterator it)
          : iter(it)
        {}
    double const&  read() const
    {
      return *iter;
    }
    bool equal(cursor const &that) const
    {
      return iter == that.iter;
    }
    void next()
    {
      ++iter;
    }
    void prev()
    {
      --iter;
    }
    std::ptrdiff_t distance_to(cursor const &that) const
    {
      return that.iter - iter;
    }
    void advance(std::ptrdiff_t n)
    {
      iter += n;
    }
  };
  cursor begin_cursor() const
  {
    return {vec.begin()};
  }
  cursor end_cursor() const
  {
    return {vec.end()};
  }
public:
  rcpp_doubles() = default;
  explicit rcpp_doubles(Rcpp::RObject sz) : vec(SEXP(sz)){}
  explicit rcpp_doubles(Rcpp::Vector<INTSXP> sz) : vec(sz) {}
};



class rcpp_ints : public ranges::view_facade<rcpp_ints>
{
  friend ranges::range_access;
  Rcpp::Vector<INTSXP> vec;
  struct cursor
  {
  private:
    Rcpp::Vector<INTSXP>::const_iterator iter;
    public:
        cursor() = default;
    cursor(Rcpp::Vector<INTSXP>::const_iterator it)
          : iter(it)
        {}
    int const&  read() const
    {
      return *iter;
    }
    bool equal(cursor const &that) const
    {
      return iter == that.iter;
    }
    void next()
    {
      ++iter;
    }
    void prev()
    {
      --iter;
    }
    std::ptrdiff_t distance_to(cursor const &that) const
    {
      return that.iter - iter;
    }
    void advance(std::ptrdiff_t n)
    {
      iter += n;
    }
  };
  cursor begin_cursor() const
  {
    return {vec.begin()};
  }
  cursor end_cursor() const
  {
    return {vec.end()};
  }
public:
  rcpp_ints() : vec(Rcpp::IntegerVector::create()) {}
  explicit rcpp_ints(Rcpp::RObject sz) : vec(SEXP(sz)){}
  explicit rcpp_ints(Rcpp::Vector<INTSXP> sz) : vec(sz) {}
};



class rcpp_raws : public ranges::view_facade<rcpp_raws>
{
  friend ranges::range_access;
  Rcpp::Vector<RAWSXP> vec;
  struct cursor
  {
  private:
    Rcpp::Vector<RAWSXP>::const_iterator iter;
    public:
        cursor() = default;
    cursor(Rcpp::Vector<RAWSXP>::const_iterator it)
          : iter(it)
        {}
    Rbyte const&  read() const
    {
      return *iter;
    }
    bool equal(cursor const &that) const
    {
      return iter == that.iter;
    }
    void next()
    {
      ++iter;
    }
    void prev()
    {
      --iter;
    }
    std::ptrdiff_t distance_to(cursor const &that) const
    {
      return that.iter - iter;
    }
    void advance(std::ptrdiff_t n)
    {
      iter += n;
    }
  };
  cursor begin_cursor() const
  {
    return {vec.begin()};
  }
  cursor end_cursor() const
  {
    return {vec.end()};
  }
public:
  rcpp_raws() : vec(Rcpp::RawVector::create()) {}
  explicit rcpp_raws(Rcpp::RObject sz) : vec(SEXP(sz)){}
  explicit rcpp_raws(Rcpp::Vector<RAWSXP> sz) : vec(sz) {}
};



class rcpp_strings : public ranges::view_facade<rcpp_strings>
{
  friend ranges::range_access;
  Rcpp::StringVector vec;
  struct cursor
  {
  private:
    SEXP vec;
    int i=0;
    //    std::string_view sv;
    public:
    cursor() = default;
    cursor(SEXP v,int it): vec(v),i(it)
    {}
    std::string_view const  read() const
    {
      auto s = STRING_ELT(vec,i);
      auto l = LENGTH( s ) ;
      const char* bgv = CHAR(s);
      std::string_view sv(bgv,l);
      return sv;
    }
    bool equal(cursor const &that) const
    {
      return (i == that.i) and (vec == that.vec);
    }
    void next()
    {
      ++i;
    }
    void prev()
    {
      --i;
    }
    std::ptrdiff_t distance_to(cursor const &that) const
    {
      return that.i - i;
    }
    void advance(std::ptrdiff_t n)
    {
      i += n;
    }
  };
  cursor begin_cursor() const
  {
    return {SEXP(vec),0};
  }
  cursor end_cursor() const
  {
    return {SEXP(vec),LENGTH(vec)};
  }
public:
  rcpp_strings() : vec(Rcpp::StringVector::create()) {}
  explicit rcpp_strings(Rcpp::RObject sz) : vec(sz){}
  explicit rcpp_strings(Rcpp::Vector<STRSXP> sz) : vec(sz) {}
};




inline bool is_compact_seq(SEXP x) {
  if(!ALTREP(x)){
    return false;
  }
  R_xlen_t info_l=::Rf_xlength(R_altrep_data1(x));
  if(info_l!=3){
    return false;
  }
  return true;
}


inline std::variant<Rcpp::NumericVector,Rcpp::IntegerVector,Rcpp::StringVector,Rcpp::LogicalVector>  get_index_variant(SEXP x){

  //  Rcpp::RObject ro(x);
  // auto attr = ro.attributeNames();
  switch (TYPEOF(x)) {
  case INTSXP:{
    return Rcpp::IntegerVector(x);
  }
  case REALSXP:{
    return Rcpp::NumericVector(x);
  }
  case STRSXP:{
    return Rcpp::StringVector(x);
  }
  case LGLSXP:{
    return Rcpp::LogicalVector(x);
  }
  default:
    Rcpp::stop("unexpected type in index_variant (should be integer, numeric, string or logical)");
  }
}



// struct LDMR :ranges::view_facade<LDMR>{
// private:
//   friend ranges::range_access;
//   Rcpp::NumericVector x;
//   struct cursor
//   {
//   private:
//     Rcpp::NumericVector::const_iterator iter;
//     Region cur;

//     public:
//         cursor() = default;
//     cursor(Rcpp::NumericVector::const_iterator it)
//           : iter(it)
//         {}
//         Region read() const
//         {
//           return Region(*iter);
//         }
//         bool equal(cursor const &that) const
//         {
//             return iter == that.iter;
//         }
//         void next()
//         {
//             ++iter;
//         }
//         void prev()
//         {
//             --iter;
//         }
//         std::ptrdiff_t distance_to(cursor const &that) const
//         {
//             return that.iter - iter;
//         }
//         void advance(std::ptrdiff_t n)
//         {
//             iter += n;
//         }
//     };
//     cursor begin_cursor() const
//     {
//         return {x.begin()};
//     }
//     cursor end_cursor() const
//     {
//         return {x.end()};
//     }
// public:
//   LDMR(Rcpp::NumericVector x_):x(x_){

//   }
// };

template<typename elem_type>  
class LDmap_v{
protected:
  Rcpp::NumericVector x;
  double* begin_p;
  size_t p;
  double* end_p;
  bool sorted;
public:
  LDmap_v(Rcpp::NumericVector x_,std::optional<bool> is_sorted=std::nullopt):
    x(x_),
    begin_p(&(*(x.begin()))),
    p(x.size()),
    end_p(begin_p+p),
    sorted(is_sorted.value_or(std::is_sorted(begin_p,end_p,[](Region a, Region b){return a<b;})))
  {
    if constexpr( std::is_same_v<elem_type,Region>){
      if(!x.inherits("ldmap_region"))
        Rcpp::stop("must inherit from ldmap_region");     
    }else{
      static_assert( std::is_same_v<elem_type,SNP>);
      if(!x.inherits("ldmap_snp"))
        Rcpp::stop("must inherit from ldmap_snp");
    }
  }
  bool is_sorted() const{
    return sorted;
  }
  template<typename T>
  int within(const T &x) const{
    if(p==1)
      return elem_type(*begin_p).contains(x) ? 1 : NA_INTEGER;
    if(sorted){
      // Find the first element of target that starts on or after x.start()
      auto candidate_l = std::lower_bound(begin_p,end_p,x,
                                          [](elem_type a, const T &b) {
                                            return(a.end_SNP()<b.start_SNP());
                                          });
      if(candidate_l==end_p)
        return NA_INTEGER;
      return elem_type(*candidate_l).contains(x) ? std::distance(begin_p,candidate_l)+1 : NA_INTEGER;
    }else{
      auto candidate_l = std::find_if(begin_p,end_p,[&](elem_type a){
                                                   return a.contains(x);
                                                 });
      if(candidate_l==end_p)
        return NA_INTEGER;
      return std::distance(begin_p,candidate_l)+1;
    }
  }
  template<typename T>
  int overlap(const T &x) const{
    
    if(p==1)
      return elem_type(*begin_p).overlap(x) ? 1 : NA_INTEGER;
    
    if(sorted){
      // Find the first element of target that starts on or after x.start()
      auto candidate_l = std::lower_bound(begin_p,end_p,x,
                                          [](elem_type a, const T &b) {
                                            return(a.end_SNP()<=b.start_SNP());
                                          });
      if(candidate_l==end_p)
        return NA_INTEGER;
      return elem_type(*candidate_l).overlap(x) ? std::distance(begin_p,candidate_l)+1 : NA_INTEGER;
    }else{
      auto candidate_l = std::find_if(begin_p,end_p,[&](elem_type a){
                                                   return a.overlap(x);
                                                 });
      if(candidate_l==end_p)
        return NA_INTEGER;
      return std::distance(begin_p,candidate_l)+1;
    }
  }

  template<typename T>
  int nearest(const T& x,int max_dist) const{
    if(p==1)
      return abs(elem_type(*begin_p).distance(x))<max_dist ? 1 : NA_INTEGER;
    
    int best_dist=std::numeric_limits<int>::max();
    if(sorted){
      auto candidate_l = std::lower_bound(begin_p,end_p,x,
                                          [](elem_type a, const T &b) {
                                            return(a.end_SNP()<b.start_SNP());
                                          });
      if(candidate_l!=end_p){
        best_dist=abs(elem_type(*candidate_l).distance(x));
      }
      if(best_dist==0){
        return std::distance(begin_p,candidate_l)+1;
      }
      if(candidate_l==begin_p){
        return best_dist<max_dist ? 1 : NA_INTEGER;
      }
      int nbdist = abs(elem_type(*(candidate_l-1)).distance(x));
      if(nbdist < best_dist)
        return nbdist < max_dist ? std::distance(begin_p,candidate_l) : NA_INTEGER;

      return best_dist < max_dist ? std::distance(begin_p,candidate_l)+1 : NA_INTEGER;
    }else{
      int best_idx=0;
      best_dist=std::abs(elem_type(*begin_p).distance(x));
      for(auto it = begin_p; it!=end_p; it++){
        auto tdist = std::abs(elem_type(*it).distance(x));
        if(tdist>best_dist)
          break;
        best_dist=tdist;
        best_idx++;
        if(best_dist==0)
          break;
      }
      return best_dist<max_dist ? best_idx : NA_INTEGER;
    }
  }

    template<typename T>
    Rcpp::IntegerVector overlaps_vec(const LDmap_v<T> &other) const{
      const size_t op = other.p;
      Rcpp::IntegerVector ret=Rcpp::no_init(op);
      if(p==1){
        elem_type rt(*begin_p);
        std::transform(other.begin_p,other.end_p,ret.begin(),[&](T x){
            return rt.overlap(x) ? 1 : NA_INTEGER;
          });
      }
      std::transform(other.begin_p,other.end_p,ret.begin(),[&](T x){
                                                             return overlap(x);
                                                           });
      return ret;
    }

   template<typename T>
   Rcpp::IntegerVector contains_vec(const LDmap_v<T> &other) const{
     const size_t op = other.p;
     Rcpp::IntegerVector ret=Rcpp::no_init(op);
     if(p==1){
       elem_type rt(*begin_p);
       std::transform(other.begin_p,other.end_p,ret.begin(),[&](T x){
                                                              return rt.contains(x) ? 1 : NA_INTEGER;
                                                            });
     }
     std::transform(other.begin_p,other.end_p,ret.begin(),[&](T x){
                                                            return this->within(x);
                                                          });
     return ret;
   }

  template<typename T>
  Rcpp::IntegerVector nearest_vec(const LDmap_v<T> &other,int max_dist=NA_INTEGER) const{
    const size_t op = other.p;
    Rcpp::IntegerVector ret=Rcpp::no_init(op);
    if(Rcpp::IntegerVector::is_na(max_dist)){
      max_dist=std::numeric_limits<int>::max()-1;
    }
    std::transform(other.begin_p,other.end_p,ret.begin(),[&](T x){
                                                           return this->nearest(x,max_dist);
                                                         });
    return ret;
  }
  
  template<typename R>
  friend class LDmap_v;
};

using LDmapSNP=LDmap_v<SNP>;
using LDmapRegion=LDmap_v<Region>;


template<typename A>
auto cycle_n(A rng,const size_t n){
  return ranges::views::common(ranges::views::take_exactly(ranges::views::cycle(rng),n));
}

