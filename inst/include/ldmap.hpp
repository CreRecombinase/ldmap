#include "R_ext/Arith.h"
#include "alleles.hpp"
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
#include <range/v3/utility/semiregular_box.hpp>
#include <Rcpp.h>
#include <variant>

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
class LDmapSNP{
  Rcpp::NumericVector x;
  double* begin_p;
  size_t p;
  double* end_p;
  bool sorted;
public:
  using elem_type=SNP;
  LDmapSNP(Rcpp::NumericVector x_,std::optional<bool> is_sorted=std::nullopt):
    x(x_),
    begin_p(&(*(x.begin()))),
    p(x.size()),
    end_p(begin_p+p),
    sorted(is_sorted.value_or(std::is_sorted(begin_p,end_p,[](SNP a, SNP b){return a<b;})))
  {
    if(!x.inherits("ldmap_snp")){
      Rcpp::stop("must inherit from ldmap_snp");
    }
  }
};



class LDmapRegion{
  Rcpp::NumericVector x;
  double* begin_p;
  size_t p;
  double* end_p;
  
  bool sorted;
public:
  using elem_type=Region;
  LDmapRegion(Rcpp::NumericVector x_,std::optional<bool> is_sorted=std::nullopt):
    x(x_),
    begin_p(&(*(x.begin()))),
    p(x.size()),
    end_p(begin_p+p),
    sorted(is_sorted.value_or(std::is_sorted(begin_p,end_p,[](Region a, Region b){return a<b;})))
  {
    if(!x.inherits("ldmap_region")){
      Rcpp::stop("must inherit from ldmap_region");
    }
  }
  bool is_sorted() const{
    return sorted;
  }
  template<typename T>
  int within(const T &x) const{
    if(p==1)
      return Region(*begin_p).contains(x) ? 1 : NA_INTEGER;
    if(sorted){
      // Find the first element of target that starts on or after x.start()
      auto candidate_l = std::lower_bound(begin_p,end_p,x,
                                          [](Region a, const T &b) {
                                            return(a.end_SNP()<b.start_SNP());
                                          });
      if(candidate_l==end_p)
        return NA_INTEGER;
      return Region(*candidate_l).contains(x) ? std::distance(begin_p,candidate_l)+1 : NA_INTEGER;
    }else{
      auto candidate_l = std::find(begin_p,end_p,[&](Region a){
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
      return Region(*begin_p).overlap(x) ? 1 : NA_INTEGER;
    
    if(sorted){
      // Find the first element of target that starts on or after x.start()
      auto candidate_l = std::lower_bound(begin_p,end_p,x,
                                          [](Region a, const T &b) {
                                            return(a.end_SNP()<b.start_SNP());
                                          });
      if(candidate_l==end_p)
        return NA_INTEGER;
      return Region(*candidate_l).contains(x) ? std::distance(begin_p,candidate_l)+1 : NA_INTEGER;
    }else{
      auto candidate_l = std::find(begin_p,end_p,[&](Region a){
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
      return abs(Region(*begin_p).distance(x))<max_dist ? 1 : NA_INTEGER;
    
    int best_dist=std::numeric_limits<int>::max();
    if(sorted){
      auto candidate_l = std::lower_bound(begin_p,end_p,x,
                                          [](Region a, const T &b) {
                                            return(a.end_SNP()<b.start_SNP());
                                          });
      if(candidate_l!=end_p){
        best_dist=abs(Region(*candidate_l).distance(x));
      }
      if(best_dist==0){
        return std::distance(begin_p,candidate_l)+1;
      }
      if(candidate_l==begin_p){
        return best_dist<max_dist ? 1 : NA_INTEGER;
      }
      int nbdist = abs(Region(*(candidate_l-1)).distance(x));
      if(nbdist < best_dist)
        return nbdist < max_dist ? std::distance(begin_p,candidate_l) : NA_INTEGER;

      return best_dist < max_dist ? std::distance(begin_p,candidate_l)+1 : NA_INTEGER;
    }else{
      int best_idx=0;
      best_dist=std::abs(Region(*begin_p).distance(x));
      for(auto it = begin_p; it!=end_p; it++){
        auto tdist = std::abs(Region(*it).distance(x));
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
  Rcpp::IntegerVector overlaps_vec(const T &other) const{
    const size_t op = other.p;
    Rcpp::IntegerVector ret=Rcpp::no_init(op);
    if(p==1){
      Region rt(*begin_p);
      std::transform(other.begin_p,other.end_p,[&](T::elem_type x){
                                                 return rt.overlap(x) ? 1 : NA_INTEGER;
                                               });
    }
    std::transform(ret.begin(),ret.end(),[&](T::elem_type x){
                                           return overlap(x);
                                         });
    return ret;
  }

  template<typename T>
  Rcpp::IntegerVector contains_vec(const T &other) const{
    const size_t op = other.p;
    Rcpp::IntegerVector ret=Rcpp::no_init(op);
    if(p==1){
      Region rt(*begin_p);
      std::transform(other.begin_p,other.end_p,ret.begin(),[&](T::elem_type x){
                                                             return rt.overlap(x) ? 1 : NA_INTEGER;
                                                           });
    }
    std::transform(other.begin_p,other.end_p,ret.begin(),[&](T::elem_type x){
                                                           return this->overlap(x);
                                                         });
    return ret;
  }

  template<typename T>
  Rcpp::IntegerVector nearest_vec(const T &other) const{
    const size_t op = other.p;
    Rcpp::IntegerVector ret=Rcpp::no_init(op);
    std::transform(other.begin_p,other.end_p,ret.begin(),[&](T::elem_type x){
                                                           return this->nearest(x);
                                                         });
    return ret;
  }
  

};


template<typename A>
auto cycle_n(A rng,const size_t n){
  return ranges::views::common(ranges::views::take_exactly(ranges::views::cycle(rng),n));
}

