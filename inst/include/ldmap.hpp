#include "alleles.hpp"
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
#include <range/v3/view/zip_with.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/algorithm/upper_bound.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/algorithm/for_each.hpp> // specific includes
#include <range/v3/view/transform.hpp>
#include <Rcpp.h>

template<typename T>
inline int range_within_range(Region x, T target){
  using namespace ranges;
  SNP sx=x.start_SNP();
  SNP esx=x.last_SNP();
  const size_t target_size = size(target);

  if(target_size==1){
    Region rt= *begin(target);
    return rt.contains(x) ? 1 : NA_INTEGER;
  }

  auto candidate_l = lower_bound(target,x,
                                 [&](Region a, Region b) {
                                   return(a.start_SNP()<=b.start_SNP());
                                 });

  auto xbc = distance(begin(target),candidate_l);

  if( xbc==target_size)
    return back(target).contains(x) ? target_size : NA_INTEGER;

  auto candidate_lr = views::drop(target,xbc-1);
  auto candidate_u = upper_bound(candidate_lr,x,
                                 [&](Region a, Region b){
                                   return(a.end_SNP()<b.end_SNP());
                                 });
  size_t xbe = xbc+distance(begin(candidate_lr),candidate_u);

  if(xbe > xbc)
    return NA_INTEGER;

  return xbc;
}



template<typename T>
inline int range_overlap_range(const Region x, T target){

  using namespace ranges;
  const size_t target_size = size(target);
  if(target_size==1){
    Region rt= *begin(target);
    return rt.overlap(x) ? 1 : NA_INTEGER;
  }

  auto bc = [](Region range_a, Region range_b){
              return(range_a<range_b);
            };

  auto xb=lower_bound(target,x,bc);
  if( xb==end(target) || !((*xb).overlap(x)))
    return NA_INTEGER;
  return std::distance(begin(target),xb)+1;
}




template<typename T>
inline int nearest_snp(const SNP x, T target){

  using namespace ranges;
  const size_t target_size = size(target);
  if(target_size == 1){
    SNP rt= *begin(target);
    return 1;
  }


  auto candidate_l = lower_bound(target,x,
                                 [&](Region a, SNP b) {
                                   return(a.start_SNP() <= b);
                                 });
  if(candidate_l == end(target))
    return target_size;

  if(candidate_l > begin(target))
    candidate_l--;

  auto cand_d = min(abs(



  auto candidate_lr = views::drop(target,xbc-1);

  auto candidate_u = lower_bound(target,x,
                                 [&](Region a, SNP b) {
                                   return(a.end_SNP() <= b);
                                 });


  auto xbc = distance(begin(target),candidate_l);

  if( xbc==target_size)
    return back(target).contains(x) ? target_size : NA_INTEGER;

  auto candidate_lr = views::drop(target,xbc-1);
  auto candidate_u = upper_bound(candidate_lr,x,
                                 [&](Region a, Region b){
                                   return(a.end_SNP()<b.end_SNP());
                                 });




  size_t xbe = xbc+distance(begin(candidate_lr),candidate_u);




  auto bc = [](SNP snp_a, SNP snp_b){
    return(snp_a<snp_b);
  };

  auto xb=lower_bound(target,x,bc);
  if(xb==begin(target))
    return 1;

  if( xb==end(target))
    return target_size;

  //{1 3 5 7} find 8 -> end
  //{1 3 5 7} find 6 -> 7
  int xbi = distance(begin(target),xb);
  auto xb_dist = std::abs(x.distance(*xb));
  auto xbp = xb-1;
  auto xbp_dist = std::abs(x.distance(*xbp));
  return xb_dist<xbp_dist ? xbi:xbi-1;
}
