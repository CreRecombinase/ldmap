#include "alleles.hpp"
#include "ldmap/genetic_map.hpp"
#include "ldmap.hpp"
#include <progress.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/algorithm/upper_bound.hpp>
#include <range/v3/action/sort.hpp>

double ConstantGeneticMap::interpolate_post(const int pos)const {
  const std::string error_mess="position"+std::to_string(pos)+" is after final position";
  if (strict) {
    Rcpp::stop(error_mess);
  }
  //  Rcpp::Rcerr << error_mess << std::endl;
  const auto it_end = std::prev(genmap.end());
  const auto it_prev = std::prev(it_end);
  return (linear_interp(it_prev, it_end, pos));
}

double ConstantGeneticMap::interpolate_prev(const int pos)const{
  const std::string error_mess="position "+std::to_string(pos)+" is before first position";
  if (strict) {
    Rcpp::stop(error_mess);
  }
  //  Rcpp::Rcerr << error_mess << std::endl;
  const auto gb=genmap.begin();
  if (gb->first == pos) {
    return (gb->second);
  }
  const auto gbp = std::next(gb);
  return (linear_interp(gb, gbp, pos));
}

double ConstantGeneticMap::linear_interp(const idmap::const_iterator prev_it,
                                          const idmap::const_iterator next_it,
                                          const int pos) const {

  const int prev_mappos = prev_it->first;
  const int mappos = next_it->first;
  const double prev_map = prev_it->second;
  const double cur_map = next_it->second;
  double frac = static_cast<double>(pos - prev_mappos) /
    static_cast<double>(mappos - prev_mappos);
  if(strict){
    if (frac < 0 || frac > 1) {
      Rcpp::Rcerr << "in position: " << pos << std::endl;
      Rcpp::Rcerr << prev_it->first << "," << prev_it->second << std::endl;
      Rcpp::Rcerr << next_it->first << "," << next_it->second << std::endl;
      Rcpp::Rcerr << "frac:" << frac << std::endl;
      Rcpp::stop("frac must be positive (and less than 1)");
    }
  }
  return (prev_map + frac * (cur_map - prev_map));
}

double ConstantGeneticMap::interpolate(const int pos)const {
  //Return first value that is greater than or equal to pos
  auto n_it = genmap.lower_bound(pos);
  if(n_it == genmap.end()){
    return (interpolate_post(pos));
  }
  if (n_it->first == pos) {
    return (n_it->second);
  }
  if(n_it == genmap.begin()){

    return (interpolate_prev(pos));
  }

  const auto prev_it = std::prev(n_it);
  return (linear_interp(prev_it, n_it, pos));
}

ConstantGeneticMap::ConstantGeneticMap(const Rcpp::IntegerVector &pos_vec,
                                       const Rcpp::NumericVector &map_vec,
                                       const bool strict_)
  :       strict(strict_),
	  genmap([](const Rcpp::IntegerVector &pos_vec,
		    const Rcpp::NumericVector &map_vec,
		    const bool strict) {
		   auto pos_begin = pos_vec.begin();
		   auto pos_end = pos_vec.end();
		   auto map_begin = map_vec.begin();
		   auto map_end = map_vec.end();
		   const size_t p = pos_end - pos_begin;
		   if (p != (map_end - map_begin)) {
		     Rcpp::stop("position and genetic map must be the same size");
		   }
		   if (p < 2) {
		     Rcpp::stop("There must be at least 2 elements in the genetic map in "
				"order to interpolate!");
		   }
		   std::map<int, double> retmap;
		   auto rb = retmap.begin();
		   auto pb = pos_begin;
		   auto mb = map_begin;
		   double orb=*mb;
		   rb = retmap.insert(rb, std::make_pair(*(pb), *(mb)));
		   bool	has_warned=false;
		   for (size_t i = 1; i < p; i++) {
		     pb++;
		     mb++;
		     if (orb >= *mb) {
                       if (!has_warned) {
                         Rcpp::Rcerr << "Genetic map is not strictly sorted at "
                                        "position: "
                                     << i << " (" << orb << ">=" << *mb
                                     << ")\n(Warning once per invocation)";
			 has_warned=true;
                       }
                       if (strict || orb > *mb) {
			 Rcpp::stop("Genetic map must be strictly sorted ");
		       }
		     } else {
		       rb = retmap.insert(rb, std::make_pair(*(pb), *(mb)));
		     }
		     orb = *mb;
		   }
		   return (retmap);
		 }(pos_vec, map_vec, strict_)) {}


//' Linear interpolation of genetic map values
//'
//' @param map  is a length `p` vector of cumulative genetic map values. `map` must be _strictly_ _sorted_
//' @param map_pos  is a length `p` vector of genome coordinates corresponding to the reference genetic map. `map_pos` must be _strictly_ _sorted_
//' @param target_pos is a vector of coordinates to interpolate
//' @param strict a boolean indicating whether to strictly interpolate
//' @param progress a boolean indicating whether to indicate progress with a progress bar
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector interpolate_genetic_map(const Rcpp::NumericVector &map,
					    const Rcpp::IntegerVector map_pos,
					    const Rcpp::IntegerVector target_pos,
					    const bool strict = true,
					    const bool progress = false) {

  const size_t p = target_pos.size();
  ConstantGeneticMap ref_map(map_pos, map, strict);
  Rcpp::NumericVector ret(p);
  Progress prog_bar(p, progress);
  std::transform(target_pos.begin(), target_pos.end(), ret.begin(),
                 [&](const int i) {
		   prog_bar.increment();
		   return (ref_map.interpolate(i)); });
  return (ret);
}



double linear_interp(const std::pair<SNP,double> prev_it,
                     const std::pair<SNP,double> next_it,
                     const SNP pos){

  const auto [prev_snp,prev_map] = prev_it;
  const auto [next_snp,next_map] = next_it;

  // const int prev_mappos = prev_it->first;
  // const int mappos = next_it->first;
  // const double prev_map = prev_it->second;
  // const double cur_map = next_it->second;


  const double frac = pos.relative_distance(prev_snp,next_snp);
  if(std::isnan(frac)){
    Rcpp::stop("cannot perform linear interpolation between chromosomes");
  }
  // double frac = static_cast<double>(pos - prev_mappos) /
  //   static_cast<double>(mappos - prev_mappos);
  // if(strict){
  //   if (frac < 0 || frac > 1) {
  //     Rcpp::Rcerr << "in position: " << pos << std::endl;
  //     Rcpp::Rcerr << prev_it->first << "," << prev_it->second << std::endl;
  //     Rcpp::Rcerr << next_it->first << "," << next_it->second << std::endl;
  //     Rcpp::Rcerr << "frac:" << frac << std::endl;
  //     Rcpp::stop("frac must be positive (and less than 1)");
  //   }

  return (prev_map + frac * (next_map - prev_map));
}



// find a snp less than or equal to and a snp greater than pos
template<typename T>
std::pair<T, T> find_pair(T zipped_b,T zipped_e, const SNP pos){

  static_assert(SNP::make_SNP<true>(1,15).distance(SNP::make_SNP<true>(1,10))==5);
  static_assert(SNP::make_SNP<true>(1,10).distance(SNP::make_SNP<true>(1,15))==(-5));
  //  static_assert(SNP::make_SNP<true>(1,10).distance(SNP::make_SNP<true>(1,20))==(-5));
  //  static_assert(SNP::make_SNP<true>(1,15).relative_distance(SNP::make_SNP<true>(1,10),SNP::make_SNP<true>(1,15)));
  auto zipped_range = ranges::make_subrange(zipped_b,zipped_e);
  if(ranges::size(zipped_range)<2){
    Rcpp::stop("zipped_range must be at least size 2 to interpolate");
  }

  auto ret_it = ranges::lower_bound(zipped_range,pos,ranges::less(),&std::pair<SNP,double>::first);

  // there is no value greater than or equal to pos in zipped_range
  // this means every value is less than pos in zipped_range
  // and it also means that you can't interpolate between two
  // values, instead we'll do "right interpolation"
  if(ret_it==std::end(zipped_range)){
    auto lp = ranges::prev(ranges::end(zipped_range));
    auto fp = ranges::prev(lp);
    return std::make_pair(fp,lp);
  }


  if(ret_it->first.distance(pos)==0)
    return std::make_pair(ret_it,next(ret_it));


  // If the first element is the smallest geq pos, then there is nothing lt
  // This means we'll do "left interpolation"
  if(ret_it == zipped_b)
    return std::make_pair(ret_it,next(ret_it));

  // The element before the smallest geq must be lt (if strictly sorted)

  auto fp = ranges::prev(ret_it);
  auto lp = ret_it;
  // Rcpp::Rcerr<<fp->first<<" vs "<<lp->first<<std::endl;
  // Rcpp::Rcerr<<static_cast<int>(fp->first.chrom())<<" vs "<<static_cast<int>(lp->first.chrom())<<std::endl;
  // Rcpp::Rcerr<<std::boolalpha<<! ((fp->first.chrom() ) == (lp->first.chrom() ))<<std::endl;
  if(! ((fp->first.chrom() ) == (lp->first.chrom() ))){
    if(pos.chrom()==fp->first.chrom()){
      if(fp > zipped_b)
        return std::make_pair(ranges::prev(fp),fp);
      else{
        //        Rcpp::Rcerr<<fp->first.chrom()<<" vs "<<lp->first.chrom()<<" looking for  "<<pos<<std::endl;
        Rcpp::stop("cannot perform (left) linear interpolation between chromosomes");
      }
    }else{
      if(ranges::next(lp) < zipped_e)
        return std::make_pair(lp,ranges::next(lp));
      else{
        // Rcpp::Rcerr<<fp->first<<" vs "<<lp->first<<std::endl;
        // Rcpp::Rcerr<<fp->first.chrom()<<" vs "<<lp->first.chrom()<<" looking for  "<<pos<<std::endl;
        Rcpp::stop("cannot perform (right) linear interpolation between chromosomes");
      }
    }
  }
  return std::make_pair(fp,lp);
}



//' Linear interpolation of genetic map values
//'
//' @param map  is a length `p` vector of cumulative genetic map values. `map` must be _strictly_ _sorted_
//' @param map_pos  is a length `p` vector of genome coordinates corresponding to the reference genetic map. `map_pos` must be _strictly_ _sorted_
//' @param target_pos is a vector of coordinates to interpolate
//' @param strict a boolean indicating whether to strictly interpolate
//' @param progress a boolean indicating whether to indicate progress with a progress bar
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector new_interpolate_genetic_map(const Rcpp::NumericVector &map,
					    const Rcpp::NumericVector map_pos,
					    const Rcpp::NumericVector target_pos,
					    const bool strict = true,
					    const bool progress = false) {


  auto make_snp = [](const double x) -> SNP{
                       return SNP::make_SNP(bit_cast<uint64_t>(x));
                     };

  const size_t p = target_pos.size();
  auto ref_pos=ranges::views::transform(rcpp_doubles(map_pos),make_snp);
  auto q_pos =  ranges::views::common(ranges::views::transform(rcpp_doubles(target_pos),make_snp));


  auto dub_map = rcpp_doubles(map);

  auto zip_v = ranges::views::zip(ref_pos,dub_map) |
    ranges::to<std::vector>() |
    ranges::actions::sort([](std::pair<SNP,double> a,std::pair<SNP,double> b){
                            return a.first < b.first;
                          });



  auto zip_b = std::begin(zip_v);
  auto zip_e = std::end(zip_v);



  return Rcpp::NumericVector::import_transform(
                                               std::begin(q_pos),
                                               std::end(q_pos),
                                               [&](SNP ix) mutable  -> double {
                                                 auto [first_p, last_p] = find_pair(zip_b,zip_e,ix);
                                                 if(first_p->first.distance(ix) == 0)
                                                   return first_p->second;
                                                 double ret_v = linear_interp(*first_p,*last_p,ix);
                                                 zip_b = first_p;
                                                 return ret_v;
                                               });



}
