#include <Rcpp.h>
#include <cstdint>
#include "Rcpp/vector/instantiation.h"
#include "ldmap.hpp"
#include <range/v3/core.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/stride.hpp>
#include <numeric>
#include <libpopcnt.h>


size_t num_bits(const int N){
  return std::ceil(static_cast<double>(N)/static_cast<double>(sizeof(double)*8));
}

//static_assert('1'==0b110001);


template< class Rng>
double pack_range(Rng inp){
  using namespace ranges;
  std::uint64_t ret_v=0;
  ret_v = accumulate(views::enumerate(inp),ret_v,[](const uint64_t a, auto b){
    const auto idx = std::get<0>(b);
    const std::uint64_t val = std::get<1>(b);
    return a|(val<<idx);
  });
  return(bit_cast<double>(ret_v));
}





Rcpp::NumericVector parse_ldmap_ht(SEXP data){

  using namespace ranges;
  const size_t str_len= ::Rf_length(data);
  const char* xb = CHAR(data);
  const char* xe = xb+str_len;
  const size_t N=(str_len+1)/2;

  auto chunk_x  = subrange(xb,xe) |
    views::stride(2) | views::transform([](char c) -> std::uint64_t{
      return c-48;
    }) |
    views::chunk(64);
  auto chunk_r = views::common(chunk_x);
  Rcpp::NumericVector ret = Rcpp::NumericVector::import_transform(std::begin(chunk_r),std::end(chunk_r),[](auto chnk) -> double{
    const double retval = pack_range(chnk);
    return retval;
  });
  ret.attr("N")=Rcpp::IntegerVector::create(N);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_ht","vctrs_vctr");
  return ret;
}


//[[Rcpp::export]]
Rcpp::List parse_hap(Rcpp::StringVector x){

  Rcpp::String itx = x[0];
  Rcpp::List retl = Rcpp::List::import_transform(x.begin(),x.end(),[&](SEXP ix)-> Rcpp::NumericVector{
    return parse_ldmap_ht(ix);
  });
  // retl.attr("N")=Rcpp::IntegerVector::create(N);
  // retl.attr("class")=Rcpp::StringVector::create("ldmap_ht","vctrs_vctr");
  return retl;
}



//' New vector of haplotypes
//'
//' @param x data vector
//' @param N sample size
//' @return vector of (packed, dense) haplotypes
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector new_ldmap_ht(Rcpp::IntegerVector x){
  static_assert(sizeof(double)*8==64);
  // = Rcpp::no_init(num_bits(N));
  const int N =x.size();
  using namespace ranges;
  const int* beg_p = x.begin();
  const int* end_p = x.end();

  auto chunk_r  = subrange(beg_p,end_p) | views::chunk(64) | views::transform([](auto chnk) -> double{
    const double retval = pack_range(chnk);
    return retval;
  });
  auto chunk_x = views::common(chunk_r);

  Rcpp::NumericVector retvec= Rcpp::NumericVector(std::begin(chunk_x),std::end(chunk_x));
  retvec.attr("N")=Rcpp::IntegerVector::create(N);
  retvec.attr("class")=Rcpp::StringVector::create("ldmap_ht","vctrs_vctr");
  return retvec;
}


constexpr int hap_2_gt(const std::uint64_t x,const int offset){
  const std::uint64_t bitma1=(1<<offset);
  return (x&bitma1) >> offset;
}



//[[Rcpp::export]]
int sum_ldmap_ht(Rcpp::NumericVector x){
  const void* d = &(x[0]);
  return popcnt(d,x.size()*sizeof(double));
}

int get_N(Rcpp::NumericVector x){
  Rcpp::IntegerVector Nv= x.hasAttribute("N") ? x.attr("N") : Rcpp::IntegerVector::create(x.size()*64);
  return Nv[0];
}

//[[Rcpp::export]]
int dot_ht(Rcpp::NumericVector x,Rcpp::NumericVector y){
  int retval=0;
  int sample_size = get_N(x);
  if((get_N(y)!=sample_size) or (x.size()!=y.size())){
    Rcpp::stop("x and y must have the same size");
  }
  retval = std::inner_product(x.begin(),x.end(),y.begin(),retval,std::plus<int>(),[](const double a,const double b) noexcept -> int{
    return popcount64(bit_cast<std::uint64_t>(a)&bit_cast<std::uint64_t>(b));
  });
  return retval;
}

//[[Rcpp::export]]
double cov_ht(Rcpp::NumericVector x,Rcpp::NumericVector y){
  double sample_size = get_N(x);
  if( (get_N(y)!=sample_size) or (x.size()!=y.size()) ){
      Rcpp::stop("x and y must have the same size");
  }
  double tret = dot_ht(x,y);
  double pop_x = sum_ldmap_ht(x);
  double pop_y = sum_ldmap_ht(y);

  double mret=(sample_size * tret - pop_x * pop_y)/(sample_size * (sample_size-1.0));
  return mret;
}


//[[Rcpp::export]]
Rcpp::NumericMatrix cov_htm(Rcpp::List x,const bool cov_2_cor=false){

  const size_t p=x.size();
  Rcpp::NumericMatrix ret(p,p);
  for(int i=0; i<p; i++){
    Rcpp::NumericVector a(x[i]);
    for(int j=i; j<p; j++){
      Rcpp::NumericVector b(x[j]);
      double cr=cov_ht(a,b);
      ret(i,j)=cr;
      ret(j,i)=cr;
    }
  }
  if(cov_2_cor){
    for(int i=0; i<p; i++){
      ret(i,i)=1.0;
    }
  }
  return ret;
}



[[gnu::pure]]
double calc_nmsum(const double m) noexcept {
  int nms=(2 * m - 1);
  double sums=0;
  for(double i=1; i<nms; i++){
    sums += (1.0 / i) ;
  }

  return sums;
}

[[gnu::pure]]
double calc_theta(const double m) noexcept {
  double nmsum=calc_nmsum(m);
  return((1/nmsum)/(2*m+1/nmsum));
}

[[gnu::pure]]
double calculate_rho(const double map_a,const double map_b, const double m, const double Ne,const double cutoff) noexcept {
  if(map_a == map_b)
    return 1.0;
  auto ddiff = 4.0 * Ne*(map_a-map_b)/100.0;
  ddiff = std::exp(ddiff / (2 * m));
  return ddiff>cutoff ? ddiff : 0.0;
}


//[[Rcpp::export]]
Rcpp::NumericMatrix ldshrink_S(Rcpp::List x,Rcpp::NumericVector map,const double m=85,const double Ne=11490.672741, const double cutoff=0.001,const bool cov_2_cor=true){

  const size_t p=map.size();
  if(x.size()!=p){
    Rcpp::stop("x and map must have equal length in `ldshrink_scale`");
  }
  const double theta = calc_theta(m);
  const double one_minus_theta_sq = ((1.0-theta)*(1.0-theta));
  const double theta_diag = (theta/2.0)*(1-(theta/2.0));
  Rcpp::NumericMatrix ret(p,p);
  for(int i=0; i<p; i++){
    Rcpp::NumericVector a(x[i]);
    const double map_a =map[i];
    for(int j=i; j<p; j++){
      Rcpp::NumericVector b(x[j]);
      const double map_b = map[j];
      const double sigma=cov_ht(a,b);
      const double shrnk =calculate_rho(map_a,map_b,m,Ne,cutoff);
      const double s = shrnk * sigma;
      double sig_hat = one_minus_theta_sq * s;
      if(i==j){
        sig_hat += theta_diag;
      }
      ret(i,j)=sig_hat;
      ret(j,i)=sig_hat;
    }
  }
  if(cov_2_cor){
    for(int i=0; i<p; i++){
      const int min_i = std::min(i+1,static_cast<int>(p));
      const double i_rat = std::sqrt(1.0/ret(i,i));
      for(int j=min_i; j<p; j++){
        const double j_rat =std::sqrt(1.0/ret(j,j));
        ret(i,j)*=i_rat*j_rat;
        ret(j,i)=ret(i,j);
      }
    }
  }

  if(cov_2_cor){
    for(int i=0; i<p; i++){
      ret(i,i)= 1.0;
    }
  }
  return ret;
}


//[[Rcpp::export]]
Rcpp::IntegerVector ht2int(Rcpp::NumericVector x){

  // static_assert(hap_2_gt(0b0001,0)==1);
  // static_assert(hap_2_gt(0b01,1)==0);
  // static_assert(hap_2_gt(0b11,1)==1);
  // static_assert(hap_2_gt(0b111,2)==1);


  const int N=get_N(x);
  const size_t p=x.size();
  Rcpp::IntegerVector ret=Rcpp::no_init(N);
  int j=0;
  const size_t full_size=N/64;
  for(int i=0; i<full_size; i++){
    const uint64_t rb=bit_cast<uint64_t>(x(i));
    for(int b=0; b<64; b++){
      ret[j++]=hap_2_gt(rb,b);
    }
  }
  if (full_size<p){
    const uint64_t rb=bit_cast<uint64_t>(x(full_size));
    int b=0;
    while(j<N){
      ret[j++]=hap_2_gt(rb,b++);
    }
  }
  return ret;
}




Rcpp::String format_bin(const double x,const int N=64){

  std::string rv;
  rv.reserve(N);
  const std::uint64_t rx=bit_cast<std::uint64_t>(x);
  const std::uint64_t bitma1=1;
  for( int i=0; i<N; i++){
    rv.push_back(((rx&(bitma1<<i))>>i) == 1 ? '1' : '0');
  }
  return Rcpp::String(rv.c_str());
}

//[[Rcpp::export]]
Rcpp::StringVector format_ht(Rcpp::NumericVector x) {
  const size_t p=x.size();
  Rcpp::IntegerVector Nv=x.attr("N");
  const int N=Nv[0];
  const int rem = (N % 64);
  const int norm_max = N/64;

  int i=0;
  std::bitset<64> ix;
  std::uint64_t tx;
  return Rcpp::StringVector::import_transform(x.begin(),x.end(),[&](const double x) -> Rcpp::String{
    if(i<norm_max){
      i++;
      tx=bit_cast<std::uint64_t>(x);
      ix=tx;
      return Rcpp::String(ix.to_string().c_str());
    }
    return format_bin(x,rem);
  });
}
