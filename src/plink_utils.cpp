#include <Rcpp.h>
#include <bitset>
#include "ldmap.hpp"
#include <string_view>





template<int N>
constexpr std::string_view extract_bits(const Rbyte x){
  const Rbyte bitma1= ((0b00000011u<<(2*N))&x)>>(2*N);
  const std::array<std::string_view,4> fmt{"0/0","./.","0/1","1/1"};
  return fmt[bitma1];
}

static_assert(extract_bits<0>(static_cast<Rbyte>(0xdc)) == "0/0");
static_assert(extract_bits<1>(static_cast<Rbyte>(0xdc)) == "1/1");
static_assert(extract_bits<2>(static_cast<Rbyte>(0xdc)) == "./.");
static_assert(extract_bits<3>(static_cast<Rbyte>(0xdc)) == "1/1");





std::string format_bitset_full(const Rbyte &x){

  return std::string(extract_bits<0>(x))+";"+
    std::string(extract_bits<1>(x))+";"+
    std::string(extract_bits<2>(x))+";"+
    std::string(extract_bits<3>(x));
}

template<int N>
std::string format_bitset(const Rbyte x){

  static_assert(N>0,"Must display at least one genotype");
  static_assert(N<=4,"Cannot contain more than 4 genotypes");

  if constexpr(N==4){
    return format_bitset_full(x);
  }else{
    if constexpr(N==1){
      return std::string(extract_bits<0>(x));
    }
    if constexpr(N==2){
      return std::string(extract_bits<0>(x))+";"+std::string(extract_bits<1>(x));
    }
    if constexpr(N==3){
      return std::string(extract_bits<0>(x))+";"+std::string(extract_bits<1>(x))+";"+std::string(extract_bits<2>(x));
    }
  }
}

//[[Rcpp::export]]
Rcpp::RawVector gt_subset(Rcpp::RawVector x, SEXP i){
  //We're assuming this was called by our overload for `[` so that

  GT_view gtv(x);
  auto index_v = get_index_variant(i);
  return std::visit(gtv,index_v);
}




// Rcpp::RawVector gt_ss(Rcpp::RawVector x, Rcpp::IntegerVector i){
//   //We're assuming this was called by our overload for `[` so that

// }



// `[.vctrs_rcrd` <-  function(x, i, ...) {
//   vec_index(x, i, ...)
// }


//[[Rcpp::export]]
Rcpp::StringVector format_strings(Rcpp::RawVector x) {

  const size_t p=x.size();
  Rcpp::IntegerVector Nv=x.attr("N");
  const int N=Nv[0];
  const int rem = (N % 4);
  const int norm_max = N/4;
  
  int i=0;
  return Rcpp::StringVector::import_transform(x.begin(),x.end(),[&](const Rbyte x){
                                                                  if(i<norm_max){
                                                                    i++;
                                                                    return std::string(format_bitset_full(x));
                                                                  }
                                                                  if(rem==1)
                                                                    return format_bitset<1>(x);
                                                                  if(rem==2)
                                                                    return format_bitset<2>(x);
                                                                  return format_bitset<3>(x);
                                                                });
}



template<int N>
constexpr double byte_2_gt(const Rbyte x){
  constexpr std::array<double, 4> plink_values{0, std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0};
  static_assert(N>=0);

  static_assert(N<4);

  if constexpr(N==0){
    constexpr Rbyte bitma1= 0b00000011u;
    return plink_values[x&bitma1];
  }
  if constexpr(N==1){
    constexpr Rbyte bitma1= 0b00001100u;
    return plink_values[(x&bitma1)>>2];
  }
  if constexpr(N==2){
    constexpr Rbyte bitma1= 0b00110000u;
    return plink_values[(x&bitma1)>>4];
  }else{
    constexpr Rbyte bitma1= 0b11000000u;
    return plink_values[(x&bitma1)>>6];
  }
}



//' Convert ldmap_gt to numeric vector
//'
//' @param x a vector of type ldmap_gt
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector gt2double(const Rcpp::RawVector x){

  Rcpp::IntegerVector Nv=x.attr("N");
  const int N=Nv[0];
  const size_t p=x.size();
  Rcpp::NumericVector ret=Rcpp::no_init(N);
  constexpr Rbyte bitma= ~0;

  static_assert(0xdc == 0b11011100);
  static_assert((0xdc&0b00001100u)>>2 == 0b00000011u);

  static_assert((bitma & 0b00000011)==3u);
  static_assert((bitma & (0b00000011<<2) )==0b00001100u);
  static_assert((0xdc & 0b00110000u)>>4 ==0b00000001u);
  static_assert(byte_2_gt<0>(0xdc)==0);
  static_assert(byte_2_gt<1>(0xdc)==2);
  static_assert(byte_2_gt<3>(0xdc)==2);
  int j=0;
  Rcpp::NumericVector::iterator it=ret.begin();
  const size_t full_size=N/4;
  for(int i=0; i<full_size; i++){
    const Rbyte rb=x(i);
    ret[j++]=byte_2_gt<0>(rb);
    ret[j++]=byte_2_gt<1>(rb);
    ret[j++]=byte_2_gt<2>(rb);
    ret[j++]=byte_2_gt<3>(rb);
  }

  if (full_size<p){
    const Rbyte rb=x(full_size);
    ret[j++]=byte_2_gt<0>(rb);
    if(j<N){
      ret[j++]=byte_2_gt<1>(rb);
      if(j<N){
        ret[j++]=byte_2_gt<2>(rb);
      }
    }
  }
  return ret;
}


GT_view::GT_view(Rcpp::RawVector x_):x(x_),p(x.size()),N(Rcpp::as<Rcpp::IntegerVector>(x.attr("N"))[0]),begin_x(x.begin()){
    if(!x.inherits("ldmap_gt"))
      Rcpp::stop("Internal error: x must be of type ldmap_gt in Internal function `GT_view`");
}



constexpr Rbyte bit_mask(const int i){
  //  static_assert(i<4 and i>=0,"i must be [0 and 4)");
  return 0b00000011u<<(2*i);
}

Rbyte GT_view::operator()(const int i) const noexcept{
  Rbyte ret_x=*(begin_x+(i/4));
  const int offset= i%4;
  return (ret_x&bit_mask(offset))>>(2*offset);
}

Rbyte pack_bytes(const std::vector<Rbyte> &rb){

  const size_t p=rb.size();
  Rbyte ret_b=0;
  int i=0;
  for(auto b:rb){
    ret_b=ret_b|(b<<(2*(i++)));
  }
  return ret_b;
}


Rcpp::RawVector GT_view::operator()(Rcpp::IntegerVector index) const{
  const size_t out_N=index.size();
  const size_t out_p=ceil(static_cast<double>(out_N)/4.0);
  std::vector<Rbyte> rb;
  rb.reserve(4);
  Rcpp::RawVector ret(out_p);
  ret.attr("N")=Rcpp::IntegerVector::create(out_N);
  ret.attr("class")=Rcpp::StringVector::create("ldmap_gt","vctrs_vctr");
  int j=0;
  static_assert(bit_mask(0)==0b00000011u);
  static_assert(bit_mask(1)==0b00001100u);
  const size_t full_size = out_p*4;
  for(int i=0; i<out_N; i++){
    const int index_i=index[i];
    const int offset= (index_i-1)%4;
    const Rbyte tx=*(begin_x+((index_i-1)/4));
    const Rbyte ret_x=(tx&bit_mask(offset))>>(2*offset);
    if(i>0 and offset==0){
      ret[j++]=pack_bytes(rb);
      rb.clear();
    }
    rb.push_back(ret_x);
  }
  if(rb.size()>0){
    ret[j++]=pack_bytes(rb);
    rb.clear();
  }
  return ret;
}


Rcpp::RawVector GT_view::operator()(Rcpp::LogicalVector i) const{
  std::vector<int> rv;
  const size_t fp=i.size();
  for(int j=0;j<p; j++){
    if(i[j % fp]){
      rv.push_back(j);
    }
  }
  return (*this)(Rcpp::IntegerVector(rv.begin(),rv.end()));
}

Rcpp::RawVector GT_view::operator()(Rcpp::StringVector i) const{
  Rcpp::stop("can't subset ldmap_gt by name");
}
Rcpp::RawVector GT_view::operator()(Rcpp::NumericVector i) const{
  return (*this)(Rcpp::as<Rcpp::IntegerVector>(i));
}
