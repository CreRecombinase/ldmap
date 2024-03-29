#include <Rcpp.h>
#include <bitset>
#include "ldmap.hpp"
#include <filesystem>
#include <string_view>
#include <fstream>
#include <range/v3/core.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/iota.hpp>

constexpr Rbyte pack_bytes(Rbyte a,Rbyte b,Rbyte c,Rbyte d){
  return a|(b<<2)|(c<<4)|(d<<6);
}

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


bool has_na_gt(const Rbyte x){
  for( int N=0; N<4; N++){
    if((((0b00000011u<<(2*N))&x)>>(2*N) )== 1)
      return true;
  }
  return false;
}

size_t num_bytes(const int N){
  return std::ceil(static_cast<double>(N)/4.0);
}

Rcpp::RawVector new_ldmap_gt(const int N,const bool init=false){
  Rcpp::RawVector retvec =Rcpp::no_init(num_bytes(N));
  if(init){
    std::fill_n(retvec.begin(),num_bytes(N),0);
  }
  retvec.attr("N")=Rcpp::IntegerVector::create(N);
  retvec.attr("class")=Rcpp::StringVector::create("ldmap_gt","vctrs_vctr");
  return retvec;
}

// static_assert(has_na_gt(0xdc));
// static_assert(!has_na_gt(0));
// static_assert(!has_na_gt(0b11111111));

//[[Rcpp::export]]
void write_plink_bed_f(Rcpp::List x,std::string file_name,bool append=false){

  auto open_type = std::ios::binary;
  const bool fe = std::filesystem::exists(file_name);
  if(append and fe){
    open_type |= std::ios::app;
  }

  std::ofstream osf(file_name, open_type);

  if(!append or !fe){
    Rbyte h[3] = {0x6c,0x1b,0x01};
    osf.write((char*) &h, sizeof(h));
  }

  const size_t p = x.size();
  std::for_each(x.begin(),x.end(),[&](SEXP tx){
    Rcpp::RawVector x(tx);
    Rbyte* xb = &x[0];
    osf.write((char*) xb, x.size());
  });
}

//[[Rcpp::export]]
Rcpp::ListOf<Rcpp::RawVector> read_plink_bed_idx(std::string file_name,Rcpp::IntegerVector id, size_t N){

  //  std::FILE* f = std::fopen(file_name.c_str(), "rb");
  if(!std::filesystem::exists(file_name))
    Rcpp::stop("file_name: "+file_name+" does not exist");

  std::ifstream f(file_name,std::ios_base::binary);
  std::vector<Rbyte> magic_b(3);
  if(magic_b.size()!=3){
    Rcpp::stop("I don't get how initialization works...");
  }
  f.read(reinterpret_cast<char*>(magic_b.data()), magic_b.size());
  std::array<Rbyte,3> magic_t = {0x6c, 0x1b, 0x01};
  if(!std::equal(begin(magic_b),end(magic_b),begin(magic_t),end(magic_t))){
    Rcpp::Rcerr<<"magic sequence should be: "<<std::hex<<0x6c<<","<<0x1b<<","<<0x01<<std::endl;
    Rcpp::Rcerr<<"magic sequence is: ";
    for(unsigned char tb : magic_b){
      Rcpp::Rcerr<<std::hex<<static_cast<int>(tb)<<",";
    }
    Rcpp::Rcerr<<std::endl;
    Rcpp::stop(file_name+" is not a valid plink bed file");
  }
  using namespace ranges;
  auto vec_sizes = num_bytes(N);
  Rcpp::List ret_l = Rcpp::List::import_transform(id.begin(),id.end(),[&](int i){
                                                                        f.seekg(3+(i-1)*vec_sizes);
                                                                        auto rd = new_ldmap_gt(N);
                                                                        Rbyte* trd = &rd[0];
                                                                        f.read(reinterpret_cast<char*>(trd),vec_sizes);
                                                                        return rd;
                                                                      });
  ret_l.attr("class")=Rcpp::StringVector::create("vctrs_list_of","vctrs_vctr");
  ret_l.attr("ptype")=new_ldmap_gt(N,true);
  return ret_l;

}


//[[Rcpp::export]]
Rcpp::ListOf<Rcpp::RawVector> read_plink_bed_l(std::string file_name,const  int p,const size_t N){


  //  std::basic_ifstream<unsigned char> ifs(file_name, std::ios_base::binary);
  if(!std::filesystem::exists(file_name))
    Rcpp::stop("file_name: "+file_name+" does not exist");
  std::FILE* f = std::fopen(file_name.c_str(), "rb");
  std::vector<Rbyte> magic_b(3);
  if(magic_b.size()!=3){
    Rcpp::stop("I don't get how initialization works...");
  }
  auto ret = std::fread(magic_b.data(), sizeof(Rbyte),magic_b.size(),f);
  std::array<Rbyte,3> magic_t = {0x6c, 0x1b, 0x01};
  if(!std::equal(begin(magic_b),end(magic_b),begin(magic_t),end(magic_t))){

    Rcpp::Rcerr<<"magic sequence should be: "<<std::hex<<0x6c<<","<<0x1b<<","<<0x01<<std::endl;
    Rcpp::Rcerr<<"magic sequence is: ";
    for(unsigned char tb : magic_b){
      Rcpp::Rcerr<<std::hex<<static_cast<int>(tb)<<",";
    }
    Rcpp::Rcerr<<std::endl;
    Rcpp::stop(file_name+" is not a valid plink bed file");
  }
  using namespace ranges;
  auto cmr=views::common(views::ints(0,p));
  auto vec_sizes = num_bytes(N);
  Rcpp::List ret_l = Rcpp::List::import_transform(cmr.begin(),cmr.end(),[&](int i){
                                                                          auto rd = new_ldmap_gt(N);
                                                                          auto ret = std::fread(&rd[0],sizeof(Rbyte),vec_sizes,f);
                                                                          return rd;
                                                                        });
  ret_l.attr("class")=Rcpp::StringVector::create("vctrs_list_of","vctrs_vctr");
  ret_l.attr("ptype")=new_ldmap_gt(N,true);
  return ret_l;
}

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



// double cov_gt(Rcpp::RawVector x,Rcpp::RawVector y){
//   double sample_size = get_N(x);
//   if( (get_N(y)!=sample_size) or (x.size()!=y.size()) ){
//       Rcpp::stop("x and y must have the same size");
//   }
//   double tret = dot_ht(x,y);
//   double pop_x = sum_ldmap_ht(x);
//   double pop_y = sum_ldmap_ht(y);

//   double mret=(sample_size * tret - pop_x * pop_y)/(sample_size * (sample_size-1.0));
//   return mret;
// }



// Rcpp::RawVector gt_subset(Rcpp::RawVector x, SEXP i){
//   //We're assuming this was called by our overload for `[` so that

//   GT_view gtv(x);
//   auto index_v = get_index_variant(i);
//   return std::visit(gtv,index_v);
// }


//[[Rcpp::export]]
double gt_af(Rcpp::RawVector x,const bool na_rm=false){
  return GT_view(x).af(na_rm);
}

//[[Rcpp::export]]
Rcpp::NumericVector gt_afs(Rcpp::List x,const bool na_rm=false){
  Rcpp::NumericVector ret = Rcpp::no_init(x.size());
  std::transform(x.begin(),x.end(),ret.begin(),[=](SEXP tx){
                                     return GT_view(tx).af(na_rm);
                                   });
  return ret;
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


int get_N(Rcpp::List x){
  if(x.hasAttribute("ptype")){
    Rcpp::RObject tx = x.attr("ptype");
    if(tx.hasAttribute("N")){
      Rcpp::IntegerVector n=tx.attr("N");
      return n[0];
    }
    else
      Rcpp::stop("I don't think x is a ldmap_gt or ldmap_ht");
  }
  Rcpp::RObject tx = x[0];
  if(tx.hasAttribute("N")){
    Rcpp::IntegerVector n=tx.attr("N");
    return n[0];
  }else
    Rcpp::stop("I don't think x is a ldmap_gt or ldmap_ht");
  return NA_INTEGER;
}


constexpr Rbyte byte_encode(int a) noexcept{
  if(a==0)
    return 0;
  if(a==1)
    return 2;
  if(a==2)
    return 3;
  return 1;
}



Rbyte* copy_int_gt(const int* b,const int *e, Rbyte* d,const int N) {

  const size_t check_n = e-b;
  if(check_n!=N){
    Rcpp::stop("check_n!=N");
  }
  const int full_bytes = (N/4);
  const int full_n = (N/4)*4;
  static_assert(pack_bytes(byte_encode(0),byte_encode(2),byte_encode(3),byte_encode(2))==0xdc);
  const int* rb=b;
  Rbyte* d_e = d+full_bytes;
  for(int i=0; i<full_bytes; i++){
    *d++=pack_bytes(byte_encode(*rb),
                    byte_encode(*(rb+1)),
                    byte_encode(*(rb+2)),
                    byte_encode(*(rb+3)));
    rb+=4;
  }
  const int N_rel = N%4;
  if(N_rel==1){
    *d++=pack_bytes(byte_encode(*rb),
                    byte_encode(5),
                    byte_encode(5),
                    byte_encode(5));
  }
  if(N_rel==2){
    *d++=pack_bytes(byte_encode(*rb),
                    byte_encode(*(rb+1)),
                    byte_encode(5),
                    byte_encode(5));
  }
  if(N_rel==3){
    *d++=pack_bytes(byte_encode(*rb),
                    byte_encode(*(rb+1)),
                    byte_encode(*(rb+2)),
                    byte_encode(5));
  }
  return d;
}


double* copy_gt_double(const Rbyte* b,const Rbyte *e,double* d,const int N) {

  const size_t nb = e-b;
  const int full_s = N/4;

  const Rbyte* rb=b;
  double* d_e = d+N;
  const Rbyte* e_full = rb + full_s;
  double* d_p=d;
  while(rb<e_full){
    *d++=byte_2_gt<0>(*rb);
    *d++=byte_2_gt<1>(*rb);
    *d++=byte_2_gt<2>(*rb);
    *d++=byte_2_gt<3>(*rb);
    rb++;
  }
  if(d<d_e)
    *d++=byte_2_gt<0>(*rb);
  if(d<d_e)
    *d++=byte_2_gt<1>(*rb);
  if(d<d_e)
    *d++=byte_2_gt<2>(*rb);
  if(d<d_e){
    Rcpp::Rcerr<<"d-d_p: "<<d-d_p<<" which is != "<<N<<std::endl;
    Rcpp::stop("bug in copy_gt_double");
  }
  return d;
}

int* copy_gt_int(const Rbyte* b,const Rbyte *e,int* d,const int N) {

  const size_t nb = e-b;
  const int full_s = N/4;

  const Rbyte* rb=b;
  int* d_e = d+N;
  const Rbyte* e_full = rb + full_s;
  int* d_p=d;
  while(rb<e_full){
    *d++=byte_2_gt<0>(*rb);
    *d++=byte_2_gt<1>(*rb);
    *d++=byte_2_gt<2>(*rb);
    *d++=byte_2_gt<3>(*rb);
    rb++;
  }
  if(d<d_e)
    *d++=byte_2_gt<0>(*rb);
  if(d<d_e)
    *d++=byte_2_gt<1>(*rb);
  if(d<d_e)
    *d++=byte_2_gt<2>(*rb);
  if(d<d_e){
    Rcpp::Rcerr<<"d-d_p: "<<d-d_p<<" which is != "<<N<<std::endl;
    Rcpp::stop("bug in copy_gt_double");
  }
  return d;
}


//' Convert list of ldmap_gt to numeric matrix
//'
//' @param x a vector of type ldmap_gt
//' @export
//[[Rcpp::export]]
SEXP gt2matrix(const Rcpp::List x){

  const size_t p = x.size();
  const size_t N = get_N(x);
  Rcpp::NumericVector x_mat= Rcpp::no_init(N*p);
  const size_t full_size=N/4;
  const size_t vec_size=num_bytes(N);
  double* dest_p = &x_mat[0];
  for(int i=0; i<p; i++){
    Rcpp::RawVector tx= x[i];
    const Rbyte* xb = &tx[0];
    const Rbyte* xe = xb+tx.size();
    dest_p = copy_gt_double(xb,xe,dest_p,N);
  }
  x_mat.attr("dim")=Rcpp::IntegerVector::create(N,p);
  return x_mat;
}



//' Convert integer vector to ldmap_gt
//'
//' @param x a vector of type integer
//' @export
//[[Rcpp::export]]
Rcpp::RawVector int2gt(const Rcpp::IntegerVector x){

  const int N=x.size();
  auto retvec = new_ldmap_gt(N);
  const auto num_full = (N/4)*4;
  const int* xb = &x[0];
  const int* xe = xb+N;
  Rbyte* rx = &retvec[0];
  copy_int_gt(xb,xe,rx,N);
  return retvec;
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


GT_view::GT_view(Rcpp::RawVector x_):x(x_),p(x.size()),N(Rcpp::as<Rcpp::IntegerVector>(x.attr("N"))[0]),begin_x(x.begin()),c_af(std::nullopt){
    if(!x.inherits("ldmap_gt"))
      Rcpp::stop("Internal error: x must be of type ldmap_gt in Internal function `GT_view`");
}

GT_view::GT_view(SEXP x_):x(x_),p(x.size()),N(Rcpp::as<Rcpp::IntegerVector>(x.attr("N"))[0]),begin_x(x.begin()),c_af(std::nullopt){
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




constexpr int num_missing(Rbyte a,Rbyte b,Rbyte c,Rbyte d){
  int ina = a==1 ? -1 : 0;
  int inb = b==1 ? -1 : 0;
  int inc = c==1 ? -1 : 0;
  int ind = d==1 ? -1 : 0;
  return ina+inb+inc+ind;
}
static_assert(pack_bytes(3, 3, 3, 3) == 0b11111111);
static_assert(pack_bytes(3, 3, 3, 3) == 255);
static_assert(pack_bytes(3, 2, 3, 3) == 0b11111011);
static_assert(pack_bytes(3,1,3,3) == 0b11110111);
static_assert(pack_bytes(1, 1, 1, 1) == 0b01010101);

template<bool NA_RM=false>
struct CT {
  constexpr CT() : popcnt2_arr(),dosage_arr(),arr(),n_missing(),re_encode() {
    std::array<int,4> dosage{0,0,1,2};
    std::array<Rbyte,4> recode{0,0,1,2};
    //    std::array<char,4> neg_encode{-1,0,0,1};
    for(Rbyte i=0; i<4; i++){
      for(Rbyte j=0; j<4; j++){
        for(Rbyte k=0; k<4; k++){
          for(Rbyte l=0; l<4; l++){
            n_missing[pack_bytes(i,j,k,l)]=num_missing(i,j,k,l);
            arr[pack_bytes(i,j,k,l)]=dosage[i]+dosage[j]+dosage[k]+dosage[l];
            re_encode[pack_bytes(i,j,k,l)]=pack_bytes(recode[i],recode[j],recode[k],recode[l]);
            popcnt2_arr[pack_bytes(i,j,k,l)]=i+j+k+l;
            dosage_arr[pack_bytes(i,j,k,l)] = pack_bytes(dosage[i],dosage[j],dosage[k],dosage[l]);
          }
        }
      }
    }
  }
  int sum_tot(const Rbyte bt,int &N) const{
    if constexpr(NA_RM){
      N+=n_missing[bt];
      return arr[bt];
    }else{
      if(n_missing[bt]!=0)
        return NA_INTEGER;
      return arr[bt];
    }
  }
  //popcount2
  // constexpr char dot(const Rbyte x, const Rbyte y) const{
  //   char z = (x|y) & 0x55;
  //   char w = ((x^y) & (0xaa-z))|z;
  //   return popcnt2_arr[w];
  // }
  char popcnt2_arr[256];
  Rbyte dosage_arr[256];
  char arr[256];
  char n_missing[256];
  int re_encode[256];
};

int GT_view::sum(const bool na_rm=false) const {
  int retv=0;
  int tN=N;
  if(na_rm){
    constexpr CT<true> ctv;
    static_assert(ctv.popcnt2_arr[pack_bytes(1,0,2,1)]==4);
  for(Rbyte xb : x)
    retv+=ctv.sum_tot(xb,tN);
  }
  else{
    constexpr CT<false> ctv;
    for(Rbyte xb : x){
      auto atb = ctv.sum_tot(xb,tN);
      if(Rcpp::IntegerVector::is_na(atb)){
        return NA_INTEGER;
      }
      retv+=atb;
    }
  }
  return retv;
}


Rcpp::RawVector GT_view::dosage(const bool na_rm=false) const {
  Rcpp::RawVector result = Rcpp::clone(x);
  if(!na_rm){
    Rcpp::Rcerr<<"na_rm is currently ignored..."<<std::endl;
  }
  constexpr CT<true> ctv;
  for(int i=0; i<p; i++){
    result[i]=ctv.dosage_arr[*(begin_x+i)];
  }
  result.attr("N")=Rcpp::IntegerVector::create(N);
  result.attr("class")=Rcpp::StringVector::create("ldmap_dosage","vctrs_vctr");
  return result;
}


double GT_view::af(const bool na_rm=false) const {
  static_assert(std::numeric_limits<Rbyte>::max()==255);
  if(!c_af.has_value()){
    int retv=0;
    int tN=N;
    if(na_rm){
      CT<true> ctv;
      for(Rbyte xb : x)
        retv+=ctv.sum_tot(xb,tN);
    }
    else{
      CT<false> ctv;
      for(Rbyte xb : x){
        auto atb = ctv.sum_tot(xb,tN);
        if(Rcpp::IntegerVector::is_na(atb)){
          *c_af =  NA_REAL;
          return NA_REAL;
        }
        retv+=atb;
      }
    }
    *c_af =static_cast<double>(retv)/static_cast<double>(2*tN);
  }
  return *c_af;
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
    if(rb.size()==4){
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
