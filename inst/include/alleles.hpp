#pragma once
// #include <bits/stdint-intn.h>
// #include <bits/stdint-uintn.h>
#include <algorithm>

#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/interval_bounds.hpp>
#include <cstdint>
#include <boost/icl/split_interval_set.hpp>
#include <optional>
#include <cstring>
#include <stdint.h>
#include <tuple>
#include <bitset>
#include <iostream>
enum class Ref : char { A = 0, C = 1, G = 2, T = 3, N = 4 };
template <class To, class From>
typename std::enable_if<
    (sizeof(To) == sizeof(From)) &&
    std::is_trivially_copyable<From>::value &&
    std::is_trivial<To>::value,
    // this implementation requires that To is trivially default constructible
    To>::type
// constexpr support needs compiler magic
inline bit_cast(const From &src) noexcept
{
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}


// Design goals:

/*
  1. Represent as much information about a (human) Single nucleotide
  polymorphism as possible using a 64 bit integer
  2. Support the case of missing data/NA values for ref and alt alleles
  3. Find the best match in a target for a query SNP, and efficiently and
  effectively communicate operations that need to be done to make a match
 */

// Match criteria;

/*
  exact equality: chromosome and position
  inexact equality: ref and alt
  ambiguity allowed in alternate alleles only (optionally)
 */


constexpr char ascii2Nuc(const char l){
  if(l <= ((1) | (2<<0) | (2<<1) | (2<<2))){
    return l;
  }
  switch(l){
  case 'A':
    return 1;
  case 'C':
    return 2 << 0;
  case 'G':
    return 2 << 1;
  case 'T':
    return 2 << 2;
  case 'M': //A or C
    return (1)|(2<<0);
  case 'R': //A or G
    return (1)|(2<<1);
  case 'W': // A or T
    return (1)|(2<<2);
  case 'S': // C or G
    return ( 2<<0 )| (2<<1);
  case 'Y': // C or T
    return (2<<0) | (2<<2);
  case 'K': //G or T
    return (2<<1) | (2<<2);
  case 'V': //A or C or G
    return (1)|(2<<0)|(2<<1);
  case 'H'://A or C or T
    return (1)|(2<<0)|(2<<2);
  case 'D': //A or G or T
    return (1)|(2<<1)|(2<<2);
  case 'B': //C or G or T
    return (2<<0)|(2<<1)|(2<<2);
  case 'N': // A or C or G or T
    return (1) | (2<<0) | (2<<1) | (2<<2);
  default:
    return 0;
  }
}

constexpr inline char ascii2Nuc(const int l){
  static_assert(static_cast<char>(65L)=='A');

  return ascii2Nuc(static_cast<char>(l));
}

constexpr inline char ascii2Nuc(const char *l){
  return ascii2Nuc(l[0]);
}
inline constexpr char operator"" _L(char l){return ascii2Nuc(l);}


struct Nuc{
  char let;
  constexpr Nuc operator^(const Nuc &l)const{
    return Nuc{static_cast<char>(let^l.let)};
  }
  constexpr bool operator==(const Nuc &l)const {
    return let==l.let;
  }
  constexpr bool match(const Nuc &l)const {
    return (let&l.let)>0;
  }
  constexpr bool operator!=(const Nuc &l)const {
    return let!=l.let;
  }
  constexpr Nuc operator|(const Nuc &l)const{
    return Nuc{static_cast<char>(let|l.let)};
  }
  constexpr bool is_ambiguous() const{
    return __builtin_popcount(let)>1;
  }
  constexpr bool is_NA() const{
    return let=='\0';
  }
  constexpr Nuc complement(const Nuc& other);
};


inline constexpr Nuc operator"" _N(char l) { return Nuc{ascii2Nuc(l)}; }


inline constexpr Nuc complement(const Nuc &l){
  switch(l.let){
  case 'A'_L:
    return 'T'_N;
  case 'C'_L:
    return 'G'_N;
  case 'G'_L:
    return 'C'_N;
  case 'T'_L:
    return 'A'_N;
  case 'M'_L: //A or C
    return 'K'_N;
  case 'R'_L: //A or G
    return 'Y'_N;
  case 'W'_L: // A or T
    return 'W'_N;
  case 'S'_L: // C or G
    return 'S'_N;
  case 'Y'_L: // C or T
    return 'R'_N;
  case 'K'_L: //G or T
    return 'M'_N;
  case 'V'_L: //A or C or G
    return 'B'_N;
  case 'H'_L://A or C or T
    return 'D'_N;
  case 'D'_L: //A or G or T
    return 'H'_N;
  case 'B'_L: //C or G or T
    return 'V'_N;
  case 'N'_L: // A or C or G or T
    return 'N'_N;
  default:
    return '0'_N;
  }
}

inline constexpr Nuc Nuc::complement(const Nuc &other){
  return(complement(other));
}

constexpr char Nuc2char(Nuc l){
  switch(l.let){
  case 'A'_L:
    return 'A';
  case 'C'_L:
    return 'C';
  case 'G'_L:
    return 'G';
  case 'T'_L:
    return 'T';
  case 'M'_L: //A or C
    return 'M';
  case 'R'_L: //A or G
    return 'R';
  case 'W'_L: // A or T
    return 'W';
  case 'S'_L: // C or G
    return 'S';
  case 'Y'_L: // C or T
    return 'Y';
  case 'K'_L: //G or T
    return 'K';
  case 'V'_L: //A or C or G
    return 'V';
  case 'H'_L://A or C or T
    return 'H';
  case 'D'_L: //A or G or T
    return 'S';
  case 'B'_L: //C or G or T
    return 'B';
  case 'N'_L: // A or C or G or T
    return 'N';
  default:
    return '\0';
  }
}

constexpr const char* Nuc2string(Nuc l){
  switch(l.let){
  case 'A'_L:
    return "A";
  case 'C'_L:
    return "C";
  case 'G'_L:
    return "G";
  case 'T'_L:
    return "T";
  case 'M'_L: //A or C
    return "M";
  case 'R'_L: //A or G
    return "R";
  case 'W'_L: // A or T
    return "W";
  case 'S'_L: // C or G
    return "S";
  case 'Y'_L: // C or T
    return "Y";
  case 'K'_L: //G or T
    return "K";
  case 'V'_L: //A or C or G
    return "V";
  case 'H'_L://A or C or T
    return "H";
  case 'D'_L: //A or G or T
    return "S";
  case 'B'_L: //C or G or T
    return "B";
  case 'N'_L: // A or C or G or T
    return "N";
  default:
    return "\0";
  }
}

struct bit_region{
  uint64_t end :29;
  uint64_t start :29;
  unsigned char chrom : 6;
} __attribute__((packed));


union bed_region{
  bit_region str;
  uint64_t dat;
  double flt;
};

struct bit_snp{
  char alt : 8; //8
  char ref : 8; //8
  uint64_t pos :43; //42
  unsigned char chrom : 5; //6
} __attribute__((packed));

union Snp{
  bit_snp str;
  std::int64_t dat;
  double flt;
};


enum class Snp_match {
                      perfect_match = 1,
                      reverse_match = 2,                  // match -> reverse
                      complement_match = 3,               // match -> complement
                      reverse_complement_match = 4,       // match -> reverse
                      ambig_match = 5,                    // match on ref and alt
                      reverse_ambig_match = 6,            // match on ref and alt (ambig)
                      complement_ambig_match = 7,         // match on ref and alt (ambig)
                      reverse_complement_ambig_match = 8, // match on ref and alt (ambig)
};

inline std::optional<Snp_match> is_perfect(const char lref, const char lalt, const char rref, const char ralt) {
    const auto match_ref = __builtin_popcount(lref & rref);
    if(match_ref>0){
      const auto match_alt = __builtin_popcount(lalt & ralt);
      if(match_alt>0){
        return match_ref+match_alt==2 ? Snp_match::perfect_match : Snp_match::ambig_match;
      }
    }
    return std::nullopt;
}

inline std::optional<Snp_match> is_reverse(const char lref, const char lalt, const char rref, const char ralt) {
    const auto match_ref = __builtin_popcount(lref & ralt);
    if(match_ref>0){
      const auto match_alt = __builtin_popcount(lalt & rref);
      if(match_alt>0){
        return match_ref+match_alt==2 ? Snp_match::reverse_match : Snp_match::reverse_ambig_match;
      }
    }
    return std::nullopt;
}

inline std::optional<Snp_match> is_complement(const char lref, const char lalt, const char rref, const char ralt) {
  const auto match_ref = __builtin_popcount(complement(Nuc{lref}).let & rref);
    if(match_ref>0){
      const auto match_alt = __builtin_popcount(complement(Nuc{lalt}).let & ralt);
      if(match_alt>0){
        return match_ref+match_alt==2 ? Snp_match::complement_match : Snp_match::complement_ambig_match;
      }
    }
    return std::nullopt;
}

inline std::optional<Snp_match> is_reverse_complement(const char lref, const char lalt, const char rref, const char ralt) {
  const auto match_ref = __builtin_popcount(complement(Nuc{lref}).let & ralt);
    if(match_ref>0){
      const auto match_alt = __builtin_popcount(complement(Nuc{lalt}).let & rref);
      if(match_alt>0){
        return match_ref+match_alt==2 ? Snp_match::reverse_complement_match : Snp_match::reverse_complement_ambig_match;
      }
    }
    return std::nullopt;
}


[[nodiscard]] inline constexpr uint64_t make_mask(const int offset, const int size){
  uint64_t retval=0;
  const uint64_t inn=1;
  for(uint64_t i=offset; i<offset+size; i++){
    retval|=inn<<i;
  }
  return retval;
}


// static_assert(make_mask(0, 1) == 1);//1
// static_assert(make_mask(0, 2) == 3);//2
// static_assert(make_mask(0, 3) == 7);//4
// static_assert(make_mask(0, 4) == 15); // 8
// static_assert(make_mask(0, 5) == 31);//16
// static_assert(make_mask(0, 5) == 31); // 32
// static_assert(make_mask(0, 6) == 63); // 64
// static_assert(make_mask(0, 7) == 127); // 128
// static_assert(make_mask(0, 8) == 255); // 256   //
// static_assert(make_mask(0, 58) == 288230376151711743); // 32




template<typename T>
[[nodiscard]] inline constexpr uint64_t set_chrom(const uint64_t input,const T chrom)noexcept{
  return (~make_mask(58,6)&input)|(static_cast<uint64_t>(chrom) << (29+29));
}

template<typename T>
[[nodiscard]] inline constexpr T get_chrom(const uint64_t input)noexcept{
  return static_cast<T>((input|make_mask(0,58))>>(29+29));
}

template<typename T>
[[nodiscard]] inline constexpr uint64_t set_end(const uint64_t input,const T end)noexcept{
  return (~make_mask(0,29)&input)|(static_cast<uint64_t>(end));
}

template<typename T>
[[nodiscard]] inline constexpr T get_end(const uint64_t input)noexcept{
  return static_cast<T>(~make_mask(29,29)&input);
}




template<typename T>
[[nodiscard]] inline constexpr uint64_t set_pos(const uint64_t input,const T pos)noexcept{
  return (~make_mask(16,43)&input)|(static_cast<uint64_t>(pos) << (16));
}

template<typename T>
[[nodiscard]] inline constexpr T get_pos(const uint64_t input)noexcept{
  // clear bits except for offset 16 length 43
  return static_cast<T>((make_mask(16,43)&input)>>(16));
}

template<typename T>
[[nodiscard]] inline constexpr uint64_t set_ref(const uint64_t input,const T ref)noexcept{
  return (~make_mask(8,8)&input)|(static_cast<uint64_t>(ref) << (8));
}


[[nodiscard]] inline constexpr uint64_t clear_alleles(const uint64_t input)noexcept{
  return ~make_mask(0,16)&input;
}

template<typename T>
[[nodiscard]] inline constexpr T get_ref(const uint64_t input)noexcept{
  // clear bits except for offset 8 length 8
  return static_cast<T>((make_mask(8,8)&input)>>(8));
}

template<typename T>
[[nodiscard]] inline constexpr uint64_t set_alt(const uint64_t input,const T alt)noexcept{
  return (~make_mask(0,8)&input)|(static_cast<uint64_t>(alt));
}

template<typename T>
[[nodiscard]] inline constexpr T get_alt(const uint64_t input)noexcept{
  // clear bits except for offset 0 length 8
  return static_cast<T>((make_mask(0,8)&input));
}

template<typename T>
[[nodiscard]] inline constexpr uint64_t set_start(const uint64_t input,const T start)noexcept{
  return (~make_mask(29,29)&input)|(static_cast<uint64_t>(start) << (29));
}

template<typename T>
[[nodiscard]] inline constexpr T get_start(const uint64_t input)noexcept{
  return static_cast<T>((~make_mask(58,6)&input)>>(29));
}

// static_assert(0x0000000F << 4 == 0x000000F0);
// static_assert(0x0000000F << 8 == 0x00000F00);
// static_assert(make_mask(0, 4) == 0x0000000F);
// static_assert(make_mask(0, 8) == 0x000000FF);

static_assert(get_start<int>(set_start(288230376151711745ull, 536870911)) == 536870911);
static_assert(get_start<int>(set_start(18446744073709551615ull, 0)) == 0);
static_assert(get_start<int>(set_start(288230376151711745ull, 24000)) == 24000);
static_assert(get_chrom<int>(set_chrom(288230376151711745ull, 24)) == 24);


// static_assert(get_pos<int>(set_pos(288230376151711745ull, 536870911)) == 536870911);
// static_assert(get_pos<int>(set_pos(18446744073709551615ull, 0)) == 0);
// static_assert(get_pos<int>(set_pos(288230376151711745ull, 24000)) == 24000);
// static_assert(get_pos<int>(set_pos(288230376151711745ull, 24)) == 24);


// static_assert(get_ref<int>(set_ref(288230376151711745ull, 11)) == 11);
// static_assert(get_ref<int>(set_ref(18446744073709551615ull, 0)) == 0);
// static_assert(get_ref<int>(set_ref(288230376151711745ull, 15)) == 15);
// static_assert(get_ref<int>(set_ref(288230376151711745ull, 24)) == 24);


static_assert(get_alt<int>(set_alt(288230376151711745ull, 11)) == 11);
static_assert(get_alt<int>(set_alt(18446744073709551615ull, 0)) == 0);
static_assert(get_alt<int>(set_alt(288230376151711745ull, 15)) == 15);
static_assert(get_alt<int>(set_alt(288230376151711745ull, 24)) == 24);


static_assert(get_end<int>(set_end(0, 10)) == 10);
static_assert(get_end<int>(set_end(288230376151711745ull, 536870911)) == 536870911);
static_assert(get_end<int>(set_end(288230376151711745ull, 240)) == 240);
static_assert(get_end<int>(set_end(18446744073709551615ull, 24000)) == 24000);
static_assert(get_end<int>(set_end(288230376151711745ull, 24000)) == 24000);
static_assert(get_end<int>(set_end(288230376151711745ull, 24)) == 24);

[[nodiscard]]
inline constexpr uint64_t make_ldmap_region(const int chrom, const unsigned int start, const unsigned int stop) noexcept{
  //take first 29 bits from stop
  // next 29 bits from start
  // next 6 bits from chrom
  uint64_t retval=0;
  return(set_end(set_start(set_chrom(retval,chrom),start),stop));
}

[[nodiscard]]
inline constexpr uint64_t make_ldmap_snp(const int chrom, const unsigned int pos,const unsigned char ref='\0', const unsigned char alt='\0') noexcept{
  //take first 29 bits from stop
  // next 29 bits from start
  // next 6 bits from chrom
  //  std::cerr<<std::endl;
  uint64_t retval=0;
  retval = set_ref(retval,ref);
  //  std::cerr<<"ref: "<<static_cast<int>(ref)<<" "<<std::bitset<8>(ref)<<" "<<std::bitset<64>(retval)<<std::endl;
  retval = set_alt(retval,alt);
  //  std::cerr<<"alt: "<<static_cast<int>(alt)<<" "<<std::bitset<8>(alt)<<" "<<std::bitset<64>(retval)<<std::endl;
  retval = set_pos(retval,pos);
  //  std::cerr<<"pos: "<<static_cast<int>(pos)<<" "<<std::bitset<43>(pos)<<" "<<std::bitset<64>(retval)<<std::endl;
  retval = set_chrom(retval,chrom);
  //  std::cerr<<"chrom: "<<chrom<<" "<<std::bitset<64>(retval)<<std::endl;
  return retval;
}

inline std::tuple<int, int, int>
get_ldmap_region(const double x) = delete;

inline std::tuple<int, int, int> get_ldmap_region(const float x) = delete;

inline std::tuple<int, int, int>
get_ldmap_region(const long double x) = delete;

[[nodiscard]]
inline std::tuple<int,int,unsigned char,unsigned char> get_ldmap_snp(const uint64_t x) noexcept{
  return {get_chrom<int>(x),get_pos<int>(x),get_ref<unsigned char>(x),get_alt<unsigned char>(x)};
}

[[nodiscard]]
inline
std::tuple<int,int,int> get_ldmap_region(const uint64_t x) noexcept{
  return {get_chrom<int>(x),get_start<int>(x),get_end<int>(x)};
}



// static_assert(get_chrom<int>(make_ldmap_region(10, 100, 1000)) == 10);
// static_assert(get_chrom<int>(make_ldmap_region(11, 120, 3456)) == 11);
// static_assert(get_start<int>(make_ldmap_region(25, 120, 3456)) == 120);
// static_assert(get_end<int>(make_ldmap_region(25, 120, 3456)) == 3456);
// static_assert(make_ldmap_region(25, 120, 3456) <
//               make_ldmap_region(25, 120, 3457));
// static_assert(make_ldmap_region(24, 120, 3456) <
//               make_ldmap_region(25, 120, 3456));
// static_assert(make_ldmap_region(24, 110, 3456) <
//               make_ldmap_region(25, 120, 3457));









class Region;

class SNP{
public:
  uint64_t snp;
  static constexpr SNP  make_SNP(const uint64_t &x) noexcept{
    return SNP(x);
  }
  static constexpr SNP  make_SNP(const float &x) = delete;

  constexpr SNP(uint64_t x) noexcept :snp(x){};
  SNP(double x) noexcept :snp(bit_cast<uint64_t>(x)){};
  SNP() noexcept{};

  static  SNP  make_SNP(const double &x) noexcept {
    return make_SNP(bit_cast<uint64_t>(x));
  }
  static constexpr SNP  make_SNP(const long double &x) = delete;


  constexpr unsigned char chrom() const noexcept{
    return get_chrom<unsigned char>(snp);
  }
  constexpr unsigned char ref() const noexcept{
    return get_ref<unsigned char>(snp);
  }
  constexpr unsigned char alt() const noexcept{
    return get_alt<unsigned char>(snp);
  }
  constexpr int pos() const noexcept {
    return get_pos<int>(snp);
  }
  constexpr int start() const noexcept {
    return get_pos<int>(snp);
  }
  constexpr int end() const noexcept{
    return start()+1;
  }

  boost::icl::discrete_interval<uint64_t> interval() const{
    return boost::icl::construct<boost::icl::discrete_interval<uint64_t>>(start(),end(),boost::icl::interval_bounds::left_open());
  }

  template<bool NA2N>
  static constexpr SNP  make_SNP(const unsigned char chrom, const uint64_t pos, const Nuc ref, const Nuc alt) noexcept{
    SNP snp(make_ldmap_snp(chrom,pos,ref.let,alt.let));
    if constexpr(NA2N){
      if(ref.is_NA()){
        snp.snp=set_ref(snp.snp,Nuc{'N'_N}.let);
      }
      if(alt.is_NA()){
        snp.snp=set_alt(snp.snp,Nuc{'N'_N}.let);
      }
    }
    return(SNP(snp.snp));
  }


  template<bool NA2N>
  static constexpr SNP  make_SNP(const unsigned char chrom, const uint64_t pos, const Nuc ref) noexcept{
    SNP tsnp(make_ldmap_snp(chrom,pos,ref.let,'\0'));
    if constexpr(NA2N){
      if(ref.is_NA()){
        tsnp.snp=set_ref(tsnp.snp,'N'_N.let);
      }
      if(Nuc{static_cast<char>(tsnp.alt())}.is_NA()){
        tsnp.snp = set_alt(tsnp.snp,'N'_N.let);
      }
    }
    return(tsnp);
  }


  template<bool NA2N>
  static constexpr  SNP  make_SNP(const unsigned char chrom, const uint64_t pos) noexcept{
    SNP snp(make_ldmap_snp(chrom,pos,'\0','\0'));
    if constexpr(NA2N){
      snp.snp=set_alt(set_ref(snp.snp,'N'_N.let),'N'_N.let);
    }
    return(SNP(snp));
  }

  constexpr bool operator!=(const SNP &other) const noexcept{
    return snp!=other.snp;
  };
  constexpr bool operator==(const SNP &other) const noexcept{
    return snp==other.snp;
  };
  constexpr bool operator<(const SNP &b) const{

    return (clear_alleles(snp)<clear_alleles(b.snp));
  }
  constexpr bool operator>(const SNP &b) const{
    return (clear_alleles(snp)>clear_alleles(b.snp));
  }
  constexpr bool operator>=(const SNP &b) const{
    return clear_alleles(snp)>=clear_alleles(b.snp);
  }
  constexpr bool operator<=(const SNP &b) const{
    return clear_alleles(snp)<=clear_alleles(b.snp);
  }
  constexpr int distance(const SNP &other) const{
    if(other.chrom()>chrom())
      return std::numeric_limits<int>::max();
    if(other.chrom()<chrom())
      return std::numeric_limits<int>::min();
    return pos()-other.pos();
  }

  constexpr bool overlap(const SNP &other) const noexcept{
    return clear_alleles(snp)==clear_alleles(other.snp);
  }

  constexpr int distance(const Region &other) const noexcept;
  constexpr bool overlap(const Region &other) const noexcept;
  constexpr bool contains(const Region &other) const noexcept;

  constexpr bool operator<(const Region &other) const noexcept;
  constexpr bool operator>(const Region &other) const noexcept;
  constexpr bool operator==(const Region &other) const noexcept;
  constexpr bool operator!=(const Region &other) const noexcept;

  constexpr SNP start_SNP() const noexcept{
    return *this;
  }
  constexpr SNP end_SNP() const noexcept{
    return SNP(clear_alleles(set_pos(snp,pos()+1)));
  }

  double to_double() const{
    return bit_cast<double>(snp);
  }

  bool is_strand_ambiguous() const{
    return (ref()== 'A'_L and alt() == 'T'_L) or
      (ref()== 'T'_L and alt() == 'A'_L) or
      (ref()== 'G'_L and alt() == 'C'_L) or
      (ref()== 'C'_L and alt() == 'G'_L);
  }

  // Assume LHS is the query and RHS is a potential target.

  std::optional<Snp_match> allele_match(const SNP &other)const{
    if(clear_alleles(snp) != clear_alleles(other.snp))
      return std::nullopt;

    const auto lref=ref();
    const auto lalt=alt();

    const auto rref = other.ref();
    const auto ralt = other.alt();

    if(auto res = is_perfect(lref,lalt,rref,ralt)){
      return res;
    }
    if(auto res = is_reverse(lref,lalt,rref,ralt)){
      return res;
    }
    if(auto res = is_complement(lref,lalt,rref,ralt)){
      return res;
    }
    if(auto res = is_reverse_complement(lref,lalt,rref,ralt)){
      return res;
    }
    return std::nullopt;
  }
  friend Region;
  friend std::ostream & operator << (std::ostream &out, const SNP &c);
};


inline std::ostream& operator<<( std::ostream& out, const SNP &lra ){

  out<<"chr"<<static_cast<int>(lra.chrom())<<":"<<lra.pos()<<"_"<<static_cast<int>(lra.ref())<<"."<<static_cast<int>(lra.alt());
  return out;
}

class Region{
  public:
  uint64_t br;

  constexpr unsigned char chrom() const noexcept{
    return get_chrom<unsigned char>(br);
  }
  constexpr int start() const noexcept {
    return get_start<int>(br);
  }
  constexpr int end() const noexcept{
    auto ret = get_end<int>(br);
    return ret;
  }
  boost::icl::discrete_interval<uint64_t> interval() const{
    return boost::icl::construct<boost::icl::discrete_interval<uint64_t>>(start(),end(),boost::icl::interval_bounds::left_open());
  }

  template<typename C,typename S>
  static constexpr Region make_Region(const C chrom,const S start, const S end) noexcept{
    return Region(make_ldmap_region(chrom,start,end));
  }
  static constexpr Region make_Region(const uint64_t x) noexcept{
    return Region(x);
  }
  static Region make_Region(const double x){
    return Region(bit_cast<uint64_t>(x));
  }
  constexpr Region(uint64_t x) noexcept :br(x){}
  Region() noexcept {}
  Region(double x) noexcept :br(bit_cast<uint64_t,double>(x)){}
  Region(int x) = delete;
  Region(float x) = delete;
  Region(long double x) = delete;


  static Region make_Region(const long double x) = delete;
  static Region make_Region(const float x) = delete;
  static Region make_Region(const int x) = delete;

  double to_double() const{
    return bit_cast<double>(br);
  }
  static constexpr Region make_Region(const SNP& a,const SNP& b) noexcept{
    if( a.chrom()!=b.chrom())
      return make_Region(0ul);
    const auto pai=a.pos();
    const auto pab=b.pos();
    auto [pa,pb] = std::minmax(pai,pab);
    return make_Region<unsigned char,int>(a.chrom(),pa,pb+1);
  }
  constexpr SNP start_SNP() const noexcept{
    return SNP::make_SNP<true>(chrom(),start());
  }
  constexpr SNP end_SNP() const noexcept{
    return SNP::make_SNP<true>(chrom(),end());
  }
  constexpr SNP last_SNP() const noexcept{
    return SNP::make_SNP<true>(chrom(),end()-1);
  }
  constexpr bool operator ==(const Region& other)const noexcept{
    return this->br==other.br;
  }
  constexpr bool operator!=(const Region& other)const noexcept{
    return this->br!=other.br;
  }
  constexpr bool starts_before(const Region& other)const noexcept{
    if(chrom() > other.chrom())
      return false;

    if(chrom() < other.chrom()){
      return true;
    }
    return start() < other.start();
  }

  constexpr bool ends_before(const Region& other)const noexcept{
    if(get_chrom<unsigned char>(this->br) > other.chrom())
      return false;

    if(get_chrom<unsigned char>(this->br) < other.chrom()){
      return true;
    }
    return this->end() < other.end();
  }

  constexpr bool operator<(const Region& other)const noexcept{
    return this->br < other.br;
  }
  constexpr bool operator>(const Region& other)const noexcept{
    return this->br > other.br;
  }

  constexpr Region operator|(const Region &other) const noexcept{

    if(get_chrom<unsigned char>(this->br)!=other.chrom()){
      return Region(0ul);
    }
    auto [sa,sb] = std::minmax(*this,other);
    if( sa.end() < sb.start()){
      return Region(0ul);
    }
    return Region(make_ldmap_region(sa.chrom(),sa.start(),sb.end()));
  }
  constexpr Region operator|=(const Region &other) noexcept{

    if(get_chrom<unsigned char>(this->br)!=other.chrom()){
      this->br=0;
      return *this;
    }
    auto [sa,sb] = std::minmax(*this,other);
    if( sa.end() < sb.start()){
      this->br=0;
      return *this;
    }
    this->br=make_ldmap_region(sa.chrom(),sa.start(),sb.end());
    return *this;
  }

  constexpr bool overlap(const Region &other) const noexcept{

    if(get_chrom<unsigned char>(this->br)!=other.chrom())
      return false;
    return start() < other.end() && other.start() < end();
  }

  constexpr int distance(const SNP &other) const noexcept{
    if(other.chrom()>chrom())
      return std::numeric_limits<int>::max();
    if(other.chrom()<chrom())
      return std::numeric_limits<int>::min();
    if (start() <= other.pos()){
      if(other.pos() < end())
        return 0;
      else
        return other.pos()-end();
    }
    return other.pos()-start();
  }


  constexpr int distance(const Region &other) const noexcept{
    if(other.chrom()>chrom())
      return std::numeric_limits<int>::max();
    if(other.chrom()<chrom())
      return std::numeric_limits<int>::min();
    if(overlap(other))
      return 0;
    if(other.start()>=end())
      return other.start()-(end()-1);
    return other.end()-start();
  }


  constexpr bool contains(const Region &other) const noexcept{

    if(get_chrom<unsigned char>(this->br)!=other.chrom())
      return false;
    return start() <= other.start() && other.end() <= end();
  }

  constexpr bool contains(const SNP &other) const noexcept{
    return overlap(other);
  }

  constexpr bool overlap(const SNP &other) const noexcept{

    if(get_chrom<unsigned char>(this->br)!=other.chrom())
      return false;
    return (this->start() <= other.pos() && other.pos() < this->end());
  }

  constexpr bool operator==(const SNP& other)const{
    return (chrom()==other.chrom()) and (this->start() <= other.pos()) and (this->end() > other.pos());
  }
  constexpr bool operator<(const SNP& other)const{
    return std::make_pair(chrom(),end()) <= std::make_pair(other.chrom(),other.pos());
  }

  constexpr bool operator>(const SNP& other)const{
    return std::make_pair(chrom(),start()) > std::make_pair(other.chrom(),other.pos());
  }

  friend SNP;
  friend std::ostream & operator << (std::ostream &out, const Region &c);

};


inline constexpr int SNP::distance(const Region &other) const noexcept{
    if(other.chrom()>chrom())
      return std::numeric_limits<int>::max();
    if(other.chrom()<chrom())
      return std::numeric_limits<int>::min();
    if(overlap(other))
      return 0;
    if(abs(pos()-other.start())>abs(pos()-other.end()))
      return pos()-other.end();
    return pos()-other.start();
}

inline constexpr bool SNP::overlap(const Region &other) const noexcept {
  if(chrom()!=other.chrom())
    return false;
  return (other.start() <= pos() && pos() < other.end());
}

inline constexpr bool SNP::contains(const Region &other) const noexcept {
  if(chrom()!=other.chrom())
    return false;
  return (other.start() == pos());
}



inline constexpr bool SNP::operator==(const Region &other) const noexcept{
  return (other.chrom()==this->chrom()) and (other.start() <= this->pos()) and (other.end() > this->pos());
}
// inline constexpr bool SNP::operator==(const  &other) const{
//   return (other.chrom()==this->chrom()) and (other.start() <= this->pos()) and (other.end() > this->pos());
// }
inline constexpr bool SNP::operator!=(const Region &other) const noexcept{
  return !(*this==other);
}
inline constexpr bool SNP::operator<(const Region &other) const noexcept{
  return (other.chrom()==this->chrom()) and (this->pos() < other.start());
}
inline constexpr bool SNP::operator>(const Region &other) const noexcept{
  return (other.chrom()==this->chrom()) and (this->pos() >= other.end());
}


inline std::ostream& operator<<( std::ostream& out, const Region &lra ){

  out<<"chr"<<static_cast<int>(lra.chrom())<<":"<<lra.start()<<"-"<<lra.end();
  return out;
}

inline Snp reverse_snp(Snp a){
  auto tmp=a.str.ref;
  a.str.ref=a.str.alt;
  a.str.alt=tmp;
  return(a);
}

inline Snp complement_snp(Snp a ){
  a.str.ref = complement(Nuc{a.str.ref}).let;
  a.str.alt = complement(Nuc{a.str.alt}).let;
  return(a);
}



class RegionSet{
  Region convex_hull;
  std::vector<Region> regions;
public:
  RegionSet() noexcept: convex_hull(0ul),regions(){};
  RegionSet(Region reg) noexcept: convex_hull(reg),regions{reg}{};
  std::optional<Region> insert(const Region reg){
    if(!convex_hull.overlap(reg))
      return std::nullopt;
    regions.push_back(reg);
    convex_hull=Region::make_Region(convex_hull.chrom(),std::min(convex_hull.start(),reg.start()),std::max(convex_hull.end(),reg.end()));
    return convex_hull;
  }
  const std::vector<Region>& get_regions() const{
    return regions;
  }
  const Region& get_region() const {
    return convex_hull;
  }
};

class RegionSets{
  std::vector<RegionSet> rsets;
  size_t count;
public:
  RegionSets():count(0),rsets(){};
  size_t insert_region(Region r){
    for(auto trs =rsets.rbegin(); trs!=rsets.rend(); trs++){
      if(auto newreg = trs->insert(r))
        return count;
    }
    rsets.emplace_back(r);
    return count++;
  }
  size_t get_count() const {
    return count;
  }
  const std::vector<RegionSet>& get_set() const{
    return rsets;
  }

};

// constexpr bool do_tests(){
//   const auto a= Region::make_Region(1,1,1892607);
//   const auto b=Region::make_Region(1,1892607,3582736);

//   static_assert(Region::make_Region(1,1892607,3582736).chrom()==Region::make_Region(1,1,1892607).chrom());
//   return a==b;
// }

static_assert(!Region::make_Region(1,1,1892607).overlap(Region::make_Region(1,1892607,3582736)));

static_assert(Region::make_Region(SNP::make_SNP<true>(1,50),SNP::make_SNP<true>(1,55)).start_SNP()==SNP::make_SNP<true>(1,50));

static_assert(make_ldmap_region(24, 536870911, 3456) <
              make_ldmap_region(25, 120, 3457));
static_assert(120 - 100 == 20);


static_assert(Region::make_Region(23, 100, 200)
              .distance(SNP(make_ldmap_snp(23, 120))) == 0);
static_assert(Region::make_Region(23, 100, 200)
                  .distance(SNP(make_ldmap_snp(23, 220))) == 20);
static_assert(Region::make_Region(23, 100, 200)
                  .distance(SNP(make_ldmap_snp(23, 90))) == -10);
static_assert(Region::make_Region(23, 100, 200)
                  .distance(SNP(make_ldmap_snp(23, 90))) == -10);

static_assert(Region::make_Region(23, 100, 200)
                  .distance(Region::make_Region(23, 100, 200)) == 0);
static_assert(Region::make_Region(23, 100, 200)
                  .distance(Region::make_Region(23, 120, 250)) == 0);
static_assert(Region::make_Region(23, 100, 200)
                  .distance(Region::make_Region(23, 90, 250)) == 0);

static_assert(Region::make_Region(23, 100, 200)
                  .distance(Region::make_Region(23, 90, 91)) == -9);

static_assert(Region::make_Region(23, 100, 200)
              .distance(Region::make_Region(23, 201, 202)) == 2);
