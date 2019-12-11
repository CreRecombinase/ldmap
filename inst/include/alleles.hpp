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
#include <tuple>
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
  polymorphism as possible using 64 bit integers
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


struct bit_range{
  uint64_t end :29;
  uint64_t start :29;
  unsigned char chrom : 6;
} __attribute__((packed));


union bed_range{
  bit_range str;
  std::int64_t dat;
  double flt;
};

struct bit_snp{
  char alt : 8;
  char ref : 8;
  uint64_t pos :43;
  unsigned char chrom : 5;
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

class Region;

class SNP{
public:
  Snp snp;

    static constexpr SNP  make_snp(const double &x) noexcept{
      return SNP{.snp={.flt=x}};
    }

  constexpr unsigned char chrom() const noexcept{
    return snp.str.chrom;
  }

  constexpr uint64_t start() const noexcept {
    return snp.str.pos;
  }

  constexpr uint64_t end() const noexcept{
    return snp.str.pos+1;
  }

  boost::icl::discrete_interval<uint64_t> interval() const{
    return boost::icl::construct<boost::icl::discrete_interval<uint64_t>>(snp.str.pos,snp.str.pos+1,boost::icl::interval_bounds::left_open());
  }

  template<bool NA2N>
  static constexpr SNP  make_snp(const unsigned char chrom, const uint64_t pos, const Nuc ref, const Nuc alt) noexcept{
    Snp snp{.str={.alt=alt.let,.ref=ref.let,.pos=pos,.chrom=chrom}};
    if constexpr(NA2N){
      if(ref.is_NA()){
        snp.str.ref='N'_L;
      }
      if(alt.is_NA()){
        snp.str.alt='N'_L;
      }
    }
    return(SNP{.snp=snp});
  }

  template<bool NA2N>
  static constexpr SNP  make_snp(const unsigned char chrom, const uint64_t pos, const Nuc ref) noexcept{
    Snp snp{.str={.alt='\0',.ref=ref.let,.pos=pos,.chrom=chrom}};
    if constexpr(NA2N){
      if(ref.is_NA()){
        snp.str.ref='N'_L;
      }
      if(Nuc{snp.str.alt}.is_NA()){
        snp.str.alt='N'_L;
      }
    }
    return(SNP{.snp=snp});
  }


  template<bool NA2N>
  static constexpr  SNP  make_snp(const unsigned char chrom, const uint64_t pos) noexcept{
    Snp snp{.str={.alt='\0',.ref='\0',.pos=pos,.chrom=chrom}};
    if constexpr(NA2N){
        snp.str.ref='N'_L;
        snp.str.alt='N'_L;
    }
    return(SNP{.snp=snp});
  }

  constexpr bool operator<(const SNP &b) const{
    return (std::tie(snp.str.chrom,snp.str.pos) < std::tie(b.snp.str.chrom,b.snp.str.pos));
  }
  constexpr bool operator>=(const SNP &b) const{
    return (std::tie(snp.str.chrom,snp.str.pos) >= std::tie(b.snp.str.chrom,b.snp.str.pos));
  }
  constexpr bool operator<=(const SNP &b) const{
    return (std::tie(snp.str.chrom,snp.str.pos) <= std::tie(b.snp.str.chrom,b.snp.str.pos));
  }

  constexpr bool operator<(const Region &other) const;
  constexpr bool operator>(const Region &other) const;
  constexpr bool operator==(const Region &other) const;
  constexpr bool operator!=(const Region &other) const;


  double to_double() const{
    return snp.flt;
  }


  bool is_strand_ambiguous() const{
    return (snp.str.ref== 'A'_L and snp.str.alt == 'T'_L) or
      (snp.str.ref== 'T'_L and snp.str.alt == 'A'_L) or
      (snp.str.ref== 'G'_L and snp.str.alt == 'C'_L) or
      (snp.str.ref== 'C'_L and snp.str.alt == 'G'_L);
  }

  // Assume LHS is the query and RHS is a potential target.

  std::optional<Snp_match> allele_match(const SNP &other)const{
    if(std::tie(this->snp.str.chrom,this->snp.str.pos) != std::tie(other.snp.str.chrom,other.snp.str.pos)){
      return std::nullopt;
    }
    const auto &lref=this->snp.str.ref;
    const auto &lalt=this->snp.str.alt;
    const auto &rref = other.snp.str.ref;
    const auto &ralt = other.snp.str.alt;

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
};


class Region{
  public:
  bed_range br;

  constexpr unsigned char chrom() const noexcept{
    return br.str.chrom;
  }
  constexpr uint64_t start() const noexcept {
    return br.str.start;
  }
  constexpr uint64_t end() const noexcept{
    return br.str.end;
  }
  boost::icl::discrete_interval<uint64_t> interval() const{
    return boost::icl::construct<boost::icl::discrete_interval<uint64_t>>(br.str.start,br.str.end,boost::icl::interval_bounds::left_open());
  }

  template<typename C,typename S>
  static constexpr Region make_Region(const C chrom,const S start, const S end) noexcept{
    return Region{.br={.str={.end=end,.start=start,.chrom=chrom}}};
  }

  static constexpr Region make_Region(const double x) noexcept{
    return Region{.br={.flt=x}};
  }

  static constexpr Region make_Region(const SNP& a,const SNP& b) noexcept{
    if( a.snp.str.chrom!=b.snp.str.chrom)
      return Region{.br={.str={.end=0,
                               .start=0,
                               .chrom=0}}};
    return Region{.br={.str={.end=b.snp.str.pos+1,
                             .start=a.snp.str.pos,
                             .chrom=a.snp.str.chrom}}};
  }
  constexpr SNP start_SNP() const noexcept{

    return SNP::make_snp<true>(chrom(),start());
  }
  constexpr SNP end_SNP() const noexcept{
    return SNP::make_snp<true>(chrom(),end());
  }
  constexpr SNP last_SNP() const noexcept{
    return SNP::make_snp<true>(chrom(),end()-1);
  }

  // constexpr Region(const unsigned char chrom, const uint64_t offset,const uint64_t size)noexcept :br({chrom,offset,size}){}
  // constexpr Region(std::int64_t data):br({.dat=data}){
  // }
  constexpr bool operator ==(const Region& other)const noexcept{
    return this->br.dat==other.br.dat;
  }
  constexpr bool operator!=(const Region& other)const noexcept{
    return this->br.dat!=other.br.dat;
  }
  // constexpr bool operator<=(const Region& other)const{
  //   return this->br.dat <= other.br.dat;
  // }
  // constexpr bool operator>=(const Region& other)const{
  //   return this->br.dat >= other.br.dat;
  // }
  constexpr bool starts_before(const Region& other)const noexcept{
    if(this->br.str.chrom > other.br.str.chrom)
      return false;

    if(this->br.str.chrom < other.br.str.chrom){
      return true;
    }
    return this->br.str.start < other.br.str.start;
  }

  constexpr bool ends_before(const Region& other)const noexcept{
    if(this->br.str.chrom > other.br.str.chrom)
      return false;

    if(this->br.str.chrom < other.br.str.chrom){
      return true;
    }
    return this->br.str.end < other.br.str.end;
  }


  constexpr bool operator<(const Region& other)const noexcept{
    return this->br.dat < other.br.dat;
  }
  constexpr bool operator>(const Region& other)const noexcept{
    return this->br.dat > other.br.dat;
  }

  constexpr Region operator|(const Region &other) const noexcept{

    if(this->br.str.chrom!=other.br.str.chrom){
      return Region{.br={.dat=0}};
    }
    auto [sa,sb] = std::minmax(*this,other);
    if( sa.br.str.end < sb.br.str.start){
      return Region{.br={.dat=0}};
    }
    return Region{.br={.str={.end=sb.br.str.end,
                             .start=sa.br.str.start,
                             .chrom=sa.br.str.chrom}}};
  }
  constexpr Region operator|=(const Region &other) noexcept{

    if(this->br.str.chrom!=other.br.str.chrom){
      this->br.dat=0;
      return *this;
    }
    auto [sa,sb] = std::minmax(*this,other);
    if( sa.br.str.end < sb.br.str.start){
      this->br.dat=0;
      return *this;
    }
    this->br.str={.end=sb.br.str.end,
                  .start=sa.br.str.start,
                  .chrom=sa.br.str.chrom};
    return *this;
  }
  constexpr bool overlap(const Region &other) const noexcept{

    if(this->br.str.chrom!=other.br.str.chrom)
      return false;

    auto [sa,sb] = std::minmax(*this,other);
    return sa.br.str.end > sb.br.str.start;
  }

  constexpr bool overlap(const SNP &other) const noexcept{

    if(this->br.str.chrom!=other.snp.str.chrom)
      return false;
    return (this->br.str.start <= other.snp.str.pos && other.snp.str.pos < this->br.str.end);
  }

  std::optional<Region> operator+(const SNP &other) {

    if(this->br.str.chrom!=other.snp.str.chrom){
      return std::nullopt;
    }
    return Region{.br={.str={.end=std::max(this->br.str.end,other.snp.str.pos+1),
                             .start=std::min(this->br.str.start,other.snp.str.pos),
                             .chrom=this->br.str.chrom}}};
  }


  constexpr bool operator==(const SNP& other)const{
    return (this->br.str.chrom==other.snp.str.chrom) and (this->br.str.start <= other.snp.str.pos) and (this->br.str.end > other.snp.str.pos);
  }
  constexpr bool operator<(const SNP& other)const{
    return std::tie(this->br.str.chrom,this->br.str.end) <= std::tie(other.snp.str.chrom,other.snp.str.pos);
  }
  // constexpr bool operator<=(const SNP& other)const{
  //   return std::tie(this->br.str.chrom,this->br.str.end) < std::tie(other.snp.str.chrom,other.snp.str.pos);
  // }
  constexpr bool operator>(const SNP& other)const{
    return std::tie(this->br.str.chrom,this->br.str.start) > std::tie(other.snp.str.chrom,other.snp.str.pos);
  }

  friend SNP;

};

inline constexpr bool SNP::operator==(const Region &other) const{
  return (other.br.str.chrom==this->snp.str.chrom) and (other.br.str.start <= this->snp.str.pos) and (other.br.str.end > this->snp.str.pos);
}
// inline constexpr bool SNP::operator==(const  &other) const{
//   return (other.br.str.chrom==this->snp.str.chrom) and (other.br.str.start <= this->snp.str.pos) and (other.br.str.end > this->snp.str.pos);
// }
inline constexpr bool SNP::operator!=(const Region &other) const{
  return !(*this==other);
}
inline constexpr bool SNP::operator<(const Region &other) const{
  return (other.br.str.chrom==this->snp.str.chrom) and (this->snp.str.pos < other.br.str.start);
}
inline constexpr bool SNP::operator>(const Region &other) const{
  return (other.br.str.chrom==this->snp.str.chrom) and (this->snp.str.pos >= other.br.str.end);
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
