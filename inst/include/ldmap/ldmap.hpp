#pragma once
#ifndef EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#endif
#include <RcppEigen.h>
#include <fmt/format.h>
enum class Ref :char{ A = 0, C = 1, G = 2, T = 3, N = 4 };

inline Ref complement(const Ref l){

  switch (l){
  case Ref::A:
    return Ref::T;
  case Ref::T:
    return Ref::A;
  case Ref::C:
    return Ref::G;
  case Ref::G:
    return Ref::C;
  case Ref::N:
    return Ref::N;
  }
}

inline Ref ascii2Ref(const char l){

  switch (l){
  case 'A':
    return Ref::A;

  case 'T':
    return Ref::T;

  case 'C':
    return Ref::C;

  case 'G':
    return Ref::G;

  default :
    return Ref::N;
  }
}



inline Ref ascii2Ref(const int l){

  switch (l){
  case 65:
    return Ref::A;

  case 67:
    return Ref::C;

  case 84:
    return Ref::T;

  case 71:
    return Ref::G;

  default :
    return Ref::N;
  }
}


inline int Ref2int(const Ref l){

  switch (l){
  case Ref::T:
    return int('T');

  case Ref::A:
    return int('A');

  case Ref::C:
    return int('C');

  case Ref::G:
    return int('G');
  case Ref::N:
  default:
    return int('N');
  }
}

inline char Ref2char(const Ref l){

  switch (l){
  case Ref::T:
    return 'T';

  case Ref::A:
    return 'A';

  case Ref::C:
    return 'C';

  case Ref::G:
    return 'G';
  case Ref::N:
  default:
    return 'N';
  }
}



inline const char* Ref2ascii(const Ref l){

  switch (l){
  case Ref::T:
    return "T";

  case Ref::A:
    return "A";

  case Ref::C:
    return "C";

  case Ref::G:
    return "G";
  case Ref::N:
  default:
    return "N";
  }
}



// namespace fmt {
// template <>
// struct formatter<Ref>:formatter<std::string> {

//   template <typename FormatContext>
//   auto format(const Ref &p, FormatContext &ctx) {
//     std::string name = "N";
//       switch (p){
//       case Ref::T:
//       name="T"; break;

//       case Ref::A:
//       name="A"; break;

//       case Ref::C:
//       name="C";break;

//       case Ref::G:
//       name="G";break;
//       case Ref::N:
//       default:
//       name="N";break;
//       }
//       return formatter<std::string>::format(name, ctx);
//   }
// };
// }





struct bit_snp{
  unsigned char chrom : 5;
  uint64_t pos :43;
  Ref ref : 8;
  Ref alt : 8;
} __attribute__((packed));

union Snp{
  bit_snp str;
  double dat;
};

enum class Snp_match {no_match = 0,
                      mismatch = 1,      // mismatch
                       one_mono_allelic=2, //mismatch
                       reverse_ref_alt=3,  //match -> reverse
                       strand_flip_match=4,//match -> complement
                       strand_flip_reverse_match=5,//match -> reverse
                      snp_match=6};


class SNP{
public:
  Snp snp;
  SNP(const unsigned char chrom, uint64_t pos, Ref ref, Ref alt):
    snp({chrom,pos,ref,alt}){}
  SNP(double dat){
    snp.dat=dat;
  }
  bool operator<(const SNP &b) const{
    return (std::tie(snp.str.chrom,snp.str.pos) < std::tie(b.snp.str.chrom,b.snp.str.pos));
  }
  bool operator>=(const SNP &b) const{
    return (std::tie(snp.str.chrom,snp.str.pos) >= std::tie(b.snp.str.chrom,b.snp.str.pos));
  }
  bool operator<=(const SNP &b) const{
    return (std::tie(snp.str.chrom,snp.str.pos) <= std::tie(b.snp.str.chrom,b.snp.str.pos));
  }
  bool is_ambiguous() const{
    return (snp.str.ref== Ref::A and snp.str.alt == Ref::T) or
      (snp.str.ref== Ref::T and snp.str.alt == Ref::A) or
      (snp.str.ref== Ref::G and snp.str.alt == Ref::C) or
      (snp.str.ref== Ref::C and snp.str.alt == Ref::G);
  }
  Snp_match allele_match(SNP other){
    
    if((other.snp.str.ref==Ref::N) or (this->snp.str.ref==Ref::N) or (this->snp.str.alt == Ref::N) or (other.snp.str.alt == Ref::N) or (this->snp.str.ref==this->snp.str.alt) or (other.snp.str.ref ==other.snp.str.alt))
    return Snp_match::one_mono_allelic;

  // #     # GOOD/TRUE: same reference/alternative alleles between datasets
  // #     (Var1 == Var3) & (Var2 == Var4) ~ TRUE,
  if(this->snp.str.ref==other.snp.str.ref and this->snp.str.alt ==other.snp.str.alt)
    return Snp_match::snp_match;

  // #     # GOOD/FALSE: reverse reference/alternative alleles
  // #     (Var1 == Var4) & (Var2 == Var3) ~ FALSE,
  if( (this->snp.str.ref == other.snp.str.alt) and (this->snp.str.alt == other.snp.str.ref))
    return Snp_match::reverse_ref_alt;

  // #     # GOOD/TRUE: same reference/alternative alleles after strand flip
  // #     (REV_ACTG[Var1] == Var3) & (REV_ACTG[Var2] == Var4) ~ TRUE,
  if( complement(this->snp.str.ref) == other.snp.str.ref and complement(this->snp.str.alt) ==other.snp.str.alt)
    return Snp_match::strand_flip_match;

  // #     # GOOD/FALSE: reverse reference/alternative alleles after strand flip
  // #     (REV_ACTG[Var1] == Var4) & (REV_ACTG[Var2] == Var3) ~ FALSE,
  if( complement(this->snp.str.ref) == other.snp.str.alt and complement(this->snp.str.alt) ==other.snp.str.alt)
    return Snp_match::strand_flip_reverse_match;

  return Snp_match::mismatch;
  }
  Rcpp::String char_rep()const{
    if(snp.str.pos==0 && snp.str.chrom==0){
      return R_NaString;
    }else{

    return "chr"+
      std::to_string(snp.str.chrom) +
      ":" +
      std::to_string(static_cast<uint64_t>(snp.str.pos)) +
      "_"+
      std::string(1,Ref2char(snp.str.ref))+
      "_"+
      std::string(1,Ref2char(snp.str.alt));

    }
  }
  
  
};








inline Snp reverse_snp(Snp a ){
  Ref tmp=a.str.ref;
  a.str.ref=a.str.alt;
  a.str.alt=tmp;
  return(a);
}

inline Snp complement_snp(Snp a ){
  a.str.ref = complement(a.str.ref);
  a.str.alt = complement(a.str.alt);
  return(a);
}


