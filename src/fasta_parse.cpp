#include <Rcpp.h>
#include "ctre.hpp"
#include "alleles.hpp"
#include <charconv>
// [[Rcpp::plugins(cpp17)]]

static constexpr ctll::fixed_string range_pattern =  R"(^(chr)?([0-9XY]+):([0-9]+)[\-_:]([0-9]+)$)";

static constexpr ctll::fixed_string  snp_pattern= R"(^(chr)?([0-9XY]+):([0-9]+)[\-_:]([ACGTMRWSYKVHDBN])[ACGTMRWSYKVHDBN]*[\-_:]([ACGTMRWSYKVHDBN])[ACGTMRWSYKVHDBN]*$)";


std::optional<Region> parse_Region(std::string_view sr){

  int chrom;
  int start;
  int stop;

  if (auto m = ctre::match<range_pattern>(sr)) {
    auto sv = m.get<2>().to_view();
    if (auto [p, ec] =std::from_chars(sv.begin(), sv.end(), chrom); ec != std::errc())
      chrom = sv=="X" ? 23 : 24;
    sv=m.get<3>().to_view();
    if (auto [p, ec] =std::from_chars(sv.begin(), sv.end(), start); ec != std::errc())
      return std::nullopt;
    sv=m.get<4>().to_view();
    if (auto [p, ec] =std::from_chars(sv.begin(), sv.end(), stop); ec != std::errc())
      return std::nullopt;
    return Region::make_Region(chrom,start,stop);
  } else {
    return std::nullopt;
  }
}




std::optional<SNP> parse_SNP(std::string_view sr){

  int chrom;
  int pos;
  char ref;
  char alt;

  if (auto m = ctre::match<snp_pattern>(sr)) {
    auto sv = m.get<2>().to_view();
    if (auto [p, ec] =std::from_chars(sv.begin(), sv.end(), chrom); ec != std::errc())
      chrom = sv=="X" ? 23 : 24;
    sv=m.get<3>().to_view();
    if (auto [p, ec] =std::from_chars(sv.begin(), sv.end(), pos); ec != std::errc())
      return std::nullopt;
    sv=m.get<4>().to_view();
    ref=ascii2Nuc(*sv.begin());
    sv=m.get<5>().to_view();
    alt=ascii2Nuc(*sv.begin());

    return SNP::make_SNP<false>(static_cast<unsigned char>(chrom),
                                static_cast<uint64_t>(pos),
                                Nuc{ref},
                                Nuc{alt});
  } else {
    return std::nullopt;
  }
}






//' Creation of new ldmap_region vector from character
//'
//' @param input a character vector
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector parse_ldmap_region(Rcpp::StringVector input){

  Rcpp::NumericVector ret=Rcpp::NumericVector::import_transform(
                                                                input.begin(),
                                                                input.end(),[](SEXP inp){
                                                                              const size_t p= LENGTH(inp);
                                                                              const char* charp=CHAR(inp);
                                                                              std::string_view sv(charp,p);
                                                                              if(auto reg = parse_Region(sv))
                                                                                return(bit_cast<double>(reg->br));
                                                                              return NA_REAL;
                                                                            });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_region","vctrs_vctr");
  return ret;
}





//' Creation of new ldmap_snp vector from character
//'
//' @param input a character vector
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector parse_ldmap_SNP(Rcpp::StringVector input){



  Rcpp::NumericVector ret=Rcpp::NumericVector::import_transform(
                                                                input.begin(),
                                                                input.end(),[](SEXP inp){
                                                                              const size_t p= LENGTH(inp);
                                                                              const char* charp=CHAR(inp);
                                                                              std::string_view sv(charp,p);
                                                                              if(auto snp=parse_SNP(sv))
                                                                                return(bit_cast<double>(snp->snp));
                                                                              return NA_REAL;
                                                                            });
  ret.attr("class")=Rcpp::StringVector::create("ldmap_snp","vctrs_vctr");
  return ret;
}
