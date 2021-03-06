#include <Rcpp.h>
#include "ctre.hpp"
#include "alleles.hpp"
#include <fstream>
#include "fmt/format.h"
#include <zlib.h>
#include <charconv>
#include <filesystem>
#include <range/v3/core.hpp>
#include <range/v3/view/split.hpp>
#include <range/v3/view/split_when.hpp>
#include <range/v3/view/counted.hpp>
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/view/common.hpp"
#include "range/v3/view/c_str.hpp"
#if HAVE_HTS
#include "htslib/bgzf.h"
#endif
#include "range/v3/view/join.hpp"
#include "range/v3/view/transform.hpp"
#include <range/v3/action/push_back.hpp>



#if HAVE_HTS
void close_wrapper( BGZF* fp){
  int ret  = bgzf_close(fp);
}
using R_BGZF = Rcpp::XPtr<BGZF, Rcpp::PreserveStorage, &close_wrapper, false>;


#endif

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


struct starts_with_header
    {
        template<typename I, typename S>
        std::pair<bool, I> operator()(I b, S) const
        {
          return {b->front() == ">", b};
        }
    };


//[[Rcpp::export]]
Rcpp::List read_fasta_file(std::string filename){

  if(!std::filesystem::exists(filename)){
    Rcpp::stop(filename+" does not exist in read_fasta_file");
  }
  using namespace ranges;
  std::ifstream sin(filename,std::ios_base::in);
  std::vector<std::string> names;
  std::vector<std::vector<Rbyte>> fastas;
  std::string line;
  int fbmax=0;
  std::getline(sin,line);
    if(line.front()!='>'){
    Rcpp::stop(" this isn't a fasta file because the header line is: '"+line+"'");
    }
    while(sin){
      names.emplace_back(line.begin()+1,line.end());
      std::vector<Rbyte> &fb = fastas.emplace_back();
      fb.reserve(fbmax);
      std::getline(sin,line);
      while(line.front() != '>' && !sin.eof()){
        fb |= ranges::action::push_back(line | ranges::view::transform([](const char x) -> Rbyte {
                                                                         return ascii2Nuc(x);
                                                                       }));
        std::getline(sin,line);
      }
    }

    Rcpp::List retl = Rcpp::List::import_transform(fastas.begin(),fastas.end(),[](const std::vector<Rbyte>& x) -> Rcpp::RawVector{
                                                                                 Rcpp::RawVector rx = Rcpp::wrap(x);
                                                                                 rx.attr("class")=Rcpp::StringVector::create("ldmap_allele","vctrs_vctr");

                                                                                 return rx;
                                                                               });
    retl.attr("names") = Rcpp::StringVector(Rcpp::wrap(names));
    return retl;

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



//[[Rcpp::export]]
SEXP open_bgzf(std::string fstr,bool read_only = true) {

  #if HAVE_HTS
  if(read_only){
    if(std::filesystem::exists(fstr)){
      auto fp = bgzf_open(fstr.c_str(),"r");
      if(std::filesystem::exists(fstr+".gzi")){
        int ret = bgzf_index_load(fp,fstr.c_str(),".gzi");
        if (ret <0)
          Rcpp::warning("error loading index...");
      }
      return R_BGZF(fp);
    }else{
      Rcpp::stop("file not found: "+fstr);
    }
  }
  return R_BGZF(bgzf_open(fstr.c_str(),"w"));
  #else
  Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif
}

//[[Rcpp::export]]
int read_bgzf(SEXP fpsexp){
#if HAVE_HTS
  R_BGZF xfp(fpsexp);
  BGZF* fp = xfp;
  auto ir = bgzf_read_block(fp);
  if(ir!=0)
    return ir;
  return(fp->block_length !=0 ? ir : -1);
  #else
  Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif
}

//[[Rcpp::export]]
int close_bgzf(SEXP fpsexp){
#if HAVE_HTS
  return 1;
#else
    Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif

}

//[[Rcpp::export]]
int num_bgzf_blocks(SEXP fpsexp){

  #if HAVE_HTS
  R_BGZF xfp(fpsexp);
  BGZF* fp = xfp;

  // auto iw = fp->is_write;
  // auto ic = fp->is_compressed;
  // auto bl = fp->block_length;
  // auto bcl = fp->block_clength;
  // auto
  int num_blocks=0;
  auto fpi = fp->idx;
  if(fpi!=nullptr){
    std::memcpy(&num_blocks,fpi,sizeof(num_blocks));
    return num_blocks;
  }else{
    auto tell_v = bgzf_tell(fp);
    int i=0;
    do{
      auto ir = bgzf_read_block(fp);
      i++;
    }while(bgzf_check_EOF(fp)==0);
    auto ret = bgzf_seek(fp,tell_v,SEEK_SET);
    if(ret!=0){
      Rcpp::warning("cannot seek back to where we started");
    }
    return i;
  }
  #else
  Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif
}

//[[Rcpp::export]]
Rcpp::StringVector format_bgzf(SEXP fpsexp){
  #if HAVE_HTS
  R_BGZF xfp(fpsexp);
  const BGZF* fp = xfp;
  // auto iw = fp->is_write;
  // auto ic = fp->is_compressed;
  // auto bl = fp->block_length;
  // auto bcl = fp->block_clength;
  // auto
  int num_blocks=0;
  auto fpi = fp->idx;
  if(fpi!=nullptr){
    std::memcpy(&num_blocks,fpi,sizeof(num_blocks));
  }
  auto ret = fmt::format("is_write: {}\nnum_blocks: {}\nblock_length: {}\nblock_clength: {}\nblock_offset: {}\nblock_address: {}\nuncompressed_address: {}\nuncompressed_block: {:p}",
                         fp->is_write, num_blocks, fp->block_length,
                         fp->block_clength, fp->block_offset,fp->block_address,fp->uncompressed_address,fp->uncompressed_block);
  return Rcpp::wrap(ret);
    #else
  Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif
}

//[[Rcpp::export]]
Rcpp::StringVector get_bgzf_data(SEXP fpsexp){
  #if HAVE_HTS
  R_BGZF xfp(fpsexp);
  BGZF* fp = xfp;
  if(fp->uncompressed_block != nullptr){
    const char* rp = static_cast<const char*>(fp->uncompressed_block);
    std::string rv(rp,fp->block_length);
    return Rcpp::wrap(rv);
  }
  return Rcpp::StringVector::create();
    #else
  Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif
}






//[[Rcpp::export]]
Rcpp::StringVector readlines_chunk_bgzf(SEXP fpsexp){
  #if HAVE_HTS
  R_BGZF xfp(fpsexp);
  BGZF* fp = xfp;
  using namespace ranges;
  if(fp->uncompressed_block != nullptr){
    const char* rp = static_cast<char*>(fp->uncompressed_block);
    std::string rv(rp,fp->block_length);
    auto cs = views::c_str(rv.c_str());
    std::string_view rpv(rv);
    auto dv = views::split(rpv,'\n') | views::common;
    std::string tstring;
    return Rcpp::StringVector::import_transform(dv.begin(),dv.end(),[&](auto el) -> Rcpp::String {
                                                                      tstring.clear();
                                                                      std::copy_n(begin(el),distance(el),std::back_inserter(tstring));
                                                                      return Rcpp::String(tstring);
                                                                    });

  }
  return Rcpp::StringVector::create();
    #else
  Rcpp::stop("htslib functions not available, please install htslib and then reinstall ldmap");
  return NILSXP;
#endif
}

// Rcpp::DataFrame readlines_df_bgzf(SEXP fpsexp,Rcpp::IntegerVector rsid){
//   R_BGZF xfp(fpsexp);
//   BGZF* fp = xfp;
//   using namespace ranges;
//   // Rcpp::DataFrame::
//   if(fp->uncompressed_block != nullptr){
//     const char* rp = static_cast<char*>(fp->uncompressed_block);
//     std::string rv(rp,fp->block_length);
//     auto cs = views::c_str(rv.c_str());
//     std::string_view rpv(rv);
//     auto dv = views::split(rpv,'\n') | views::common;
//     std::string tstring;
//     return Rcpp::StringVector::import_transform(dv.begin(),dv.end(),[&](auto el) -> Rcpp::String {
//                                                                       tstring.clear();
//                                                                       std::copy_n(begin(el),distance(el),std::back_inserter(tstring));
//                                                                       return Rcpp::String(tstring);
//                                                                     });

//   }
//   return Rcpp::StringVector::create();
// }
