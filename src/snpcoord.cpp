#include <ldmap/ldmap.hpp>
#include <ldmap/genetic_map.hpp>
#include <progress.hpp>
//[[Rcpp::export]]
bool sorted_snp_df(const Rcpp::IntegerVector chr,  const Rcpp::IntegerVector pos){


  const size_t p=chr.size();
    //std::vector<std::pair<int,int>> r(p);
  std::pair<int,int> old_p{chr[0],pos[0]};
  std::pair<int,int> new_p;
  for(int i=1; i<p;i++){
    new_p={chr[i],pos[i]};
    if(new_p<old_p){
      return(false);
    }
    old_p=new_p;
  }
  return(true);
}


struct SNPpos{
  std::pair<int,int> c_p;
  SNPpos(const int chrom,const int pos):c_p(std::make_pair(chrom,pos)){};
};

struct GRange{
  std::pair<SNPpos,SNPpos> range;
  GRange(const int chrom,const int start,const int end):range(std::make_pair(SNPpos(chrom,start),SNPpos(chrom,end))){
    if(end<start){
      Rcpp::Rcerr<<"For range: "<<chrom<<":"<<start<<"-"<<end<<std::endl;
      Rcpp::stop("Invalid range! end<start");
    }};
};
bool operator<(const SNPpos p, const GRange g){
  return(p.c_p<g.range.first.c_p);
}
bool operator<=(const SNPpos p, const GRange g){
  return(p.c_p<=g.range.first.c_p);
}
bool operator>(const SNPpos p, const GRange g){
  return(p.c_p>=g.range.second.c_p);
}
// bool operator>=(const SNPpos p, const GRange g){
//   return(p.c_p>=g.range.second.c_p);
// }
bool operator==(const SNPpos p, const GRange g){
  return((p.c_p>=g.range.first.c_p) && (p.c_p<g.range.second.c_p));
}

// 
// struct pairhash {
// public:
//   std::size_t operator()(const std::pair<uint32_t, uint32_t> &x) const
//   {
//     return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
//   }
// };

//' Assign snps to regions of the genome, breaking up regions based on number of SNPs
// '@param ld_chr vector of chromosomes	of per region
// '@param ld_start vector of start positions for each region
// '@param ld_stop vector of end positions for each region
// '@return returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
//' @export
//[[Rcpp::export]]
Rcpp::StringVector set_ld_region(const Rcpp::IntegerVector ld_chr,
				    const Rcpp::IntegerVector ld_start,
				    const Rcpp::IntegerVector ld_stop,
				    const Rcpp::IntegerVector ld_region_id,
				    const Rcpp::IntegerVector chr,
				    const Rcpp::IntegerVector pos,
				    uint32_t max_size=0,
				    int min_size=1,
				    const bool assign_all=true){

  const size_t ld_size=ld_chr.size();

  if(!std::is_sorted(ld_chr.begin(),ld_chr.end())){
    Rcpp::stop("break regions must be sorted by chromosome!");
  }

  const size_t snp_size = chr.size();
  Rcpp::StringVector ret_region(snp_size);
  using pair_i = uint64_t;
  std::vector< pair_i > tret_region(snp_size);
  using ref_p =std::pair<std::vector<pair_i>::iterator,std::vector<pair_i>::iterator >;
  std::unordered_map<uint32_t,
                     int> ref_map(ld_size*2);
  if(!std::is_sorted(chr.begin(),chr.end())){
    Rcpp::stop("SNPs must be sorted by chromosome!");
  }
  size_t idx1=0;
  size_t idx2=0;
  GRange gr(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
  max_size = max_size == 0 ? snp_size + 1 : max_size;
  uint32_t grc=0;
  
  Progress p(snp_size, true);
  
  while(idx1<snp_size){
    if (Progress::check_abort() ){
      Rcpp::stop("aborted");
    }
    const SNPpos s_pos(chr[idx1],pos[idx1]);
    if(s_pos==gr){
      //point in interval
      grc++;
      pair_i rel= (ld_region_id[idx2] << 16) | (grc/max_size);
      tret_region[idx1]=rel;
      ref_map[rel]++;
      idx1++;
      p.increment();
    }else{
      if(s_pos<gr){
        // point before interval
        if(assign_all){
          Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
          Rcpp::Rcerr<<"Can't map to "<<ld_chr[idx2]<<":"<<ld_start[idx2]<<"-"<<ld_stop[idx2]<<" (ahead)"<<std::endl;
          if(idx2>0){
            Rcpp::Rcerr<<"Can't map to "<<ld_chr[idx2-1]<<":"<<ld_start[idx2-1]<<"-"<<ld_stop[idx2-1]<<" (behind)"<<std::endl;
          }
          Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
        }else{
          tret_region[idx1]=(NA_INTEGER << 16) | NA_INTEGER;
        }
        idx1++;
        p.increment();
      }
      if(s_pos>gr){
        //point after interval
        if((idx2+1)<ld_size){
          grc=0;
          idx2++;
          gr=GRange(ld_chr[idx2],ld_start[idx2],ld_stop[idx2]);
        }else{
          if(assign_all){
            Rcpp::Rcerr<<"For SNP: "<<chr[idx1]<<":"<<pos[idx1]<<std::endl;
            Rcpp::stop("Can't map SNP to region! Set assign_all to false to ignore");
          }else{
            tret_region[idx1]=(NA_INTEGER << 16) | NA_INTEGER;
          }
          idx1++;
          p.increment();
        }
      }
    }
  }

  std::transform(tret_region.begin(),tret_region.end(),ret_region.begin(),[&ref_map,min_size](pair_i eli){
    int high16 = eli >> 16;
    int low16 = eli & 0xFFFF;
    low16 = (ref_map[eli] < min_size) ? std::max(0,low16-1) : low16;
      return((std::to_string(high16)+"."+std::to_string(low16))); 
  });
  return(ret_region);
}

