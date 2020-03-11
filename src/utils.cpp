#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <libpopcnt.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <bitset>

inline char flip(const char t) noexcept{
  switch(t){
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return '0';
  }
};



int allele_check(const char* query,const char* alt){
  const char q_ref=query[0];
  const char q_alt=query[2];
  const char t_ref=alt[0];
  const char t_alt=alt[2];
  // Alt and ref being equal is a no-go
  if((q_ref==q_alt) or (t_ref==t_alt)){
    return 0;
  }

  //2 scenarios if first letters match
  if(q_ref==t_ref){
    //matches in both places
    if(q_alt==t_alt){
      return 1;
    }
    //matches in only one place
    return 0;
  }
  //this also falls under only matches in one place
  if(q_alt==t_alt){
    return 0;
  }
  //"classic" flip
  if(q_ref==t_alt){
    if(q_alt==t_ref){
      return -1;
    }
    return 0;
  }
  //strand flipped ref matches ref
  if((flip(q_ref)==t_ref)){
    //strand flipped alt mathces too, this means flip
    if(flip(q_alt) == t_alt){
      return -1;
    }

  }
  //Last chance for a match
  if(flip(q_ref)==t_alt){
    if(flip(t_ref)==q_alt){
      return 1;
    }
  }
  //Neither match
  return 0;

}



//' Determine whether 2 alleles are compatible
//' @param query_ref_alt is a ref/alt pair
//' @param target_ref_alt is another ref/alt pair
//' @return raw matrix
//' @export
//[[Rcpp::export]]
Rcpp::RawMatrix snp2raw(Rcpp::IntegerMatrix input_matrix){
  const int p = input_matrix.cols();
  const size_t N = input_matrix.rows();
  const size_t newsize = static_cast<size_t>(std::ceil(static_cast<double>(N)/ 8.0f));
  const size_t rem = N % 8;
  const size_t newsize_b = N-rem;
  // Rcpp::Rcout<<newsize_b<<std::endl;
  // Rcpp::Rcout<<rem<<std::endl;

  Rcpp::RawMatrix retmat(newsize,p);

  // RcppParallel::RMatrix<char> prm(retmat);
  // RcppParallel::RMatrix<int>  inpm(input_matrix);
  if(rem!=0){
    for(int i=0; i<p; i++){
      for(int j=0; j<newsize_b; j+=8){
	retmat(j/8,i)=
	  input_matrix(j+0,i)|
	  (input_matrix(j+1,i)<<1)|
	  (input_matrix(j+2,i)<<2)|
	  (input_matrix(j+3,i)<<3)|
	  (input_matrix(j+4,i)<<4)|
	  (input_matrix(j+5,i)<<5)|
	  (input_matrix(j+6,i)<<6)|
	  (input_matrix(j+7,i)<<7);
      }

      Rbyte btmp = 0;
      for(int k=0; k<rem; k++){
        btmp|=	(input_matrix(newsize_b+k,i)<<k);
        // Rcpp::Rcout<<i<<"\t"<<k<<"\t"<<static_cast<int>(btmp)<<std::endl;
      }
      retmat(newsize-1,i)=btmp;

    }
  }else{
        for(size_t i=0; i<p; i++){

	  for(int j=0; j<newsize_b; j+=8){
	    Rbyte btmp = 0;
	    btmp=
	      input_matrix(j+0,i)|
	      (input_matrix(j+1,i)<<1)|
	      (input_matrix(j+2,i)<<2)|
	      (input_matrix(j+3,i)<<3)|
	      (input_matrix(j+4,i)<<4)|
	      (input_matrix(j+5,i)<<5)|
	      (input_matrix(j+6,i)<<6)|
	      (input_matrix(j+7,i)<<7);
	    // Rcpp::Rcout<<i<<"\t"<<j<<"\t"<<static_cast<int>(btmp)<<std::endl;
	    retmat(j/8,i)=btmp;
	  }

	}
  }

  return(retmat);
}


//' @export
//[[Rcpp::export]]
Rcpp::NumericVector popcnt_v(Rcpp::RawMatrix X, double sample_size=0){
  
  if(sample_size==0){
    sample_size=X.rows()*8.0;
  }
  const size_t p = X.cols();
  std::vector<uint64_t> popc(p);
  using bel = Rcpp::RawMatrix::elem_type;
  const size_t n_el=X.rows();
  for(int i=0; i<p; i++){
    bel* Xi=&X(0,i);
    if(popc[i]==0){
      popc[i]=popcnt(Xi,n_el);
    }
  }
  
  return(Rcpp::wrap(popc));
}


//' @export
//[[Rcpp::export]]
Rcpp::NumericMatrix covbin(Rcpp::RawMatrix X, double sample_size=0){
  if(sample_size==0){
    sample_size=X.rows()*8.0;
  }
  const size_t p = X.cols();
  Rcpp::NumericMatrix S(p,p);
  std::vector<uint64_t> popc(p);
  const size_t n_el=X.rows();
  using bel = Rcpp::RawMatrix::elem_type;
  // Rcpp::Rcout<<sample_size;
  for(int i=0; i<p; i++){
    bel* Xi=&X(0,i);
    if(popc[i]==0){
      popc[i]=popcnt(Xi,n_el);
    }
    for(int j=i; j<p; j++){
      bel* Xj=&X(0,j);
      if(popc[j]==0){
        popc[j]=popcnt(Xj,n_el);
      }

      auto tret = std::inner_product(Xi,Xi+n_el,Xj,0,std::plus<bel>(),[](const bel a,const bel b) noexcept{
        return __builtin_popcount(a&b);
      });
      // Rcpp::Rcout<<i<<"\t"<<j<<"\t"<<popc[i]<<"\t"<<popc[j]<<"\t"<<tret<<std::endl;
      auto mret=(sample_size*tret-popc[i]*popc[j])/(sample_size*(sample_size-1));
      S(i,j)=mret;
      S(j,i)=mret;
    }
  }
  return(S);
}






//' Determine whether 2 alleles are compatible
//' @param query_ref_alt is a ref/alt pair
//' @param target_ref_alt is another ref/alt pair
//' @return returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
//' @export
//[[Rcpp::export]]
Rcpp::StringVector strand_flip(Rcpp::StringVector ref_alt,bool reverse=false){
  Rcpp::StringVector retvec(ref_alt.size());


  if(!reverse){
    std::transform(ref_alt.begin(),ref_alt.end(),retvec.begin(),[](const char* query)  {
      if(strlen( query)!=3 ){
        Rcpp::stop("alleles must all be length 3");
      }
      const char ret[]={flip(query[0]),',',flip(query[2])};
      return(Rcpp::String(ret));
    });
  }else{
    std::transform(ref_alt.begin(),ref_alt.end(),retvec.begin(),[](const char* query)  {
      if(strlen( query)!=3 ){
        Rcpp::stop("alleles must all be length 3");
      }
      const char ret[]={flip(query[2]),',',flip(query[0])};
      return(Rcpp::String(ret));
    });
  }


  return(retvec);
}


class Reference_Snps{
  const Rcpp::IntegerVector chrom;
  const Rcpp::IntegerVector pos;
  const Rcpp::IntegerVector chunk;
  const size_t p;
  Rcpp::IntegerVector::const_iterator chrom_b;
  Rcpp::IntegerVector::const_iterator chrom_e;
  Rcpp::IntegerVector::const_iterator pos_b;
  Rcpp::IntegerVector::const_iterator pos_e;

  Rcpp::IntegerVector::const_iterator chunk_b;
  Rcpp::IntegerVector::const_iterator chunk_e;

  bool use_chunk;
public:
  Reference_Snps(const Rcpp::IntegerVector &ref_chrom,
                 const Rcpp::IntegerVector &ref_pos,
                 const Rcpp::IntegerVector &query_chunk):chrom(ref_chrom),
                 pos(ref_pos),
                 chunk(query_chunk),
                 p(chrom.size()),
                 use_chunk(query_chunk.size()==p){
    chrom_b = chrom.begin();
    chrom_e = std::upper_bound(chrom_b,chrom.end(),*chrom_b);
    pos_b = pos.begin();
    pos_e = pos.begin()+(chrom_e-chrom_b);
    chunk_b = chunk.begin();
    chunk_e = std::upper_bound(chunk_b,chunk_b+(chrom_e-chrom_b),*chunk_b);
  }
  int index(int r_chunk, int r_chrom, int r_pos) noexcept{

    if((r_chunk != *chunk_b)||(r_chrom != *chrom_b)){
      auto t_ret_p =std::equal_range(chrom_b,chrom.end(),r_chrom);
      chrom_b = t_ret_p.first;
      chrom_e = t_ret_p.second;
      t_ret_p =std::equal_range(chunk.begin()+(chrom_b-chrom.begin()),chunk.begin()+(chrom_e-chrom.begin()),r_chunk);
      chunk_b = t_ret_p.first;
      chunk_e = t_ret_p.second;
      pos_b = pos.begin()+(chunk_b-chunk.begin());
      pos_e = pos.begin()+(chunk_e-chunk.begin());
    }
    auto trp = std::lower_bound(pos_b,pos_e,r_pos);
    if(trp!=pos_e){
      pos_b = trp;
      if(*trp==r_pos){
        return(trp-pos.begin()+1);
      }
    }
    return(NA_INTEGER);
  }
  int index(int r_chrom, int r_pos) noexcept{

    if((r_chrom!=*chrom_b)){
      auto t_ret_p =std::equal_range(chrom_e,chrom.end(),r_chrom);
      chrom_b = t_ret_p.first;
      chrom_e = t_ret_p.second;
      pos_b = pos.begin()+(chrom_b-chunk.begin());
      pos_e = pos.begin()+(chrom_e-chunk.begin());
    }
    auto trp = std::lower_bound(pos_b,pos_e,r_pos);
    if(trp!=pos_e){
      pos_b = trp;
      if(*trp==r_pos){
        return(trp-pos.begin()+1);
      }
    }
    return(NA_INTEGER);
  }
};


//' Find query SNP in a list of reference snps
//' @param query_chrom query chromosome
//' @param query_pos query position
//' @param ref_chrom reference chromosome
//' @param ref_pos reference position
//' @param query_chunk region assignment for query (optional)
//' @param ref_chunk region assignment for reference (optional)
//' @return returns a vector with the position of the match, or NA if no match is found.
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector find_alleles(Rcpp::IntegerVector query_chrom,
                                 Rcpp::IntegerVector query_pos,
                                 Rcpp::IntegerVector ref_chrom,
                                 Rcpp::IntegerVector ref_pos,
                                 Rcpp::IntegerVector query_chunk=Rcpp::IntegerVector::create(),
                                 Rcpp::IntegerVector ref_chunk=Rcpp::IntegerVector::create()
                                  ){
  const size_t q_p = query_chrom.size();

  Progress p(0, false); // we need an instance, should be improved in next version

  if(query_chrom.size()!=query_pos.size()){
    Rcpp::stop("query chrom and pos must be the same size");
  }
  if(ref_chrom.size()!=ref_pos.size()){
    Rcpp::stop("ref chrom and pos must be the same size");
  }
  const bool qchunk =(query_chunk.size()>0);
  const bool rchunk =(ref_chunk.size()>0);
  const  bool use_chunk = qchunk && rchunk;
  if(use_chunk){
    if( query_chunk.size()!=query_chrom.size()){
      Rcpp::stop("query chrom and chunk must be the same size if query chunk is specified");
    }
    if(ref_chunk.size()!=ref_chrom.size()){
      Rcpp::stop("ref chrom and pos must be the same size");
    }
  }


  Reference_Snps refsnps(ref_chrom,ref_pos,ref_chunk);



  // auto q_chunk_e = query_chrom.begin();

  Rcpp::IntegerVector ret_p(q_p);
  if(use_chunk){
    for(size_t i=0; i<q_p; i++){
      if (Progress::check_abort() )
        return ret_p;
      ret_p[i]=refsnps.index(query_chunk[i],query_chrom[i],query_pos[i]);
    }

  }else{
    for(size_t i=0; i<q_p; i++){
      if (Progress::check_abort() )
        return ret_p;
      ret_p[i]=refsnps.index(query_chrom[i],query_pos[i]);
    }
  }

  return ret_p;
}



//' Determine whether 2
//' @param query_ref_alt is a ref/alt pair
//' @param target_ref_alt is another ref/alt pair
//' @return returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
//' Note that in the case that the reference and/or alternate allele are not single characters, they will only be checked for equality
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector flip_alleles(Rcpp::StringVector query_ref_alt,Rcpp::StringVector target_ref_alt ){

  if(query_ref_alt.size()!=target_ref_alt.size()){
    Rcpp::stop("query ref_alt and target ref_alt must be the same size!");
  }
  Rcpp::IntegerVector retvec(target_ref_alt.size());
  std::transform(query_ref_alt.begin(),query_ref_alt.end(),target_ref_alt.begin(),retvec.begin(),[](const char* query,const char* alt) {
												   if(strlen( query)!=3 || strlen(alt)!=3 ){
												     if(strcmp(query,alt)!=0){
												       return(0);
												     }else{
												       return(1);
												     }
												     Rcpp::stop("alleles must all be length 3");
												   }
												   return(allele_check(query,alt));
												 });

  return(retvec);
}
