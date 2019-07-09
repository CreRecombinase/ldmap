// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// interpolate_genetic_map
Rcpp::NumericVector interpolate_genetic_map(const Rcpp::NumericVector& map, const Rcpp::IntegerVector map_pos, const Rcpp::IntegerVector target_pos, const bool strict, const bool progress);
RcppExport SEXP _ldmap_interpolate_genetic_map(SEXP mapSEXP, SEXP map_posSEXP, SEXP target_posSEXP, SEXP strictSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type map_pos(map_posSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type target_pos(target_posSEXP);
    Rcpp::traits::input_parameter< const bool >::type strict(strictSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(interpolate_genetic_map(map, map_pos, target_pos, strict, progress));
    return rcpp_result_gen;
END_RCPP
}
// sorted_snp_df
bool sorted_snp_df(const Rcpp::IntegerVector chr, const Rcpp::IntegerVector pos);
RcppExport SEXP _ldmap_sorted_snp_df(SEXP chrSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(sorted_snp_df(chr, pos));
    return rcpp_result_gen;
END_RCPP
}
// set_ld_region
Rcpp::IntegerVector set_ld_region(const Rcpp::IntegerVector ld_chr, const Rcpp::IntegerVector ld_start, const Rcpp::IntegerVector ld_stop, const Rcpp::IntegerVector ld_region_id, const Rcpp::IntegerVector chr, const Rcpp::IntegerVector pos, const bool assign_all);
RcppExport SEXP _ldmap_set_ld_region(SEXP ld_chrSEXP, SEXP ld_startSEXP, SEXP ld_stopSEXP, SEXP ld_region_idSEXP, SEXP chrSEXP, SEXP posSEXP, SEXP assign_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_chr(ld_chrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_start(ld_startSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_stop(ld_stopSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_region_id(ld_region_idSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< const bool >::type assign_all(assign_allSEXP);
    rcpp_result_gen = Rcpp::wrap(set_ld_region(ld_chr, ld_start, ld_stop, ld_region_id, chr, pos, assign_all));
    return rcpp_result_gen;
END_RCPP
}
// snp2raw
Rcpp::RawMatrix snp2raw(Rcpp::IntegerMatrix input_matrix);
RcppExport SEXP _ldmap_snp2raw(SEXP input_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type input_matrix(input_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(snp2raw(input_matrix));
    return rcpp_result_gen;
END_RCPP
}
// popcnt_v
Rcpp::NumericVector popcnt_v(Rcpp::RawMatrix X, double sample_size);
RcppExport SEXP _ldmap_popcnt_v(SEXP XSEXP, SEXP sample_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type sample_size(sample_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(popcnt_v(X, sample_size));
    return rcpp_result_gen;
END_RCPP
}
// covbin
Rcpp::NumericMatrix covbin(Rcpp::RawMatrix X, double sample_size);
RcppExport SEXP _ldmap_covbin(SEXP XSEXP, SEXP sample_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type sample_size(sample_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(covbin(X, sample_size));
    return rcpp_result_gen;
END_RCPP
}
// strand_flip
Rcpp::StringVector strand_flip(Rcpp::StringVector ref_alt, bool reverse);
RcppExport SEXP _ldmap_strand_flip(SEXP ref_altSEXP, SEXP reverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type ref_alt(ref_altSEXP);
    Rcpp::traits::input_parameter< bool >::type reverse(reverseSEXP);
    rcpp_result_gen = Rcpp::wrap(strand_flip(ref_alt, reverse));
    return rcpp_result_gen;
END_RCPP
}
// find_alleles
Rcpp::IntegerVector find_alleles(Rcpp::IntegerVector query_chrom, Rcpp::IntegerVector query_pos, Rcpp::IntegerVector ref_chrom, Rcpp::IntegerVector ref_pos, Rcpp::IntegerVector query_chunk, Rcpp::IntegerVector ref_chunk);
RcppExport SEXP _ldmap_find_alleles(SEXP query_chromSEXP, SEXP query_posSEXP, SEXP ref_chromSEXP, SEXP ref_posSEXP, SEXP query_chunkSEXP, SEXP ref_chunkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type query_chrom(query_chromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type query_pos(query_posSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_chrom(ref_chromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_pos(ref_posSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type query_chunk(query_chunkSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_chunk(ref_chunkSEXP);
    rcpp_result_gen = Rcpp::wrap(find_alleles(query_chrom, query_pos, ref_chrom, ref_pos, query_chunk, ref_chunk));
    return rcpp_result_gen;
END_RCPP
}
// flip_alleles
Rcpp::IntegerVector flip_alleles(Rcpp::StringVector query_ref_alt, Rcpp::StringVector target_ref_alt);
RcppExport SEXP _ldmap_flip_alleles(SEXP query_ref_altSEXP, SEXP target_ref_altSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type query_ref_alt(query_ref_altSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type target_ref_alt(target_ref_altSEXP);
    rcpp_result_gen = Rcpp::wrap(flip_alleles(query_ref_alt, target_ref_alt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ldmap_interpolate_genetic_map", (DL_FUNC) &_ldmap_interpolate_genetic_map, 5},
    {"_ldmap_sorted_snp_df", (DL_FUNC) &_ldmap_sorted_snp_df, 2},
    {"_ldmap_set_ld_region", (DL_FUNC) &_ldmap_set_ld_region, 7},
    {"_ldmap_snp2raw", (DL_FUNC) &_ldmap_snp2raw, 1},
    {"_ldmap_popcnt_v", (DL_FUNC) &_ldmap_popcnt_v, 2},
    {"_ldmap_covbin", (DL_FUNC) &_ldmap_covbin, 2},
    {"_ldmap_strand_flip", (DL_FUNC) &_ldmap_strand_flip, 2},
    {"_ldmap_find_alleles", (DL_FUNC) &_ldmap_find_alleles, 6},
    {"_ldmap_flip_alleles", (DL_FUNC) &_ldmap_flip_alleles, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ldmap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
