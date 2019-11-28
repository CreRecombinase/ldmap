// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

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
// parse_ldmap_range
Rcpp::NumericVector parse_ldmap_range(Rcpp::StringVector input);
RcppExport SEXP _ldmap_parse_ldmap_range(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_ldmap_range(input));
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
Rcpp::StringVector set_ld_region(const Rcpp::IntegerVector ld_chr, const Rcpp::IntegerVector ld_start, const Rcpp::IntegerVector ld_stop, const Rcpp::IntegerVector ld_region_id, const Rcpp::IntegerVector chr, const Rcpp::IntegerVector pos, uint32_t max_size, int min_size, const bool assign_all);
RcppExport SEXP _ldmap_set_ld_region(SEXP ld_chrSEXP, SEXP ld_startSEXP, SEXP ld_stopSEXP, SEXP ld_region_idSEXP, SEXP chrSEXP, SEXP posSEXP, SEXP max_sizeSEXP, SEXP min_sizeSEXP, SEXP assign_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_chr(ld_chrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_start(ld_startSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_stop(ld_stopSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type ld_region_id(ld_region_idSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< uint32_t >::type max_size(max_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type min_size(min_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type assign_all(assign_allSEXP);
    rcpp_result_gen = Rcpp::wrap(set_ld_region(ld_chr, ld_start, ld_stop, ld_region_id, chr, pos, max_size, min_size, assign_all));
    return rcpp_result_gen;
END_RCPP
}
// new_ldmap_range
Rcpp::NumericVector new_ldmap_range(Rcpp::IntegerVector chrom, Rcpp::IntegerVector start, Rcpp::IntegerVector end);
RcppExport SEXP _ldmap_new_ldmap_range(SEXP chromSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type chrom(chromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(new_ldmap_range(chrom, start, end));
    return rcpp_result_gen;
END_RCPP
}
// snp_in_range
Rcpp::IntegerVector snp_in_range(Rcpp::NumericVector ldmap_snp, Rcpp::NumericVector ldmap_range);
RcppExport SEXP _ldmap_snp_in_range(SEXP ldmap_snpSEXP, SEXP ldmap_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_range(ldmap_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(snp_in_range(ldmap_snp, ldmap_range));
    return rcpp_result_gen;
END_RCPP
}
// snp_in_ranges
Rcpp::List snp_in_ranges(Rcpp::NumericVector ldmap_snp, Rcpp::ListOf<Rcpp::NumericVector> ldmap_ranges);
RcppExport SEXP _ldmap_snp_in_ranges(SEXP ldmap_snpSEXP, SEXP ldmap_rangesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<Rcpp::NumericVector> >::type ldmap_ranges(ldmap_rangesSEXP);
    rcpp_result_gen = Rcpp::wrap(snp_in_ranges(ldmap_snp, ldmap_ranges));
    return rcpp_result_gen;
END_RCPP
}
// format_ldmap_range
Rcpp::StringVector format_ldmap_range(Rcpp::NumericVector x);
RcppExport SEXP _ldmap_format_ldmap_range(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(format_ldmap_range(x));
    return rcpp_result_gen;
END_RCPP
}
// match_ranges_snps
Rcpp::List match_ranges_snps(Rcpp::DataFrame df, Rcpp::NumericVector ldmap_range, const std::string snp_col);
RcppExport SEXP _ldmap_match_ranges_snps(SEXP dfSEXP, SEXP ldmap_rangeSEXP, SEXP snp_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_range(ldmap_rangeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type snp_col(snp_colSEXP);
    rcpp_result_gen = Rcpp::wrap(match_ranges_snps(df, ldmap_range, snp_col));
    return rcpp_result_gen;
END_RCPP
}
// window_ldmap_range
Rcpp::NumericVector window_ldmap_range(Rcpp::NumericVector ldmap_snp, Rcpp::NumericVector cm, const double window);
RcppExport SEXP _ldmap_window_ldmap_range(SEXP ldmap_snpSEXP, SEXP cmSEXP, SEXP windowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cm(cmSEXP);
    Rcpp::traits::input_parameter< const double >::type window(windowSEXP);
    rcpp_result_gen = Rcpp::wrap(window_ldmap_range(ldmap_snp, cm, window));
    return rcpp_result_gen;
END_RCPP
}
// split_ldmap_range_overlap
Rcpp::NumericVector split_ldmap_range_overlap(Rcpp::NumericVector x);
RcppExport SEXP _ldmap_split_ldmap_range_overlap(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(split_ldmap_range_overlap(x));
    return rcpp_result_gen;
END_RCPP
}
// merge_ldmap_ranges
Rcpp::NumericVector merge_ldmap_ranges(Rcpp::NumericVector x, Rcpp::NumericVector y);
RcppExport SEXP _ldmap_merge_ldmap_ranges(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(merge_ldmap_ranges(x, y));
    return rcpp_result_gen;
END_RCPP
}
// ldmap_range_2_data_frame
SEXP ldmap_range_2_data_frame(Rcpp::NumericVector ldmap_range);
RcppExport SEXP _ldmap_ldmap_range_2_data_frame(SEXP ldmap_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_range(ldmap_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(ldmap_range_2_data_frame(ldmap_range));
    return rcpp_result_gen;
END_RCPP
}
// sample_interval
Rcpp::IntegerVector sample_interval(Rcpp::IntegerVector n, Rcpp::IntegerVector begin, Rcpp::IntegerVector end, const bool replace);
RcppExport SEXP _ldmap_sample_interval(SEXP nSEXP, SEXP beginSEXP, SEXP endSEXP, SEXP replaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type begin(beginSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type end(endSEXP);
    Rcpp::traits::input_parameter< const bool >::type replace(replaceSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_interval(n, begin, end, replace));
    return rcpp_result_gen;
END_RCPP
}
// new_ldmap_snp
Rcpp::NumericVector new_ldmap_snp(Rcpp::IntegerVector chrom, Rcpp::NumericVector pos, Rcpp::RObject ref, Rcpp::RObject alt, const bool NA2N);
RcppExport SEXP _ldmap_new_ldmap_snp(SEXP chromSEXP, SEXP posSEXP, SEXP refSEXP, SEXP altSEXP, SEXP NA2NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type chrom(chromSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type ref(refSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type alt(altSEXP);
    Rcpp::traits::input_parameter< const bool >::type NA2N(NA2NSEXP);
    rcpp_result_gen = Rcpp::wrap(new_ldmap_snp(chrom, pos, ref, alt, NA2N));
    return rcpp_result_gen;
END_RCPP
}
// is_strand_ambiguous
Rcpp::LogicalVector is_strand_ambiguous(Rcpp::NumericVector x);
RcppExport SEXP _ldmap_is_strand_ambiguous(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(is_strand_ambiguous(x));
    return rcpp_result_gen;
END_RCPP
}
// new_ldmap_allele
Rcpp::RawVector new_ldmap_allele(Rcpp::RObject allele);
RcppExport SEXP _ldmap_new_ldmap_allele(SEXP alleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type allele(alleleSEXP);
    rcpp_result_gen = Rcpp::wrap(new_ldmap_allele(allele));
    return rcpp_result_gen;
END_RCPP
}
// order_snps
Rcpp::IntegerVector order_snps(Rcpp::NumericVector struct_vec);
RcppExport SEXP _ldmap_order_snps(SEXP struct_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type struct_vec(struct_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(order_snps(struct_vec));
    return rcpp_result_gen;
END_RCPP
}
// rank_snps
Rcpp::IntegerVector rank_snps(Rcpp::NumericVector struct_vec);
RcppExport SEXP _ldmap_rank_snps(SEXP struct_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type struct_vec(struct_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(rank_snps(struct_vec));
    return rcpp_result_gen;
END_RCPP
}
// chromosomes
Rcpp::IntegerVector chromosomes(Rcpp::NumericVector struct_vec);
RcppExport SEXP _ldmap_chromosomes(SEXP struct_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type struct_vec(struct_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(chromosomes(struct_vec));
    return rcpp_result_gen;
END_RCPP
}
// starts
Rcpp::IntegerVector starts(Rcpp::NumericVector ldmap_range);
RcppExport SEXP _ldmap_starts(SEXP ldmap_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_range(ldmap_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(starts(ldmap_range));
    return rcpp_result_gen;
END_RCPP
}
// ends
Rcpp::IntegerVector ends(Rcpp::NumericVector ldmap_range);
RcppExport SEXP _ldmap_ends(SEXP ldmap_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_range(ldmap_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(ends(ldmap_range));
    return rcpp_result_gen;
END_RCPP
}
// positions
Rcpp::NumericVector positions(Rcpp::NumericVector ldmap_snp);
RcppExport SEXP _ldmap_positions(SEXP ldmap_snpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    rcpp_result_gen = Rcpp::wrap(positions(ldmap_snp));
    return rcpp_result_gen;
END_RCPP
}
// ref_alleles
SEXP ref_alleles(Rcpp::NumericVector ldmap_snp, const bool as_ldmap_allele);
RcppExport SEXP _ldmap_ref_alleles(SEXP ldmap_snpSEXP, SEXP as_ldmap_alleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    Rcpp::traits::input_parameter< const bool >::type as_ldmap_allele(as_ldmap_alleleSEXP);
    rcpp_result_gen = Rcpp::wrap(ref_alleles(ldmap_snp, as_ldmap_allele));
    return rcpp_result_gen;
END_RCPP
}
// format_ldmap_allele
Rcpp::StringVector format_ldmap_allele(Rcpp::RawVector x);
RcppExport SEXP _ldmap_format_ldmap_allele(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(format_ldmap_allele(x));
    return rcpp_result_gen;
END_RCPP
}
// as_integer_ldmap_allele
Rcpp::IntegerVector as_integer_ldmap_allele(Rcpp::RawVector x);
RcppExport SEXP _ldmap_as_integer_ldmap_allele(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(as_integer_ldmap_allele(x));
    return rcpp_result_gen;
END_RCPP
}
// as_integer_ldmap_range
Rcpp::NumericVector as_integer_ldmap_range(Rcpp::NumericVector x);
RcppExport SEXP _ldmap_as_integer_ldmap_range(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(as_integer_ldmap_range(x));
    return rcpp_result_gen;
END_RCPP
}
// as_integer_ldmap_snp
Rcpp::NumericVector as_integer_ldmap_snp(Rcpp::NumericVector x);
RcppExport SEXP _ldmap_as_integer_ldmap_snp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(as_integer_ldmap_snp(x));
    return rcpp_result_gen;
END_RCPP
}
// alt_alleles
SEXP alt_alleles(Rcpp::NumericVector ldmap_snp, const bool as_ldmap_allele);
RcppExport SEXP _ldmap_alt_alleles(SEXP ldmap_snpSEXP, SEXP as_ldmap_alleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    Rcpp::traits::input_parameter< const bool >::type as_ldmap_allele(as_ldmap_alleleSEXP);
    rcpp_result_gen = Rcpp::wrap(alt_alleles(ldmap_snp, as_ldmap_allele));
    return rcpp_result_gen;
END_RCPP
}
// ldmap_snp_2_dataframe
SEXP ldmap_snp_2_dataframe(Rcpp::NumericVector ldmap_snp, bool alleles_to_character);
RcppExport SEXP _ldmap_ldmap_snp_2_dataframe(SEXP ldmap_snpSEXP, SEXP alleles_to_characterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ldmap_snp(ldmap_snpSEXP);
    Rcpp::traits::input_parameter< bool >::type alleles_to_character(alleles_to_characterSEXP);
    rcpp_result_gen = Rcpp::wrap(ldmap_snp_2_dataframe(ldmap_snp, alleles_to_character));
    return rcpp_result_gen;
END_RCPP
}
// join_snp
Rcpp::List join_snp(Rcpp::NumericVector query, Rcpp::NumericVector reference, Rcpp::IntegerVector rsid);
RcppExport SEXP _ldmap_join_snp(SEXP querySEXP, SEXP referenceSEXP, SEXP rsidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type query(querySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type reference(referenceSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type rsid(rsidSEXP);
    rcpp_result_gen = Rcpp::wrap(join_snp(query, reference, rsid));
    return rcpp_result_gen;
END_RCPP
}
// extract_alt
Rcpp::StringVector extract_alt(Rcpp::StringVector ref, Rcpp::StringVector alleles_as_ambig);
RcppExport SEXP _ldmap_extract_alt(SEXP refSEXP, SEXP alleles_as_ambigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type alleles_as_ambig(alleles_as_ambigSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_alt(ref, alleles_as_ambig));
    return rcpp_result_gen;
END_RCPP
}
// fast_str2int
Rcpp::IntegerVector fast_str2int(Rcpp::StringVector input, int offset, const std::string prefix, const int na_val);
RcppExport SEXP _ldmap_fast_str2int(SEXP inputSEXP, SEXP offsetSEXP, SEXP prefixSEXP, SEXP na_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type input(inputSEXP);
    Rcpp::traits::input_parameter< int >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const std::string >::type prefix(prefixSEXP);
    Rcpp::traits::input_parameter< const int >::type na_val(na_valSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_str2int(input, offset, prefix, na_val));
    return rcpp_result_gen;
END_RCPP
}
// make_ambig
Rcpp::StringVector make_ambig(Rcpp::StringVector A1, Rcpp::StringVector A2);
RcppExport SEXP _ldmap_make_ambig(SEXP A1SEXP, SEXP A2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type A2(A2SEXP);
    rcpp_result_gen = Rcpp::wrap(make_ambig(A1, A2));
    return rcpp_result_gen;
END_RCPP
}
// format_ldmap_snp
Rcpp::StringVector format_ldmap_snp(Rcpp::NumericVector x);
RcppExport SEXP _ldmap_format_ldmap_snp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(format_ldmap_snp(x));
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
    {"_ldmap_parse_ldmap_range", (DL_FUNC) &_ldmap_parse_ldmap_range, 1},
    {"_ldmap_sorted_snp_df", (DL_FUNC) &_ldmap_sorted_snp_df, 2},
    {"_ldmap_set_ld_region", (DL_FUNC) &_ldmap_set_ld_region, 9},
    {"_ldmap_new_ldmap_range", (DL_FUNC) &_ldmap_new_ldmap_range, 3},
    {"_ldmap_snp_in_range", (DL_FUNC) &_ldmap_snp_in_range, 2},
    {"_ldmap_snp_in_ranges", (DL_FUNC) &_ldmap_snp_in_ranges, 2},
    {"_ldmap_format_ldmap_range", (DL_FUNC) &_ldmap_format_ldmap_range, 1},
    {"_ldmap_match_ranges_snps", (DL_FUNC) &_ldmap_match_ranges_snps, 3},
    {"_ldmap_window_ldmap_range", (DL_FUNC) &_ldmap_window_ldmap_range, 3},
    {"_ldmap_split_ldmap_range_overlap", (DL_FUNC) &_ldmap_split_ldmap_range_overlap, 1},
    {"_ldmap_merge_ldmap_ranges", (DL_FUNC) &_ldmap_merge_ldmap_ranges, 2},
    {"_ldmap_ldmap_range_2_data_frame", (DL_FUNC) &_ldmap_ldmap_range_2_data_frame, 1},
    {"_ldmap_sample_interval", (DL_FUNC) &_ldmap_sample_interval, 4},
    {"_ldmap_new_ldmap_snp", (DL_FUNC) &_ldmap_new_ldmap_snp, 5},
    {"_ldmap_is_strand_ambiguous", (DL_FUNC) &_ldmap_is_strand_ambiguous, 1},
    {"_ldmap_new_ldmap_allele", (DL_FUNC) &_ldmap_new_ldmap_allele, 1},
    {"_ldmap_order_snps", (DL_FUNC) &_ldmap_order_snps, 1},
    {"_ldmap_rank_snps", (DL_FUNC) &_ldmap_rank_snps, 1},
    {"_ldmap_chromosomes", (DL_FUNC) &_ldmap_chromosomes, 1},
    {"_ldmap_starts", (DL_FUNC) &_ldmap_starts, 1},
    {"_ldmap_ends", (DL_FUNC) &_ldmap_ends, 1},
    {"_ldmap_positions", (DL_FUNC) &_ldmap_positions, 1},
    {"_ldmap_ref_alleles", (DL_FUNC) &_ldmap_ref_alleles, 2},
    {"_ldmap_format_ldmap_allele", (DL_FUNC) &_ldmap_format_ldmap_allele, 1},
    {"_ldmap_as_integer_ldmap_allele", (DL_FUNC) &_ldmap_as_integer_ldmap_allele, 1},
    {"_ldmap_as_integer_ldmap_range", (DL_FUNC) &_ldmap_as_integer_ldmap_range, 1},
    {"_ldmap_as_integer_ldmap_snp", (DL_FUNC) &_ldmap_as_integer_ldmap_snp, 1},
    {"_ldmap_alt_alleles", (DL_FUNC) &_ldmap_alt_alleles, 2},
    {"_ldmap_ldmap_snp_2_dataframe", (DL_FUNC) &_ldmap_ldmap_snp_2_dataframe, 2},
    {"_ldmap_join_snp", (DL_FUNC) &_ldmap_join_snp, 3},
    {"_ldmap_extract_alt", (DL_FUNC) &_ldmap_extract_alt, 2},
    {"_ldmap_fast_str2int", (DL_FUNC) &_ldmap_fast_str2int, 4},
    {"_ldmap_make_ambig", (DL_FUNC) &_ldmap_make_ambig, 2},
    {"_ldmap_format_ldmap_snp", (DL_FUNC) &_ldmap_format_ldmap_snp, 1},
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
