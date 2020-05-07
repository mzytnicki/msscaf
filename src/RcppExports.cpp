// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// parseBamFileCpp
DataFrame parseBamFileCpp(String fileName, int binSize);
RcppExport SEXP _tenxchecker_parseBamFileCpp(SEXP fileNameSEXP, SEXP binSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int >::type binSize(binSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(parseBamFileCpp(fileName, binSize));
    return rcpp_result_gen;
END_RCPP
}
// parseHicCpp
DataFrame parseHicCpp(std::string& fname, int resolution);
RcppExport SEXP _tenxchecker_parseHicCpp(SEXP fnameSEXP, SEXP resolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< int >::type resolution(resolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(parseHicCpp(fname, resolution));
    return rcpp_result_gen;
END_RCPP
}
// parsePafCpp
DataFrame parsePafCpp(std::string& fname, int resolution, int minAlnLen, int minCount, int minNCells);
RcppExport SEXP _tenxchecker_parsePafCpp(SEXP fnameSEXP, SEXP resolutionSEXP, SEXP minAlnLenSEXP, SEXP minCountSEXP, SEXP minNCellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< int >::type resolution(resolutionSEXP);
    Rcpp::traits::input_parameter< int >::type minAlnLen(minAlnLenSEXP);
    Rcpp::traits::input_parameter< int >::type minCount(minCountSEXP);
    Rcpp::traits::input_parameter< int >::type minNCells(minNCellsSEXP);
    rcpp_result_gen = Rcpp::wrap(parsePafCpp(fname, resolution, minAlnLen, minCount, minNCells));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tenxchecker_parseBamFileCpp", (DL_FUNC) &_tenxchecker_parseBamFileCpp, 2},
    {"_tenxchecker_parseHicCpp", (DL_FUNC) &_tenxchecker_parseHicCpp, 2},
    {"_tenxchecker_parsePafCpp", (DL_FUNC) &_tenxchecker_parsePafCpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_tenxchecker(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
