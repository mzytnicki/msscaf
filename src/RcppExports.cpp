// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// parseBamFileCpp
DataFrame parseBamFileCpp(String fileName, int32_t binSize);
RcppExport SEXP _msscaf_parseBamFileCpp(SEXP fileNameSEXP, SEXP binSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int32_t >::type binSize(binSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(parseBamFileCpp(fileName, binSize));
    return rcpp_result_gen;
END_RCPP
}
// computeMeanTrianglesCpp
DataFrame computeMeanTrianglesCpp(DataFrame& data, int distance, int metaSize, IntegerVector sizesIn, DataFrame outliers);
RcppExport SEXP _msscaf_computeMeanTrianglesCpp(SEXP dataSEXP, SEXP distanceSEXP, SEXP metaSizeSEXP, SEXP sizesInSEXP, SEXP outliersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type distance(distanceSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type outliers(outliersSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMeanTrianglesCpp(data, distance, metaSize, sizesIn, outliers));
    return rcpp_result_gen;
END_RCPP
}
// removeNearEqualBreaksCpp
DataFrame removeNearEqualBreaksCpp(DataFrame& breaks, int distance);
RcppExport SEXP _msscaf_removeNearEqualBreaksCpp(SEXP breaksSEXP, SEXP distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< int >::type distance(distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(removeNearEqualBreaksCpp(breaks, distance));
    return rcpp_result_gen;
END_RCPP
}
// computeCornerSize
int computeCornerSize(int size1, int size2, int distance);
RcppExport SEXP _msscaf_computeCornerSize(SEXP size1SEXP, SEXP size2SEXP, SEXP distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type size1(size1SEXP);
    Rcpp::traits::input_parameter< int >::type size2(size2SEXP);
    Rcpp::traits::input_parameter< int >::type distance(distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeCornerSize(size1, size2, distance));
    return rcpp_result_gen;
END_RCPP
}
// computeOtherSize
int computeOtherSize(int size1, int size2, int distance);
RcppExport SEXP _msscaf_computeOtherSize(SEXP size1SEXP, SEXP size2SEXP, SEXP distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type size1(size1SEXP);
    Rcpp::traits::input_parameter< int >::type size2(size2SEXP);
    Rcpp::traits::input_parameter< int >::type distance(distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeOtherSize(size1, size2, distance));
    return rcpp_result_gen;
END_RCPP
}
// filterCornersCpp
DataFrame filterCornersCpp(DataFrame data, IntegerVector sizesIn, int cornerSize, int metaSize);
RcppExport SEXP _msscaf_filterCornersCpp(SEXP dataSEXP, SEXP sizesInSEXP, SEXP cornerSizeSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< int >::type cornerSize(cornerSizeSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(filterCornersCpp(data, sizesIn, cornerSize, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// sumCornerCpp
DataFrame sumCornerCpp(DataFrame interactions, DataFrame outliers, IntegerVector sizesIn, int cornerSize, int metaSize);
RcppExport SEXP _msscaf_sumCornerCpp(SEXP interactionsSEXP, SEXP outliersSEXP, SEXP sizesInSEXP, SEXP cornerSizeSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type outliers(outliersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< int >::type cornerSize(cornerSizeSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sumCornerCpp(interactions, outliers, sizesIn, cornerSize, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// extractCornersCpp
DataFrame extractCornersCpp(DataFrame interactions, DataFrame selectedCorners, DataFrame outliers, IntegerVector sizesIn, int minCornerSize, int maxCornerSize, int metaSize);
RcppExport SEXP _msscaf_extractCornersCpp(SEXP interactionsSEXP, SEXP selectedCornersSEXP, SEXP outliersSEXP, SEXP sizesInSEXP, SEXP minCornerSizeSEXP, SEXP maxCornerSizeSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type selectedCorners(selectedCornersSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type outliers(outliersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< int >::type minCornerSize(minCornerSizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxCornerSize(maxCornerSizeSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(extractCornersCpp(interactions, selectedCorners, outliers, sizesIn, minCornerSize, maxCornerSize, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// classifyCornerPointsCpp
DataFrame classifyCornerPointsCpp(DataFrame interactions, int size1, int size2, int metaSize, int maxDistance);
RcppExport SEXP _msscaf_classifyCornerPointsCpp(SEXP interactionsSEXP, SEXP size1SEXP, SEXP size2SEXP, SEXP metaSizeSEXP, SEXP maxDistanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< int >::type size1(size1SEXP);
    Rcpp::traits::input_parameter< int >::type size2(size2SEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    rcpp_result_gen = Rcpp::wrap(classifyCornerPointsCpp(interactions, size1, size2, metaSize, maxDistance));
    return rcpp_result_gen;
END_RCPP
}
// extractCornersFullCpp
DataFrame extractCornersFullCpp(DataFrame interactions, DataFrame selectedCorners, DataFrame outliers, IntegerVector sizesIn, int cornerSize, int metaSize);
RcppExport SEXP _msscaf_extractCornersFullCpp(SEXP interactionsSEXP, SEXP selectedCornersSEXP, SEXP outliersSEXP, SEXP sizesInSEXP, SEXP cornerSizeSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type selectedCorners(selectedCornersSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type outliers(outliersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< int >::type cornerSize(cornerSizeSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(extractCornersFullCpp(interactions, selectedCorners, outliers, sizesIn, cornerSize, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// computeCornerDifferenceOffsetCpp
double computeCornerDifferenceOffsetCpp(int offset, DataFrame corner, DataFrame background, int maxDistance);
RcppExport SEXP _msscaf_computeCornerDifferenceOffsetCpp(SEXP offsetSEXP, SEXP cornerSEXP, SEXP backgroundSEXP, SEXP maxDistanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type corner(cornerSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type background(backgroundSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeCornerDifferenceOffsetCpp(offset, corner, background, maxDistance));
    return rcpp_result_gen;
END_RCPP
}
// computeCornerDifferenceBothOffsetCpp
double computeCornerDifferenceBothOffsetCpp(int offset, DataFrame corner, DataFrame background, int maxDistance);
RcppExport SEXP _msscaf_computeCornerDifferenceBothOffsetCpp(SEXP offsetSEXP, SEXP cornerSEXP, SEXP backgroundSEXP, SEXP maxDistanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type corner(cornerSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type background(backgroundSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeCornerDifferenceBothOffsetCpp(offset, corner, background, maxDistance));
    return rcpp_result_gen;
END_RCPP
}
// estimateDistanceCountCpp
DataFrame estimateDistanceCountCpp(DataFrame& data, DataFrame& outliers, IntegerVector& sizesIn, int distance, int metaSize, int nOutputElements);
RcppExport SEXP _msscaf_estimateDistanceCountCpp(SEXP dataSEXP, SEXP outliersSEXP, SEXP sizesInSEXP, SEXP distanceSEXP, SEXP metaSizeSEXP, SEXP nOutputElementsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< DataFrame& >::type outliers(outliersSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< int >::type distance(distanceSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    Rcpp::traits::input_parameter< int >::type nOutputElements(nOutputElementsSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateDistanceCountCpp(data, outliers, sizesIn, distance, metaSize, nOutputElements));
    return rcpp_result_gen;
END_RCPP
}
// sampleTriangles
DataFrame sampleTriangles(DataFrame& data, DataFrame& outliers, IntegerVector& sizesIn, int distance, int metaSize, int nSamples);
RcppExport SEXP _msscaf_sampleTriangles(SEXP dataSEXP, SEXP outliersSEXP, SEXP sizesInSEXP, SEXP distanceSEXP, SEXP metaSizeSEXP, SEXP nSamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< DataFrame& >::type outliers(outliersSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizesIn(sizesInSEXP);
    Rcpp::traits::input_parameter< int >::type distance(distanceSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    Rcpp::traits::input_parameter< int >::type nSamples(nSamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleTriangles(data, outliers, sizesIn, distance, metaSize, nSamples));
    return rcpp_result_gen;
END_RCPP
}
// estimateMetaSizeCpp
int estimateMetaSizeCpp(std::vector < double >& rowAvg, int maxDistance, int nMeta, int minCount);
RcppExport SEXP _msscaf_estimateMetaSizeCpp(SEXP rowAvgSEXP, SEXP maxDistanceSEXP, SEXP nMetaSEXP, SEXP minCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector < double >& >::type rowAvg(rowAvgSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    Rcpp::traits::input_parameter< int >::type nMeta(nMetaSEXP);
    Rcpp::traits::input_parameter< int >::type minCount(minCountSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateMetaSizeCpp(rowAvg, maxDistance, nMeta, minCount));
    return rcpp_result_gen;
END_RCPP
}
// estimateMoleculeSizeCpp
int estimateMoleculeSizeCpp(std::vector < double >& metaSums, int maxDistance, int minCount, int metaSize);
RcppExport SEXP _msscaf_estimateMoleculeSizeCpp(SEXP metaSumsSEXP, SEXP maxDistanceSEXP, SEXP minCountSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector < double >& >::type metaSums(metaSumsSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    Rcpp::traits::input_parameter< int >::type minCount(minCountSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateMoleculeSizeCpp(metaSums, maxDistance, minCount, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// estimateMaxMoleculeSizeCpp
int estimateMaxMoleculeSizeCpp(std::vector < double >& metaSums, int minDistance, int maxDistance, int minCount, int metaSize);
RcppExport SEXP _msscaf_estimateMaxMoleculeSizeCpp(SEXP metaSumsSEXP, SEXP minDistanceSEXP, SEXP maxDistanceSEXP, SEXP minCountSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector < double >& >::type metaSums(metaSumsSEXP);
    Rcpp::traits::input_parameter< int >::type minDistance(minDistanceSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    Rcpp::traits::input_parameter< int >::type minCount(minCountSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateMaxMoleculeSizeCpp(metaSums, minDistance, maxDistance, minCount, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// estimateMetaBinsMoleculeSizeCpp
List estimateMetaBinsMoleculeSizeCpp(DataFrame& data, IntegerVector& sizes, int minCount, int nMeta, bool moleculeSize);
RcppExport SEXP _msscaf_estimateMetaBinsMoleculeSizeCpp(SEXP dataSEXP, SEXP sizesSEXP, SEXP minCountSEXP, SEXP nMetaSEXP, SEXP moleculeSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< int >::type minCount(minCountSEXP);
    Rcpp::traits::input_parameter< int >::type nMeta(nMetaSEXP);
    Rcpp::traits::input_parameter< bool >::type moleculeSize(moleculeSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateMetaBinsMoleculeSizeCpp(data, sizes, minCount, nMeta, moleculeSize));
    return rcpp_result_gen;
END_RCPP
}
// keepScaffoldsCpp
DataFrame keepScaffoldsCpp(DataFrame& data, CharacterVector keptRefs);
RcppExport SEXP _msscaf_keepScaffoldsCpp(SEXP dataSEXP, SEXP keptRefsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type keptRefs(keptRefsSEXP);
    rcpp_result_gen = Rcpp::wrap(keepScaffoldsCpp(data, keptRefs));
    return rcpp_result_gen;
END_RCPP
}
// keepScaffoldsPairsCpp
DataFrame keepScaffoldsPairsCpp(DataFrame& data, DataFrame& keptRefs);
RcppExport SEXP _msscaf_keepScaffoldsPairsCpp(SEXP dataSEXP, SEXP keptRefsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< DataFrame& >::type keptRefs(keptRefsSEXP);
    rcpp_result_gen = Rcpp::wrap(keepScaffoldsPairsCpp(data, keptRefs));
    return rcpp_result_gen;
END_RCPP
}
// extractLines
DataFrame extractLines(DataFrame matrix, DataFrame lines, int maxDistance);
RcppExport SEXP _msscaf_extractLines(SEXP matrixSEXP, SEXP linesSEXP, SEXP maxDistanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lines(linesSEXP);
    Rcpp::traits::input_parameter< int >::type maxDistance(maxDistanceSEXP);
    rcpp_result_gen = Rcpp::wrap(extractLines(matrix, lines, maxDistance));
    return rcpp_result_gen;
END_RCPP
}
// splitChromosomeCpp
void splitChromosomeCpp(DataFrame& data, int prevRef, int newRef, int shiftedRef, long splitPoint, bool firstPart);
RcppExport SEXP _msscaf_splitChromosomeCpp(SEXP dataSEXP, SEXP prevRefSEXP, SEXP newRefSEXP, SEXP shiftedRefSEXP, SEXP splitPointSEXP, SEXP firstPartSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type prevRef(prevRefSEXP);
    Rcpp::traits::input_parameter< int >::type newRef(newRefSEXP);
    Rcpp::traits::input_parameter< int >::type shiftedRef(shiftedRefSEXP);
    Rcpp::traits::input_parameter< long >::type splitPoint(splitPointSEXP);
    Rcpp::traits::input_parameter< bool >::type firstPart(firstPartSEXP);
    splitChromosomeCpp(data, prevRef, newRef, shiftedRef, splitPoint, firstPart);
    return R_NilValue;
END_RCPP
}
// computeRefSizesCpp
IntegerVector computeRefSizesCpp(DataFrame& data);
RcppExport SEXP _msscaf_computeRefSizesCpp(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(computeRefSizesCpp(data));
    return rcpp_result_gen;
END_RCPP
}
// computeNRrows
int computeNRrows(DataFrame& data, IntegerVector& sizes);
RcppExport SEXP _msscaf_computeNRrows(SEXP dataSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizes(sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeNRrows(data, sizes));
    return rcpp_result_gen;
END_RCPP
}
// computeSymmetricColSum
DataFrame computeSymmetricColSum(DataFrame& data, IntegerVector& sizes);
RcppExport SEXP _msscaf_computeSymmetricColSum(SEXP dataSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizes(sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeSymmetricColSum(data, sizes));
    return rcpp_result_gen;
END_RCPP
}
// computeSymmetricColSumMeta
DataFrame computeSymmetricColSumMeta(DataFrame& data, IntegerVector& sizes, int metaSize);
RcppExport SEXP _msscaf_computeSymmetricColSumMeta(SEXP dataSEXP, SEXP sizesSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(computeSymmetricColSumMeta(data, sizes, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// removeLowCountRowsCpp
List removeLowCountRowsCpp(DataFrame& data, IntegerVector& sizes, int threshold);
RcppExport SEXP _msscaf_removeLowCountRowsCpp(SEXP dataSEXP, SEXP sizesSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< int >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(removeLowCountRowsCpp(data, sizes, threshold));
    return rcpp_result_gen;
END_RCPP
}
// normalizeHighCountRowsCpp
void normalizeHighCountRowsCpp(DataFrame& data, IntegerVector& sizes);
RcppExport SEXP _msscaf_normalizeHighCountRowsCpp(SEXP dataSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type sizes(sizesSEXP);
    normalizeHighCountRowsCpp(data, sizes);
    return R_NilValue;
END_RCPP
}
// parseHicCpp
DataFrame parseHicCpp(std::string& fname, int resolution);
RcppExport SEXP _msscaf_parseHicCpp(SEXP fnameSEXP, SEXP resolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< int >::type resolution(resolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(parseHicCpp(fname, resolution));
    return rcpp_result_gen;
END_RCPP
}
// makeSymmetricRefCpp
DataFrame makeSymmetricRefCpp(DataFrame& data);
RcppExport SEXP _msscaf_makeSymmetricRefCpp(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(makeSymmetricRefCpp(data));
    return rcpp_result_gen;
END_RCPP
}
// removeOutliersRefCpp
DataFrame removeOutliersRefCpp(DataFrame& data, DataFrame& outliersTibble, int size, int minLim, int maxLim);
RcppExport SEXP _msscaf_removeOutliersRefCpp(SEXP dataSEXP, SEXP outliersTibbleSEXP, SEXP sizeSEXP, SEXP minLimSEXP, SEXP maxLimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< DataFrame& >::type outliersTibble(outliersTibbleSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type minLim(minLimSEXP);
    Rcpp::traits::input_parameter< int >::type maxLim(maxLimSEXP);
    rcpp_result_gen = Rcpp::wrap(removeOutliersRefCpp(data, outliersTibble, size, minLim, maxLim));
    return rcpp_result_gen;
END_RCPP
}
// removeOutliersCpp
DataFrame removeOutliersCpp(DataFrame& data, DataFrame& outliersTibble, IntegerVector sizes);
RcppExport SEXP _msscaf_removeOutliersCpp(SEXP dataSEXP, SEXP outliersTibbleSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< DataFrame& >::type outliersTibble(outliersTibbleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizes(sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(removeOutliersCpp(data, outliersTibble, sizes));
    return rcpp_result_gen;
END_RCPP
}
// parsePafCpp
DataFrame parsePafCpp(std::string& fname, uint32_t resolution, int minAlnLen, int minCount, int minNCells);
RcppExport SEXP _msscaf_parsePafCpp(SEXP fnameSEXP, SEXP resolutionSEXP, SEXP minAlnLenSEXP, SEXP minCountSEXP, SEXP minNCellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< uint32_t >::type resolution(resolutionSEXP);
    Rcpp::traits::input_parameter< int >::type minAlnLen(minAlnLenSEXP);
    Rcpp::traits::input_parameter< int >::type minCount(minCountSEXP);
    Rcpp::traits::input_parameter< int >::type minNCells(minNCellsSEXP);
    rcpp_result_gen = Rcpp::wrap(parsePafCpp(fname, resolution, minAlnLen, minCount, minNCells));
    return rcpp_result_gen;
END_RCPP
}
// getRefOrders
List getRefOrders(DataFrame& joins, IntegerVector sizes);
RcppExport SEXP _msscaf_getRefOrders(SEXP joinsSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type joins(joinsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizes(sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(getRefOrders(joins, sizes));
    return rcpp_result_gen;
END_RCPP
}
// scaffoldSizes
IntegerVector scaffoldSizes(List orders, IntegerVector orderIds, IntegerVector sizes);
RcppExport SEXP _msscaf_scaffoldSizes(SEXP ordersSEXP, SEXP orderIdsSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type orders(ordersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type orderIds(orderIdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizes(sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(scaffoldSizes(orders, orderIds, sizes));
    return rcpp_result_gen;
END_RCPP
}
// scaffoldContigs
CharacterVector scaffoldContigs(CharacterVector contigs, List orders, IntegerVector sizes, int binSize);
RcppExport SEXP _msscaf_scaffoldContigs(SEXP contigsSEXP, SEXP ordersSEXP, SEXP sizesSEXP, SEXP binSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type contigs(contigsSEXP);
    Rcpp::traits::input_parameter< List >::type orders(ordersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< int >::type binSize(binSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(scaffoldContigs(contigs, orders, sizes, binSize));
    return rcpp_result_gen;
END_RCPP
}
// scaffoldCounts
List scaffoldCounts(DataFrame matrices, DataFrame outlierBins, List groups, IntegerVector scaffoldRefs, IntegerVector sizes, int metaSize);
RcppExport SEXP _msscaf_scaffoldCounts(SEXP matricesSEXP, SEXP outlierBinsSEXP, SEXP groupsSEXP, SEXP scaffoldRefsSEXP, SEXP sizesSEXP, SEXP metaSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type matrices(matricesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type outlierBins(outlierBinsSEXP);
    Rcpp::traits::input_parameter< List >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type scaffoldRefs(scaffoldRefsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< int >::type metaSize(metaSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(scaffoldCounts(matrices, outlierBins, groups, scaffoldRefs, sizes, metaSize));
    return rcpp_result_gen;
END_RCPP
}
// splitCpp
S4 splitCpp(S4 object);
RcppExport SEXP _msscaf_splitCpp(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(splitCpp(object));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_msscaf_parseBamFileCpp", (DL_FUNC) &_msscaf_parseBamFileCpp, 2},
    {"_msscaf_computeMeanTrianglesCpp", (DL_FUNC) &_msscaf_computeMeanTrianglesCpp, 5},
    {"_msscaf_removeNearEqualBreaksCpp", (DL_FUNC) &_msscaf_removeNearEqualBreaksCpp, 2},
    {"_msscaf_computeCornerSize", (DL_FUNC) &_msscaf_computeCornerSize, 3},
    {"_msscaf_computeOtherSize", (DL_FUNC) &_msscaf_computeOtherSize, 3},
    {"_msscaf_filterCornersCpp", (DL_FUNC) &_msscaf_filterCornersCpp, 4},
    {"_msscaf_sumCornerCpp", (DL_FUNC) &_msscaf_sumCornerCpp, 5},
    {"_msscaf_extractCornersCpp", (DL_FUNC) &_msscaf_extractCornersCpp, 7},
    {"_msscaf_classifyCornerPointsCpp", (DL_FUNC) &_msscaf_classifyCornerPointsCpp, 5},
    {"_msscaf_extractCornersFullCpp", (DL_FUNC) &_msscaf_extractCornersFullCpp, 6},
    {"_msscaf_computeCornerDifferenceOffsetCpp", (DL_FUNC) &_msscaf_computeCornerDifferenceOffsetCpp, 4},
    {"_msscaf_computeCornerDifferenceBothOffsetCpp", (DL_FUNC) &_msscaf_computeCornerDifferenceBothOffsetCpp, 4},
    {"_msscaf_estimateDistanceCountCpp", (DL_FUNC) &_msscaf_estimateDistanceCountCpp, 6},
    {"_msscaf_sampleTriangles", (DL_FUNC) &_msscaf_sampleTriangles, 6},
    {"_msscaf_estimateMetaSizeCpp", (DL_FUNC) &_msscaf_estimateMetaSizeCpp, 4},
    {"_msscaf_estimateMoleculeSizeCpp", (DL_FUNC) &_msscaf_estimateMoleculeSizeCpp, 4},
    {"_msscaf_estimateMaxMoleculeSizeCpp", (DL_FUNC) &_msscaf_estimateMaxMoleculeSizeCpp, 5},
    {"_msscaf_estimateMetaBinsMoleculeSizeCpp", (DL_FUNC) &_msscaf_estimateMetaBinsMoleculeSizeCpp, 5},
    {"_msscaf_keepScaffoldsCpp", (DL_FUNC) &_msscaf_keepScaffoldsCpp, 2},
    {"_msscaf_keepScaffoldsPairsCpp", (DL_FUNC) &_msscaf_keepScaffoldsPairsCpp, 2},
    {"_msscaf_extractLines", (DL_FUNC) &_msscaf_extractLines, 3},
    {"_msscaf_splitChromosomeCpp", (DL_FUNC) &_msscaf_splitChromosomeCpp, 6},
    {"_msscaf_computeRefSizesCpp", (DL_FUNC) &_msscaf_computeRefSizesCpp, 1},
    {"_msscaf_computeNRrows", (DL_FUNC) &_msscaf_computeNRrows, 2},
    {"_msscaf_computeSymmetricColSum", (DL_FUNC) &_msscaf_computeSymmetricColSum, 2},
    {"_msscaf_computeSymmetricColSumMeta", (DL_FUNC) &_msscaf_computeSymmetricColSumMeta, 3},
    {"_msscaf_removeLowCountRowsCpp", (DL_FUNC) &_msscaf_removeLowCountRowsCpp, 3},
    {"_msscaf_normalizeHighCountRowsCpp", (DL_FUNC) &_msscaf_normalizeHighCountRowsCpp, 2},
    {"_msscaf_parseHicCpp", (DL_FUNC) &_msscaf_parseHicCpp, 2},
    {"_msscaf_makeSymmetricRefCpp", (DL_FUNC) &_msscaf_makeSymmetricRefCpp, 1},
    {"_msscaf_removeOutliersRefCpp", (DL_FUNC) &_msscaf_removeOutliersRefCpp, 5},
    {"_msscaf_removeOutliersCpp", (DL_FUNC) &_msscaf_removeOutliersCpp, 3},
    {"_msscaf_parsePafCpp", (DL_FUNC) &_msscaf_parsePafCpp, 5},
    {"_msscaf_getRefOrders", (DL_FUNC) &_msscaf_getRefOrders, 2},
    {"_msscaf_scaffoldSizes", (DL_FUNC) &_msscaf_scaffoldSizes, 3},
    {"_msscaf_scaffoldContigs", (DL_FUNC) &_msscaf_scaffoldContigs, 4},
    {"_msscaf_scaffoldCounts", (DL_FUNC) &_msscaf_scaffoldCounts, 6},
    {"_msscaf_splitCpp", (DL_FUNC) &_msscaf_splitCpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_msscaf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
