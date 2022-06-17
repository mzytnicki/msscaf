#include <vector>
#include <cmath>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "sharedFunctions.h"
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// Compute the difference between two (distance, count) distribution
// Supposes that the observed and background counts are sorted by distance, with no missing value
// Reads nElements from each distribution (with possible offsets), and compute the distance
double cornerDifferenceCpp (DataFrame observedDistribution, DataFrame backgroundDistribution, int nElements, int observedOffset = 0, int backgroundOffset = 0) {
    NumericVector observedCount   = observedDistribution["count"];
    NumericVector backgroundCount = backgroundDistribution["count"];
    int nObserved = 0;
    double s = 0;
    for (int i = 0; i < nElements; ++i) {
        // Skip NAs
        if (! NumericVector::is_na(observedCount[i])) {
            double a = observedCount[i + observedOffset] - backgroundCount[i + backgroundOffset];
            s += a * a;
            ++nObserved;
        }
    }
    if (nObserved == 0) return -1;
    return (sqrt(s) / nObserved);
}

// Apply an offset to the background distribution, and compare them.
// [[Rcpp::export]]
double computeCornerDifferenceOffsetCpp (int offset, DataFrame corner, DataFrame background, int maxDistance) {
    return cornerDifferenceCpp(corner, background, maxDistance - offset, 0, offset);
}

// Apply an offset to the corner and background distribution, and compare them.
// [[Rcpp::export]]
double computeCornerDifferenceBothOffsetCpp (int offset, DataFrame corner, DataFrame background, int maxDistance) {
    return cornerDifferenceCpp(corner, background, maxDistance - offset, 0, 0);
}

// Compute a 2-column tibble, with distance and count
// Sample of handful of lines
// [[Rcpp::export]]
DataFrame estimateDistanceCountCpp (DataFrame &data, DataFrame &outliers, IntegerVector &sizesIn, int distance, int metaSize, int nOutputElements) {
    bool  useMetaBins    = (metaSize > 1);
    IntegerVector sizes  = clone(sizesIn);
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    IntegerVector outlierRefs  = outliers["ref"];
    IntegerVector outlierBins  = outliers["bin"];
    std::vector < long int > offsets (sizes.size());
    std::vector < std::vector < double > > fullMatrix;
    std::vector < bool > outlierBinsBool;
    if (useMetaBins) {
        for (int i = 0; i < sizes.size(); ++i) {
            sizes[i] = sizes[i] / metaSize;
        }
    }
    long int nElements = computeOffsets(sizes, offsets);
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    if (nOutputElements > nElements) {
        nOutputElements = nElements;
    }
    IntegerVector sampledLines = sample(nElements, nOutputElements, false, R_NilValue, false);
    std::vector < bool > sampledLinesBool (nElements, false);
    std::vector < int > sampleIndex (nElements, 0);
    for (int i = 0; i < sampledLines.size(); ++i) {
        int sampledLine = sampledLines[i];
        if (sampledLine >= nElements) stop("Problem while sampling.");
        sampleIndex[sampledLine] = i;
        sampledLinesBool[sampledLine] = true;
    }
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    long int n = refs1.size();
    IntegerVector outputCounts    (nOutputElements * (distance + 1), 0.0);
    IntegerVector outputDistances (nOutputElements * (distance + 1), 0.0);
    int cpt = 0;
    for (long int i = 0; i < nOutputElements; ++i) {
        for (int d = 0; d <= distance; ++d) {
            outputDistances[cpt++] = d;
        }
    }
    if (useMetaBins) {
        Progress progress (n, true);
        for (long int i = 0; i < n; ++i) {
            int ref1 = refs1[i];
            int ref2 = refs2[i];
            if (ref1 == ref2) {
                int offset = offsets[ref1 - 1];
                int bin1 = bins1[i];
                int bin2 = bins2[i];
                int d    = bin1 - bin2;
                bin1 /= metaSize;
                bin2 /= metaSize;
                if ((! outlierBinsBool[offset + bin1]) && (! outlierBinsBool[offset + bin2]) && (sampledLinesBool[offset + bin1])) {
                    if (d < 0) stop("Bin1 should not be less than bin2 in 'setFullMatrix'.\n");
                    d /= metaSize;
                    if ((counts[i] > 0) && (d <= distance)) {
                        outputCounts[sampleIndex[offset + bin1] * (distance + 1) + d] += counts[i];
                    }
                }
            }
            progress.increment();
        }
    }
    else {
        Progress progress (n, true);
        for (long int i = 0; i < n; ++i) {
            int ref1 = refs1[i];
            int ref2 = refs2[i];
            if (ref1 == ref2) {
                assert(ref1 <= static_cast<int>(offsets.size()));
                int offset = offsets[ref1 - 1];
                int bin1 = bins1[i];
                int bin2 = bins2[i];
                assert(offset + bin1 < static_cast<int>(outlierBinsBool.size()));
                assert(offset + bin1 < static_cast<int>(sampledLinesBool.size()));
                if ((! outlierBinsBool[offset + bin1]) && (! outlierBinsBool[offset + bin2]) && (sampledLinesBool[offset + bin1])) {
                    int d = bin1 - bin2;
                    assert(d >= 0);
                    if ((counts[i] > 0) && (d <= distance)) {
                        if (sampleIndex[offset + bin1] * (distance + 1) + d >= outputCounts.size()) stop("Count is out of range 'estimateDistanceCountCpp'.\n");
                        outputCounts[sampleIndex[offset + bin1] * (distance + 1) + d] = counts[i];
                    }
                }
            }
            progress.increment();
        }
    }
/*
    int nRefs = sizes.size();
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    long int      n      = bins1.size();
    IntegerVector outlierRefs  = outliers["ref"];
    IntegerVector outlierBins  = outliers["bin"];
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsets(sizes, offsets);
    std::vector < std::vector < double > > fullMatrix;
    std::vector < bool > outlierBinsBool;
    IntegerVector outputDistance (n * (distance + 1));
    IntegerVector outputCount (n * (distance + 1));
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    setFullMatrix(refs1, refs2, bins1, bins2, counts, offsets, nElements, outlierBinsBool, distance, fullMatrix);
    long int cpt = 0;
    for (int refId = 0; refId < nRefs; ++refId) {
        int offset = offsets[refId];
        for (int bin1 = 0; bin1 <= sizes[refId]; ++bin1) {
            if (! outlierBinsBool[offset + bin1]) {
                for (int d = 0; (d <= distance) && (d <= bin1); ++d) {
                    int bin2 = bin1 - d;
                    if (! outlierBinsBool[offset + bin2]) {
                        outputDistance[cpt] = d;
                        outputCount[cpt] = fullMatrix[bin1][d];
                        ++cpt;
                    }
                }
            }
        }
    }
    outputDistance.erase(cpt, n * (distance + 1));
    outputCount.erase(cpt, n * (distance + 1));
*/
    return Rcpp::DataFrame::create(_["distance"] = outputDistances, _["count"] = outputCounts);
}


// Find random positions
// Extract triangles near them
//    bin1 ->
//   +---------------------------
// b |+     
// i | +   | 
// n |  +  |X 
// 2 |   + |XX 
//   |    +|XXX 
// | | bin +----
// v |      +                    
//   |       +                   
// [[Rcpp::export]]
DataFrame sampleTriangles (DataFrame &data, DataFrame &outliers, IntegerVector &sizesIn, int distance, int metaSize, int nSamples) {
    bool  useMetaBins    = (metaSize > 1);
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    IntegerVector outlierRefs  = outliers["ref"];
    IntegerVector outlierBins  = outliers["bin"];
    IntegerVector sizes  = clone(sizesIn);
    CharacterVector refNames = sizes.names();
    std::vector < long int > offsets (sizes.size());
    std::vector < std::vector < double > > fullMatrix;
    std::vector < bool > outlierBinsBool;
    if (useMetaBins) {
        for (int i = 0; i < sizes.size(); ++i) {
            sizes[i] = sizes[i] / metaSize;
        }
    }
    long int nElements   = computeOffsets(sizes, offsets);
    int      nMaxSamples = std::max<int>(nElements / (distance + 1), 1);
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    if (nSamples > nMaxSamples) {
        nSamples = nMaxSamples;
    }
    // Sample few elements, which should not overlap
    IntegerVector sampledBins = sample(nMaxSamples, nSamples, false, R_NilValue, false);
    for (int i = 0; i < sampledBins.size(); ++i) {
        sampledBins[i] *= distance + 1;
    }
    std::vector < int > sampledIndices (nElements + 1, -1);
    for (int i = 0; i < sampledBins.size(); ++i) {
        int sampledBin = sampledBins[i];
        // May have a problem if the sampled line is near the end.
        for (int d = 0; (d <= distance) && (sampledBin + d <= nElements); ++d) {
            if (sampledIndices[sampledBin + d] != -1) {
                Rcerr << "Problem in the sampling: data overlap.\n";
            }
            sampledIndices[sampledBin + d] = i;
        }
    }
    long int n = refs1.size();
    int nElementsPerSample = (distance + 1) * (distance + 2) / 2;
    long int nOutput = nSamples * nElementsPerSample;
    IntegerVector outputIndex     (nOutput, 0);
    IntegerVector outputCounts    (nOutput, 0);
    IntegerVector outputDistances (nOutput, 0);
    long int cpt = 0;
    for (long int i = 0; i < nSamples; ++i) {
        for (long int j = 0; j <= distance; ++j) {
            for (long int k = 0; k <= j; ++k) {
                outputIndex[cpt] = i;
                outputDistances[cpt] = j;
                ++cpt;
            }
        }
    }
    if (cpt != nOutput) Rcerr << "Problem while filling output matrix.\n";
    if (useMetaBins) {
        Progress progress (n, true);
        for (long int i = 0; i < n; ++i) {
            int ref1 = refs1[i];
            int ref2 = refs2[i];
            if (ref1 == ref2) {
                int bin1 = bins1[i] / metaSize;
                int bin2 = bins2[i] / metaSize;
                int offset   = offsets[ref1 - 1];
                int metaBin1 = offset + bins1[i] / metaSize;
                int metaBin2 = offset + bins2[i] / metaSize;
                if ((! outlierBinsBool[metaBin1]) && (! outlierBinsBool[metaBin2])) {
                    int sampledIndex = sampledIndices[metaBin1];
                    if (sampledIndex >= 0) {
                        int bin = sampledBins[sampledIndex];
//Rcout << "Reading " << ref1 << " " << bin1 << "-" << bin2 << " -> " << metaBin1 << "-" << metaBin2 << " -> " << bin << ", " << sampleId;
                        bin1 = metaBin1 - bin;
                        bin2  = bin - metaBin2;
                        int d = bin1 + bin2;
                        if ((bin2 >= 0) && (d <= distance)) {
                            int index = d * (d + 1) / 2 + bin1;
// Rcout << " -> " << bin1 << "-" << bin2 << " (" << index << ") -> " << (bin1 + bin2) << ", " << counts[i] << " -- " << outputCounts[sampleId * nElementsPerSample + index] << "\n";
// if (bin1 + bin2 != outputDistances[sampleId * nElementsPerSample + index]) Rcerr << "Problem 1 in estimate corner variance\n";
// if (sampleId != outputIndex[sampleId * nElementsPerSample + index]) Rcerr << "Problem 2 in estimate corner variance\n";
                            outputCounts[sampledIndex * nElementsPerSample + index] += counts[i];
                        }
                    }
                }
            }
            progress.increment();
        }
    }
    else {
        Progress progress (n, true);
        for (long int i = 0; i < n; ++i) {
            int ref1 = refs1[i];
            int ref2 = refs2[i];
            if (ref1 == ref2) {
                int offset = offsets[ref1 - 1];
                int bin1 = offset + bins1[i];
                int bin2 = offset + bins2[i];
                if ((! outlierBinsBool[bin1]) && (! outlierBinsBool[bin2])) {
                    int sampledIndex = sampledIndices[bin1];
                    // if ((sampledIndex >= 0) && (sampledIndex == sampledIndices[bin2])) {
                    if (sampledIndex >= 0) {
                        int bin = sampledBins[sampledIndex];
                        bin1 -= bin;
                        bin2  = bin - bin2;
                        int d = bin1 + bin2;
                        if ((bin2 >= 0) && (d <= distance)) {
                            int index = sampledIndex * nElementsPerSample + d * (d + 1) / 2 + bin1;
                            assert(outputCounts[index] == 0);
                            assert(outputDistances[index] == d);
                            outputCounts[index] = counts[i];
                        }
                    }
                }
            }
            progress.increment();
        }
    }
/*
    int nRefs = sizes.size();
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    long int      n      = bins1.size();
    IntegerVector outlierRefs  = outliers["ref"];
    IntegerVector outlierBins  = outliers["bin"];
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsets(sizes, offsets);
    std::vector < std::vector < double > > fullMatrix;
    std::vector < bool > outlierBinsBool;
    IntegerVector outputDistance (n * (distance + 1));
    IntegerVector outputCount (n * (distance + 1));
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    setFullMatrix(refs1, refs2, bins1, bins2, counts, offsets, nElements, outlierBinsBool, distance, fullMatrix);
    long int cpt = 0;
    for (int refId = 0; refId < nRefs; ++refId) {
        int offset = offsets[refId];
        for (int bin1 = 0; bin1 <= sizes[refId]; ++bin1) {
            if (! outlierBinsBool[offset + bin1]) {
                for (int d = 0; (d <= distance) && (d <= bin1); ++d) {
                    int bin2 = bin1 - d;
                    if (! outlierBinsBool[offset + bin2]) {
                        outputDistance[cpt] = d;
                        outputCount[cpt] = fullMatrix[bin1][d];
                        ++cpt;
                    }
                }
            }
        }
    }
    outputDistance.erase(cpt, n * (distance + 1));
    outputCount.erase(cpt, n * (distance + 1));
*/
    return Rcpp::DataFrame::create(_["index"] = outputIndex, _["distance"] = outputDistances, _["count"] = outputCounts);
}


// Check whether bins should be aggregated
// There should at least nMeta meta-bins with sum count minCount.
// The meta bins are packed horizontally and vertically
// [[Rcpp::export]]
int estimateMetaSizeCpp (std::vector < double > &rowAvg, int maxDistance, int nMeta, int minCount) {
    int currentD = 1;
    int firstD = currentD;
    int metaSize = 1;
//Rcerr << "max size: " << maxDistance << ", n meta: " << nMeta << ", min count: " << minCount << "\n";
    for (int metaId = 0; metaId < nMeta; ++metaId) {
        double sumAverages = 0.0;
        firstD = currentD;
//Rcerr << "meta " << metaId << "-> " << currentD << "\n";
        for (int d = 0; d < metaSize; ++d, ++currentD) {
            sumAverages += rowAvg[currentD];
        }
//Rcerr << "  sum: " << sumAverages << " * " << metaSize << " = " << (sumAverages * metaSize) << "\n";
        for (; sumAverages * metaSize <= minCount; ++currentD, ++metaSize) {
            for (int d = 0; d < metaId; ++d, ++firstD, ++currentD) {
                if (currentD == maxDistance) return 0;
                sumAverages += rowAvg[currentD] - rowAvg[firstD];
            }
            if (currentD == maxDistance) return 0;
            sumAverages += rowAvg[currentD];
//Rcerr << "  d: " << currentD << " - " << firstD << " = " << metaSize << ": " << sumAverages << " * " << metaSize << " = " << (sumAverages * metaSize) << "\n";
        }
//Rcerr << "  metaSize: " << metaSize << "\n";
    }
    return metaSize;
}


// Check the length of the signal
// [[Rcpp::export]]
int estimateMoleculeSizeCpp (std::vector < double > &metaSums, int maxDistance, int minCount, int metaSize) {
    for (int metaD = 1; metaD < maxDistance / metaSize; ++metaD) {
//Rcerr << metaD << " -> " << (metaSums[metaD] * metaSize) << " / " << minCount << "\n";
        if (metaSums[metaD] * metaSize < minCount) {
            return metaD - 1;
        }
    }
    return 0;
}

// Check the max length of the signal
// [[Rcpp::export]]
int estimateMaxMoleculeSizeCpp (std::vector < double > &metaSums, int minDistance, int maxDistance, int minCount, int metaSize) {
    for (int metaD = minDistance; metaD < maxDistance / metaSize - 1; ++metaD) {
        if (metaSums[metaD + 1] >= 0.99 * metaSums[metaD]) {
            return metaD;
        }
    }
    return minDistance;
}

// Compute meta bins size and molecule sizes
// [[Rcpp::export]]
List estimateMetaBinsMoleculeSizeCpp (DataFrame &data, IntegerVector &sizes, int minCount, int nMeta, bool moleculeSize) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsets(sizes, offsets);
    std::vector < bool > outlierBinsBool (nElements, false);
    std::vector < long int > rowSums;
    std::vector < int > nRowSums;
    int maxDistance = max(sizes);
    Rcout << "\tComputing distance/count matrix.\n";
    computeDistanceCounts(refs1, refs2, bins1, bins2, counts, offsets, outlierBinsBool, maxDistance, sizes, rowSums, nRowSums);
    std::vector < double > rowAvg (rowSums.size());
    for (unsigned int i = 0; i < rowSums.size(); ++i) {
        if (nRowSums[i] == 0) {
            rowAvg[i] = 0.0;
        }
        else {
            rowAvg[i] = static_cast < double > (rowSums[i]) / nRowSums[i];
        }
    }
    Rcout << "\tGrouping bins to meta-bins.\n";
    int metaSize = estimateMetaSizeCpp(rowAvg, maxDistance, nMeta, minCount);
    int minLinkRange = nMeta;
    int maxLinkRange = minLinkRange;
    Rcout << "\tEstimating molecule size.\n";
    std::vector < double > metaSums(rowAvg.size());
    // Compute average count per distance
    for (int d = 1; d <= maxDistance; ++d) {
        int metaD = 1 + (d - 1) / metaSize;
        metaSums[metaD] += rowAvg[d];
    }
    if ((! moleculeSize) && (metaSize == 1)) {
        minLinkRange = estimateMoleculeSizeCpp(metaSums, maxDistance, minCount, metaSize);
    }
    maxLinkRange = estimateMaxMoleculeSizeCpp(metaSums, minLinkRange, maxDistance, minCount, metaSize);
    return List::create(_["metaSize"] = metaSize, _["minLinkRange"] = minLinkRange, _["maxLinkRange"] = maxLinkRange);
}
