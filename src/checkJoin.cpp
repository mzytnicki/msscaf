#include <cstdlib>
#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "sharedFunctions.h"
                                                                                                                                                                                                                   
// [[Rcpp::plugins(cpp11)]]
                                                                                                                                                                                                                   
using namespace Rcpp;

enum CornerType { BB, EB, BE, EE, interior };
const int nCornerType { 4 };
const int nCornerIntType { 5 };

int getCornerDistanceCpp(int bin1, int bin2, int center1, int center2, int metaSize) {
    return std::abs((bin1 / metaSize) - (center1 / metaSize)) + std::abs((bin2 / metaSize) - (center2 / metaSize));
}

bool isInCorner(int bin1, int bin2, int center1, int center2, int metaSize, int maxDistance) {
    return getCornerDistanceCpp(bin1, bin2, center1, center2, metaSize) <= maxDistance;
}

// Compute corner size
// [[Rcpp::export]]
int computeCornerSize(int size1, int size2, int distance) {
    int cornerSize = distance * (distance + 1) / 2;
    size1 = distance - (size1 / 2);
    size2 = distance - (size2 / 2);
    if (size1 > 0) {
        cornerSize = cornerSize - (size1 * (size1 + 1) / 2);
    }
    if (size2 > 0) {
        cornerSize = cornerSize - (size2 * (size2 + 1) / 2);
    }
    return cornerSize;
}


// Compute non-corner size
// [[Rcpp::export]]
int computeOtherSize(int size1, int size2, int distance) {
    return size1 * size2 - computeCornerSize(size1, size2, distance);
}


// Read interaction matrix,
//   stores the number of non-zero for each corner,
//   output the non-near-empty corners.
// [[Rcpp::export]]
DataFrame filterCornersCpp (DataFrame data, IntegerVector sizesIn, int cornerSize, int metaSize) {
    IntegerVector sizes  = clone(sizesIn);
    for (int i = 0; i < sizes.size(); ++i) {
        sizes[i] = sizes[i] / metaSize;
    }
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    CharacterVector refs = refs1.attr("levels");
    int nRefs = sizes.size();
    std::vector < std::vector < std::array < int, nCornerType > > > cornerCounts (nRefs + 1); // indexes start with 1 in R
    std::vector < int > cornerRefs1;
    std::vector < int > cornerRefs2;
    std::vector < int > cornerTypes;
    int minCount = cornerSize * (cornerSize + 1) / 4; // This will be our threshold
    cornerSize *= metaSize;
    for (int refId = 0; refId < nRefs + 1; ++refId) {
        cornerCounts[refId] = std::vector < std::array < int, nCornerType > > (refId);
    }
    long long int nData = refs1.size();
    for (long long int dataId = 0; dataId < nData; ++dataId) {
        int ref1  = refs1[dataId];
        int ref2  = refs2[dataId];
        int bin1  = bins1[dataId];
        int bin2  = bins2[dataId];
        int count = counts[dataId];
        int size1 = sizes[ref1-1];
        int size2 = sizes[ref2-1];
        if ((ref1 != ref2) && (count > 0)) {
            assert(ref1 <= nRefs);
            assert(ref2 < ref1);
            if (isInCorner(bin1, bin2, 0, 0, metaSize, cornerSize)) {
                ++cornerCounts[ref1][ref2][BB];
            }
            if (isInCorner(bin1, bin2, size1, 0, metaSize, cornerSize)) {
                ++cornerCounts[ref1][ref2][EB];
            }
            if (isInCorner(bin1, bin2, 0, size2, metaSize, cornerSize)) {
                ++cornerCounts[ref1][ref2][BE];
            }
            if (isInCorner(bin1, bin2, size1, size2, metaSize, cornerSize)) {
                ++cornerCounts[ref1][ref2][EE];
            }
        }
    }
    for (int ref1 = 1; ref1 <= nRefs; ++ref1) {
        for (int ref2 = 1; ref2 < ref1; ++ref2) {
            for (int cornerType = 0; cornerType < nCornerType; ++cornerType) {
                if (cornerCounts[ref1][ref2][cornerType] >= minCount) {
                    cornerRefs1.push_back(ref1);
                    cornerRefs2.push_back(ref2);
                    cornerTypes.push_back(cornerType);
                }
            }
        }
    }
    CharacterVector cornerTypesFactor  = {"BB", "EB", "BE", "EE", "interior"};
    IntegerVector   cornerTypesR = wrap(cornerTypes);
    cornerTypesR = cornerTypesR + 1; // factors start with 1
    cornerTypesR.attr("class")   = "factor";
    cornerTypesR.attr("levels")  = cornerTypesFactor;
    IntegerVector   cornerRefs1R = wrap(cornerRefs1);
    cornerRefs1R.attr("class")   = "factor";
    cornerRefs1R.attr("levels")  = refs;
    IntegerVector   cornerRefs2R = wrap(cornerRefs2);
    cornerRefs2R.attr("class")   = "factor";
    cornerRefs2R.attr("levels")  = refs;
    return DataFrame::create(_["ref1"] = cornerRefs1R, _["ref2"] = cornerRefs2R, _["corner"] = cornerTypesR);
}

bool isInStrictCorner(int bin1, int bin2, int center1, int center2, int size1, int size2, int metaSize, int maxDistance) {
    return ((getCornerDistanceCpp(bin1, bin2, center1, center2, metaSize) <= maxDistance) && (std::abs((bin1 / metaSize) - (center1 / metaSize)) < size1 / 2) && (std::abs((bin2 / metaSize) - (center2 / metaSize)) < size2 / 2));
}

CornerType classifyCornerPoint(int bin1, int bin2, int size1, int size2, int metaSize, int maxDistance) {
    if (isInStrictCorner(bin1, bin2, 0, 0, size1, size2, metaSize, maxDistance)) {
        return BB;
    }
    if (isInStrictCorner(bin1, bin2, size1, 0, size1, size2, metaSize, maxDistance)) {
        return EB;
    }
    if (isInStrictCorner(bin1, bin2, 0, size2, size1, size2, metaSize, maxDistance)) {
        return BE;
    }
    if (isInStrictCorner(bin1, bin2, size1, size2, size1, size2, metaSize, maxDistance)) {
        return EE;
    }
    return interior;
}

int getCornerDistanceCpp(int bin1, int bin2, CornerType ct, int size1, int size2, int metaSize) {
    if (ct == BB) {
        return getCornerDistanceCpp(bin1, bin2, 0, 0, metaSize);
    }
    if (ct == BE) {
        return getCornerDistanceCpp(bin1, bin2, 0, size2, metaSize);
    }
    if (ct == EB) {
        return getCornerDistanceCpp(bin1, bin2, size1, 0, metaSize);
    }
    if (ct == EE) {
        return getCornerDistanceCpp(bin1, bin2, size1, size2, metaSize);
    }
    return 0;
}

// Sum corners
//   output the corners.
// [[Rcpp::export]]                                                                                                                                                                                                
DataFrame sumCornerCpp (DataFrame interactions, DataFrame outliers, IntegerVector sizesIn, int cornerSize, int metaSize) {
    IntegerVector sizes  = clone(sizesIn);
    for (int i = 0; i < sizes.size(); ++i) {
        sizes[i] = sizes[i] / metaSize;
    }
    IntegerVector refs1  = interactions["ref1"];
    IntegerVector refs2  = interactions["ref2"];
    IntegerVector bins1  = interactions["bin1"];
    IntegerVector bins2  = interactions["bin2"];
    IntegerVector counts = interactions["count"];
    CharacterVector refs = refs1.attr("levels");
    assert(refs.size() == sizes.size());
    long long int nInteractions = refs1.size();
    int nRefs = sizes.size();
    int nOutput = nRefs * (nRefs - 1) / 2 * nCornerType;
    std::vector < int > outputRefs1       (nOutput);
    std::vector < int > outputRefs2       (nOutput);
    std::vector < int > outputCornerTypes (nOutput);
    std::vector < int > outputSums        (nOutput, 0);
    std::vector < int > outputCounts      (nOutput, 0);
    int cpt = 0;
    IntegerVector outlierRefs = outliers["ref"];
    IntegerVector outlierBins = outliers["bin"];
    std::vector < long int > offsets (nRefs);
    std::vector < bool > outlierBinsBool;
    long int nElements = computeOffsets(sizes, offsets);
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    for (int ref1 = 0; ref1 < nRefs; ++ref1) {
        for (int ref2 = 0; ref2 < ref1; ++ref2) {
            for (int ct = 0; ct < nCornerType; ++ct) {
                outputRefs1[cpt]       = ref1 + 1; // factors start with 1 in R
                outputRefs2[cpt]       = ref2 + 1;
                outputCornerTypes[cpt] = ct   + 1;
                ++cpt;
            }
        }
    }
    assert(cpt == nOutput);
    for (long long int interactionId = 0; interactionId < nInteractions; ++interactionId) {
        int ref1  = refs1[interactionId] - 1;
        int ref2  = refs2[interactionId] - 1;
        int bin1  = bins1[interactionId] / metaSize;
        int bin2  = bins2[interactionId] / metaSize;
        int count = counts[interactionId];
        int size1 = sizes[ref1];
        int size2 = sizes[ref2];
        if (ref1 != ref2) {
            CornerType ct = classifyCornerPoint(bin1, bin2, size1, size2, 1, cornerSize); // We already have metabins
            if (ct < nCornerType) {
                int offset1 = offsets[ref1];
                int offset2 = offsets[ref2];
                assert(offset1 + bin1 < static_cast<int>(outlierBinsBool.size()));
                assert(offset2 + bin2 < static_cast<int>(outlierBinsBool.size()));
                if ((! outlierBinsBool[offset1 + bin1]) && (! outlierBinsBool[offset2 + bin2])) {
                    int index = (ref1 * (ref1 - 1) / 2 + ref2) * nCornerType + ct;
                    assert(outputRefs1[index]       == ref1 + 1);
                    assert(outputRefs2[index]       == ref2 + 1);
                    assert(outputCornerTypes[index] == ct   + 1);
                    outputSums[index] += count;
                    ++outputCounts[index];
                }
            }
        }
    }
    CharacterVector cornerTypesFactor  = {"BB", "EB", "BE", "EE", "interior"};
    IntegerVector   outputCornerTypesR = wrap(outputCornerTypes);
    outputCornerTypesR.attr("class")   = "factor";
    outputCornerTypesR.attr("levels")  = cornerTypesFactor;
    IntegerVector   outputRefs1R = wrap(outputRefs1);
    outputRefs1R.attr("class")   = "factor";
    outputRefs1R.attr("levels")  = refs;
    IntegerVector   outputRefs2R = wrap(outputRefs2);
    outputRefs2R.attr("class")   = "factor";
    outputRefs2R.attr("levels")  = refs;
    return DataFrame::create(_["ref1"] = outputRefs1R, _["ref2"] = outputRefs2R, _["corner"] = outputCornerTypesR, _["count"] = wrap(outputSums), _["n"] = wrap(outputCounts));
}

// Read interaction matrix and a set of pairs of corners
//   extract the corner, and the "rest" (with negative counts)
// [[Rcpp::export]]
DataFrame extractCornersCpp (DataFrame interactions, DataFrame selectedCorners, DataFrame outliers, IntegerVector sizesIn, int minCornerSize, int maxCornerSize, int metaSize) {
    IntegerVector sizes  = clone(sizesIn);
    for (int i = 0; i < sizes.size(); ++i) {
        sizes[i] = sizes[i] / metaSize;
    }
    IntegerVector refs1  = interactions["ref1"];
    IntegerVector refs2  = interactions["ref2"];
    IntegerVector bins1  = interactions["bin1"];
    IntegerVector bins2  = interactions["bin2"];
    IntegerVector counts = interactions["count"];
    CharacterVector refs = refs1.attr("levels");
    long long int nInteractions = refs1.size();
    int nRefs = sizes.size();
    IntegerVector selectedRefs1 = selectedCorners["ref1"];
    IntegerVector selectedRefs2 = selectedCorners["ref2"];
    IntegerVector selectedEnds  = selectedCorners["corner"];
    int nSelectedCorners = selectedRefs1.size();
    std::vector < std::vector < std::vector < CornerType > > > selectedCornersVector (nRefs + 1);
    for (int refId = 0; refId <= nRefs; ++refId) {
        selectedCornersVector[refId] = std::vector < std::vector < CornerType > > (refId + 1);
    }
    IntegerVector outlierRefs = outliers["ref"];
    IntegerVector outlierBins = outliers["bin"];
    std::vector < long int > offsets (nRefs);
    std::vector < bool > outlierBinsBool;
    long int nElements = computeOffsets(sizes, offsets);
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    std::vector < int > outputRefs1;
    std::vector < int > outputRefs2;
    std::vector < int > outputDistances;
    std::vector < int > outputCornerTypes;
    std::vector < int > outputCounts;
    for (int selectedCornerId = 0; selectedCornerId < nSelectedCorners; ++selectedCornerId) {
        selectedCornersVector[selectedRefs1[selectedCornerId]][selectedRefs2[selectedCornerId]].push_back(static_cast < CornerType > (selectedEnds[selectedCornerId] - 1));
    }
    for (long long int interactionId = 0; interactionId < nInteractions; ++interactionId) {
        int ref1  = refs1[interactionId];
        int ref2  = refs2[interactionId];
        int bin1  = bins1[interactionId] / metaSize;
        int bin2  = bins2[interactionId] / metaSize;
        int count = counts[interactionId];
        int size1 = sizes[ref1-1];
        int size2 = sizes[ref2-1];
        assert(bin1 <= size1);
        assert(bin2 <= size2);
        int offset1 = offsets[ref1 - 1];
        int offset2 = offsets[ref2 - 1];
        assert(offset1 + bin1 < static_cast<int>(outlierBinsBool.size()));
        assert(offset2 + bin2 < static_cast<int>(outlierBinsBool.size()));
        if ((! outlierBinsBool[offset1 + bin1]) && (! outlierBinsBool[offset2 + bin2])) {
            for (CornerType corner: selectedCornersVector[ref1][ref2]) {
                int distance = getCornerDistanceCpp(bin1, bin2, corner, size1, size2, 1);        // Bins are already metabins
                if (distance <= minCornerSize) {
                    outputRefs1.push_back(ref1);
                    outputRefs2.push_back(ref2);
                    outputDistances.push_back(distance);
                    outputCornerTypes.push_back(corner);
                    outputCounts.push_back(count);
                }
                else if (distance >= maxCornerSize) {
                    outputRefs1.push_back(ref1);
                    outputRefs2.push_back(ref2);
                    outputDistances.push_back(distance);
                    outputCornerTypes.push_back(corner);
                    outputCounts.push_back(- count);
                }
            }
        }
    }
    CharacterVector cornerTypesFactor  = {"BB", "EB", "BE", "EE", "interior"};
    IntegerVector   outputCornerTypesR = wrap(outputCornerTypes);
    outputCornerTypesR = outputCornerTypesR + 1; // factors start with 1
    outputCornerTypesR.attr("class")   = "factor";
    outputCornerTypesR.attr("levels")  = cornerTypesFactor;
    IntegerVector   outputRefs1R = wrap(outputRefs1);
    outputRefs1R.attr("class")   = "factor";
    outputRefs1R.attr("levels")  = refs;
    IntegerVector   outputRefs2R = wrap(outputRefs2);
    outputRefs2R.attr("class")   = "factor";
    outputRefs2R.attr("levels")  = refs;
    return DataFrame::create(_["ref1"] = outputRefs1R, _["ref2"] = outputRefs2R, _["distance"] = wrap(outputDistances), _["corner"] = outputCornerTypesR, _["count"] = wrap(outputCounts));
}

// Read interaction matrix,
//   classify points into corners,
//   output a (count, corner) tibble.
// [[Rcpp::export]]
DataFrame classifyCornerPointsCpp (DataFrame interactions, int size1, int size2, int metaSize, int maxDistance) {
    IntegerVector bins1  = interactions["bin1"];
    IntegerVector bins2  = interactions["bin2"];
    IntegerVector counts = interactions["count"];
    std::vector < int > countsCpp = as < std::vector < int > > (counts);
    std::vector < int > nCounts (nCornerIntType + 1, 0);
    long long int nInteractions = bins1.size();
    std::vector < int > cornerTypes (nInteractions);
    int cornerSize   = computeCornerSize(size1, size2, maxDistance);
    int interiorSize = size1 * size2 - 4 * cornerSize;
    for (long long int interactionId = 0; interactionId < nInteractions; ++interactionId) {
        int type = classifyCornerPoint(bins1[interactionId], bins2[interactionId], size1, size2, metaSize, maxDistance);
        cornerTypes[interactionId] = type;
	++nCounts[type];
    }
    for (int type = 0; type < nCornerIntType; ++type) {
        int size = (type == nCornerType)? interiorSize: cornerSize;
	int sizeDiff = size - nCounts[type];
	assert(sizeDiff >= 0);
        if (sizeDiff > 0) {
	    std::vector < int > fillerType  (sizeDiff, type);
	    std::vector < int > fillerCount (sizeDiff, 0);
	    cornerTypes.insert(cornerTypes.end(), fillerType.begin(), fillerType.end());
	    countsCpp.insert(countsCpp.end(), fillerCount.begin(), fillerCount.end());
	}
    }
    CharacterVector cornerTypesFactor = {"BB", "EB", "BE", "EE", "interior"};
    IntegerVector   cornerTypesR      = wrap(cornerTypes);
    cornerTypesR                 = cornerTypesR + 1; // factors start with 1
    cornerTypesR.attr("class")   = "factor";
    cornerTypesR.attr("levels")  = cornerTypesFactor;
    counts = wrap(countsCpp);
    return DataFrame::create(_["type"] = cornerTypesR, _["count"] = counts);
}

int cornerTypeBoolToBin (bool after1, bool after2) {
    int cornerType = 0;
    if (after1) {
        cornerType = 1;
    }
    if (after2) {
        cornerType += 2;
    }
    return cornerType;
}

// Read interaction matrix, a set of pairs of references, and a set of corners
//   output the full corners: the full list of distance/count per corner
// [[Rcpp::export]]                                                                                                                                                                                                
DataFrame extractCornersFullCpp (DataFrame interactions, DataFrame selectedCorners, DataFrame outliers, IntegerVector sizesIn, int cornerSize, int metaSize) {
    IntegerVector sizes  = clone(sizesIn);
    for (int i = 0; i < sizes.size(); ++i) {
        sizes[i] = sizes[i] / metaSize;
    }
    IntegerVector refs1  = interactions["ref1"];
    IntegerVector refs2  = interactions["ref2"];
    IntegerVector bins1  = interactions["bin1"];
    IntegerVector bins2  = interactions["bin2"];
    IntegerVector counts = interactions["count"];
    CharacterVector refs = refs1.attr("levels");
    long long int nInteractions = refs1.size();
    int nRefs = sizes.size();
    IntegerVector selectedRefs1   = selectedCorners["ref1"];
    IntegerVector selectedRefs2   = selectedCorners["ref2"];
    LogicalVector selectedAfters1 = selectedCorners["after1"];
    LogicalVector selectedAfters2 = selectedCorners["after2"];
    int nSelectedCorners = selectedRefs1.size();
    IntegerVector outlierRefs = outliers["ref"];
    IntegerVector outlierBins = outliers["bin"];
    std::vector < std::vector < std::array < int, nCornerType > > > selectedCornersBool (nRefs + 1);
    for (int refId = 0; refId <= nRefs; ++refId) {
        selectedCornersBool[refId] = std::vector < std::array < int, nCornerType > > (refId + 1, std::array < int, nCornerType > ( { -1, -1, -1, -1 } ));
    }
    std::vector < long int > offsets (nRefs);
    std::vector < bool > outlierBinsBool;
    long int nElements = computeOffsets(sizes, offsets);
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    int multiplier = (cornerSize + 1) * (cornerSize + 2) / 2;
    int nOutputElements = nSelectedCorners * multiplier;
    std::vector < int >  outputIndex     (nOutputElements);
    std::vector < int >  outputRefs1     (nOutputElements);
    std::vector < int >  outputRefs2     (nOutputElements);
    std::vector < bool > outputAfters1   (nOutputElements);
    std::vector < bool > outputAfters2   (nOutputElements);
    std::vector < int >  outputDistances (nOutputElements);
    std::vector < int >  outputCounts    (nOutputElements, 0);
    for (int selectedRefId = 0; selectedRefId < nSelectedCorners; ++selectedRefId) {
        selectedCornersBool[selectedRefs1[selectedRefId]][selectedRefs2[selectedRefId]][cornerTypeBoolToBin(selectedAfters1[selectedRefId], selectedAfters2[selectedRefId])] = selectedRefId;
    }
    // Prepare output columns
    long int cpt = 0;
    for (long int i = 0; i < nSelectedCorners; ++i) {
        for (long int j = 0; j <= cornerSize; ++j) {
            for (long int k = 0; k <= j; ++k) {
                outputIndex[cpt]     = i;
                outputRefs1[cpt]     = selectedRefs1[i];
                outputRefs2[cpt]     = selectedRefs2[i];
                outputAfters1[cpt]   = selectedAfters1[i];
                outputAfters2[cpt]   = selectedAfters2[i];
                outputDistances[cpt] = j;
                ++cpt;
            }
        }
    }
    // Set NAs for outliers
    for (long int i = 0; i < nSelectedCorners; ++i) {
        int  ref1    = selectedRefs1[i];
        int  ref2    = selectedRefs2[i];
        bool after1  = selectedAfters1[i];
        bool after2  = selectedAfters2[i];
        int  size1   = sizes[ref1 - 1];
        int  size2   = sizes[ref2 - 1];
        int  offset1 = offsets[ref1 - 1];
        int  offset2 = offsets[ref2 - 1];
        for (int b1 = 0; b1 <= std::min<int>(cornerSize, size1); ++b1) {
            for (int b2 = 0; (b2 <= size2) && (b1 + b2 <= cornerSize); ++b2) {
                int distance = b1 + b2;
                int bin1     = (after1)? size1 - b1: b1;
                int bin2     = (after2)? size2 - b2: b2;
                assert(offset1 + bin1 < static_cast<int>(outlierBinsBool.size()));
                assert(offset2 + bin2 < static_cast<int>(outlierBinsBool.size()));
                if (outlierBinsBool[offset1 + bin1] || outlierBinsBool[offset2 + bin2]) {
                    int index = i * multiplier + distance * (distance + 1) / 2 + b1;
                    assert(index < static_cast<int>(outputCounts.size()));
                    outputCounts[index] = NA_INTEGER;
                }
            }
        }
    }
    // Set to -1 cells that extend the range of the corner in small scaffolds
    for (long int index = 0; index < nSelectedCorners; ++index) {
        for (int bin1 = sizes[selectedRefs1[index] - 1] / 2; bin1 <= cornerSize; ++bin1) {
            for (int bin2 = 0; (bin2 <= bin1) && (bin1 + bin2 <= cornerSize); ++bin2) {
                int distance = bin1 + bin2;
                int cellIndex = index * multiplier + distance * (distance + 1) / 2 + bin1;
                assert(cellIndex < static_cast<int>(outputCounts.size()));
                outputCounts[cellIndex] = -1;
            }
        }
        for (int bin2 = sizes[selectedRefs2[index] - 1] / 2; bin2 <= cornerSize; ++bin2) {
            for (int bin1 = bin2; bin1 + bin2 <= cornerSize; ++bin1) {
                int distance = bin1 + bin2;
                int cellIndex = index * multiplier + distance * (distance + 1) / 2 + bin1;
                outputCounts[cellIndex] = -1;
            }
        }
    }
    for (long long int interactionId = 0; interactionId < nInteractions; ++interactionId) {
        int ref1  = refs1[interactionId];
        int ref2  = refs2[interactionId];
        int bin1  = bins1[interactionId] / metaSize;
        int bin2  = bins2[interactionId] / metaSize;
        int count = counts[interactionId];
        int size1 = sizes[ref1-1];
        int size2 = sizes[ref2-1];
        if ((! outlierBinsBool[offsets[ref1-1] + bin1]) && (! outlierBinsBool[offsets[ref2-1] + bin2])) {
            CornerType ct = classifyCornerPoint(bin1, bin2, size1, size2, 1, cornerSize); // we already have meta-bins
            if (ct < 4) {
                int index = selectedCornersBool[ref1][ref2][ct];
                if (index >= 0) {
                    int distance = getCornerDistanceCpp(bin1, bin2, ct, size1, size2, 1); // we already have meta-bins
                    if ((ct & 1) == 1) {
                        bin1 = size1 - bin1;
                    }
                    index = index * multiplier + distance * (distance + 1) / 2 + bin1;
                    outputCounts[index] += count;
                }
            }
        }
    }
    IntegerVector outputIndexR     = wrap(outputIndex);
    IntegerVector outputRefs1R     = wrap(outputRefs1);
    IntegerVector outputRefs2R     = wrap(outputRefs2);
    LogicalVector outputAfters1R   = wrap(outputAfters1);
    LogicalVector outputAfters2R   = wrap(outputAfters2);
    IntegerVector outputDistancesR = wrap(outputDistances);
    IntegerVector outputCountsR    = wrap(outputCounts);
    outputRefs1R.attr("class")  = "factor";
    outputRefs1R.attr("levels") = refs;
    outputRefs2R.attr("class")  = "factor";
    outputRefs2R.attr("levels") = refs;
    return DataFrame::create(_["index"]    = outputIndexR,
                             _["ref1"]     = outputRefs1R,
                             _["ref2"]     = outputRefs2R,
                             _["after1"]   = outputAfters1R,
                             _["after2"]   = outputAfters2R,
                             _["distance"] = outputDistancesR,
                             _["count"]    = outputCountsR);
}
