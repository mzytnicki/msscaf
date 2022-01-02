#include <cstdlib>
#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
                                                                                                                                                                                                                   
// [[Rcpp::plugins(cpp11)]]
                                                                                                                                                                                                                   
using namespace Rcpp;

enum CornerType { BB, EB, BE, EE, interior };
const int nCornerType { 4 };
const int nCornerIntType { 5 };

bool isInCorner(int bin1, int bin2, int center1, int center2, int maxDistance) {
    return (std::abs(bin1 - center1) + std::abs(bin2 - center2) < maxDistance);
}


// Compute corner size
// [[Rcpp::export]]
int computeCornerSize(int size1, int size2, int maxDistance) {
    int cornerSize = maxDistance * (maxDistance + 1) / 2;
    size1 = maxDistance - (size1 / 2);
    size2 = maxDistance - (size2 / 2);
    if (size1 > 0) {
        cornerSize = cornerSize - (size1 * (size1 + 1) / 2);
    }
    if (size2 > 0) {
        cornerSize = cornerSize - (size2 * (size2 + 1) / 2);
    }
    return cornerSize;
}


// Read interaction matrix,
//   stores the number of non-zero for each corner,
//   output the non-near-empty corners.
// [[Rcpp::export]]
DataFrame filterCornersCpp (DataFrame data, IntegerVector sizes, int cornerSize) {
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
            if (ref1 > nRefs) {
                Rcerr << "Error #1 in filter corner: ref #" << ref1 << " >= " << nRefs << " refs.\n";
                stop("Aborting.");
            }
            if (ref2 >= ref1) {
                Rcerr << "Error #2 in filter corner: first ref #" << ref1 << " <= second ref " << ref2 << ".\n";
                stop("Aborting.");
            }
            if (isInCorner(bin1, bin2, 0, 0, cornerSize)) {
                ++cornerCounts[ref1][ref2][BB];
            }
            if (isInCorner(bin1, bin2, size1, 0, cornerSize)) {
                ++cornerCounts[ref1][ref2][EB];
            }
            if (isInCorner(bin1, bin2, 0, size2, cornerSize)) {
                ++cornerCounts[ref1][ref2][BE];
            }
            if (isInCorner(bin1, bin2, size1, size2, cornerSize)) {
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


bool isInStrictCorner(int bin1, int bin2, int center1, int center2, int size1, int size2, int maxDistance) {
    return ((std::abs(bin1 - center1) + std::abs(bin2 - center2) < maxDistance) && (std::abs(bin1 - center1) < size1 / 2) && (std::abs(bin2 - center2) < size2 / 2));
}

CornerType classifyCornerPoint(int bin1, int bin2, int size1, int size2, int maxDistance) {
    if (isInStrictCorner(bin1, bin2, 0, 0, size1, size2, maxDistance)) {
        return BB;
    }
    if (isInStrictCorner(bin1, bin2, size1, 0, size1, size2, maxDistance)) {
        return EB;
    }
    if (isInStrictCorner(bin1, bin2, 0, size2, size1, size2, maxDistance)) {
        return BE;
    }
    if (isInStrictCorner(bin1, bin2, size1, size2, size1, size2, maxDistance)) {
        return EE;
    }
    return interior;
}

int getCornerDistanceCpp(int bin1, int bin2, int center1, int center2) {
    return std::abs(bin1 - center1) + std::abs(bin2 - center2);
}

int getCornerDistanceCpp(int bin1, int bin2, CornerType ct, int size1, int size2) {
    if (ct == BB) {
        return getCornerDistanceCpp(bin1, bin2, 0, 0);
    }
    if (ct == BE) {
        return getCornerDistanceCpp(bin1, bin2, 0, size2);
    }
    if (ct == EB) {
        return getCornerDistanceCpp(bin1, bin2, size1, 0);
    }
    if (ct == EE) {
        return getCornerDistanceCpp(bin1, bin2, size1, size2);
    }
    return 0;
}

// Read interaction matrix and a set of pairs of references
//   output the corners.
// [[Rcpp::export]]                                                                                                                                                                                                
DataFrame extractCornersCpp (DataFrame interactions, DataFrame selectedRefs, IntegerVector sizes, int cornerSize) {
    IntegerVector refs1  = interactions["ref1"];
    IntegerVector refs2  = interactions["ref2"];
    IntegerVector bins1  = interactions["bin1"];
    IntegerVector bins2  = interactions["bin2"];
    IntegerVector counts = interactions["count"];
    CharacterVector refs = refs1.attr("levels");
    long long int nInteractions = refs1.size();
    int nRefs = sizes.size();
    IntegerVector selectedRefs1 = selectedRefs["ref1"];
    IntegerVector selectedRefs2 = selectedRefs["ref2"];
    int nSelectedRefs = selectedRefs1.size();
    std::vector < std::vector < bool > > selectedRefsBool (nRefs + 1);
    for (int refId = 0; refId <= nRefs; ++refId) {
        selectedRefsBool[refId] = std::vector < bool > (refId + 1, false);
    }
    std::vector < int > outputRefs1;
    std::vector < int > outputRefs2;
    std::vector < int > outputDistances;
    std::vector < int > outputCornerTypes;
    std::vector < int > outputCounts;
    for (int selectedRefId = 0; selectedRefId < nSelectedRefs; ++selectedRefId) {
        selectedRefsBool[selectedRefs1[selectedRefId]][selectedRefs2[selectedRefId]] = true;
    }
    for (long long int interactionId = 0; interactionId < nInteractions; ++interactionId) {
        int ref1  = refs1[interactionId];
        int ref2  = refs2[interactionId];
        int bin1  = bins1[interactionId];
        int bin2  = bins2[interactionId];
        int count = counts[interactionId];
        int size1 = sizes[ref1-1];
        int size2 = sizes[ref2-1];
        if (selectedRefsBool[ref1][ref2]) {
            CornerType ct = classifyCornerPoint(bin1, bin2, size1, size2, cornerSize);
            int        distance = getCornerDistanceCpp(bin1, bin2, ct, size1, size2);
            outputRefs1.push_back(ref1);
            outputRefs2.push_back(ref2);
            outputDistances.push_back(distance);
            outputCornerTypes.push_back(ct);
            outputCounts.push_back(count);
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
    return DataFrame::create(_["ref1"] = outputRefs1R, _["ref2"] = outputRefs2R, _["distance"] = wrap(outputDistances), _["type"] = outputCornerTypesR, _["count"] = wrap(outputCounts));
}

// Read interaction matrix,
//   classify points into corners,
//   output a (count, corner) tibble.
// [[Rcpp::export]]
DataFrame classifyCornerPointsCpp (DataFrame interactions, int size1, int size2, int maxDistance) {
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
        int type = classifyCornerPoint(bins1[interactionId], bins2[interactionId], size1, size2, maxDistance);
        cornerTypes[interactionId] = type;
	++nCounts[type];
    }
    for (int type = 0; type < nCornerIntType; ++type) {
        int size = (type == nCornerType)? interiorSize: cornerSize;
	int sizeDiff = size - nCounts[type];
	if (sizeDiff < 0) Rcerr << "Error while counting corner sizes: " << size << " vs " << nCounts[type] << "\n";
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
