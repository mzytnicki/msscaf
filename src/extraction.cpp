#include <vector>
#include <Rcpp.h>
#include "sharedFunctions.h"
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame keepScaffoldsCpp (DataFrame &data, CharacterVector keptRefs) {
    IntegerVector   refs1    = data["ref1"];
    IntegerVector   refs2    = data["ref2"];
    IntegerVector   bins1    = data["bin1"];
    IntegerVector   bins2    = data["bin2"];
    IntegerVector   counts   = data["count"];
    CharacterVector refNames = refs1.attr("levels");
    int             nRefs    = refNames.size();
    int             nKept    = keptRefs.size();
    long int        nRowsIn  = refs1.size();
    long int        nRowsOut = 0;
    std::vector < bool > keptIds (nRefs + 1, false); // Levels start with 1 in R
    // Recompute factors for the refs
    std::vector < int > translateRefIds (nRefs + 1);
    int refId = 1;                     // factors start with 1
    for (int i = 0; i < nRefs; ++i) {
        for (int j = 0; j < nKept; ++j) {
            std::string refName = as< std::string >(refNames[i]);
            std::string keptRef = as< std::string >(keptRefs[j]);
            if (caseInsensitiveEqual(refName, keptRef)) { // HiC data may change refs to upper case
                keptIds[i+1] = true; // Levels start with 1 in R
                translateRefIds[i+1] = refId;
                ++refId;
                break;
            }
        }
    }
    if (refId - 1 != nKept) {
        Rcerr << "Problem during transfer: I found " << (refId-1) << " refs, but was told to keep " << nKept << ".\n";
        stop("Stopping here.");
    }
    // In-place update of the data
    for (long int i = 0; i < nRowsIn; ++i) {
        if (keptIds[refs1[i]] && keptIds[refs2[i]]) {
            //Rcerr << i << ", " << refs1[i] << " is now " << nRowsOut << ", " << translateRefIds[refs1[i]] << ".\n";
            refs1[nRowsOut]  = translateRefIds[refs1[i]];
            refs2[nRowsOut]  = translateRefIds[refs2[i]];
            if (nRowsOut != i) {
                bins1[nRowsOut]  = bins1[i];
                bins2[nRowsOut]  = bins2[i];
                counts[nRowsOut] = counts[i];
            }
            ++nRowsOut;
        }
    }
    refs1.erase(nRowsOut, nRowsIn);
    refs2.erase(nRowsOut, nRowsIn);
    bins1.erase(nRowsOut, nRowsIn);
    bins2.erase(nRowsOut, nRowsIn);
    counts.erase(nRowsOut, nRowsIn);
    if (refs1.size() != nRowsOut) {
        Rcerr << "Problem while resizing vectors.\n";
        stop("Stopping here.");
    }
    refs1.attr("class") = "factor";
    refs2.attr("class") = "factor";
    refs1.attr("levels") = keptRefs;
    refs2.attr("levels") = keptRefs;
    return Rcpp::DataFrame::create(_["ref1"]= refs1, _["bin1"]= bins1, _["ref2"]= refs2, _["bin2"]= bins2, _["count"]= counts);
}

// Extract matrix lines such that (ref1, ref2) match given data.
// [[Rcpp::export]]
DataFrame keepScaffoldsPairsCpp (DataFrame &data, DataFrame &keptRefs) {
    IntegerVector   refs1     = data["ref1"];
    IntegerVector   refs2     = data["ref2"];
    IntegerVector   bins1     = data["bin1"];
    IntegerVector   bins2     = data["bin2"];
    IntegerVector   counts    = data["count"];
    CharacterVector refNames  = refs1.attr("levels");
    int             nRefs     = refNames.size();
    IntegerVector   keptRefs1 = keptRefs["ref1"];
    IntegerVector   keptRefs2 = keptRefs["ref2"];
    int             nKeptRefs = keptRefs1.size();
    std::vector < int > newRefs1;
    std::vector < int > newRefs2;
    std::vector < int > newBins1;
    std::vector < int > newBins2;
    std::vector < int > newCounts;
    long long int   nDataIn   = refs1.size();
    long long int   nDataOut  = 0;
    std::vector < std::vector < bool > > keptRefsBool (nRefs + 1);
    for (int refId = 0; refId <= nRefs; ++refId) {
        keptRefsBool[refId] = std::vector < bool > (refId + 1, false);
    }
    for (int keptRefId = 0; keptRefId < nKeptRefs; ++keptRefId) {
        if (keptRefs1[keptRefId] >= static_cast < int > (keptRefsBool.size())) {
            Rcerr << "Problem #1 in keepScaffoldsPairsCpp: " << keptRefs1[keptRefId] << " >= " << keptRefsBool.size() << ".\n";
        }
        if (keptRefs2[keptRefId] >= static_cast < int > (keptRefsBool[keptRefs1[keptRefId]].size())) {
            Rcerr << "Problem #2 in keepScaffoldsPairsCpp: " << keptRefs2[keptRefId] << " >= " << keptRefsBool[keptRefs1[keptRefId]].size() << ".\n";
        }
        keptRefsBool[keptRefs1[keptRefId]][keptRefs2[keptRefId]] = true;
    }
    for (long long int dataId = 0; dataId < nDataIn; ++dataId) {
        if (keptRefsBool[refs1[dataId]][refs2[dataId]]) {
	    ++nDataOut;
        }
    }
    newRefs1.reserve(nDataOut);
    newRefs2.reserve(nDataOut);
    newBins1.reserve(nDataOut);
    newBins2.reserve(nDataOut);
    newCounts.reserve(nDataOut);
    for (long long int dataId = 0; dataId < nDataIn; ++dataId) {
        if (keptRefsBool[refs1[dataId]][refs2[dataId]]) {
            newRefs1.push_back(refs1[dataId]);
            newRefs2.push_back(refs2[dataId]);
            newBins1.push_back(bins1[dataId]);
            newBins2.push_back(bins2[dataId]);
            newCounts.push_back(counts[dataId]);
        }
    }
    IntegerVector newRefs1R = wrap(newRefs1);
    IntegerVector newRefs2R = wrap(newRefs2);
    IntegerVector newBins1R = wrap(newBins1);
    IntegerVector newBins2R = wrap(newBins2);
    IntegerVector newCountsR = wrap(newCounts);
    newRefs1R.attr("class") = "factor";
    newRefs2R.attr("class") = "factor";
    newRefs1R.attr("levels") = refNames;
    newRefs2R.attr("levels") = refNames;
    return Rcpp::DataFrame::create(_["ref1"] = newRefs1R, _["bin1"] = newBins1R, _["ref2"] = newRefs1R, _["bin2"] = newBins2R, _["count"] = newCountsR);
}

uint64_t combineInts(int i1, int i2) {
    if (i1 < 0) stop("Problem in 'combineInts': first data is " + std::to_string(i1) + ".");
    if (i2 < 0) stop("Problem in 'combineInts': second data is " + std::to_string(i2) + ".");
    return (static_cast < uint64_t > (i1) << 32) | static_cast < uint64_t > (i2);
}

// Extract matrix lines such that (ref, bin) match given data.
// Convert them into (ref, bin, distance, count)
// [[Rcpp::export]]
DataFrame extractLines (DataFrame matrix, DataFrame lines, int maxDistance) {
    IntegerVector   refs1         = matrix["ref1"];
    IntegerVector   refs2         = matrix["ref2"];
    IntegerVector   bins1         = matrix["bin1"];
    IntegerVector   bins2         = matrix["bin2"];
    IntegerVector   counts        = matrix["count"];
    CharacterVector refNames      = refs1.attr("levels");
    IntegerVector   selectedRefs  = lines["ref"];
    IntegerVector   selectedBins  = lines["bin"];
    long int        nMatrixLines  = refs1.size();
    long int        nLines        = selectedRefs.size();
    long int        nOutputLines  = nLines * (maxDistance + 1);
    std::vector < std::vector < int > > lineDistance (nLines, std::vector < int > (maxDistance + 1, 0));
    std::unordered_map < uint64_t, size_t > keptLines;
    IntegerVector   outputRefs      (nOutputLines, 0);
    IntegerVector   outputBins      (nOutputLines, 0);
    IntegerVector   outputDistances (nOutputLines, 0);
    IntegerVector   outputCounts    (nOutputLines, 0);
    // Store lines
    if (is_true(any(selectedBins < 0))) stop("Error in 'extractLines', got negative bins.");
    keptLines.reserve(nLines);
    for (long int lineId = 0; lineId < nLines; ++lineId) {
        keptLines[combineInts(selectedRefs[lineId], selectedBins[lineId])] = lineId;
    }
    // Read matrix
    for (long int i = 0; i < nMatrixLines; ++i) {
        // Select matching lines
        if (refs1[i] == refs2[i]) {
            int distance = bins1[i] - bins2[i];
            if (distance < 0) stop("Error in 'extractLines', got negative distances: ref: " + std::to_string(refs1[i]) + ", bin1: " + std::to_string(bins1[i]) + ", bin2: " + std::to_string(bins2[i]) + ".");
            if (distance <= maxDistance) {
                std::unordered_map < uint64_t, size_t >::iterator it = keptLines.find(combineInts(refs1[i], bins1[i]));
                if (it != keptLines.end()) {
                    lineDistance[it->second][distance] = counts[i];
                }
            }
        }
    }
    // Create output data frame
    unsigned int i = 0;
    for (long int lineId = 0; lineId < nLines; ++lineId) {
        for (int distance = 0; distance <= maxDistance; ++distance) {
            outputRefs[i]      = selectedRefs[lineId];
            outputBins[i]      = selectedBins[lineId];
            outputDistances[i] = distance;
            outputCounts[i]    = lineDistance[lineId][distance];
            ++i;
        }
    }
    outputRefs.attr("class") = "factor";
    outputRefs.attr("levels") = refNames;
    return Rcpp::DataFrame::create(_["ref"] = outputRefs, _["bin"] = outputBins, _["distance"] = outputDistances, _["count"] = outputCounts);
}
