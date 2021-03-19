#include <vector>
#include <string>
#include <algorithm>
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector computeRefSizesCpp (DataFrame &data) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    CharacterVector refNames = refs1.attr("levels");
    std::vector < int > sizes(refNames.size());
    long long int nElements = refs1.size();
    for (long long int i = 0; i < nElements; ++i) {
        sizes[refs1[i]] = std::max<int>(sizes[refs1[i]], bins1[i]);
        sizes[refs2[i]] = std::max<int>(sizes[refs2[i]], bins2[i]);
    }
    IntegerVector sizesR;
    for (size_t i = 0; i < sizes.size(); ++i) {
        sizesR.push_back(sizes[i], as < std::string > (refNames[i]));
    }
    return sizesR;
}

// [[Rcpp::export]]
int computeNRrows (DataFrame &data, IntegerVector &sizes) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    int nElements = 0;
    std::vector < int > offsets (sizes.size());
    std::vector < bool > seen;
    int nRows = 0;
    for (int i = 0; i < sizes.size(); ++i) {
        offsets[i] = nElements;
        nElements += sizes[i] + 1;
    }
    seen.resize(nElements, false);
    for (int i = 0; i < refs1.size(); ++i) {
        int bin1 = offsets[refs1[i] - 1] + bins1[i]; // factors (such as 'refs') start with 1 in R
        if (! seen[bin1]) {
          ++nRows;
          seen[bin1] = true;
        }
        int bin2 = offsets[refs2[i] - 1] + bins2[i]; // factors (such as 'refs') start with 1 in R
        if (! seen[bin2]) {
          ++nRows;
          seen[bin2] = true;
        }
    }
    return nRows;
}

int computeOffsets (DataFrame &data, IntegerVector &sizes, std::vector < int > &offsets) {
    int nElements = 0;
    for (int i = 0; i < sizes.size(); ++i) {
        offsets[i] = nElements;
        nElements += sizes[i] + 1;
    }
    return nElements;
}

IntegerVector computeSymmetricColSum (DataFrame &data, IntegerVector &sizes, std::vector < int > &offsets, int nElements) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    std::vector < int > sums; // Col sums of all refs will be concatenated here.
    sums.resize(nElements, 0);
    for (int i = 0; i < refs1.size(); ++i) {
        sums[offsets[refs1[i] - 1] + bins1[i]] += counts[i]; // factors (such as 'refs') start with 1 in R
        sums[offsets[refs2[i] - 1] + bins2[i]] += counts[i];
    }
    return wrap(sums);
}

// [[Rcpp::export]]
IntegerVector computeSymmetricColSum (DataFrame &data, IntegerVector &sizes) {
    std::vector < int > offsets (sizes.size());
    int nElements = computeOffsets(data, sizes, offsets);
    return computeSymmetricColSum(data, sizes, offsets, nElements);
}

// [[Rcpp::export]]
DataFrame removeLowCountRowsCpp (DataFrame &data, IntegerVector &sizes, int threshold) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    CharacterVector refNames = refs1.attr("levels");
    std::vector < int > offsets (sizes.size());
    int nElements = computeOffsets(data, sizes, offsets);
    IntegerVector colSums = computeSymmetricColSum(data, sizes, offsets, nElements);
    std::vector < bool > passedThreshold (colSums.size());
    for (int i = 0; i < colSums.size(); ++i) {
        passedThreshold[i] = (colSums[i] >= threshold);
    }
    int nPassed = std::count(passedThreshold.begin(), passedThreshold.end(), true);
    Rcerr << "\tKeeping " << nPassed << "/" << nElements << " rows.\n";
    std::vector < int > refs1C;
    std::vector < int > refs2C;
    std::vector < int > bins1C;
    std::vector < int > bins2C;
    std::vector < int > countsC;
    for (int i = 0; i < refs1.size(); ++i) {
        int bin1 = offsets[refs1[i] - 1] + bins1[i]; // factors (such as 'refs') start with 1 in R
        int bin2 = offsets[refs2[i] - 1] + bins2[i]; // factors (such as 'refs') start with 1 in R
        if (passedThreshold[bin1] && passedThreshold[bin2]) {
            refs1C.push_back(refs1[i]);
            refs2C.push_back(refs2[i]);
            bins1C.push_back(bins1[i]);
            bins2C.push_back(bins2[i]);
            countsC.push_back(counts[i]);
        }
    }
    refs1  = wrap(refs1C);
    refs2  = wrap(refs2C);
    bins1  = wrap(bins1C);
    bins2  = wrap(bins2C);
    counts = wrap(countsC);
    refs1.attr("class") = "factor";
    refs2.attr("class") = "factor";
    refs1.attr("levels") = refNames;
    refs2.attr("levels") = refNames;
    DataFrame outputData = Rcpp::DataFrame::create(_["ref1"]= refs1, _["bin1"]= bins1, _["ref2"]= refs2, _["bin2"]= bins2, _["count"]= counts);
    return outputData;
}

bool caseInsensitiveEqual(std::string &s1, std::string &s2) {
    if (s1.size() != s2.size()) {
        return false;
    }
    for (size_t i = 0; i < s1.size(); ++i) {
        if (tolower(s1[i]) != tolower(s2[i])) {
            return false;
        }
    }
    return true;
}

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


// [[Rcpp::export]]
DataFrame computeMeanTrianglesCpp (DataFrame &data) {
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    long int      size   = std::max < int > (max(bins1), max(bins2)) + 1;
    long int      n      = bins1.size();
    std::vector < int > sumBin1 (size);
    std::vector < int > sumBin2 (size);
    std::vector < int > nBin1   (size);
    std::vector < int > nBin2   (size);
    IntegerVector bin           (size);
    NumericVector meanCount     (size);
    IntegerVector nCells        (size);
    for (long int i = 0; i < n; ++i) {
        sumBin1[bins1[i]] += counts[i];
        sumBin2[bins2[i]] += counts[i];
        ++nBin1[bins1[i]];
        ++nBin2[bins2[i]];
    }
    meanCount[0] = sumBin2[0];
    nCells[0]    = nBin2[0];
    for (long int i = 1; i < size; ++i) {
        meanCount[i] = meanCount[i-1] + sumBin2[i] - sumBin1[i-1];
        nCells[i]     = nCells[i-1]   + nBin2[i]   - nBin1[i-1];
    }
    for (long int i = 0; i < size; ++i) {
        bin[i] = i;
        if (nCells[i] > 0) {
            meanCount[i] /= nCells[i];
        }
    }
    return Rcpp::DataFrame::create(_["bin"] = bin, _["meanCount"] = meanCount, _["nCells"] = nCells);
}
