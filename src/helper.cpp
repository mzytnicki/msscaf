#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <Rcpp.h>
#include "sharedFunctions.h"
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


// [[Rcpp::export]]
void splitChromosomeCpp (DataFrame &data, int prevRef, int newRef, int shiftedRef, long splitPoint, bool firstPart) {
    auto f = (firstPart)? std::function<bool(long)>([splitPoint] (long x) { return (x >= splitPoint); }): std::function<bool(long)>([splitPoint] (long x) { return (x <= splitPoint); });
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    long long int nElements = refs1.size();
    long offset = splitPoint - 1;
    for (long long int i = 0; i < nElements; ++i) {
        if ((refs1[i] == prevRef) && (f(bins1[i]))) {
            refs1[i] = newRef;
        }
        if ((refs2[i] == prevRef) && (f(bins2[i]))) {
            refs2[i] = newRef;
        }
        if (refs1[i] == shiftedRef) {
            bins1[i] -= offset;
        }
        if (refs2[i] == shiftedRef) {
            bins2[i] -= offset;
        }
        if ((refs1[i] == refs2[i]) && (bins1[i] < bins2[i])) stop("Error in 'splitChromosomeCpp', got negative distances. Prev ref: " + std::to_string(prevRef) + ", new ref: " + std::to_string(newRef) + ", shifted ref: " + std::to_string(shiftedRef) + ", split point: " + std::to_string(splitPoint) + ", first part: " + std::to_string(firstPart) + ".");
    }
}


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

void computeSymmetricColSum (DataFrame &data, IntegerVector &sizes, std::vector < long int > &offsets, int nElements, std::vector < int > &sums) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    sums = std::vector < int > (nElements, 0); // Col sums of all refs will be concatenated here.
    for (int i = 0; i < refs1.size(); ++i) {
        sums[offsets[refs1[i] - 1] + bins1[i]] += counts[i]; // factors (such as 'refs') start with 1 in R
        sums[offsets[refs2[i] - 1] + bins2[i]] += counts[i];
    }
}

void computeSymmetricColSumMeta (DataFrame &data, IntegerVector &sizes, std::vector < long int > &offsets, int nElements, int metaSize, std::vector < int > &sums) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    sums = std::vector < int > (nElements, 0); // Col sums of all refs will be concatenated here.
    for (int i = 0; i < refs1.size(); ++i) {
        sums[offsets[refs1[i] - 1] + (bins1[i] / metaSize)] += counts[i]; // factors (such as 'refs') start with 1 in R
        sums[offsets[refs2[i] - 1] + (bins2[i] / metaSize)] += counts[i];
    }
}

/*
DataFrame computeSymmetricColSum (DataFrame &data, IntegerVector &sizes, std::vector < int > &offsets, int nElements) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    std::vector < int > sums (nElements, 0); // Col sums of all refs will be concatenated here.
    std::vector < int > refs (nElements); // Output reference ids
    std::vector < int > bins (nElements); // Output reference bins
    int cpt = 0;
    for (int ref = 0; ref < sizes.size(); ++ref) {
        for (int bin = 0; bin < sizes[ref]; ++bin) {
            refs[cpt] = ref + 1;
            bins[cpt] = bin;
            ++cpt;
        }
    }
    for (int i = 0; i < refs1.size(); ++i) {
        sums[offsets[refs1[i] - 1] + bins1[i]] += counts[i]; // factors (such as 'refs') start with 1 in R
        sums[offsets[refs2[i] - 1] + bins2[i]] += counts[i];
    }
    return Rcpp::DataFrame::create(_["ref"] = wrap(refs), _["bin"] = wrap(bins), _["sum"] = wrap(sums));
}
*/

// [[Rcpp::export]]
DataFrame computeSymmetricColSum (DataFrame &data, IntegerVector &sizes) {
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsets(sizes, offsets);
    std::vector < int > sums;
    IntegerVector refs; // Output reference ids
    IntegerVector bins; // Output reference bins
    computeRefsBins(sizes, nElements, refs, bins);
    computeSymmetricColSum(data, sizes, offsets, nElements, sums);
    return Rcpp::DataFrame::create(_["ref"] = refs, _["bin"] = bins, _["sum"] = wrap(sums));
}

// [[Rcpp::export]]
DataFrame computeSymmetricColSumMeta (DataFrame &data, IntegerVector &sizes, int metaSize) {
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsetsMeta(sizes, metaSize, offsets);
    std::vector < int > sums;
    IntegerVector refs; // Output reference ids
    IntegerVector bins; // Output reference bins
    computeRefsBinsMeta(sizes, nElements, metaSize, true, refs, bins);
    computeSymmetricColSumMeta(data, sizes, offsets, nElements, metaSize, sums);
    return Rcpp::DataFrame::create(_["ref"] = refs, _["bin"] = bins, _["sum"] = wrap(sums));
}

// [[Rcpp::export]]
List removeLowCountRowsCpp (DataFrame &data, IntegerVector &sizes, int threshold) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    CharacterVector refNames = refs1.attr("levels");
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsets(sizes, offsets);
    std::vector < int > colSums;
    computeSymmetricColSum(data, sizes, offsets, nElements, colSums);
    std::vector < bool > passedThreshold (colSums.size());
    for (size_t i = 0; i < colSums.size(); ++i) {
        passedThreshold[i] = (colSums[i] >= threshold);
    }
    int nPassed = std::count(passedThreshold.begin(), passedThreshold.end(), true);
    int nFailed = colSums.size() - nPassed;
    Rcerr << "\tKeeping " << nPassed << "/" << nElements << " rows.\n";
    std::vector < int > refs1C;
    std::vector < int > refs2C;
    std::vector < int > bins1C;
    std::vector < int > bins2C;
    std::vector < int > countsC;
    for (int i = 0; i < refs1.size(); ++i) {
        long int bin1 = offsets[refs1[i] - 1] + bins1[i]; // factors (such as 'refs') start with 1 in R
        long int bin2 = offsets[refs2[i] - 1] + bins2[i]; // factors (such as 'refs') start with 1 in R
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
    DataFrame outputMatrix = Rcpp::DataFrame::create(_["ref1"] = refs1, _["bin1"] = bins1, _["ref2"] = refs2, _["bin2"] = bins2, _["count"] = counts);
    // Compute the set of low count rows
    IntegerVector refs (nFailed);
    IntegerVector bins (nFailed);
    int id = 0;
    for (int refId = 0; refId < sizes.size(); ++refId) {
        for (int bin = 0; bin <= sizes[refId]; ++bin) {
            if (! passedThreshold[offsets[refId] + bin]) {
                if (id >= refs.size()) {
                    Rcerr << "Problem during removeLowCountRows: I found " << nFailed << " rows, but I am reading the #" << id << ".\n";
                    stop("Stopping here.");
                }
                refs[id] = refId;
                bins[id] = bin;
                ++id;
            }
        }
    }
    DataFrame outputRows = Rcpp::DataFrame::create(_["ref"] = refs, _["bin"] = bins);
    List outputData = List::create(_["matrix"] = outputMatrix, _["rows"] = outputRows);
    return outputData;
}

// [[Rcpp::export]]
void normalizeHighCountRowsCpp (DataFrame &data, IntegerVector &sizes) {
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    std::vector < long int > offsets (sizes.size());
    long int nElements = computeOffsets(sizes, offsets);
    std::vector < int > colSums;
    computeSymmetricColSum(data, sizes, offsets, nElements, colSums);
    NumericVector colSumsR = wrap(colSums);
    int nCols = colSums.size();
    double meanColSum = mean(colSumsR);
    double sdColSum   = sd(colSumsR);
    NumericVector normFactors (nCols, 1.0);
    double threshold = meanColSum + 3 * sdColSum;
    int nColChanged = 0;
    for (int i = 0; i < nCols; ++i) {
        if (colSums[i] > threshold) {
            normFactors[i] = meanColSum / colSums[i];
            ++nColChanged;
        }
    }
    for (int i = 0; i < refs1.size(); ++i) {
        long int bin1 = offsets[refs1[i] - 1] + bins1[i]; // factors (such as 'refs') start with 1 in R
        long int bin2 = offsets[refs2[i] - 1] + bins2[i]; // factors (such as 'refs') start with 1 in R
        counts[i] = counts[i] * sqrt(normFactors[bin1] * normFactors[bin2]);
    }
    Rcerr << "\t\tModified " << nColChanged << " / " << nCols << "\n";
}
