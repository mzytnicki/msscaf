#include <vector>                                                                                                                                                                                                  
#include <Rcpp.h>                                                                                                                                                                                                  
#include "sharedFunctions.h"
// [[Rcpp::plugins(cpp11)]]                                                                                                                                                                                        
                                                                                                                                                                                                                   
using namespace Rcpp;                                                                                                                                                                                              

// [[Rcpp::export]]
DataFrame makeSymmetricRefCpp (DataFrame &data) {
    IntegerVector bins1R = data["bin1"];
    IntegerVector bins2R = data["bin2"];
    IntegerVector countsR = data["count"];
    std::vector < int > bins1  = as < std::vector < int > > (bins1R);
    std::vector < int > bins2  = as < std::vector < int > > (bins2R);
    std::vector < int > counts = as < std::vector < int > > (countsR);
    unsigned long int n = bins1.size();
    for (unsigned long i = 0; i < n; ++i) {
        int bin1  = bins1[i];
        int bin2  = bins2[i];
        int count = counts[i];
        if (bin1 != bin2) {
            bins1.push_back(bin2);
            bins2.push_back(bin1);
            counts.push_back(count);
        }
    }
    return Rcpp::DataFrame::create(_["bin1"] = wrap(bins1), _["bin2"] = wrap(bins2), _["count"] = wrap(counts));
}

// Set all the outlier bins to NA
// Possibly overwrite current data, and create new cells
// Data: a symmetric tibble: bin1, bin2, count
// Outliers: a 1-row tibble: bin
// Size: size of the reference
// MinLim: min. bin to consider
// MaxLim: max. bin to consider
// [[Rcpp::export]]
DataFrame removeOutliersRefCpp (DataFrame &data, DataFrame &outliersTibble, int size, int minLim = -1, int maxLim = -1) {
    IntegerVector bins1R   = data["bin1"];
    IntegerVector bins2R   = data["bin2"];
    IntegerVector countsR  = data["count"];
    IntegerVector outliers = outliersTibble["bin"];
    unsigned long int n = bins1R.size();
    if (outliers.size() == 0) return data;
    std::vector < bool > outliersBool (size + 1, false);
    if (minLim == -1) minLim = 0;
    if (maxLim == -1) minLim = size;
    for (int outlier: outliers) {
        outliersBool[outlier] = true;
    }
    std::vector < int > bins1;
    std::vector < int > bins2;
    std::vector < int > counts;
    for (int outlier: outliers) {
        if ((outlier >= minLim) && (outlier <= maxLim)) {
            for (int bin = minLim; bin <= maxLim; ++bin) {
                bins1.push_back(outlier);
                bins2.push_back(bin);
                counts.push_back(NA_INTEGER);
                if (bin != outlier) {
                    bins1.push_back(bin);
                    bins2.push_back(outlier);
                    counts.push_back(NA_INTEGER);
                }
            }
        }
    }
    for (unsigned long i = 0; i < n; ++i) {
        int bin1  = bins1R[i];
        int bin2  = bins2R[i];
        int count = countsR[i];
        if ((bin1 >= minLim) && (bin1 <= maxLim) && (bin2 >= minLim) && (bin2 <= maxLim)) {
            if ((! outliersBool[bin1]) && (! outliersBool[bin2])) {
                bins1.push_back(bin1);
                bins2.push_back(bin2);
                counts.push_back(count);
            }
        }
    }
    return Rcpp::DataFrame::create(_["bin1"] = wrap(bins1), _["bin2"] = wrap(bins2), _["count"] = wrap(counts));
}

// Set all the outlier bins to NA
// Data already are metabins
// Possibly overwrite current data, and create new cells
// Data: a (symmetric) tibble: ref1, re2, bin1, bin2, count
// Outliers: a tibble: ref, bin
// Size: sizes of the reference
// [[Rcpp::export]]
DataFrame removeOutliersCpp (DataFrame &data, DataFrame &outliersTibble, IntegerVector sizes) {
    IntegerVector refs1R   = data["ref1"];
    IntegerVector refs2R   = data["ref2"];
    IntegerVector bins1R   = data["bin1"];
    IntegerVector bins2R   = data["bin2"];
    IntegerVector countsR  = data["count"];
    IntegerVector outlierRefs = outliersTibble["ref"];
    IntegerVector outlierBins = outliersTibble["bin"];
    CharacterVector refNames = refs1R.attr("levels");
    assert(refNames.size() == sizes.size()); 
    std::vector < int > refs1;
    std::vector < int > refs2;
    std::vector < int > bins1;
    std::vector < int > bins2;
    std::vector < int > counts;
    unsigned long int nInput = bins1R.size();
    int nOutliers = outlierRefs.size();
    if (nOutliers == 0) return data;
    int genomeSize = sum(sizes);
    refs1.reserve(2 * nOutliers * genomeSize + nInput);
    refs2.reserve(2 * nOutliers * genomeSize + nInput);
    bins1.reserve(2 * nOutliers * genomeSize + nInput);
    bins2.reserve(2 * nOutliers * genomeSize + nInput);
    counts.reserve(2 * nOutliers * genomeSize + nInput);
    std::vector < std::vector < bool > > outliersBool (sizes.size() + 1); // factors start with 1
    for (unsigned int i = 0; i < sizes.size(); ++i) {
        assert(outliersBool.size() >= i);
        assert(sizes.size() > i);
        outliersBool[i + 1] = std::vector < bool > (sizes[i] + 1, false);
    }
    for (int i = 0; i < nOutliers; ++i) {
        assert(outlierRefs[i] < static_cast<int>(outliersBool.size()));
        assert(outlierBins[i] < static_cast<int>(outliersBool[outlierRefs[i]].size()));
        outliersBool[outlierRefs[i]][outlierBins[i]] = true;
    }
    std::vector < int > allRefs;
    std::vector < int > allBins;
    for (int ref = 1; ref <= sizes.size(); ++ref) {
        for (int bin = 0; bin <= sizes[ref - 1]; ++bin) {
            allRefs.push_back(ref);
            allBins.push_back(bin);
        }
    }
    unsigned int nBins = allBins.size();
    for (int i = 0; i < outlierRefs.size(); ++i) {
        int outlierRef = outlierRefs[i];
        int outlierBin = outlierBins[i];
        refs1.insert(refs1.end(), nBins, outlierRef);
        bins1.insert(bins1.end(), nBins, outlierBin);
        refs2.insert(refs2.end(), allRefs.begin(), allRefs.end());
        bins2.insert(bins2.end(), allBins.begin(), allBins.end());
        counts.insert(counts.end(), nBins, NA_INTEGER);
        refs1.insert(refs1.end(), allRefs.begin(), allRefs.end());
        bins1.insert(bins1.end(), allBins.begin(), allBins.end());
        refs2.insert(refs2.end(), nBins, outlierRef);
        bins2.insert(bins2.end(), nBins, outlierBin);
        counts.insert(counts.end(), nBins, NA_INTEGER);
/*
        for (int ref = 1; ref <= sizes.size(); ++ref) {
            for (int bin = 0; bin <= sizes[ref - 1]; ++bin) {
                refs1.push_back(ref);
                refs2.push_back(outlierRef);
                bins1.push_back(bin);
                bins2.push_back(outlierBin);
                counts.push_back(NA_INTEGER);
                if ((outlierRef != ref) || (outlierBin != bin)) {
                    refs1.push_back(outlierRef);
                    refs2.push_back(ref);
                    bins1.push_back(outlierBin);
                    bins2.push_back(bin);
                    counts.push_back(NA_INTEGER);
                }
            }
        }
*/
    }
    for (unsigned long i = 0; i < nInput; ++i) {
        int ref1  = refs1R[i];
        int ref2  = refs2R[i];
        int bin1  = bins1R[i];
        int bin2  = bins2R[i];
        int count = countsR[i];
        assert(ref1 < static_cast<int>(outliersBool.size()));
        assert(bin1 < static_cast<int>(outliersBool[ref1].size()));
        assert(ref2 < static_cast<int>(outliersBool.size()));
        assert(bin2 < static_cast<int>(outliersBool[ref2].size()));
        if ((! outliersBool[ref1][bin1]) && (! outliersBool[ref2][bin2])) {
            refs1.push_back(ref1);
            refs2.push_back(ref2);
            bins1.push_back(bin1);
            bins2.push_back(bin2);
            counts.push_back(count);
        }
    }
    IntegerVector outputRefs1R = wrap(refs1);
    IntegerVector outputRefs2R = wrap(refs2);
    outputRefs1R.attr("class") = "factor";
    outputRefs2R.attr("class") = "factor";
    outputRefs1R.attr("levels") = refNames;
    outputRefs2R.attr("levels") = refNames;
    return Rcpp::DataFrame::create(_["ref1"] = outputRefs1R, _["ref2"] = outputRefs2R, _["bin1"] = wrap(bins1), _["bin2"] = wrap(bins2), _["count"] = wrap(counts));
}
