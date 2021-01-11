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
    Rcerr << "Keeping " << nPassed << "/" << nElements << " rows.\n";
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

// [[Rcpp::export]]
DataFrame removeSmallScaffoldsCpp (DataFrame &data, int nRefs, int threshold) {
    IntegerVector   tmp = data["ref1"];
    CharacterVector refNames = tmp.attr("levels");
    std::vector < int > refs1  = as < std::vector < int > > (data["ref1"]);
    std::vector < int > refs2  = as < std::vector < int > > (data["ref2"]);
    std::vector < int > bins1  = as < std::vector < int > > (data["bin1"]);
    std::vector < int > bins2  = as < std::vector < int > > (data["bin2"]);
    std::vector < int > counts = as < std::vector < int > > (data["count"]);
    std::vector < int > refMax (nRefs + 1, 0);               // factors (such as 'refs') start with 1 in R
    std::vector < bool > passedThreshold (nRefs + 1, false);
    int nElements, nPassed, previousNPassed = 0;
    Rcerr << "Removing small scaffolds (currently: " << nRefs << ").\n";
    // Do it several times, since deleting a ref can trigger another deletion
    while (true) {
        nElements = 0;
        refMax.assign(nRefs + 1, 0); // factors (such as 'refs') start with 1 in R
        for (size_t i = 0; i < refs1.size(); ++i) {
            refMax[refs1[i]] = std::max<int>(bins1[i], refMax[refs1[i]]);
            refMax[refs2[i]] = std::max<int>(bins2[i], refMax[refs2[i]]);
        }
        for (int i = 1; i <= nRefs; ++i) { // factors start with 1
            passedThreshold[i] = (refMax[i] >= threshold);
        }
        nPassed = std::count(passedThreshold.begin(), passedThreshold.end(), true);
        if ((nPassed == 0) || (nPassed == previousNPassed)) {
            break;
        }
        // In-place update of the data
        for (size_t i = 0; i < refs1.size(); ++i) {
            if (passedThreshold[refs1[i]] && passedThreshold[refs2[i]]) {
                refs1[nElements]  = refs1[i];
                refs2[nElements]  = refs2[i];
                bins1[nElements]  = bins1[i];
                bins2[nElements]  = bins2[i];
                counts[nElements] = counts[i];
                ++nElements;
            }
        }
        refs1.resize(nElements);
        refs2.resize(nElements);
        bins1.resize(nElements);
        bins2.resize(nElements);
        counts.resize(nElements);
        previousNPassed = nPassed;
    }
    // Recompute factors for the refs
    CharacterVector refNamesPassed (nPassed);
    std::vector < int > translateRefIds (nRefs + 1);
    int refId = 1;                     // factors start with 1
    for (int i = 1; i <= nRefs; ++i) { // factors start with 1
        if (passedThreshold[i]) {
            // CharacterVector start with 0
            refNamesPassed[refId - 1] = refNames[i - 1];
            translateRefIds[i] = refId;
            ++refId;
        }
    }
    for (size_t i = 0; i < refs1.size(); ++i) {
        refs1[i] = translateRefIds[refs1[i]];
        refs2[i] = translateRefIds[refs2[i]];
    }
    IntegerVector refs1R  = wrap(refs1);
    IntegerVector refs2R  = wrap(refs2);
    IntegerVector bins1R  = wrap(bins1);
    IntegerVector bins2R  = wrap(bins2);
    IntegerVector countsR = wrap(counts);
    refs1R.attr("class") = "factor";
    refs2R.attr("class") = "factor";
    refs1R.attr("levels") = refNamesPassed;
    refs2R.attr("levels") = refNamesPassed;
    DataFrame outputData = Rcpp::DataFrame::create(_["ref1"]= refs1R, _["bin1"]= bins1R, _["ref2"]= refs2R, _["bin2"]= bins2R, _["count"]= countsR);
    Rcerr << "Keeping " << nPassed << " big references.\n";
    return outputData;
}
