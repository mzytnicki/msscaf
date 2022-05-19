#include <vector>                                                                                                                                                                                                  
#include <Rcpp.h>                                                                                                                                                                                                  
#include "sharedFunctions.h"
// [[Rcpp::plugins(cpp11)]]                                                                                                                                                                                        
                                                                                                                                                                                                                   
using namespace Rcpp;                                                                                                                                                                                              

/* Matrix is sparse symmetric, represented as list of triplets: (ref1, ref2, bin1, bin2, c),
 *   with bin1 (i) >= bin2 (j).
 * Distance is the maximum link range (here, 4), values further away from
 *   diagonal are discarded: i <= j + d.  'X' cells are discarded.
 * Triangle are 'a' cells.  They are computed with row/col sums.
 *   tri[i] = tri[i-1] + rowSum[i] - colSum[i-1]
 *   When the diagonal is excluded:
 *    - the rowSum/colSums do not start from diagonal, but with an offset of 2
 *    - tri[0] = 0
 *    - tri[i] = tri[i-1] + rowSum[i-1] - colSum[i]
 * colSums are sumBins1, rowSum are sumBins2
 * Squares are 'a' and 'b' cells.  They are computed with diag sums.
 *   squ[i] = squ[i-1] + diag[2i+d] + diag[2i+d-1] - diag[2i-d] - diag[2i-d+1]
 *
 *      0               i             (bin1)
 *    |   |   |   |   |   |   |   |   |   | 
 *   -+---+---+---+---+---+---+---+---+---+-
 * 0  |   |   |   |   | X | X | X | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   |   | X | X | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   |   | a | X | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   |   | a | a | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 * j  |   |   |   |   | O |   |   |   | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   |   |   |   |   |   | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   |   |   |   |   |   | 
 * b -+---+---+---+---+---+---+---+---+---+-
 * i  |   |   |   |   |   |   |   |   |   | 
 * n -+---+---+---+---+---+---+---+---+---+-
 * 2  |   |   |   |   |   |   |   |   |   | 
 *
 * Recursion scheme:
 *      0               i                   
 *    |   |   |   |   |   |   |   |   |   | 
 *   -+---+---+---+---+---+---+---+---+---+-
 * 0  |   |   |   |   | a | X | X | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   | a | c | X | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   |   | a | c | c | X | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *    |   |   |   | 1 |   | b | b | b | X | 
 *   -+---+---+---+---+---+---+---+---+---+-
 * j  |   |   |   |   | 2 |   |   |   |   | 
 *   -+---+---+---+---+---+---+---+---+---+-
 *                         \   \   \   \
 *          diag. distance: 0   1   2   3
 *                           \   \   \   \
 * a's are discarded
 * b's are kept
 * c's are added
 * Note: the row/col sums should have min distance of 2
 */

// [[Rcpp::export]]
DataFrame computeMeanTrianglesCpp (DataFrame &data, int distance, int metaSize, IntegerVector sizesIn, DataFrame outliers) {
    IntegerVector sizes = clone(sizesIn);
    int nRefs = sizes.size();
    const double  pseudoCount = 1.1;
    bool  useMetaBins    = (metaSize > 1);
    IntegerVector refs1  = data["ref1"];
    IntegerVector refs2  = data["ref2"];
    IntegerVector bins1  = data["bin1"];
    IntegerVector bins2  = data["bin2"];
    IntegerVector counts = data["count"];
    IntegerVector outlierRefs  = outliers["ref"];
    IntegerVector outlierBins  = outliers["bin"];
    std::vector < long int > offsets (nRefs);
    std::vector < std::vector < double > > fullMatrix;
    std::vector < double > normFactors  (distance + 1, 0); // normalization factors
    std::vector < int >    nNormFactors (distance + 1, 0); // number of normalization factors
    IntegerVector refs;
    IntegerVector bins;
    std::vector < bool > outlierBinsBool;
    if (useMetaBins) {
        for (int i = 0; i < sizes.size(); ++i) {
            sizes[i] = sizes[i] / metaSize;
        }
    }
    long int nElements = computeOffsets(sizes, offsets);
    setOutlierBinsToBool(outlierRefs, outlierBins, offsets, nElements, outlierBinsBool);
    if (useMetaBins) {
        setFullMatrixMeta(refs1, refs2, bins1, bins2, counts, offsets, nElements, outlierBinsBool, distance, metaSize, fullMatrix);
    }
    else {
        setFullMatrix(refs1, refs2, bins1, bins2, counts, offsets, nElements, outlierBinsBool, distance, fullMatrix);
    }
    std::vector < double > sumBin1 (nElements);            // row sum
    std::vector < double > sumBin2 (nElements);            // col sum
    std::vector < int >    nBin1   (nElements);            // # row cells
    std::vector < int >    nBin2   (nElements);            // # col cells
    NumericVector triCounts     (nElements, 0.0);          // sum of the triangles
    IntegerVector nTriCells     (nElements, 0);            // # cells of the triangles
    NumericVector meanCounts    (nElements, 0.0);          // output mean counts
    computeRefsBins(sizes, nElements, refs, bins);
    if (useMetaBins) {
        int cpt = 0;
        for (int ref = 0; ref < sizes.size(); ++ref) {
            int bin;
            for (bin = 0; bin <= sizes[ref]; ++bin, ++cpt) {
                if (cpt >= bins.size()) stop("Errors while converting meta-bins.");
                bins[cpt] *= metaSize;
            }
            bins[cpt-1] = ((bin - 1) * metaSize + sizes[ref]) / 2;
        }
    }
    // Take log
    // Normalize w.r.t. distance to diagonal,
    // Find normalization factor.
    // Caveat: some values should be not used (distance < bin).
    for (int refId = 0; refId < nRefs; ++refId) {
        int offset = offsets[refId];
        for (int bin = 0; bin <= sizes[refId]; ++bin) {
            if ((useMetaBins) || (! outlierBinsBool[offset + bin])) {
                for (int d = 0; (d <= distance) && (d <= bin); ++d) {
                    if ((useMetaBins) || (! outlierBinsBool[offset + bin - d])) {
                        fullMatrix[offset + bin][d] = log10(fullMatrix[offset + bin][d] + pseudoCount);
                        normFactors[d] += fullMatrix[offset + bin][d];
                        ++nNormFactors[d];
                    }
                }
                // if ((nNormFactors[i] > 0) && (normFactors[i] == 0)) stop("Problem of null normalization factor (d = " + std::to_string(i) + ", # counts = " + std::to_string(nNormFactors[i]) + ") in 'computeMeanTrianglesCpp'\n");
                // Rcerr << d << ": " << normFactors[d] << "\n";
            }
        }
    }
    for (int d = 0; d <= distance; ++d) {
        normFactors[d] /= nNormFactors[d];
//Rcerr << "norm fact @ " << d << ": " << normFactors[d] << "\n";
    }
/*
for (int i = 0; i < size; ++i) {
    Rcerr << i << ":";
    for (int d = 0; d <= distance; ++d) {
       Rcerr << " " << fullMatrix[i][d] << " (" << log10(fullMatrix[i][d] / normFactors[d]) << ")  ";
    }
    Rcerr << "\n";
}
*/
/*
for (int d = 0; d <= distance; ++d) {
    double s = 0.0;
    for (int bin = 0; bin < size && d <= bin; ++bin) {
        s += fullMatrix[bin][d] - normFactors[d];
    }
    Rcerr << d << ":" << s << "\n";
}
*/
    // for (int d = 0; d <= distance; ++d) { double s = 0; for (int i = 0; i < size; ++i) s += fullMatrix[i][d]; if (abs(s) < 0.01) stop("Problem during first normalization in 'computeMeanTrianglesCpp'.\n"); }
    for (int refId = 0; refId < nRefs; ++refId) {
        int offset = offsets[refId];
        for (int bin1 = 0; bin1 <= sizes[refId]; ++bin1) {
            if ((useMetaBins) || (! outlierBinsBool[offset + bin1])) {
                int startD = (useMetaBins)? 0: 2;
                for (int d = startD; (d <= distance) && (d <= bin1); ++d) {
                    int bin2 = bin1 - d;
                    if ((useMetaBins) || (! outlierBinsBool[offset + bin2])) {
                        double count = fullMatrix[offset + bin1][d];
                        if (count < 0) stop("Problem of negative count in 'computeMeanTrianglesCpp'\n");
                        // if (normFactors[bins1Red[i] - bins2Red[i]] == 0) stop("Problem of null denominator (i = " + std::to_string(bins1Red[i]) + ", j = " + std::to_string(bins2Red[i]) + ", value = " + std::to_string(counts[i]) +  ", # counts = " + std::to_string(nNormFactors[bins1Red[i] - bins2Red[i]]) + ") in 'computeMeanTrianglesCpp'\n");
                        double norm  = normFactors[d];
                        if ((norm == 0.0) && (count != 0.0)) stop("Problem during normalization in 'computeMeanTrianglesCpp'\n");
                        //Rcerr << bin1 << "-" << bin2 << " (" << d << "): " << count << "-" << norm << " -> " << ((norm == 0.0)? 0.0: log10(count / norm)) << "\n";
                        count -= norm;
                        sumBin1[offset + bin1] += count;
                        sumBin2[offset + bin2] += count;
                        ++nBin1[offset + bin1];
                        ++nBin2[offset + bin2];
                        //Rcerr << "sumbin2 " << d << " " << bin2 << " " << sumBin2[bin2] << "\n";
                    }
                }
            }
        }
    }
    //Rcerr << "Sum bins1: " << std::accumulate(sumBin1.begin(), sumBin1.end(), 0.0) << "\n";
    if (useMetaBins) {
        for (int refId = 0; refId < nRefs; ++refId) {
            int offset         = offsets[refId];
            triCounts[offset]  = sumBin2[offset];
            nTriCells[offset]  = nBin2[offset];
            meanCounts[offset] = (nTriCells[offset] == 0)? 0: triCounts[offset] / nTriCells[offset];
            for (int bin = 1; bin <= sizes[refId]; ++bin) {
                //Rcerr << "Bin:" << i << ": " << sumBin1[i] << ", " << sumBin2[i] << "\n";
                triCounts[offset + bin]  = triCounts[offset + bin - 1] + sumBin2[offset + bin] - sumBin1[offset + bin - 1];
                nTriCells[offset + bin]  = nTriCells[offset + bin - 1] + nBin2[offset + bin]   - nBin1[offset + bin - 1];
                meanCounts[offset + bin] = (nTriCells[offset + bin] == 0)? 0: triCounts[offset + bin] / nTriCells[offset + bin];
            }
        }
    }
    else {
        for (int refId = 0; refId < nRefs; ++refId) {
            int    offset       = offsets[refId];
            double currentCount = 0.0;
            double currentN     = 0.0;
            for (int bin = 1; bin <= sizes[refId]; ++bin) {
                currentCount += sumBin2[offset + bin - 1] - sumBin1[offset + bin];
                currentN     += nBin2[offset + bin - 1]   - nBin1[offset + bin];
                if ((useMetaBins) || (! outlierBinsBool[offset + bin])) {
                    //Rcerr << "Bin:" << i << ": " << sumBin1[i] << ", " << sumBin2[i] << "\n";
                    triCounts[offset + bin] = currentCount;
                    nTriCells[offset + bin] = currentN;
                    meanCounts[offset + bin] = (currentN == 0)? 0: currentCount / currentN;
                }
            }
        }
    }
    return Rcpp::DataFrame::create(_["ref"] = refs, _["bin"] = bins, _["fcMeanCount"] = meanCounts, _["nCells"] = nTriCells);
}

bool isClose (int refId1, int bin1, int refId2, int bin2, int distance) {
    return ((refId1 == refId2) && (abs(bin1 - bin2) <= distance));
}

bool isOneClose (int refId, int bin, IntegerVector &refs, IntegerVector &bins, int distance) {
    for (int i = 0; i < refs.size(); ++i) {
        if (isClose(refId, bin, refs[i], bins[i], distance)) {
            return true;
        }
    }
    return false;
}


// Remove close breaks (closer than given distance).
// Suppose that the data is sorted by adj. p-value.
// [[Rcpp::export]]
DataFrame removeNearEqualBreaksCpp (DataFrame &breaks, int distance) {
    IntegerVector refs  = breaks["ref"];
    IntegerVector bins  = breaks["bin"];
    //NumericVector fcs   = breaks["fcMeanCount"];
    //IntegerVector ncs   = breaks["nCells"];
    NumericVector pvals = breaks["pvalue"];
    //NumericVector padjs = breaks["padj"];
    CharacterVector refNames = refs.attr("levels");
    IntegerVector keptRefs;
    IntegerVector keptBins;
    //NumericVector keptFcs;
    //IntegerVector keptNcs;
    NumericVector keptPvals;
    //NumericVector keptPadjs;
    long int      n = refs.size();
    for (int i = 0; i < n; ++i) {
        if (! isOneClose(refs[i], bins[i], keptRefs, keptBins, distance)) {
            keptRefs.push_back(refs[i]);
            keptBins.push_back(bins[i]);
            //keptFcs.push_back(fcs[i]);
            //keptNcs.push_back(ncs[i]);
            keptPvals.push_back(pvals[i]);
            //keptPadjs.push_back(padjs[i]);
        }
    }
    keptRefs.attr("class") = "factor";
    keptRefs.attr("levels") = refNames;
    return Rcpp::DataFrame::create(_["ref"] = keptRefs, _["bin"] = keptBins, _["pvalue"] = keptPvals);
    //return Rcpp::DataFrame::create(_["ref"] = keptRefs, _["bin"] = keptBins, _["fcMeanCount"] = keptFcs, _["nCells"] = keptNcs, _["pvalue"] = keptPvals, _["padj"] = keptPadjs);
}
