#include "sharedFunctions.h"
#include <progress.hpp>
#include <progress_bar.hpp>

void computeRefsBins (IntegerVector &sizes, int nElements, IntegerVector &refs, IntegerVector &bins) {
    refs = IntegerVector(nElements);
    bins = IntegerVector(nElements);
    int cpt = 0;
    for (int ref = 0; ref < sizes.size(); ++ref) {
        for (int bin = 0; bin <= sizes[ref]; ++bin) {
            if (cpt >= nElements) stop("Error while computing bin sums: got more than " + std::to_string(nElements) + " elements + (ref: " + std::to_string(ref) + "/" + std::to_string(sizes.size()) + ", bin: " + std::to_string(bin) + "/" + std::to_string(sizes[ref]) + ").");
            refs[cpt] = ref + 1;
            bins[cpt] = bin;
            ++cpt;
        }
    }
    if (cpt != nElements) stop("Error while computing bin sums: got " + std::to_string(cpt) + " vs " + std::to_string(nElements));
    CharacterVector refNames = sizes.names();
    refs.attr("class") = "factor";
    refs.attr("levels") = refNames;
}

void computeRefsBinsMeta (IntegerVector &sizes, int nElements, int metaSize, bool metaBins, IntegerVector &refs, IntegerVector &bins) {
    refs = IntegerVector(nElements);
    bins = IntegerVector(nElements);
    int cpt = 0;
    for (int ref = 0; ref < sizes.size(); ++ref) {
        int bin;
        for (bin = 0; bin <= sizes[ref] / metaSize; ++bin) {
            if (cpt >= nElements) stop("Error while computing meta bin sums: got more than " + std::to_string(nElements) + " elements + (ref: " + std::to_string(ref) + "/" + std::to_string(sizes.size()) + ", bin: " + std::to_string(bin) + "/" + std::to_string(sizes[ref]) + ").");
            refs[cpt] = ref + 1;
            if (metaBins) {
                bins[cpt] = bin;
            }
            else {
                bins[cpt] = bin * metaSize + metaSize / 2;
            }
            ++cpt;
        }
        if (! metaBins) {
            bins[cpt-1] = ((bin - 1) * metaSize + sizes[ref]) / 2;
        }
    }
    if (cpt != nElements) stop("Error while computing meta bin sums: got " + std::to_string(cpt) + " vs " + std::to_string(nElements));
    CharacterVector refNames = sizes.names();
    refs.attr("class") = "factor";
    refs.attr("levels") = refNames;
}

long int computeOffsets (IntegerVector &sizes, std::vector < long int > &offsets) {
    int nElements = 0;
    for (int i = 0; i < sizes.size(); ++i) {
        offsets[i] = nElements;
        nElements += sizes[i] + 1;
    }
    return nElements;
}

long int computeOffsetsMeta (IntegerVector &sizes, int metaSize, std::vector < long int > &offsets) {
    int nElements = 0;
    for (int i = 0; i < sizes.size(); ++i) {
        offsets[i] = nElements;
        nElements += (sizes[i] / metaSize) + 1;
    }
    return nElements;
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

void setOutlierBinsToBool (IntegerVector &refs, IntegerVector &bins, std::vector < long int > &offsets, long int genomeSize, std::vector < bool > &outlierBinsBool) {
    outlierBinsBool = std::vector < bool > (genomeSize, false);
    int nOutlierBins = refs.size();
    assert(refs.size() == bins.size());
    for (long int outlierId = 0; outlierId < nOutlierBins; ++outlierId) {
        assert(refs[outlierId] <= static_cast<int>(offsets.size()));
        assert(offsets[refs[outlierId] - 1] + bins[outlierId] < genomeSize);
        outlierBinsBool[offsets[refs[outlierId] - 1] + bins[outlierId]] = true;
    }
}

void setFullMatrix (IntegerVector &refs1, IntegerVector &refs2, IntegerVector &bins1, IntegerVector &bins2, IntegerVector &counts, std::vector < long int > &offsets, long int genomeSize, std::vector < bool > &outlierBinsBool, int distance, std::vector < std::vector < double > > &fullMatrix) {
    long int n = refs1.size();
    fullMatrix = std::vector < std::vector < double > > (genomeSize, std::vector < double > (distance + 1, 0.0)); // 1st axis: bins; 2nd: distance
    for (long int i = 0; i < n; ++i) {
        int ref1 = refs1[i];
        int ref2 = refs2[i];
        if (ref1 == ref2) {
            int offset = offsets[ref1 - 1];
            if ((! outlierBinsBool[offset + bins1[i]]) && (! outlierBinsBool[offset + bins2[i]])) {
                int d = bins1[i] - bins2[i];
                if (d < 0) stop("Bin1 should not be less than bin2 in 'setFullMatrix'.\n");
                if ((counts[i] > 0) && (d <= distance)) {
                    if (counts[i] < 1) stop("Counts should be greater than 1 in 'setFullMatrix'.\n");
                    fullMatrix[offset + bins1[i]][d] = counts[i];
                }
            }
        }
    }
}

void setFullMatrixMeta (IntegerVector &refs1, IntegerVector &refs2, IntegerVector &bins1, IntegerVector &bins2, IntegerVector &counts, std::vector < long int > &offsets, long int genomeSize, std::vector < bool > &outlierBinsBool, int distance, int metaSize, std::vector < std::vector < double > > &fullMatrix) {
    long int n = refs1.size();
    fullMatrix = std::vector < std::vector < double > > (genomeSize, std::vector < double > (distance + 1, 0.0)); // 1st axis: bins; 2nd: distance
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
            if ((! outlierBinsBool[offset + bin1]) && (! outlierBinsBool[offset + bin2])) {
                if (d < 0) stop("Bin1 should not be less than bin2 in 'setFullMatrix'.\n");
                --d; // discard diagonal counts
                if (d >= 0) {
                    d /= metaSize;
                    if ((counts[i] > 0) && (d <= distance)) {
                        fullMatrix[offset + bin1][d] += counts[i];
                    }
                }
            }
        }
    }
}

void computeDistanceCounts (IntegerVector &refs1, IntegerVector &refs2, IntegerVector &bins1, IntegerVector &bins2, IntegerVector &counts, std::vector < long int > &offsets, std::vector < bool > &outlierBinsBool, int distance, IntegerVector &sizes, std::vector < long int > &sums, std::vector < int > &nCounts) {
    long int n = refs1.size();
    sums.resize(distance + 1);
    nCounts.resize(distance + 1);
    Progress progress1(n, true);
    for (long int i = 0; i < n; ++i) {
        int ref1 = refs1[i];
        int ref2 = refs2[i];
        if (ref1 == ref2) {
            int offset = offsets[ref1 - 1];
            if ((! outlierBinsBool[offset + bins1[i]]) && (! outlierBinsBool[offset + bins2[i]])) {
                int d = bins1[i] - bins2[i];
                if (d < 0) stop("Bin1 should not be less than bin2 in 'computeMD'.\n");
                if (d <= distance) {
                    sums[d] += counts[i];
                }
            }
        }
        progress1.increment();
    }
    Progress progress2(outlierBinsBool.size(), true);
    for (int ref = 0; ref < sizes.size(); ++ref) {
        int offset = offsets[ref];
        for (int bin1 = 0; bin1 <= sizes[ref]; ++bin1) {
            if (! outlierBinsBool[offset + bin1]) {
                for (int d = 0; (d <= distance) && (d <= bin1); ++d) {
                    int bin2 = bin1 - d;
                    if (! outlierBinsBool[offset + bin2]) {
                        ++nCounts[d];
                    }
                }
            }
            progress2.increment();
        }
    }
}
