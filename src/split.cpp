#include <vector>
#include <utility>
#include <string>
#include <sstream>
#include <algorithm>
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

using PositionConvertor = std::vector < std::vector < std::pair < int, int > > >;

// Find the largest sub-contig for each contig obtained after split.
// The split dataframe contains 2 columns: ref, and bin.
// The split dataframe is ordered by increasing ref, bin.
// The output vector is the index of the largest sub-contig, for each contig
void findLargestSubContigs (DataFrame splits, IntegerVector sizes, std::vector < int > &largestSubContigs) {
    IntegerVector refs = splits["ref"];
    IntegerVector bins = splits["bin"];
    int nRefs          = sizes.size();
    int nSplits        = refs.size();
    largestSubContigs = std::vector < int > (nRefs + 1, 0); // factors start with 1
    int prevRef     = -1;
    int prevSplit   = 0;
    int prevMaxSize = 0;
    int subContig   = 0;
    if (nSplits == 0) return;
    for (int splitId = 0; splitId < nSplits; ++splitId) {
        int ref   = refs[splitId];
        int split = bins[splitId];
        int size;
        if ((ref != prevRef) && (prevRef != -1)) {
            size = sizes[prevRef - 1] - prevSplit;
            if (size > prevMaxSize) {
                largestSubContigs[prevRef] = subContig;
            }
            prevSplit   = 0;
            prevMaxSize = 0;
            subContig   = 0;
        }
        size = split - prevSplit;
        if (size > prevMaxSize) {
            prevMaxSize = size;
            largestSubContigs[ref] = subContig;
        }
        prevSplit = split;
        prevRef   = ref;
        ++subContig;
    }
    int size = sizes[prevRef - 1] - prevSplit;
    if (size > prevMaxSize) {
        prevSplit = subContig;
    }
    size = sizes[prevRef - 1] - prevSplit;
    if (size > prevMaxSize) {
        largestSubContigs[prevRef] = subContig;
    }
}

// Update the convertor between two bounds.
void updateConvertor (PositionConvertor &convertor, int ref, int largestSubContig, int &newRef, int subContig, int prevSplit, int split) {
    // discard previous split (split will be removed in the next call)
    if (prevSplit != 0) {
        convertor[ref][prevSplit] = std::make_pair(-1, -1);
        ++prevSplit;
    }
    if (subContig == largestSubContig) {
        for (int bin = prevSplit; bin < split; ++bin) {
            convertor[ref][bin] = std::make_pair(ref, bin - prevSplit);
        }
        return;
    }
    for (int bin = prevSplit; bin < split; ++bin) {
        convertor[ref][bin] = std::make_pair(newRef, bin - prevSplit);
    }
    ++newRef;
}

// Convert splits to a list of (refId, bin) -> (newRefId, newBin)
// The split dataframe contains 3 columns: ref, newRef, and bin.
// The split dataframe is ordered by increasing ref, bin.
void setConvertor (DataFrame splits, IntegerVector sizes, PositionConvertor &convertor, std::vector < int > &largestSubContigs) {
    CharacterVector refNames  = sizes.names();
    IntegerVector refs = splits["ref"];
    IntegerVector bins = splits["bin"];
    int nRefs          = sizes.size();
    int nSplits        = refs.size();
    int newRef         = nRefs + 1; // factors start with 1
    int prevRef        = -1;
    int prevSplit      = 0;
    int subContig      = 0;
    if (nSplits == 0) {
        return;
    }
    convertor.resize(nRefs + 1); // factors start with 1
    for (int refId = 0; refId < nRefs; ++refId) {
        convertor[refId + 1].resize(sizes[refId] + 1); // size is the last element
        for (int bin = 0; bin <= sizes[refId]; ++bin) {
            convertor[refId + 1][bin] = std::make_pair(refId + 1, bin);
        }
    }
    for (int splitId = 0; splitId < nSplits; ++splitId) {
        int ref   = refs[splitId];
        int split = bins[splitId];
        if ((ref != prevRef) && (prevRef != -1)) {
            updateConvertor(convertor, prevRef, largestSubContigs[prevRef], newRef, subContig, prevSplit, sizes[prevRef - 1] + 1);
            prevSplit   = 0;
            subContig   = 0;
        }
        updateConvertor(convertor, ref, largestSubContigs[ref], newRef, subContig, prevSplit, split);
        prevSplit = split;
        prevRef   = ref;
        ++subContig;
    }
    updateConvertor(convertor, prevRef, largestSubContigs[prevRef], newRef, subContig, prevSplit, sizes[prevRef - 1] + 1);
}

// Update the count matrices, given the splits
DataFrame splitCountMatrices (String name, DataFrame matrices, List splits, PositionConvertor &convertor, CharacterVector newSequences) {
    IntegerVector refs1  = matrices["ref1"];
    IntegerVector refs2  = matrices["ref2"];
    IntegerVector bins1  = matrices["bin1"];
    IntegerVector bins2  = matrices["bin2"];
    IntegerVector counts = matrices["count"];
    std::string nameCpp  = name;
    long long int nCounts = refs1.size();
    std::vector < int > refs1Cpp;
    std::vector < int > refs2Cpp;
    std::vector < int > bins1Cpp;
    std::vector < int > bins2Cpp;
    std::vector < int > countsCpp;
    refs1Cpp.reserve(nCounts);
    refs2Cpp.reserve(nCounts);
    bins1Cpp.reserve(nCounts);
    bins2Cpp.reserve(nCounts);
    countsCpp.reserve(nCounts);
    Rcout << "\tDataset '" << nameCpp << "'.\n";
    Progress progress (nCounts, true);
    for (long long countId = 0; countId < nCounts; ++countId) {
        assert(refs1[countId] < static_cast < int > (convertor.size()));
        assert(bins1[countId] < static_cast < int > (convertor[refs1[countId]].size()));
        assert(refs2[countId] < static_cast < int > (convertor.size()));
        assert(bins2[countId] < static_cast < int > (convertor[refs2[countId]].size()));
        std::pair <int, int> p1 = convertor[refs1[countId]][bins1[countId]];
        std::pair <int, int> p2 = convertor[refs2[countId]][bins2[countId]];
        if ((p1.first >= 0) && (p2.first >= 0)) {
            if ((p1.first < p2.first) || ((p1.first == p2.first) && (p1.second < p2.second))) {
                refs2Cpp.push_back(p1.first);
                bins2Cpp.push_back(p1.second);
                refs1Cpp.push_back(p2.first);
                bins1Cpp.push_back(p2.second);
            }
            else {
                refs1Cpp.push_back(p1.first);
                bins1Cpp.push_back(p1.second);
                refs2Cpp.push_back(p2.first);
                bins2Cpp.push_back(p2.second);
            }
            countsCpp.push_back(counts[countId]);
            assert(refs1Cpp.back() >= refs2Cpp.back());
        }
        progress.increment();
    }
    refs1 = refs1Cpp;
    bins1 = bins1Cpp;
    refs2 = refs2Cpp;
    bins2 = bins2Cpp;
    counts = countsCpp;
    refs1.attr("class") = "factor";
    refs1.attr("levels") = newSequences;
    refs2.attr("class") = "factor";
    refs2.attr("levels") = newSequences;
    DataFrame output = DataFrame::create(_["ref1"] = refs1, _["bin1"] = bins1, _["ref2"] = refs2, _["bin2"] = bins2, _["count"] = counts);
    return output;
}


// Update outlier bins, given the splits
DataFrame splitOutlierBins (DataFrame outlierBins, PositionConvertor &convertor, CharacterVector newSequences, int metaSize, IntegerVector newSizes) {
    IntegerVector refs  = outlierBins["ref"];
    IntegerVector bins  = outlierBins["bin"];
    long long int nBins = refs.size();
    std::vector < int > refsCpp;
    std::vector < int > binsCpp;
    refsCpp.reserve(nBins);
    binsCpp.reserve(nBins);
    for (long binId = 0; binId < nBins; ++binId) {
        assert(refs[binId] < static_cast < int > (convertor.size()));
        assert(bins[binId] * metaSize < static_cast < int > (convertor[refs[binId]].size()));
        std::pair <int, int> p = convertor[refs[binId]][bins[binId] * metaSize];
        if (p.first >= 0) {
            refsCpp.push_back(p.first);
            binsCpp.push_back(p.second / metaSize);
if (p.second > newSizes[p.first-1]) Rcerr << refs[binId] << ":" << (bins[binId] * metaSize) << " -> " << p.first << ":" << p.second << "  (" << metaSize << ")" << "\n";
            assert(p.first <= newSequences.size());
        }
    }
    refs = refsCpp;
    bins = binsCpp;
    refs.attr("class") = "factor";
    refs.attr("levels") = newSequences;
    DataFrame output = DataFrame::create(_["ref"] = refs, _["bin"] = bins);
    return output;
}

void splitSequence(CharacterVector &contigs, CharacterVector &oldContigs, CharacterVector &newContigs, int ref, int largestSubContig, int subContig, int prevSplit, int split, int binSize) {
    // Exclude the split points
    if (prevSplit != 0) ++prevSplit;
    size_t length = (split == -1)? std::string::npos: (split - prevSplit) * binSize;
    prevSplit *= binSize;
    std::string contig    = as < std::string > (contigs[ref - 1]); // Factors start with 1
    std::string newContig = contig.substr(prevSplit, length);
    if (subContig == largestSubContig) {
        oldContigs[ref - 1] = newContig;
    }
    else {
        newContigs.push_back(newContig);
    }
}

// Split the contigs sequences
CharacterVector splitSequences(CharacterVector contigs, DataFrame splits, IntegerVector sizes, int binSize, std::vector < int > &largestSubContigs) {
    CharacterVector oldContigs = clone(contigs);
    CharacterVector newContigs;
    IntegerVector refs = splits["ref"];
    IntegerVector bins = splits["bin"];
    int nSplits        = refs.size();
    int subContig      = 0;
    int prevRef        = -1;
    int prevSplit      = 0;
    Rcout << "Splitting sequences.\n";
    Progress progress (nSplits, true);
    for (int splitId = 0; splitId < nSplits; ++splitId) {
        int ref   = refs[splitId];
        int split = bins[splitId];
        if ((ref != prevRef) && (prevRef != -1)) {
            splitSequence(contigs, oldContigs, newContigs, prevRef, largestSubContigs[prevRef], subContig, prevSplit, -1, binSize);
            subContig = 0;
            prevSplit = 0;
        }
        splitSequence(contigs, oldContigs, newContigs, ref, largestSubContigs[ref], subContig, prevSplit, split, binSize);
        ++subContig;
        prevRef = ref;
        progress.increment();
        prevSplit = split;
    }
    splitSequence(contigs, oldContigs, newContigs, prevRef, largestSubContigs[prevRef], subContig, prevSplit, -1, binSize);
    for (int i = 0; i < newContigs.size(); ++i) {
        String newContig = newContigs[i];
        oldContigs.push_back(newContig, "new_ref_" + std::to_string(i + 1));
    }
    if (oldContigs.size() != contigs.size() + nSplits) {
        Rcerr << "Error while splitting sequences.\n";
    }
    return oldContigs;
}

void updateSize (IntegerVector &sizes, int ref, int subContig, int largestSubContig, int prevSplit, int split, int binSize) {
    // Exclude the split points
    if (prevSplit != 0) ++prevSplit;
    --split;
    int size = split - prevSplit;
    if (largestSubContig == subContig) {
        sizes[ref - 1] = size; // factors start with 1
    }
    else {
        sizes.push_back(size);
    }
}

// Update sizes
IntegerVector updateSizes (DataFrame splits, IntegerVector sizes, std::vector < int > &largestSubContigs, int binSize) {
    IntegerVector refs     = splits["ref"];
    IntegerVector bins     = splits["bin"];
    IntegerVector newSizes = clone(sizes);
    int nSplits            = refs.size();
    int prevRef            = -1;
    int prevSplit          = 0;
    int subContig          = 0;
    for (int splitId = 0; splitId < nSplits; ++splitId) {
        int ref   = refs[splitId];
        int split = bins[splitId];
        if ((ref != prevRef) && (prevRef != -1)) {
            updateSize(newSizes, prevRef, subContig, largestSubContigs[prevRef], prevSplit, sizes[prevRef - 1] + 1, binSize);
            prevSplit   = 0;
            subContig   = 0;
        }
        updateSize(newSizes, ref, subContig, largestSubContigs[ref], prevSplit, split, binSize);
        prevSplit = split;
        prevRef   = ref;
        ++subContig;
    }
    updateSize(newSizes, prevRef, subContig, largestSubContigs[prevRef], prevSplit, sizes[prevRef - 1] + 1, binSize);
    if (newSizes.size() != sizes.size() + nSplits) {
        Rcerr << "Error while splitting sequences: expected " << sizes.size() << " + " << nSplits << ", got " << newSizes.size() << "\n";
    }
    return newSizes;
}


// Split the contigs
// [[Rcpp::export]]
S4 splitCpp(S4 object) {
    DataFrame       splits      = wrap(object.slot("breaks"));
    List            data        = object.slot("data");
    CharacterVector chromosomes = object.slot("sequences");
    IntegerVector   sizes       = object.slot("sizes");
    int             binSize     = object.slot("binSize");
    PositionConvertor convertor;
    std::vector < int > largestSubContigs;
    findLargestSubContigs(splits, sizes, largestSubContigs);
    setConvertor(splits, sizes, convertor, largestSubContigs);
    CharacterVector newSequences = splitSequences(chromosomes, splits, sizes, binSize, largestSubContigs);
    IntegerVector   newSizes     = updateSizes(splits, sizes, largestSubContigs, binSize);
    newSizes.names()             = newSequences.names();
    Rcout << "Splitting count matrices.\n";
    for (int i = 0; i < data.size(); ++i) {
        S4 object = data[i];
        String name                      = wrap(object.slot("name"));
        DataFrame matrices               = wrap(object.slot("interactionMatrix"));
        DataFrame outlierBins            = wrap(object.slot("outlierBins"));
        S4        parameters             = wrap(object.slot("parameters"));
        int       metaSize               = parameters.slot("metaSize");
        DataFrame newMatrices            = splitCountMatrices(name, matrices, splits, convertor, newSequences.names());
        DataFrame newOutlierBins         = splitOutlierBins(outlierBins, convertor, newSequences.names(), metaSize, newSizes);
        object.slot("interactionMatrix") = newMatrices;
        object.slot("outlierBins")       = newOutlierBins;
        data[i]                          = object;
    }
    object.slot("chromosomes") = newSequences.names();
    object.slot("data")        = data;
    object.slot("sequences")   = newSequences;
    object.slot("sizes")       = newSizes;
    return object;
}
