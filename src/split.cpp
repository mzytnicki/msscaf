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
    int newRef         = nRefs;
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
            convertor[refId+1][bin] = std::make_pair(refId+1, bin);
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
DataFrame splitCountMatrices (String name, DataFrame matrices, List splits, IntegerVector sizes, PositionConvertor &convertor, CharacterVector newSequences) {
    IntegerVector refs1  = matrices["ref1"];
    IntegerVector refs2  = matrices["ref2"];
    IntegerVector bins1  = matrices["bin1"];
    IntegerVector bins2  = matrices["bin2"];
    IntegerVector counts = matrices["count"];
    long long int nCounts = refs1.size();
    Rcout << "\tDataset '" << wrap(name) << "'.\n";
    Progress progress (nCounts, true);
    for (long long countId = 0; countId < nCounts; ++countId) {
        std::pair <int, int> p;
        p = convertor[refs1[countId]][bins1[countId]];
        refs1[countId] = p.first;
        bins1[countId] = p.second;
        p = convertor[refs2[countId]][bins2[countId]];
        refs2[countId] = p.first;
        bins2[countId] = p.second;
        progress.increment();
        if (refs1[countId] < refs2[countId]) {
            std::swap < int > (refs1[countId], refs2[countId]);
            std::swap < int > (bins1[countId], bins2[countId]);
        }
        else if ((refs1[countId] == refs2[countId]) && (bins1[countId] < bins2[countId])) {
            std::swap < int > (bins1[countId], bins2[countId]);
        }
    }
    refs1.attr("class") = "factor";
    refs1.attr("levels") = newSequences;
    refs2.attr("class") = "factor";
    refs2.attr("levels") = newSequences;
    DataFrame output = DataFrame::create(_["ref1"] = refs1, _["bin1"] = bins1, _["ref2"] = refs2, _["bin2"] = bins2, _["count"] = counts);
    return output;
}

void splitSequence(CharacterVector &contigs, CharacterVector &oldContigs, CharacterVector &newContigs, int ref, int largestSubContig, int subContig, int prevSplit, int split) {
    size_t lastPos = (split == -1)? std::string::npos: split;
    std::string contig    = as < std::string > (contigs[ref - 1]); // Factors start with 1
    std::string newContig = contig.substr(prevSplit, lastPos);
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
            splitSequence(contigs, oldContigs, newContigs, prevRef, largestSubContigs[prevRef], subContig, prevSplit, -1);
            subContig = 0;
            prevSplit = 0;
        }
        splitSequence(contigs, oldContigs, newContigs, ref, largestSubContigs[ref], subContig, prevSplit, split);
        ++subContig;
        prevRef = ref;
        progress.increment();
    }
    splitSequence(contigs, oldContigs, newContigs, prevRef, largestSubContigs[prevRef], subContig, prevSplit, -1);
    for (int i = 0; i < newContigs.size(); ++i) {
        String newContig = newContigs[i];
        oldContigs.push_back(newContig, "new_ref_" + std::to_string(i + 1));
    }
    if (oldContigs.size() != contigs.size() + nSplits) {
        Rcerr << "Error while splitting sequences.\n";
    }
    return oldContigs;
}

void updateSize (IntegerVector &sizes, int ref, int subContig, int largestSubContig, int prevSplit, int split) {
    int size = split - prevSplit;
    if (largestSubContig == subContig) {
        sizes[ref - 1] = size; // factors start with 1
    }
    else {
        sizes.push_back(size);
    }
}

// Update sizes
IntegerVector updateSizes (DataFrame splits, IntegerVector sizes, std::vector < int > &largestSubContigs) {
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
            updateSize(newSizes, prevRef, subContig, largestSubContigs[prevRef], prevSplit, sizes[prevRef - 1]);
            prevSplit   = 0;
            subContig   = 0;
        }
        updateSize(newSizes, ref, subContig, largestSubContigs[ref], prevSplit, split);
        prevSplit = split;
        prevRef   = ref;
        ++subContig;
    }
    updateSize(newSizes, prevRef, subContig, largestSubContigs[prevRef], prevSplit, sizes[prevRef - 1]);
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
    IntegerVector   newSizes     = updateSizes(splits, sizes, largestSubContigs);
    newSizes.names()             = newSequences.names();
    for (int i = 0; i < data.size(); ++i) {
        S4 object = data[i];
        Rcout << "Splitting count matrices.\n"; 
        String name                      = wrap(object.slot("name"));
        DataFrame matrices               = wrap(object.slot("interactionMatrix"));
        DataFrame newMatrices            = splitCountMatrices (name, matrices, splits, sizes, convertor, newSequences.names());
        object.slot("interactionMatrix") = newMatrices;
        data[i]                          = object;
    }
    object.slot("data")      = data;
    object.slot("sequences") = newSequences;
    object.slot("sizes")     = newSizes;
    return object;
}
