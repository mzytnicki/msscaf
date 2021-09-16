#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


// Transform the isNext boolean to a vector offset
inline int boolToNext(bool isNext) {
    return (isNext? 0: 1);
}

// [[Rcpp::export]]
List getRefOrders (DataFrame &joins, List sizes) {
    CharacterVector refNames = sizes.names();
    IntegerVector references = joins["ref1"];
    IntegerVector otherRefs  = joins["ref2"];
    LogicalVector afters1    = joins["after1"];
    LogicalVector afters2    = joins["after2"];
    int nJoins = references.size();
    if (nJoins == 0) {
        return List();
    }
    int nReferences = std::max(*std::max_element(references.begin(), references.end()), *std::max_element(otherRefs.begin(), otherRefs.end())); // Starts with 1
    std::vector < std::vector < int > > refNext(nReferences + 1, std::vector < int > (2, 0)); // Indices: 0 => after, 1 => before
    std::vector < int > refIdToGroupId(nReferences + 1, -1);
    std::vector < std::vector < int > > groupIdToRefIds;
    // Create groups
    for (int joinId = 0; joinId < nJoins; ++joinId) {
        int referenceId  = references[joinId];
        int otherRefId  = otherRefs[joinId];
        bool after1     = afters1[joinId];
        bool after2     = afters2[joinId];
        // Check whether the side of the main reference is not used
        if (refNext[referenceId][boolToNext(after1)] == 0) {
            // Check whether the side of the other reference is not used
            if (refNext[otherRefId][boolToNext(after2)] == 0) {
                bool isCollinear = (after1 != after2);
                // Merge group
                int groupIdRef   = refIdToGroupId[referenceId];
                int groupIdOther = refIdToGroupId[otherRefId];
                // No one is in a group
                if ((groupIdRef == -1) && (groupIdOther == -1)) {
                    int newGroupId = groupIdToRefIds.size();
                    refIdToGroupId[referenceId] = newGroupId;
                    refIdToGroupId[otherRefId]  = newGroupId;
                    groupIdToRefIds.push_back({referenceId, otherRefId});
                }
                // Main ref is in a group
                else if (groupIdOther == -1) {
                    refIdToGroupId[otherRefId] = groupIdRef;
                    groupIdToRefIds[groupIdRef].push_back(otherRefId);
                }
                // Other ref is in a group
                else if (groupIdRef == -1) {
                    refIdToGroupId[referenceId] = groupIdOther;
                    groupIdToRefIds[groupIdOther].push_back(referenceId);
                }
                // Both refs are in a group: merge them
                else {
                    if (groupIdRef == groupIdOther) {
                        Rcpp::stop("Problem while merging groups.");
                    }
                    for (int otherRefId: groupIdToRefIds[groupIdOther]) {
                        refIdToGroupId[otherRefId] = groupIdRef;
                    }
                    groupIdToRefIds[groupIdRef].insert(groupIdToRefIds[groupIdRef].end(), groupIdToRefIds[otherRefId].begin(), groupIdToRefIds[otherRefId].end());
                    groupIdToRefIds[otherRefId].clear();
                }
                refNext[referenceId][boolToNext(after1)] = otherRefId * (isCollinear? 1: -1);
                refNext[otherRefId][boolToNext(after2)] = referenceId * (isCollinear? 1: -1);
            }
        }
    }
    // Get the (main) largest element of the group
    int nGroups = groupIdToRefIds.size();
    std::vector < int > groupLargestSizes (nGroups, 0);
    std::vector < int > groupLargestRefs  (nGroups, 0);
    for (int groupId = 0; groupId < nGroups; ++groupId) {
        for (int refId: groupIdToRefIds[groupId]) {
            if (sizes[refId-1] > groupLargestSizes[groupId]) {
                groupLargestSizes[groupId] = sizes[refId-1];
                groupLargestRefs[groupId]  = refId;
            }
        }
    }
    // Follow the threads
    List groups;
    for (int groupId = 0; groupId < nGroups; ++groupId) {
        int firstRefId = groupLargestRefs[groupId];
        if (firstRefId != 0) {
            // Go forward from the main element of each group
            std::vector < int > forwardGroup;
            std::vector < int > backwardGroup;
            int nextRefId = firstRefId;
            bool isCollinear = true;
            while (nextRefId != 0) {
                forwardGroup.push_back(nextRefId * (isCollinear? 1: -1));
                if (nextRefId < 0) {
                    nextRefId = - nextRefId;
                    isCollinear = ! isCollinear;
                }
                nextRefId = refNext[nextRefId][boolToNext(isCollinear)];
            }
            // Go backward (but skip main element)
            isCollinear = true;
            nextRefId = firstRefId;
            nextRefId = refNext[nextRefId][boolToNext(!isCollinear)];
            while (nextRefId != 0) {
                backwardGroup.push_back(nextRefId * (isCollinear? 1: -1));
                if (nextRefId < 0) {
                    nextRefId = - nextRefId;
                    isCollinear = ! isCollinear;
                }
                nextRefId = refNext[nextRefId][boolToNext(!isCollinear)];
            }
            // Merge lists
            std::reverse(backwardGroup.begin(), backwardGroup.end());
            backwardGroup.insert(backwardGroup.end(), forwardGroup.begin(), forwardGroup.end());
            IntegerVector group = wrap(backwardGroup);
            std::string refName = as<std::string>(refNames[firstRefId-1]);
            groups.push_back(group, refName);
        }
    }
    return groups;
}

char complement(char c) {
    switch (std::toupper(c)) {
        case 'A': return 'T';
        case 'B': return 'V';
        case 'C': return 'G';
        case 'D': return 'H';
        case 'G': return 'C';
        case 'H': return 'D';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'N': return 'N';
        case 'R': return 'Y';
        case 'S': return 'S';
        case 'T': case 'U': return 'A';
        case 'V': return 'B';
        case 'W': return 'W';
        case 'Y': return 'R';
    }
    Rcpp::stop("Problem in complement.");
    return ' ';
}

void reverseComplement(std::string &s) {
    std::transform(begin(s), end(s), begin(s), complement);
}

// Scaffold the contigs into scaffolds
// orders is a list of list. Each list contains the (1-based) id of a contig, * (-1) is the contig is reversed.
// [[Rcpp::export]]
CharacterVector scaffoldContigs(CharacterVector contigs, List orders, List sizes, int binSize) {
    CharacterVector scaffolds;
    int nScaffolds = orders.size();
    Progress p (nScaffolds, true);
    for (int scaffoldId = 0; scaffoldId < nScaffolds; ++scaffoldId) {
        IntegerVector order = orders[scaffoldId];
        std::ostringstream scaffold;
        std::string filler;
        for (int i = 0; i < order.size(); ++i) {
            int  contigId = order[i];
            int  contigSize;
            if (! filler.empty()) {
                scaffold << filler;
            }
            if (contigId < 0) {
                contigId = - contigId;
                std::string contig = Rcpp::as<std::string>(contigs[contigId - 1]); // contig id are factors (start with 1)
                reverseComplement(contig);
                scaffold << contig;
                contigSize = contig.size();
            }
            else {
                std::string contig = Rcpp::as<std::string>(contigs[contigId - 1]); // contig id are factors (start with 1)
                scaffold << contig;
                contigSize = contig.size();
            }
            int nBins = sizes[contigId-1];
            int nMissingNts = (nBins+1) * binSize - contigSize;
            filler = std::string(nMissingNts, 'N');
        }
        scaffolds.push_back(scaffold.str());
        p.increment();
    }
    return scaffolds;
}

// Scaffold the count matrices
// orders is a list of list. Each list contains the (1-based) id of a contig, * (-1) is the contig is reversed.
// [[Rcpp::export]]
DataFrame scaffoldCounts(DataFrame matrices, List groups, IntegerVector scaffoldRefs, List sizes) {
    int nScaffolds = groups.size();
    if (scaffoldRefs.size() != nScaffolds) {
        Rcpp::stop("Problem in group sizes.");
    }
    int nRefs = sizes.size();
    std::vector < int > newRefs (nRefs + 1); // Indexes start with 1 in R
    std::vector < int > offsets (nRefs + 1, 0);
    std::vector < bool > forwards (nRefs + 1, true);
    for (int refId = 0; refId <= nRefs; ++refId) {
        newRefs[refId] = refId;
    }
    for (int scaffoldId = 0; scaffoldId < nScaffolds; ++scaffoldId) {
        int size = 0;
        IntegerVector scaffold = groups[scaffoldId];
        for (int i = 0; i < scaffold.size(); ++i) {
            int contigId = scaffold[i];
            bool forward = true;
            if (contigId < 0) {
                forward = false;
                contigId = - contigId;
            }
            newRefs[contigId] = scaffoldRefs[scaffoldId];
            offsets[contigId] = size;
            forwards[contigId] = forward;
            size += sizes[contigId-1];
        }
    }
    IntegerVector refs1  = matrices["ref1"];
    IntegerVector refs2  = matrices["ref2"];
    IntegerVector bins1  = matrices["bin1"];
    IntegerVector bins2  = matrices["bin2"];
    IntegerVector counts = matrices["count"];
    long long int nCounts = refs1.size();
    Progress p (nCounts, true);
    for (long long countId = 0; countId < nCounts; ++countId) {
        bins1[countId] = forwards[refs1[countId]]? bins1[countId] + offsets[refs1[countId]]: sizes[refs1[countId]-1] - bins1[countId] + offsets[refs1[countId]];
        bins2[countId] = forwards[refs2[countId]]? bins2[countId] + offsets[refs2[countId]]: sizes[refs2[countId]-1] - bins2[countId] + offsets[refs2[countId]];
        refs1[countId] = newRefs[refs1[countId]];
        refs2[countId] = newRefs[refs2[countId]];
        p.increment();
    }
    DataFrame output = DataFrame::create(_["ref1"] = refs1, _["bin1"] = bins1, _["ref2"] = refs2, _["bin2"] = bins2, _["count"] = counts);
    return output;
}
