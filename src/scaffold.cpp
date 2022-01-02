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

// Group the contigs into scaffolds
// Return a list, such that:
//    the names are the ids of the longest contigs for each scaffold,
//    the values are the ids of the contigs,
//    the values are ordered (first element is left of second element)
//    the opposite of the ids of the contigs are used, if the contigs is reversed
// [[Rcpp::export]]
List getRefOrders (DataFrame &joins, List sizes) {
    CharacterVector refNames = sizes.names();
    IntegerVector references = joins["ref1"];
    IntegerVector otherRefs  = joins["ref2"];
    LogicalVector afters1    = joins["after1"];
    LogicalVector afters2    = joins["after2"];
    int nJoins               = references.size();
    int nReferences          = sizes.size();
    std::vector < bool > seen (nReferences + 1, false);
    std::vector < std::vector < int > > refNext(nReferences + 1, std::vector < int > (2, 0)); // Indices: 0 => after, 1 => before
    std::vector < int > refIdToGroupId(nReferences + 1, -1);
    std::vector < std::vector < int > > groupIdToRefIds;
    if (nJoins == 0) {
        return List();
    }
    // Create groups
    for (int joinId = 0; joinId < nJoins; ++joinId) {
        int referenceId = references[joinId];
        int otherRefId  = otherRefs[joinId];
        bool after1     = afters1[joinId];
        bool after2     = afters2[joinId];
        seen[referenceId] = true;
        seen[otherRefId]  = true;
        // Check whether the side of the main reference is not used
        if (refNext[referenceId][boolToNext(after1)] == 0) {
            // Check whether the side of the other reference is not used
            if (refNext[otherRefId][boolToNext(after2)] == 0) {
                bool isCollinear = (after1 != after2);
                // Merge group
                int groupIdRef   = refIdToGroupId[referenceId];
                int groupIdOther = refIdToGroupId[otherRefId];
// Rcout << "reading " << referenceId << " and " << otherRefId << ", " << isCollinear << ", " << groupIdRef << ", " << groupIdOther << "\n";
                // No one is in a group
                if ((groupIdRef == -1) && (groupIdOther == -1)) {
                    int newGroupId = groupIdToRefIds.size();
                    refIdToGroupId[referenceId] = newGroupId;
                    refIdToGroupId[otherRefId]  = newGroupId;
                    groupIdToRefIds.push_back({referenceId, otherRefId});
// Rcout << "  case A: " << newGroupId << "\n";
                }
                // Main ref is in a group
                else if (groupIdOther == -1) {
                    refIdToGroupId[otherRefId] = groupIdRef;
                    if (groupIdRef >= static_cast<int>(groupIdToRefIds.size())) {
                        Rcerr << "Error #1 in getRefOrder: looking for a group #" << groupIdRef << " in a vector of size " << groupIdToRefIds.size() << "\n";
                    }
                    groupIdToRefIds[groupIdRef].push_back(otherRefId);
// Rcout << "  case B: " << groupIdRef << "\n";
                }
                // Other ref is in a group
                else if (groupIdRef == -1) {
// Rcout << "  case C\n";
                    refIdToGroupId[referenceId] = groupIdOther;
                    if (groupIdOther >= static_cast<int>(groupIdToRefIds.size())) {
                        Rcerr << "Error #2 in getRefOrder: looking for a group #" << groupIdOther << " in a vector of size " << groupIdToRefIds.size() << "\n";
                    }
                    groupIdToRefIds[groupIdOther].push_back(referenceId);
// Rcout << "  case B: " << groupIdOther << "\n";
                }
                // Both refs are in a group: merge them
                // Do not merge contigs from the same group, in order to avoid cycles
                else if (groupIdRef != groupIdOther) {
                    if (groupIdRef >= static_cast<int>(groupIdToRefIds.size())) {
                        Rcerr << "Error #3 in getRefOrder: looking for a group #" << groupIdRef << " in a vector of size " << groupIdToRefIds.size() << "\n";
                    }
                    if (groupIdOther >= static_cast<int>(groupIdToRefIds.size())) {
                        Rcerr << "Error #4 in getRefOrder: looking for a group #" << groupIdOther << " in a vector of size " << groupIdToRefIds.size() << "\n";
                    }
                    for (int otherRefId: groupIdToRefIds[groupIdOther]) {
                        refIdToGroupId[otherRefId] = groupIdRef;
                    }
                    groupIdToRefIds[groupIdRef].insert(groupIdToRefIds[groupIdRef].end(), groupIdToRefIds[groupIdOther].begin(), groupIdToRefIds[groupIdOther].end());
                    groupIdToRefIds[groupIdOther].clear();
// Rcout << "  case D " << groupIdRef << " with " << groupIdToRefIds[groupIdRef].size() << "\n";
                }
                refNext[referenceId][boolToNext(after1)] = otherRefId * (isCollinear? 1: -1);
                refNext[otherRefId][boolToNext(after2)] = referenceId * (isCollinear? 1: -1);
            }
        }
    }
    // Get the (main) largest element of the group
    int nGroups = groupIdToRefIds.size();
// Rcout << nGroups << "\n";
    std::vector < int > groupLargestSizes (nGroups, 0);
    std::vector < int > groupLargestRefs  (nGroups, 0);
    for (int groupId = 0; groupId < nGroups; ++groupId) {
// Rcout << "  group #" << groupId << ":";
        for (int refId: groupIdToRefIds[groupId]) {
// Rcout << "   " << refNext[refId][0] << " " << refId << " " << refNext[refId][1];
            if (sizes[refId-1] > groupLargestSizes[groupId]) {
                groupLargestSizes[groupId] = sizes[refId-1];
                groupLargestRefs[groupId]  = refId;
            }
        }
// Rcout << "\n";
    }
//Rcpp::stop("stop");
//Rcout << "Step 1\n";
    // Follow the threads
    List groups;
    for (int groupId = 0; groupId < nGroups; ++groupId) {
        std::vector < bool > refSeen (nReferences + 1, false);
        int firstRefId = groupLargestRefs[groupId];
// Rcout << "Start with " << firstRefId << "\n";
        if (firstRefId != 0) {
            // Go forward from the main element of each group
            std::vector < int > forwardGroup;
            std::vector < int > backwardGroup;
            int nextRefId = firstRefId;
            bool isCollinear = true;
// Rcout << "Go forward: " << nextRefId << " " << isCollinear << "   ";
            while ((nextRefId != 0) && (! refSeen[nextRefId])) {
                forwardGroup.push_back(nextRefId * (isCollinear? 1: -1));
                refSeen[nextRefId] = true;
                nextRefId = refNext[nextRefId][boolToNext(isCollinear)];
                if (nextRefId < 0) {
                    nextRefId = - nextRefId;
                    isCollinear = ! isCollinear;
                }
// Rcout << nextRefId << " " << isCollinear << "   ";
            }
            // Go backward (but skip main element)
            isCollinear = true;
            nextRefId = firstRefId;
            nextRefId = refNext[nextRefId][boolToNext(!isCollinear)];
            if (nextRefId < 0) {
                nextRefId = - nextRefId;
                isCollinear = false;
            }
// Rcout << "\nGo backward: " << nextRefId << " " << isCollinear << "   ";
            while ((nextRefId != 0) && (! refSeen[nextRefId])) {
                backwardGroup.push_back(nextRefId * (isCollinear? 1: -1));
                refSeen[nextRefId] = true;
                nextRefId = refNext[nextRefId][boolToNext(!isCollinear)];
                if (nextRefId < 0) {
                    nextRefId = - nextRefId;
                    isCollinear = false;
                }
// Rcout << nextRefId << " " << isCollinear << "   ";
            }
// Rcout << "\n";
            // Merge lists
            std::reverse(backwardGroup.begin(), backwardGroup.end());
            backwardGroup.insert(backwardGroup.end(), forwardGroup.begin(), forwardGroup.end());
            IntegerVector group = wrap(backwardGroup);
            std::string refName = as<std::string>(refNames[firstRefId-1]);
            groups.push_back(group, refName);
        }
    }
    // Add unmerged contigs
    for (int refId = 1; refId <= nReferences; ++refId) {
        if (! seen[refId]) {
            IntegerVector group = { refId };
            std::string refName = as<std::string>(refNames[refId - 1]);
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

// Compute the size of the scaffolds after scaffolding
// Outputs the new sizes for the "main" contigs.
// [[Rcpp::export]]
List scaffoldSizes(List orders, IntegerVector orderIds, List sizes) {
    std::vector<int> sizesVector (sizes.size());
    std::vector<int> newSizesVector (orders.size(), 0);
    List             newSizes;
    CharacterVector  orderNames = orders.names();
    // Transform list to vector (to do arithmetics on the sizes)
    for (int sizeId = 0; sizeId < sizes.size(); ++sizeId) {
        // size is actually the last element
        sizesVector[sizeId] = as<int>(sizes[sizeId]) + 1;
    }
    // Set main contig sizes
    for (int orderId = 0; orderId < orders.size(); ++orderId) {
        IntegerVector contigIds = orders[orderId];
        for (int contigId: contigIds) {
            if (contigId < 0) contigId = - contigId;
            newSizesVector[orderId] += sizesVector[contigId - 1];
        }
    }
    // Convert back to list
    for (int orderId = 0; orderId < orders.size(); ++orderId) {
        newSizes.push_back(newSizesVector[orderId] - 1, as<std::string>(orderNames[orderId]));
    }
    return newSizes;
}

// Scaffold the contigs into scaffolds
// orders is a list of list. Each list contains the (1-based) id of a contig, * (-1) is the contig is reversed.
// Do not preserve order of the contigs
// [[Rcpp::export]]
CharacterVector scaffoldContigs(CharacterVector contigs, List orders, List sizes, int binSize) {
    CharacterVector scaffolds;
    CharacterVector chromosomeNames = contigs.names();
    int nScaffolds = orders.size();
// Rcout << "# contigs: " << contigs.size() << "\n";
    Progress p (nScaffolds, true);
    for (int scaffoldId = 0; scaffoldId < nScaffolds; ++scaffoldId) {
        IntegerVector order = orders[scaffoldId];
        std::ostringstream scaffold;
        std::string filler;
        for (int i = 0; i < order.size(); ++i) {
            int  contigId = order[i];
            int  contigSize;
// Rcout << "  contigs id: " << contigId << "\n";
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
// Rcout << "Contig size #" << contigId << ": " << contigSize << " vs " << as<int>(sizes[contigId-1]) << "\n";
            int nBins = sizes[contigId-1];
            int nMissingNts = (nBins+1) * binSize - contigSize;
// Rcout << "# missing nts: " << nMissingNts << "\n";
            if (nMissingNts < 0) {
                CharacterVector refNames = contigs.names();
                std::string refName = Rcpp::as<std::string>(refNames[contigId - 1]);
                Rcerr << "Problem while creating filler for " << refName << ".  Got contig with " << nBins << " bins (bin size " << binSize << "), but the contig size is " << contigSize << "\n";
            }
            filler = std::string(nMissingNts, 'N');
        }
// Rcout << "  scaffold id: " << scaffoldId << "\n";
        scaffolds.push_back(scaffold.str(), as<std::string>(as<CharacterVector>(orders.names())[scaffoldId]));
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
