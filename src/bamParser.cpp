#include <cstdio>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <thread>
#include <atomic>
#include <mutex>
#include "sam.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

inline uint64_t getKey (const uint32_t i, const uint32_t j) {
    return (static_cast<uint64_t>(i) << 32) | j;
}

inline uint32_t getFirstKey (const uint64_t k) {
    return static_cast<uint32_t>(k >> 32);
}

inline uint32_t getSecondKey (const uint64_t k) {
    return static_cast<uint32_t>(k & 0xFFFFFFFF);
} 

struct position_t {
    int      chrId;
    uint32_t pos;
    position_t (): chrId(0), pos(0) {}
    position_t (int c, uint32_t p): chrId(c), pos(p) {}
    position_t (const position_t &p): chrId(p.chrId), pos(p.pos) {}
};

bool operator<(const position_t &lhs, const position_t &rhs) {
    if (lhs.chrId < rhs.chrId) return true;
    if (lhs.chrId > rhs.chrId) return false;
    return (lhs.pos < rhs.pos);
}

bool operator==(const position_t &lhs, const position_t &rhs) {
    return ((lhs.chrId == rhs.chrId) && (lhs.pos == rhs.pos));
}

bool operator!=(const position_t &lhs, const position_t &rhs) {
    return ((lhs.chrId != rhs.chrId) || (lhs.pos != rhs.pos));
}

std::ostream& operator<<(std::ostream& os, position_t const &position) {
    return os << "(" << &position << ") " << position.chrId << ":" << position.pos;
}


using positions_t      = std::vector < position_t >;
using umiMap_t         = std::unordered_map<std::string, positions_t >;
using matrixCounts_t   = std::vector < unsigned long >;
using validMatrices_t  = std::vector < bool >;
using validMatrixIds_t = std::vector < size_t >;
using matrix_t         = std::unordered_map < uint64_t, unsigned long >;
using matrices_t       = std::vector < matrix_t >;

void dump(matrices_t          &matrices,
          int                  chrId1, 
          std::vector < int > &stdChrs1,
          std::vector < int > &stdChrs2,
          std::vector < int > &stdBins1,
          std::vector < int > &stdBins2,
          std::vector < int > &stdCounts,
          unsigned int         minCount) {
    for (size_t chrId2 = 0; chrId2 < matrices.size(); ++chrId2) {
        matrix_t &matrix = matrices[chrId2];
        for (auto it = matrix.begin(); it != matrix.end(); ) {
            if (it->second < minCount) {
                it = matrix.erase(it);
            }
            else {
                it++;
            }
        }
        size_t nElements = matrix.size();
        std::vector < int > tmpChrs1(nElements, chrId1);
        std::vector < int > tmpChrs2(nElements, chrId2);
        std::vector < int > tmpBins1, tmpBins2, tmpCounts;
        tmpBins1.reserve(nElements);
        tmpBins2.reserve(nElements);
        tmpCounts.reserve(nElements);
        for (auto it = matrix.begin(); it != matrix.end(); ++it) {
            tmpBins1.push_back(getFirstKey(it->first));
            tmpBins2.push_back(getSecondKey(it->first));
            tmpCounts.push_back(it->second);
        }
        stdChrs1.insert(stdChrs1.end(), tmpChrs1.begin(), tmpChrs1.end());
        stdChrs2.insert(stdChrs2.end(), tmpChrs2.begin(), tmpChrs2.end());
        stdBins1.insert(stdBins1.end(), tmpBins1.begin(), tmpBins1.end());
        stdBins2.insert(stdBins2.end(), tmpBins2.begin(), tmpBins2.end());
        stdCounts.insert(stdCounts.end(), tmpCounts.begin(), tmpCounts.end());
    }
    matrices.clear();
}

// [[Rcpp::export]]
DataFrame parseBamFileCpp(String fileName, int32_t binSize) {
    unsigned int minCount = 2;
    Rcout << "Reading " << fileName.get_cstring() << "\n";
    samFile   *inputFile  = hts_open(fileName.get_cstring(), "r");
    DataFrame emptyDataFrame;
    if (inputFile == NULL) {
        Rcerr << "File does not exist.  Exiting.\n";
        return emptyDataFrame;
    }
    bam_hdr_t *header     = sam_hdr_read(inputFile);
    bam1_t    *alignment  = bam_init1();
    int        nChrs      = header->n_targets;
    Rcout << "Found " << nChrs << " references\n";
    CharacterVector chrs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        chrs[i] = header->target_name[i];
    }
    umiMap_t umiMap;
    Rcout << "Reading BAM file.\n";
    unsigned int nAlignments;
    int previousChrId = nChrs;
    std::vector < int > stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts;
    matrices_t matrices;
    for (nAlignments = 0; sam_read1(inputFile, header, alignment) > 0; ++nAlignments) {
        uint32_t pos    = alignment->core.pos + 1;
        uint32_t posBin = pos / binSize;
        int chrId       = alignment->core.tid;
        if (chrId > nChrs) stop("Found an out-of-bound chr ID.");
        uint8_t *umi = bam_aux_get(alignment, "BX");
        if (chrId != previousChrId) {
            dump(matrices, previousChrId, stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts, minCount);
            matrices.resize(chrId + 1);
            previousChrId = chrId;
        }
        if ((chrId != -1) && (umi != NULL)) {
            std::string umiString((char *) ++umi);
            position_t position(chrId, posBin);
            if ((umiMap[umiString].empty()) || (umiMap[umiString].back() != position)) {
                for (position_t &previousPosition: umiMap[umiString]) {
                    if (previousPosition.chrId > chrId) {
                        Rcpp::stop("Input BAM file is not sorted: ref #" +
                            std::to_string(previousPosition.chrId) +
                            " (" + chrs[previousPosition.chrId] +
                            ") is seen before ref #" + std::to_string(chrId) +
                            " (" + chrs[chrId] + ").");
                    }
                    ++matrices[previousPosition.chrId][getKey(posBin, previousPosition.pos)];
                }
                umiMap[umiString].push_back(position);
            }
        }
        if ((nAlignments > 0) && (nAlignments % 10000000 == 0)) {
            Rcout << "\tRead " << nAlignments << " reads.\n";
        }
    }
    Rcout << "\tRead " << nAlignments << " reads.\n";
    dump(matrices, previousChrId, stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts, minCount);
    bam_destroy1(alignment);
    sam_close(inputFile);
    Rcout << "Read parsing done with " << chrs.size() << " references and " << nAlignments << " alignments.  Placing to buckets.\n";
    umiMap.clear();
    IntegerVector chrs1, chrs2;
    IntegerVector bins1, bins2, counts;
    Rcout << "Exported " << stdChrs1.size() << " points.\n";
    chrs1 = stdChrs1;
    chrs2 = stdChrs2;
    bins1 = stdBins1;
    bins2 = stdBins2;
    counts = stdCounts;
    chrs1 = chrs1 + 1; // factors in R start with 1
    chrs2 = chrs2 + 1;
    chrs1.attr("class") = "factor";
    chrs2.attr("class") = "factor";
    chrs1.attr("levels") = chrs;
    chrs2.attr("levels") = chrs;
    return DataFrame::create(_["ref1"]  = chrs1,
                             _["bin1"]  = bins1,
                             _["ref2"]  = chrs2,
                             _["bin2"]  = bins2,
                             _["count"] = counts);
}
