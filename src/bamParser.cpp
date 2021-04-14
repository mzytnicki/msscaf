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
//#include "gperftools/profiler.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//const long maxRecords = 1000000000;

int getMedian(std::vector<int> &v) {
    if (v.size() == 0) return 0;
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}                                                                                                                                                                                                                  

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
//using matrices_t       = std::unordered_map < uint64_t, matrix_t >;
using matrices_t       = std::vector < matrix_t >;
using chrIndices_t     = std::vector < std::pair < size_t, size_t > >;
using indicesChr_t     = std::vector < std::vector < size_t > >;

/*
int getMoleculeSize (std::vector < positions_t > &positionss) {
    std::vector < int > sizes;
    Rcout << "Computing molecule sizes.\n";
    Progress progress(positionss.size(), true);
    for (auto &positions: positionss) {
        int maxSize = 0;
        int chrId   = -1;
        int start   = 0;
        int end     = 0;
        for (auto &position: positions) {
            if (position.chrId == chrId) {
                end = position.pos;
            }
            else {
                maxSize = std::max < int > (maxSize, end - start + 1);
                chrId   = position.chrId;
                start   = position.pos;
                end     = position.pos;
            }
        }
        maxSize = std::max < int > (maxSize, end - start + 1);
        sizes.push_back(maxSize);
        progress.increment();
    }
    return getMedian(sizes);
}
*/

int getMoleculeSize (umiMap_t &umiMap) {
    std::vector < int > sizes;
    Rcout << "Computing molecule sizes.\n";
    Progress progress(umiMap.size(), true);
    for (auto &it: umiMap) {
        auto &positions = it.second;
        int maxSize = 0;
        int chrId   = -1;
        int start   = 0;
        int end     = 0;
        for (auto &position: positions) {
            if (position.chrId == chrId) {
                end = position.pos;
            }
            else {
                maxSize = std::max < int > (maxSize, end - start + 1);
                chrId   = position.chrId;
                start   = position.pos;
                end     = position.pos;
            }
        }
        maxSize = std::max < int > (maxSize, end - start + 1);
        sizes.push_back(maxSize);
        progress.increment();
    }
    return getMedian(sizes);
}

void updateMatrixCounts(matrixCounts_t &matrixCounts, positions_t &positions, indicesChr_t &indicesChr) {
    for (size_t positionId1 = 0; positionId1 < positions.size(); ++positionId1) {
        position_t &position1 = positions[positionId1];
        for (size_t positionId2 = 0; positionId2 <= positionId1; ++positionId2) {
            position_t &position2 = positions[positionId2];
            ++matrixCounts[indicesChr[position1.chrId][position2.chrId]];
        }
    }
}

void updateMatrices(matrices_t &matrices, validMatrices_t &validMatrices, validMatrices_t &validMatricesPos1, positions_t &positions, indicesChr_t &indicesChr) {
    for (size_t positionId1 = 0; positionId1 < positions.size(); ++positionId1) {
        position_t &position1 = positions[positionId1];
        if (validMatricesPos1[position1.chrId]) {
            //if (position1.start > position1.end) Rcerr << "Problem in 'updateMatrices': [" << positionId1 << "] " << position1 << "\n";
            for (size_t positionId2 = 0; positionId2 <= positionId1; ++positionId2) {
                position_t &position2 = positions[positionId2];
                //if (position1.chrId < position2.chrId) Rcerr << "Problem in 'updateMatrices': [" << positionId1 << "] " << position1 << ", [" << positionId2 << "] " << position2 << "\n";
                uint64_t chrId = indicesChr[position1.chrId][position2.chrId];
                if (validMatrices[chrId]) {
                    ++matrices[chrId][getKey(position1.pos, position2.pos)];
                }
            }
        }
    }
}

/*
void dump(matrices_t          &matrices,
          chrIndices_t        &chrIndices,
          std::vector < int > &stdChrs1,
          std::vector < int > &stdChrs2,
          std::vector < int > &stdBins1,
          std::vector < int > &stdBins2,
          std::vector < int > &stdCounts,
          unsigned int         minCount) {
    Rcout << "\tExporting matrices.\n";
    Progress progress(matrices.size(), true);
    size_t chrId1, chrId2;
    for (auto &it: matrices) {
        size_t    chrId  = it.first;
        matrix_t &matrix = it.second;
        for (auto it = matrix.begin(); it != matrix.end(); ) {
            if (it->second < minCount) {
                it = matrix.erase(it);
            }
            else {
                it++;
            }
        }
        size_t nElements = matrix.size();
        std::tie(chrId1, chrId2) = chrIndices[chrId];
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
        progress.increment();
    }
    matrices.clear();
}
*/

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

/*
// [[Rcpp::export]]
List parseBamFileCpp(String fileName, int32_t binSize) {
    unsigned int minCount = 2;
    unsigned int minCountPerMatrix = 100;
    unsigned int nMatricesPerBatch = 100;
    //ProfilerStart("/tmp/profile.out");
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
    std::vector <int> sizes(header->n_targets);
    CharacterVector chrs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        chrs[i] = header->target_name[i];
        sizes[i] = header->target_len[i] / binSize + 1;
    }
    umiMap_t umiMap;
    Rcout << "Reading BAM file.\n";
    unsigned int nAlignments;
    for (nAlignments = 0; sam_read1(inputFile, header, alignment) > 0; ++nAlignments) {
        uint32_t pos    = alignment->core.pos + 1;
        uint32_t posBin = pos / binSize;
        int chrId       = alignment->core.tid;
        if (chrId > nChrs) Rcerr << "Found an out-of-bound chr ID.\n";
        uint8_t *umi = bam_aux_get(alignment, "BX");
        if ((chrId != -1) && (umi != NULL)) {
            std::string umiString((char *) ++umi);
            umiMap[umiString].emplace_back(chrId, posBin);
        }
        if ((nAlignments > 0) && (nAlignments % 10000000 == 0)) {
            Rcout << "\tRead " << nAlignments << " reads.\n";
        }
    }
    bam_destroy1(alignment);
    sam_close(inputFile);
    Rcout << "Read parsing done with " << sizes.size() << " references and " << nAlignments << " alignments.  Placing to buckets.\n";
    chrIndices_t chrIndices (nChrs * (nChrs + 1) / 2);
    indicesChr_t indicesChr (nChrs);
    size_t nMatrices = 0;
    for (int i = 0; i < nChrs; ++i) {
        indicesChr[i] = std::vector < size_t > (i+1);
        for (int j = 0; j <= i; ++j) {
            chrIndices[nMatrices] = {i, j};
            indicesChr[i][j] = nMatrices;
            ++nMatrices;
        }
    }
    Progress progress(umiMap.size(), true);
    std::vector < positions_t > positionss;
    positionss.reserve(umiMap.size());
    for (auto &mapElement: umiMap) {
        positions_t &positions = mapElement.second;
        std::sort(positions.begin(), positions.end());
        auto it = std::unique(positions.begin(), positions.end());
        positions.resize(std::distance(positions.begin(), it));
        positionss.push_back(positions);
        progress.increment();
    }
    umiMap.clear();
    int moleculeSize = getMoleculeSize(positionss);
    Rcout << "Finding non-sparse matrices.\n";
    Progress progress1(positionss.size(), true);
    matrixCounts_t matrixCounts (nMatrices);
    for (auto &positions: positionss) {
        updateMatrixCounts(matrixCounts, positions, indicesChr);
        progress1.increment();
    }
    validMatrices_t validMatrices (nMatrices);
    for (size_t i = 0; i < matrixCounts.size(); ++i) {
        validMatrices[i] = (matrixCounts[i] >= minCountPerMatrix);
    }
    size_t nValidMatrices = std::count(validMatrices.begin(), validMatrices.end(), true);
    validMatrixIds_t validMatrixIds;
    validMatrixIds.reserve(nValidMatrices);
    for (size_t i = 0; i < matrixCounts.size(); ++i) {
        if (validMatrices[i]) {
            validMatrixIds.push_back(i);
        }
    }
    Rcout << "Found " << nValidMatrices << "/" << nMatrices << " valid matrices.  Filling them.\n";
    std::vector < int > stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts;
    matrices_t matrices (nMatrices);
    for (unsigned int batchId = 0; batchId < (nValidMatrices / nMatricesPerBatch) + 1; ++batchId) {
        Rcout << "Batch " << (batchId + 1) << " / " << ((nValidMatrices / nMatricesPerBatch) + 1) << "\n";
        auto itStart = validMatrixIds.begin() + batchId * nMatricesPerBatch;
        auto itEnd   = ((batchId+1) * nMatricesPerBatch >= validMatrixIds.size())? validMatrixIds.end(): validMatrixIds.begin() + (batchId+1) * nMatricesPerBatch;
        validMatrixIds_t validMatrixIdsBatch (itStart, itEnd);
        validMatrices_t validMatricesBatch     (nMatrices, false);
        validMatrices_t validMatricesBatchPos1 (nChrs,     false);
        for (size_t id: validMatrixIdsBatch) {
            int id1, id2;
            validMatricesBatch[id] = true;
            std::tie(id1, id2) = chrIndices[id];
            validMatricesBatchPos1[id1] = true;
        }
        Rcout << "\tFilling matrices.\n";
        Progress progress2(positionss.size(), true);
        for (auto &positions: positionss) {
            updateMatrices(matrices, validMatricesBatch, validMatricesBatchPos1, positions, indicesChr);
            progress2.increment();
        }
        dump(matrices, chrIndices, stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts, minCount);
    }
    positionss.clear();
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
    IntegerVector chrSizes;
    for (size_t i = 0; i < sizes.size(); ++i) {
        chrSizes.push_back(sizes[i], as<std::string>(chrs[i]));
    }
    DataFrame outputDataFrame = DataFrame::create(_["ref1"]  = chrs1,
                                                  _["bin1"]  = bins1,
                                                  _["ref2"]  = chrs2,
                                                  _["bin2"]  = bins2,
                                                  _["count"] = counts);
    return List::create(_["data"]  = outputDataFrame,
                        _["size"]  = moleculeSize,
                        _["sizes"] = chrSizes);
}
*/


// [[Rcpp::export]]
List parseBamFileCpp(String fileName, int32_t binSize) {
    unsigned int minCount = 2;
    //ProfilerStart("/tmp/profile.out");
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
    std::vector <int> sizes(header->n_targets);
    CharacterVector chrs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        chrs[i] = header->target_name[i];
        sizes[i] = header->target_len[i] / binSize + 1;
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
        if (chrId > nChrs) Rcerr << "Found an out-of-bound chr ID.\n";
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
                        Rcpp::stop("Input BAM file is not sorted.");
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
    dump(matrices, previousChrId, stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts, minCount);
    bam_destroy1(alignment);
    sam_close(inputFile);
    Rcout << "Read parsing done with " << sizes.size() << " references and " << nAlignments << " alignments.  Placing to buckets.\n";
    chrIndices_t chrIndices (nChrs * (nChrs + 1) / 2);
    indicesChr_t indicesChr (nChrs);
    size_t nMatrices = 0;
    for (int i = 0; i < nChrs; ++i) {
        indicesChr[i] = std::vector < size_t > (i+1);
        for (int j = 0; j <= i; ++j) {
            chrIndices[nMatrices] = {i, j};
            indicesChr[i][j] = nMatrices;
            ++nMatrices;
        }
    }
    int moleculeSize = getMoleculeSize(umiMap);
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
    IntegerVector chrSizes;
    for (size_t i = 0; i < sizes.size(); ++i) {
        chrSizes.push_back(sizes[i], as<std::string>(chrs[i]));
    }
    DataFrame outputDataFrame = DataFrame::create(_["ref1"]  = chrs1,
                                                  _["bin1"]  = bins1,
                                                  _["ref2"]  = chrs2,
                                                  _["bin2"]  = bins2,
                                                  _["count"] = counts);
    return List::create(_["data"]  = outputDataFrame,
                        _["size"]  = moleculeSize,
                        _["sizes"] = chrSizes);
}
