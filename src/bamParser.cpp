#include <cstdio>
#include <string>
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

std::ostream& operator<<(std::ostream& os, position_t const &position) {
    return os << "(" << &position << ")" << position.chrId << ":" << position.pos;
}


/*
class SparseMatrix {
public:
    size_t nRows;
    size_t nCols;
private:
    std::unordered_map<size_t, unsigned int> matrix;
    size_t translateCoordinates (size_t r, size_t c) {
        return nRows * r + c;
    }
public:
    void setSizes (size_t nr, size_t nc) {
        nRows = nr;
        nCols = nc;
    }
    SparseMatrix (size_t nr = 0, size_t nc = 0) {
        setSizes(nr, nc);
    }
    void addElement (size_t r, size_t c) {
        ++matrix[translateCoordinates(r, c)];
    }
    size_t getNElements () {
        return matrix.size();
    }
    void clear () {
        matrix.clear();
    }
    void writeToFile (std::ofstream &outFile) {
        size_t size = matrix.size();
        outFile.write(reinterpret_cast < char const * > (&size), sizeof(size_t));
        for (auto &cell: matrix) {
            outFile.write(reinterpret_cast < char const * > (&cell.first), sizeof(size_t));
            outFile.write(reinterpret_cast < char const * > (&cell.second), sizeof(unsigned int));
        }
        matrix.clear();
    }
    void loadFromFile (std::ifstream &inFile) {
        size_t size;
        size_t key;
        unsigned int value;
        inFile.read(reinterpret_cast < char * > (&size), sizeof(size_t));
        //printf("  reading %zu elements\n", size);
        matrix.reserve(size);
        for (size_t i = 0; i < size; ++i) {
            inFile.read(reinterpret_cast < char * > (&key),   sizeof(size_t));
            inFile.read(reinterpret_cast < char * > (&value), sizeof(unsigned int));
            matrix[key] = value;
        }
    }
    void merge (SparseMatrix &otherMatrix) {
        for (auto &cell: otherMatrix.matrix) {
            matrix[cell.first] += cell.second;
        }
    }
    unsigned int convertToPlain (std::vector < int > &bin1, std::vector < int > &bin2, std::vector < int > &count, unsigned int minCount) {
        if (nRows == 0) {
            printf("Error, #rows = 0\n");
        }
        unsigned int nAdded = 0;
        std::vector < int > tmpBin1, tmpBin2, tmpCount;
        tmpBin1.reserve(matrix.size());
        tmpBin2.reserve(matrix.size());
        tmpCount.reserve(matrix.size());
        for (auto &cell: matrix) {
            if (cell.second >= minCount) {
                tmpBin1.push_back(cell.first / nRows);
                tmpBin2.push_back(cell.first % nRows);
                tmpCount.push_back(cell.second);
                ++nAdded;
            }
        }
        bin1.insert(bin1.end(), tmpBin1.begin(), tmpBin1.end());
        bin2.insert(bin2.end(), tmpBin2.begin(), tmpBin2.end());
        count.insert(count.end(), tmpCount.begin(), tmpCount.end());
        clear();
        return nAdded;
    }
};


class SparseMatrices {
public:
    size_t nMatrices;
private:
    std::vector < SparseMatrix > matrices;
public:
    SparseMatrices (std::vector < int > &sizes): nMatrices(sizes.size() * (sizes.size() + 1) / 2), matrices(nMatrices) {
        size_t k = 0;
        for (size_t i = 0; i < sizes.size(); ++i) {
            for (size_t j = 0; j <= i; ++j) {
                matrices[k].setSizes(sizes[i], sizes[j]);
                ++k;
            } 
        } 
    }
    void addElement (size_t i, size_t r, size_t c) {
        matrices[i].addElement(r, c);
    }
    void clear () {
        for (auto &matrix: matrices) {
            matrix.clear();
        }
    }
    void writeToFile (std::string &fileName) {
        std::ofstream fileOut(fileName, std::ofstream::binary);
        for (auto &matrix: matrices) {
            matrix.writeToFile(fileOut);
        }
        fileOut.close();
    }
    void loadFromFile (std::string &fileName) {
        std::ifstream fileIn(fileName, std::ifstream::binary);
        //printf("load from file: %zu matrices\n", matrices.size());
        for (auto &matrix: matrices) {
            matrix.loadFromFile(fileIn);
        }
        fileIn.close();
    }
    void merge (SparseMatrices &otherMatrices) {
        for (size_t i = 0; i < matrices.size(); ++i) {
            matrices[i].merge(otherMatrices.matrices[i]);
        }
	otherMatrices.clear();
    }
    void convertToPlain (std::vector < int > &chrs1, std::vector < int > &chrs2, std::vector < int > &bins1, std::vector < int > &bins2, std::vector < int > &counts, std::vector < std::pair < size_t, size_t > > &translation, unsigned int minCount) {
        Progress progress(matrices.size(), true);
        std::vector < int > tmpChrs1, tmpChrs2;
        unsigned int nAdded;
        for (size_t i = 0; i < matrices.size(); ++i) {
            //printf("Converting matrix #%zu\n", i);
            int chr1, chr2;
            std::tie(chr1, chr2) = translation[i];
            nAdded = matrices[i].convertToPlain(bins1, bins2, counts, minCount);
            // R factors start with 1
            std::vector < int > tmpChrs1(nAdded, chr1 + 1);
            std::vector < int > tmpChrs2(nAdded, chr2 + 1);
            chrs1.insert(chrs1.end(), tmpChrs1.begin(), tmpChrs1.end());
            chrs2.insert(chrs2.end(), tmpChrs2.begin(), tmpChrs2.end());
            progress.increment();
        }
    }
};


struct contact_t {
    int chrId;
    int bin1;
    int bin2;
    contact_t (int c = -1, int b1 = 0, int b2 = 0): chrId(c), bin1(b1), bin2(b2) {}
    friend bool operator<(const contact_t &l, const contact_t &r) {
        return std::tie(l.chrId, l.bin1, l.bin2) < std::tie(r.chrId, r.bin1, r.bin2);
    }
};

struct heap_data_t {
    int fileId;
    int chrId;
    int bin1;
    int bin2;
    heap_data_t (int f = -1 , int c = -1, int b1 = 0, int b2 = 0): fileId(f), chrId(c), bin1(b1), bin2(b2) {}
    friend bool operator<(const heap_data_t &l, const heap_data_t &r) {
        return std::tie(l.chrId, l.bin1, l.bin2) > std::tie(r.chrId, r.bin1, r.bin2);
    }
    friend bool operator>(const heap_data_t &l, const heap_data_t &r) {
        return std::tie(l.chrId, l.bin1, l.bin2) < std::tie(r.chrId, r.bin1, r.bin2);
    }
};
*/

using positions_t     = std::vector < position_t >;
using umiMap_t        = std::unordered_map<std::string, positions_t >;
using matrixCounts_t  = std::vector < unsigned long >;
using validMatrices_t = std::vector < bool >;
using matrix_t        = std::unordered_map < uint64_t, unsigned long >;
using matrices_t      = std::vector < matrix_t >;
using chrIndices_t    = std::vector < std::pair < size_t, size_t > >;
using indicesChr_t    = std::vector < std::vector < size_t > >;

void updateMatrixCounts(matrixCounts_t &matrixCounts, positions_t &positions, indicesChr_t &indicesChr) {
    for (size_t positionId1 = 0; positionId1 < positions.size(); ++positionId1) {
        position_t &position1 = positions[positionId1];
        for (size_t positionId2 = 0; positionId2 <= positionId1; ++positionId2) {
            position_t &position2 = positions[positionId2];
            ++matrixCounts[indicesChr[position1.chrId][position2.chrId]];
        }
    }
}

void updateMatrices(matrices_t &matrices, positions_t &positions, validMatrices_t &validMatrices, indicesChr_t &indicesChr) {
    for (size_t positionId1 = 0; positionId1 < positions.size(); ++positionId1) {
        position_t &position1 = positions[positionId1];
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
    positions.clear();
}

void dump(matrix_t            &matrix,
          size_t               chrId,
          chrIndices_t        &chrIndices,
          std::vector < int > &stdChrs1,
          std::vector < int > &stdChrs2,
          std::vector < int > &stdBins1,
          std::vector < int > &stdBins2,
          std::vector < int > &stdCounts,
          unsigned int         minCount) {
    size_t chrId1, chrId2;
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
}

/*
std::mutex mtx;

void sortPositions (stackedPositions_t &stackedPositions, std::atomic_size_t &globalIndex, int threadId) {
    Rcout << "Starting thread #" << threadId << std::endl;
    size_t index;
    while ((index = globalIndex++) < stackedPositions.size()) {
        mtx.lock();
        Rcout << "thread " << threadId << ", index: " << index << std::endl;
        mtx.unlock();
        // Rcout << "Matrix " << index << "/" << stackedPositions.size() << " (" << stackedPositions[index].size() << " elements), thread #" << threadId << std::endl;
        //std::sort(stackedPositions[index].begin(), stackedPositions[index].end());
    }
    Rcout << "Over with this thread." << std::endl;
}

void writeContacts (std::vector < std::string > &tmpFileNames, std::vector < contact_t > &contacts) {
    std::string tmpFileName = tmpnam(NULL);
    std::ofstream tmpFile(tmpFileName, std::ofstream::binary);
    size_t size = contacts.size();
    tmpFile.write(reinterpret_cast < char const * > (&size), sizeof(size_t));
    tmpFile.write(reinterpret_cast < char const * > (contacts.data()), size * sizeof(contact_t));
    tmpFile.close();
    tmpFileNames.push_back(tmpFileName);
    contacts.clear();
}

void sortContacts (std::string &inputFileName, std::string &outputFileName) {
    std::vector < contact_t > contacts;
    size_t size;
    std::ifstream inputFile(inputFileName, std::ofstream::binary);
    inputFile.read(reinterpret_cast < char * > (&size), sizeof(size_t));
    std::cout << "Sorting '" << inputFileName << "' with " << size << " elements\n";
    contacts.resize(size);
    inputFile.read(reinterpret_cast <char * > (contacts.data()), size * sizeof(contact_t));
    inputFile.close();
    std::sort(contacts.begin(), contacts.end());
    std::ofstream outputFile(outputFileName, std::ofstream::binary);
    outputFile.write(reinterpret_cast < char const * > (contacts.data()), size * sizeof(contact_t));
    outputFile.close();
    if (remove(inputFileName.c_str()) != 0) {
        std::cerr << "Error, cannot remove temporary file '" << inputFileName << "'\n";
    }
}
*/

// [[Rcpp::export]]
List parseBamFileCpp(String fileName, int32_t binSize, int nThreads) {
    unsigned int minCount = 2;
    unsigned int minCountPerMatrix = 100;
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
        int chrId    = alignment->core.tid;
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
    Rcout << "Found " <<
        std::count(validMatrices.begin(), validMatrices.end(), true) <<
            "/" << nMatrices << " valid matrices.  Filling matrices.\n";
    Progress progress2(positionss.size(), true);
    matrices_t matrices (nMatrices);
    for (auto &positions: positionss) {
        updateMatrices(matrices, positions, validMatrices, indicesChr);
        progress2.increment();
    }
    positionss.clear();
    Rcout << "Exporting.\n";
    std::vector < int > stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts;
    IntegerVector chrs1, chrs2;
    IntegerVector bins1, bins2, counts;
    for (size_t i = 0; i < matrices.size(); ++i) {
        dump(matrices[i], i, chrIndices, stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts, minCount);
    }
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
    chrSizes = sizes;
    DataFrame outputDataFrame = DataFrame::create(_["ref1"]  = chrs1,
                                                  _["bin1"]  = bins1,
                                                  _["ref2"]  = chrs2,
                                                  _["bin2"]  = bins2,
                                                  _["count"] = counts);
    return List::create(_["data"] = outputDataFrame,
                        _["sizes"] = chrSizes);
    /*
    SparseMatrices sparseMatrices (sizes);
    std::vector < std::string > tmpFileNames;
    size_t nRecords = 0;
    size_t sumRecords = 0;
    for (auto &mapElement: umiMap) {
        auto &positions = mapElement.second;
        int nPositions = positions.size();
        for (int posId1 = 0; posId1 < nPositions; ++posId1) {
            for (int posId2 = 0; posId2 < posId1; ++posId2) {
                ref1 = positions[posId1].chrId;
                ref2 = positions[posId2].chrId;
                // Do not store "self" maps
                //if (ref1 != ref2) {
                    pos1 = positions[posId1].pos;
                    pos2 = positions[posId2].pos;
                    if (ref1 < ref2) {
                        std::swap(ref1, ref2);
                        std::swap(pos1, pos2);
                    }
                    else if ((ref1 == ref2) && (pos1 < pos2)) {
                        std::swap(pos1, pos2);
                    }
                    if (ref1 >= nRefs) {
                        Rcerr << "Error: first id reference exceeds size (" << ref1 << " vs " << nRefs << ")." << std::endl;
                    }
                    if (ref2 > ref1) {
                        Rcerr << "Error: second id reference exceeds first ref (" << ref2 << " vs " << ref1 << " vs" << nRefs << ")." << std::endl;
                    }
                    if (pos1 >= sizes[ref1]) {
                        Rcerr << "Error: row id exceeds size (" << pos1 << " vs " << sizes[ref1] << ")." << std::endl;
                    }
                    if (pos2 >= sizes[ref2]) {
                        Rcerr << "Error: col id exceeds size (" << pos2 << " vs " << sizes[ref2] << ")." << std::endl;
                    }
                    sparseMatrices.addElement(indicesChr[ref1][ref2], pos1, pos2);
                    ++nRecords;
                //}
            }
        }
        if (nRecords >= maxRecords) {
            tmpFileNames.emplace_back(tmpnam(NULL));
            sparseMatrices.writeToFile(tmpFileNames.back());
            sumRecords += nRecords;
            nRecords = 0;
        }
        progress1.increment();
        if (Progress::check_abort()) {
            return emptyDataFrame;
        }
    }
    umiMap.clear();
    if (nRecords > 0) {
        tmpFileNames.emplace_back(tmpnam(NULL));
        sparseMatrices.writeToFile(tmpFileNames.back());
        sumRecords += nRecords;
    }
    #ifdef _OPENMP
    if (nThreads > 0) {
        omp_set_num_threads(nThreads);
    }
    #endif
    Rcout << tmpFileNames.size() << " buckets with " << nRecords << " records.  Merging them.\n";
    if (tmpFileNames.empty()){
        return emptyDataFrame;
    }
    //Progress progress2(tmpFileNames.size(), true);
    //#pragma omp parallel for
    //#pragma omp parallel for default(none) private(tmpFileNames, tmpFileNames2) shared(progress2) schedule(dynamic)
    //#pragma omp parallel for default(none) shared(tmpFileNames, tmpFileNames2, progress2) schedule(dynamic)
    //#pragma omp parallel for default(none) shared(tmpFileNames, tmpFileNames2, sizes) schedule(dynamic)
    Rcout << "Merging buckets (" << tmpFileNames.size() << " buckets).\n";
    Progress progress2(tmpFileNames.size(), true);
    sparseMatrices.loadFromFile(tmpFileNames.front());
    progress2.increment();
    for (size_t i = 1; i < tmpFileNames.size(); ++i) {
        SparseMatrices otherSparseMatrices (sizes);
        otherSparseMatrices.loadFromFile(tmpFileNames[i]);
        sparseMatrices.merge(otherSparseMatrices);
        if (std::remove(tmpFileNames[i].c_str()) != 0) {
            Rcout << "Error deleting file " << tmpFileNames[i] << "\n";
        }
        progress2.increment();
    }
    Rcout << "\nMerging done, exporting to R.\n";
    std::vector < int > ref1Vector;
    std::vector < int > ref2Vector;
    std::vector < int > pos1Vector;
    std::vector < int > pos2Vector;
    std::vector < int > countVector;
    sparseMatrices.convertToPlain(ref1Vector, ref2Vector, pos1Vector, pos2Vector, countVector, chrIndices, minCount);
    Rcout << "\n" << ref1Vector.size() << " points found.\n";
    sparseMatrices.clear();
    IntegerVector ref1VectorR (ref1Vector.begin(), ref1Vector.end());
    ref1Vector.clear();
    IntegerVector ref2VectorR (ref2Vector.begin(), ref2Vector.end());
    ref2Vector.clear();
    IntegerVector pos1VectorR (pos1Vector.begin(), pos1Vector.end());
    pos1Vector.clear();
    IntegerVector pos2VectorR (pos2Vector.begin(), pos2Vector.end());
    pos2Vector.clear();
    IntegerVector countVectorR (countVector.begin(), countVector.end());
    countVector.clear();
    //Rcout << "Counting done...\n";
    ref1VectorR.attr("class") = "factor";
    ref1VectorR.attr("levels") = refs;
    ref2VectorR.attr("class") = "factor";
    ref2VectorR.attr("levels") = refs;
    DataFrame outputData = DataFrame::create(_["ref1"] = ref1VectorR, _["bin1"] = pos1VectorR, _["ref2"] = ref2VectorR, _["bin2"] = pos2VectorR, _["count"] = countVectorR);
    //ProfilerStop();
    IntegerVector chrSizes;
    chrSizes.assign(sizes.begin(), sizes.end());
    chrSizes.names() = refs;
    return List::create(_["data"] = outputData,
                        _["sizes"] = chrSizes);
     */
}
