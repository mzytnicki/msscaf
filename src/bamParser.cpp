#include <string>
#include <vector>
#include <unordered_map>
#include <valarray>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <thread>
#include <atomic>
#include <mutex>
#include "sam.h"

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

struct position_t {
    int     chrId;
    int32_t pos;
};

typedef std::unordered_map<std::string, std::vector<position_t>> umiMap_t;
typedef std::vector < std::vector < std::pair < int, int > > > stackedPositions_t;

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

// struct StackedPositions : public RcppParallel::Worker {
//     // destination matrix
//     stackedPositions_t &stackedPositions;
//     std::vector < std::pair < size_t, size_t > > chrIndices;
//     
//     StackedPositions(stackedPositions_t &sp, size_t nc): stackedPositions(sp), chrIndices(nc * (nc + 1) / 2) {
//         size_t p = 0;
//         for (size_t i = 0; i < nc; ++i) {
//             for (size_t j = 0; j <= i; ++j) {
//                 chrIndices[p++] = {i, j};
//             }
//         }
//     }
//     
//     // take the square root of the range of elements requested
//     void operator()(std::size_t begin, std::size_t end) {
//         for (size_t i = begin; i < end; ++i) {
//             std::sort(stackedPositions[chrIndices[i].first][chrIndices[i].second].begin(), stackedPositions[chrIndices[i].first][chrIndices[i].second].end());
//         }
//     }
// };

class SparseMatrix {
public:
    size_t n_rows;
    size_t n_cols;
private:
    std::unordered_map<size_t, unsigned int> matrix;
    std::unordered_map<size_t, unsigned int>::const_iterator matrixIterator;
    size_t translateCoordinates (size_t r, size_t c) {
        return n_rows * r + c;
    }
public:
    SparseMatrix (size_t nr = 0, size_t nc = 0) {
        n_rows = nr;
        n_cols = nc;
    }
    void addElement (size_t r, size_t c) {
        ++matrix[translateCoordinates(r, c)];
    }
    size_t getNElements () {
        return matrix.size();
    }
    void startIterator () {
        matrixIterator = matrix.begin();
    }
    bool iteratorOver () {
        return (matrixIterator == matrix.end());
    }
    void incIterator () {
        ++matrixIterator;
    }
    size_t getRow () {
        return matrixIterator->first / n_rows;
    }
    size_t getCol () {
        return matrixIterator->first % n_rows;
    }
    unsigned int getValue () {
        return matrixIterator->second;
    }
};

class SimpleMatrix {
public:
    size_t n_rows;
    size_t n_cols;
    size_t n_nonzero;
private:
    std::valarray<unsigned int> matrix;
    size_t currentRow;
    size_t currentCol;
    size_t currentPos;
    size_t translateCoordinates (size_t r, size_t c) {
        return n_rows * (n_rows - 1) / 2 - (n_rows - c) * (n_rows - c - 1) / 2 + r;
    }
public:
    SimpleMatrix (size_t n = 0): n_rows(n), n_cols(n), n_nonzero(0), matrix(static_cast<unsigned int>(0), n * (n + 1) / 2) { }
    void addElement (size_t r, size_t c) {
        size_t pos = translateCoordinates(r, c);
        if (matrix[pos] == 0) {
          ++n_nonzero;
        }
        ++matrix[pos];
    }
    bool iteratorOver () {
        return (currentCol >= n_cols);
    }
    void _incIterator () {
        ++currentPos;
        ++currentRow;
        if (currentRow > currentCol) {
            currentRow = 0;
            ++currentCol;
        }
    }
    void incIterator () {
        do {
            _incIterator();
        } while (matrix[currentPos] == 0);
    }
    void startIterator () {
        currentPos = 0;
        currentRow = 0;
        currentCol = 0;
        if (matrix[currentPos] == 0) {
            incIterator();
        }
    }
    size_t getRow () {
        return currentRow;
    }
    size_t getCol () {
        return currentCol;
    }
    unsigned int getValue () {
        return matrix[currentPos];
    }
};

unsigned int getChr (unsigned int position, std::vector <unsigned int> &offsets) {
    unsigned int minId = 0;
    unsigned int maxId = offsets.size() - 1;
    unsigned int midId;
    while (true) {
        midId = (minId + maxId) / 2;
        if (position < offsets[midId]) {
            maxId = midId - 1;
        }
        if (position > offsets[midId+1]) {
            minId = midId + 1;
        }
        return midId;
    }
    Rcerr << "Error!  Cannot find offset of " << position << ".\n";
    return midId;
}


// [[Rcpp::export]]
DataFrame parseBamFileCpp(String fileName, int binSize, int nThreads) {
    Rcout << "Reading " << fileName.get_cstring() << "\n";
    samFile   *inputFile  = hts_open(fileName.get_cstring(), "r");
    if (inputFile == NULL) {
        Rcerr << "File does not exist.  Exiting.\n";
        DataFrame outputData;
        return outputData;
    }
    bam_hdr_t *header     = sam_hdr_read(inputFile);
    bam1_t    *alignment  = bam_init1();
    int        nRefs      = header->n_targets;
    Rcout << "Found " << nRefs << " references\n";
    std::vector <int> sizes(header->n_targets);
    // std::vector <unsigned int> offsets(header->n_targets+1);
    CharacterVector refs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        refs[i] = header->target_name[i];
        sizes[i] = header->target_len[i] / binSize + 1;
    }
    // offsets[0];
    // for (int i = 0; i < header->n_targets; ++i) {
    //     offsets[i+1] = offsets[i] + (sizes[i] / binSize) + 1;
    // }
    // unsigned int nBins = offsets[header->n_targets];
    // std::vector<unsigned int> bins2Chr(nBins);
    // for (int i = 0; i < header->n_targets; ++i) {
    //     for (int j = offsets[i]; j < offsets[i+1]; ++j) {
    //         bins2Chr[j] = i;
    //     }
    // }
    //std::vector<std::vector<arma::sp_mat>> matrices(nRefs);
    //std::vector<std::vector<SparseMatrix>> matrices(nRefs);
    //SparseMatrix matrix(nBins, nBins);
    //arma::sp_mat matrix(nBins, nBins);
    umiMap_t   umiMap;
    /*
    for (int idRef1 = 0; idRef1 < nRefs; ++idRef1) {
        //matrices[idRef1] = std::vector<SparseMatrix>(idRef1+1);
        matrices[idRef1] = std::vector<arma::sp_mat>(idRef1+1);
        for (int idRef2 = 0; idRef2 <= idRef1; ++idRef2) {
            //Rcout << "\tMatrix [" << idRef1 << ", " << idRef2 << "] has size (" << (header->target_len[idRef1] / binSize) << ", " << (header->target_len[idRef2] / binSize) << ")\n";
            // Adding one more row and col, because division is truncated, and position "n * binSize" should go to bin "n + 1"
            //matrices[idRef1][idRef2] = SparseMatrix((header->target_len[idRef1] / binSize) + 1, (header->target_len[idRef2] / binSize) + 1);
            matrices[idRef1][idRef2] = arma::sp_mat((header->target_len[idRef1] / binSize) + 1, (header->target_len[idRef2] / binSize) + 1);
        }
    }
    */
    Rcout << "Reading BAM file.\n";
    for (unsigned int cpt = 0; sam_read1(inputFile, header, alignment) > 0; ++cpt) {
        int32_t pos    = alignment->core.pos + 1;
        int32_t posBin = pos / binSize;
        int chrId    = alignment->core.tid;
        uint8_t *umi = bam_aux_get(alignment, "BX");
        if ((chrId != -1) && (umi != NULL)) {
            std::string umiString((char *) ++umi);
            umiMap[umiString].push_back({chrId, posBin});
        }
        if (cpt % 10000000 == 0) {
            Rcout << "\tReading read #" << cpt << "\n";
        }
    }
    bam_destroy1(alignment);
    sam_close(inputFile);
    Rcout << "Read parsing done.  Placing to buckets.\n";
    //Rcout << "Read parsing done.  Filling matrices (" << nBins << " bins).\n";
    int ref1, ref2;
    int pos1, pos2;
    //SimpleMatrix matrix(nBins);
    stackedPositions_t stackedPositions(nRefs * (nRefs + 1) / 2);
    std::vector < std::pair < size_t, size_t > > chrIndices (nRefs * (nRefs + 1) / 2);
    std::vector < std::vector < size_t > > indicesChr (nRefs);
    size_t p = 0;
    for (int i = 0; i < nRefs; ++i) {
        indicesChr[i] = std::vector < size_t > (i+1);
        for (int j = 0; j <= i; ++j) {
            chrIndices[p] = {i, j};
            indicesChr[i][j] = p;
            ++p;
        }
    }
    Progress progress1(umiMap.size(), true);
    while (! umiMap.empty()) {
        auto &mapElement = *umiMap.begin();
        auto &positions = mapElement.second;
        int nPositions = positions.size();
        for (int posId1 = 0; posId1 < nPositions; ++posId1) {
            for (int posId2 = 0; posId2 < posId1; ++posId2) {
                ref1 = positions[posId1].chrId;
                ref2 = positions[posId2].chrId;
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
                stackedPositions[indicesChr[ref1][ref2]].emplace_back(pos1, pos2);
            }
        }
        umiMap.erase(umiMap.begin());
        progress1.increment();
    }
    Rcout << "Bucketting done.  Sorting them with " << nThreads << " threads.\n";
    // std::atomic_size_t i (0);
    // std::vector < std::thread > threads;
    // for (int threadId = 0; threadId < nThreads - 1; ++threadId) {
    //     threads.push_back(std::thread(sortPositions, std::ref(stackedPositions), std::ref(i), threadId));
    // }
    // sortPositions(stackedPositions, i, nThreads - 1);
    // for (auto &thread: threads) {
    //     thread.join();
    // }
    // for (size_t i = 0; i < chrIndices.size(); ++i) {
    //     std::sort(stackedPositions[chrIndices[i].first][chrIndices[i].second].begin(), stackedPositions[chrIndices[i].first][chrIndices[i].second].end());
    // }
    // StackedPositions stackedPositionsWorker(stackedPositions, nRefs);
    // parallelFor(0, nRefs * (nRefs + 1) / 2, stackedPositionsWorker);
    // //Progress progress2(nRefs * (nRefs+1) / 2, true);
    // for (int ref1Id = 0; ref1Id < nRefs; ++ref1Id) {
    //     for (int ref2Id = 0; ref2Id <= ref1Id; ++ref2Id) {
    //         Rcout << "\tSorting " << ref1Id << " vs " << ref2Id << ": " << stackedPositions[ref1Id][ref2Id].size() << " elements." << std::endl;
    //         std::sort(stackedPositions[ref1Id][ref2Id].begin(), stackedPositions[ref1Id][ref2Id].end());
    //         //std::sort(unsortedCounts[ref1Id][ref2Id].begin(), unsortedCounts[ref1Id][ref2Id].end(), [] (const std::pair < int, int > & a, const std::pair <int, int > & b) { return ((a.first < b.first) || ((a.first == b.first) && (a.second < b.second))); }));
    //         //progress2.increment();
    //     }
    // }
    // Rcout << "Sorting done.  Getting counts.\n";
    // Progress progress3(nRefs * (nRefs+1) / 2, true);
    // std::vector<int> ref1Vector, ref2Vector, pos1Vector, pos2Vector, countVector;
    // for (size_t ref = 0; ref < stackedPositions.size(); ++ref) {
    //     int ref1Id, ref2Id;
    //     std::tie(ref1Id, ref2Id) = chrIndices[ref];
    //     std::vector < int > theseRef1;
    //     std::vector < int > theseRef2;
    //     std::vector < int > thesePos1;
    //     std::vector < int > thesePos2;
    //     std::vector < int > theseCounts;
    //     std::pair < int, int > prev = { 0, 0 };
    //     int prevCount = 0;
    //     for (auto &p: stackedPositions[ref]) {
    //         if (p == prev) {
    //             ++prevCount;
    //         }
    //         else if (prevCount != 0) {
    //             thesePos1.push_back(prev.first);
    //             thesePos2.push_back(prev.second);
    //             theseCounts.push_back(prevCount);
    //             prevCount = 1;
    //             prev      = p;
    //         }
    //     }
    //     if (prevCount != 0) {
    //         thesePos1.push_back(prev.first);
    //         thesePos2.push_back(prev.second);
    //         theseCounts.push_back(prevCount);
    //     }
    //     theseRef1 = std::vector < int > (thesePos1.size(), ref1Id + 1); // Factors in R start with 1
    //     theseRef2 = std::vector < int > (thesePos1.size(), ref2Id + 1);
    //     ref1Vector.insert(ref1Vector.end(), theseRef1.begin(), theseRef1.end());
    //     ref2Vector.insert(ref2Vector.end(), theseRef2.begin(), theseRef2.end());
    //     pos1Vector.insert(pos1Vector.end(), thesePos1.begin(), thesePos1.end());
    //     pos2Vector.insert(pos2Vector.end(), thesePos2.begin(), thesePos2.end());
    //     countVector.insert(countVector.end(), theseCounts.begin(), theseCounts.end());
    //     progress3.increment();
    // }
    std::vector<int> ref1Vector, ref2Vector, pos1Vector, pos2Vector, countVector;
    //Progress progress2(nRefs * (nRefs+1) / 2, true);
    for (int ref = 0; ref < nRefs * (nRefs+1) / 2; ++ref) {
        if (! stackedPositions[ref].empty()) {
            std::vector<int> thisRef1Vector, thisRef2Vector, thisPos1Vector, thisPos2Vector, thisCountVector;
            int nNonZeros = 0;
            int ref1, ref2;
            std::tie(ref1, ref2) = chrIndices[ref];
            int size1 = sizes[ref1];
            int size2 = sizes[ref2];
            std::cout << "Ref: " << ref1 << "/" << ref2 << ": " << stackedPositions[ref].size() << " elements.\n";
            std::vector < std::vector < int > > matrix (size1, std::vector < int > (size2, 0));
            for (auto &positions: stackedPositions[ref]) {
                int pos1, pos2;
                std::tie(pos1, pos2) = positions;
                if (matrix[pos1][pos2] == 0) {
                    ++nNonZeros;
                }
                ++matrix[pos1][pos2];
            }
            std::cout << "\t" << nNonZeros << " non zero elements.\n";
            thisPos1Vector.reserve(nNonZeros);
            thisPos2Vector.reserve(nNonZeros);
            thisCountVector.reserve(nNonZeros);
            thisRef1Vector.insert(thisRef1Vector.begin(), nNonZeros, ref1 + 1);
            thisRef2Vector.insert(thisRef2Vector.begin(), nNonZeros, ref2 + 1);
            for (int pos1 = 0; pos1 < size1; ++pos1) {
                for (int pos2 = 0; pos2 < size2; ++pos2) {
                    if (matrix[pos1][pos2] != 0) {
                        thisPos1Vector.push_back(pos1);
                        thisPos2Vector.push_back(pos2);
                        thisCountVector.push_back(matrix[pos1][pos2]);
                    }
                }
            }
            ref1Vector.insert(ref1Vector.end(), thisRef1Vector.begin(), thisRef1Vector.end());
            ref2Vector.insert(ref2Vector.end(), thisRef2Vector.begin(), thisRef2Vector.end());
            pos1Vector.insert(pos1Vector.end(), thisPos1Vector.begin(), thisPos1Vector.end());
            pos2Vector.insert(pos2Vector.end(), thisPos2Vector.begin(), thisPos2Vector.end());
            countVector.insert(countVector.end(), thisCountVector.begin(), thisCountVector.end());
        }
        //progress2.increment();
    }
    //Rcout << "Counting done...\n";
    IntegerVector ref1VectorR(ref1Vector.begin(), ref1Vector.end());
    IntegerVector ref2VectorR(ref2Vector.begin(), ref2Vector.end());
    IntegerVector pos1VectorR(pos1Vector.begin(), pos1Vector.end());
    IntegerVector pos2VectorR(pos2Vector.begin(), pos2Vector.end());
    IntegerVector countVectorR(countVector.begin(), countVector.end());
    ref1VectorR.attr("class") = "factor";
    ref1VectorR.attr("levels") = refs;
    ref2VectorR.attr("class") = "factor";
    ref2VectorR.attr("levels") = refs;
    DataFrame outputData = DataFrame::create(_["ref1"] = ref1VectorR, _["bin1"] = pos1VectorR, _["ref2"] = ref2VectorR, _["bin2"] = pos2VectorR, _["count"] = countVectorR);
    return outputData;
}
