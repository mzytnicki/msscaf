#include <string>
#include <vector>
#include <unordered_map>
#include <valarray>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "sam.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

struct position_t {
    int     chrId;
    int32_t pos;
};

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

typedef std::unordered_map<std::string, std::vector<position_t>> umiMap_t;

// [[Rcpp::export]]
DataFrame parseBamFileCpp(String fileName, int binSize) {
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
    std::vector <unsigned int> sizes(header->n_targets);
    std::vector <unsigned int> offsets(header->n_targets+1);
    CharacterVector refs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        refs[i] = header->target_name[i];
        sizes[i] = header->target_len[i];
    }
    offsets[0];
    for (int i = 0; i < header->n_targets; ++i) {
        offsets[i+1] = offsets[i] + (sizes[i] / binSize) + 1;
    }
    unsigned int nBins = offsets[header->n_targets];
    std::vector<unsigned int> bins2Chr(nBins);
    for (int i = 0; i < header->n_targets; ++i) {
        for (int j = offsets[i]; j < offsets[i+1]; ++j) {
            bins2Chr[j] = i;
        }
    }
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
    Rcout << "Created matrices.\n";
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
            Rcout << "Reading read #" << cpt << "\n";
        }
    }
    bam_destroy1(alignment);
    sam_close(inputFile);
    Rcout << "Read parsing done.  Filling matrices (" << nBins << " bins).\n";
    unsigned int ref1, ref2;
    unsigned int pos1, pos2;
    SimpleMatrix matrix(nBins);
    Progress progress1(umiMap.size(), true);
    for (auto &mapElement: umiMap) {
        auto &positions = mapElement.second;
        size_t nPositions = positions.size();
        for (unsigned int posId1 = 0; posId1 < nPositions; ++posId1) {
            for (unsigned int posId2 = 0; posId2 < posId1; ++posId2) {
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
                    Rcerr << "Error: first id reference exceeds size (" << ref1 << " vs " << nRefs << ").\n";
                }
                if (ref2 >= nRefs) {
                    Rcerr << "Error: second id reference exceeds size (" << ref2 << " vs " << nRefs << ").\n";
                }
                if (pos1 >= sizes[ref1]) {
                    Rcerr << "Error: row id exceeds size (" << pos1 << " vs " << sizes[ref1] << ").\n";
                }
                if (pos2 >= sizes[ref2]) {
                    Rcerr << "Error: col id exceeds size (" << pos2 << " vs " << sizes[ref2] << ").\n";
                }
                matrix.addElement(offsets[ref1] + pos1, offsets[ref2] + pos2);
                //++matrix(offsets[ref1] + pos1, offsets[ref2] + pos2);
            }
        }
        progress1.increment();
    }
    umiMap.clear();
    Rcout << "Matrix filling done, transforming data to sparse matrices.\n";
    std::vector<int> ref1Vector, ref2Vector, pos1Vector, pos2Vector, countVector;
    int nElements = matrix.n_nonzero;
    ref1Vector.reserve(nElements);
    ref2Vector.reserve(nElements);
    pos1Vector.reserve(nElements);
    pos2Vector.reserve(nElements);
    countVector.reserve(nElements);
    Progress progress2(nElements, true);
    //for (arma::sp_mat::const_iterator matrixIt = matrix.begin(); matrixIt != matrix.end(); ++matrixIt) {
    for (matrix.startIterator(); ! matrix.iteratorOver(); matrix.incIterator()) {
        // Rcout << "\t\t[" << matrixIt.row() << ", " << matrixIt.col() << "]: " << (*matrixIt) << "\n";
        //unsigned int pos = matrixIt.row();
        unsigned int pos = matrix.getRow();
        //unsigned int chrId = getChr(pos, offsets);
        unsigned int chrId = bins2Chr[pos];
        //pos1Vector.push_back(matrixIt.row() - pos);
        pos1Vector.push_back(pos - offsets[chrId]);
        ref1Vector.push_back(chrId + 1); // these are R factors, which are 1-based
        //pos = matrixIt.col();
        pos = matrix.getCol();
        //chrId = getChr(pos, offsets);
        chrId = bins2Chr[pos];
        pos2Vector.push_back(pos - offsets[chrId]);
        ref2Vector.push_back(chrId + 1); // these are R factors, which are 1-based
        //countVector.push_back(*matrixIt);
        countVector.push_back(matrix.getValue());
        if (ref1Vector.size() >= nElements) {
            Rcerr << "\nProblem while filling the sparse matrix.\n";
        }
        progress2.increment();
    }
    Rcout << "Sparse matrix filling done.\n";
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
