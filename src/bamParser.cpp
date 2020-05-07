#include <string>
#include <vector>
#include <unordered_map>
#include <RcppArmadillo.h>
#include "sam.h"

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

struct position_t {
    int     chrId;
    int32_t pos;
};

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
    CharacterVector refs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        refs[i] = header->target_name[i];
    }
    std::vector<std::vector<arma::sp_mat>> matrices(nRefs);
    umiMap_t   umiMap;
    for (int idRef1 = 0; idRef1 < nRefs; ++idRef1) {
        matrices[idRef1] = std::vector<arma::sp_mat>(idRef1+1);
        for (int idRef2 = 0; idRef2 <= idRef1; ++idRef2) {
            //Rcout << "\tMatrix [" << idRef1 << ", " << idRef2 << "] has size (" << (header->target_len[idRef1] / binSize) << ", " << (header->target_len[idRef2] / binSize) << ")\n";
            // Adding one more row and col, because division is truncated, and position "n * binSize" should go to bin "n + 1"
            matrices[idRef1][idRef2] = arma::sp_mat((header->target_len[idRef1] / binSize) + 1, (header->target_len[idRef2] / binSize) + 1);
        }
    }
    Rcout << "Created matrices.\n";
    for (unsigned int cpt = 0; sam_read1(inputFile, header, alignment) > 0; ++cpt) {
        int32_t pos    = alignment->core.pos + 1;
        int32_t posBin = pos / binSize;
        int chrId    = alignment->core.tid;
        uint8_t *umi = bam_aux_get(alignment, "BX");
        if (umi != NULL) {
            std::string umiString((char *) ++umi);
            umiMap[umiString].push_back({chrId, posBin});
        }
        if (cpt % 1000000 == 0) {
            Rcout << "Reading read #" << cpt << "\n";
        }
    }
    bam_destroy1(alignment);
    sam_close(inputFile);
    Rcout << "Read parsing done.\n";
    unsigned int ref1, ref2;
    unsigned int pos1, pos2;
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
                ++matrices[ref1][ref2](pos1, pos2);
            }
        }
    }
    Rcout << "Matrix filling done.\n";
    std::vector<int> ref1Vector, ref2Vector, pos1Vector, pos2Vector, countVector;
    for (int idRef1 = 0; idRef1 < nRefs; ++idRef1) {
        for (int idRef2 = 0; idRef2 <= idRef1; ++idRef2) {
            // Rcout << "\tScanning matrix (" << idRef1 << ", " << idRef2 << ")\n";
            arma::sp_mat &matrix = matrices[idRef1][idRef2];
            int nElements = matrix.n_nonzero;
            std::vector<int> thisPos1Vector, thisPos2Vector, thisCountVector;
            thisPos1Vector.reserve(nElements);
            thisPos2Vector.reserve(nElements);
            thisCountVector.reserve(nElements);
            for (arma::sp_mat::const_iterator matrixIt = matrix.begin(); matrixIt != matrix.end(); ++matrixIt) {
                //Rcout << "\t\t[" << matrixIt.row() << ", " << matrixIt.col() << "]: " << (*matrixIt) << "\n";
                thisPos1Vector.push_back(matrixIt.row());
                thisPos2Vector.push_back(matrixIt.col());
                thisCountVector.push_back(*matrixIt);
            }
            std::vector<int> thisRef1Vector(nElements, idRef1 + 1); // these are R factors, which are 1-based
            std::vector<int> thisRef2Vector(nElements, idRef2 + 1);
            ref1Vector.insert(ref1Vector.end(), thisRef1Vector.begin(), thisRef1Vector.end());
            ref2Vector.insert(ref2Vector.end(), thisRef2Vector.begin(), thisRef2Vector.end());
            pos1Vector.insert(pos1Vector.end(), thisPos1Vector.begin(), thisPos1Vector.end());
            pos2Vector.insert(pos2Vector.end(), thisPos2Vector.begin(), thisPos2Vector.end());
            countVector.insert(countVector.end(), thisCountVector.begin(), thisCountVector.end());
        }
    }
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