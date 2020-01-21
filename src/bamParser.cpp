#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "sam.h"

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
    bam_hdr_t *header     = sam_hdr_read(inputFile);
    bam1_t    *alignment  = bam_init1();
    int        nRefs      = header->n_targets;
    Rcout << "Found " << nRefs << " references\n";
    CharacterVector refs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        refs[i] = header->target_name[i];
    }
    std::vector<std::vector<IntegerMatrix>> matrices(nRefs);
    umiMap_t   umiMap;
    for (int idRef1 = 0; idRef1 < nRefs; ++idRef1) {
        matrices[idRef1] = std::vector<IntegerMatrix>(idRef1+1);
        for (int idRef2 = 0; idRef2 <= idRef1; ++idRef2) {
            Rcout << "Matrix [" << idRef1 << ", " << idRef2 << "] has size (" << (header->target_len[idRef1] / binSize) << ", " << (header->target_len[idRef2] / binSize) << ")\n";
            matrices[idRef1][idRef2] = IntegerMatrix(header->target_len[idRef1] / binSize, header->target_len[idRef2] / binSize);
        }
    }
    while (sam_read1(inputFile, header, alignment) > 0) {
        int32_t pos    = alignment->core.pos + 1;
        int32_t posBin = pos / binSize;
        int chrId    = alignment->core.tid;
        uint8_t *umi = bam_aux_get(alignment, "BX");
        if (umi != NULL) {
            std::string umiString((char *) ++umi);
            umiMap[umiString].push_back({chrId, posBin});
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
    IntegerVector ref1Vector, ref2Vector, pos1Vector, pos2Vector, countVector;
    for (int idRef1 = 0; idRef1 < nRefs; ++idRef1) {
        for (int idRef2 = 0; idRef2 <= idRef1; ++idRef2) {
            for (unsigned int posId1 = 0; posId1 < header->target_len[idRef1] / binSize; ++posId1) {
                for (unsigned int posId2 = 0; posId2 < header->target_len[idRef2] / binSize; ++posId2) {
                    if (matrices[idRef1][idRef2](posId1, posId2) != 0) {
                        ref1Vector.push_back(idRef1+1);
                        ref2Vector.push_back(idRef2+1);
                        pos1Vector.push_back(posId1);
                        pos2Vector.push_back(posId2);
                        countVector.push_back(matrices[idRef1][idRef2](posId1, posId2));
                    }
                }
            }
        }
    }
    ref1Vector.attr("class") = "factor";
    ref1Vector.attr("levels") = refs;
    ref2Vector.attr("class") = "factor";
    ref2Vector.attr("levels") = refs;
    DataFrame outputData = DataFrame::create(_["ref1"] = ref1Vector, _["bin1"] = pos1Vector, _["ref2"] = ref2Vector, _["bin2"] = pos2Vector, _["count"] = countVector);
    return outputData;
}