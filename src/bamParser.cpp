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

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const long maxRecords = 100000000;

struct position_t {
    int     chrId;
    int32_t pos;
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

typedef std::unordered_map<std::string, std::vector<position_t>> umiMap_t;
typedef std::vector < contact_t > stackedPositions_t;

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
    tmpFile.write(reinterpret_cast < char const * > (&size), sizeof(size));
    tmpFile.write(reinterpret_cast < char const * > (contacts.data()), size * sizeof(contact_t));
    tmpFile.close();
    tmpFileNames.push_back(tmpFileName);
    contacts.clear();
}

void sortContacts (std::string &inputFileName, std::string &outputFileName) {
    std::vector < contact_t > contacts;
    size_t size;
    std::ifstream inputFile(inputFileName, std::ofstream::binary);
    inputFile.read(reinterpret_cast < char * > (&size), sizeof(size));
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

// [[Rcpp::export]]
DataFrame parseBamFileCpp(String fileName, int binSize, int nThreads) {
    Rcout << "Reading " << fileName.get_cstring() << "\n";
    samFile   *inputFile  = hts_open(fileName.get_cstring(), "r");
    DataFrame emptyDataFrame;
    if (inputFile == NULL) {
        Rcerr << "File does not exist.  Exiting.\n";
        return emptyDataFrame;
    }
    bam_hdr_t *header     = sam_hdr_read(inputFile);
    bam1_t    *alignment  = bam_init1();
    int        nRefs      = header->n_targets;
    Rcout << "Found " << nRefs << " references\n";
    std::vector <int> sizes(header->n_targets);
    CharacterVector refs(header->n_targets);
    for (int i = 0; i < header->n_targets; ++i) {
        refs[i] = header->target_name[i];
        sizes[i] = header->target_len[i] / binSize + 1;
    }
    umiMap_t umiMap;
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
    stackedPositions_t stackedPositions;
    stackedPositions.reserve(2 * maxRecords);
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
    std::vector < size_t > bucketSizes(nRefs * (nRefs + 1) / 2, 0);
    Progress progress1(umiMap.size(), true);
    std::vector < std::string > tmpFileNames;
    std::vector < std::string > tmpFileNames2;
    size_t nRecords = 0;
    for (auto &mapElement: umiMap) {
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
                stackedPositions.emplace_back(indicesChr[ref1][ref2], pos1, pos2);
                ++nRecords;
            }
        }
        if (nRecords >= maxRecords) {
            writeContacts(tmpFileNames, stackedPositions);
            stackedPositions.reserve(2 * maxRecords);
            nRecords = 0;
        }
        progress1.increment();
        if (Progress::check_abort()) {
            return emptyDataFrame;
        }
    }
    umiMap.clear();
    if (nRecords > 0) {
        writeContacts(tmpFileNames, stackedPositions);
    }
    tmpFileNames2.reserve(tmpFileNames.size());
    for (size_t i = 0; i < tmpFileNames.size(); ++i) {
        tmpFileNames2.push_back(tmpnam(NULL));
    }
    #ifdef _OPENMP
    if (nThreads > 0) {
        omp_set_num_threads(nThreads);
    }
    #endif
    Rcout << tmpFileNames.size() << " buckets with " << nThreads << " thread(s).  Sorting them.\n";
    //Progress progress2(tmpFileNames.size(), true);
    //#pragma omp parallel for
    //#pragma omp parallel for default(none) private(tmpFileNames, tmpFileNames2) shared(progress2) schedule(dynamic)
    //#pragma omp parallel for default(none) shared(tmpFileNames, tmpFileNames2, progress2) schedule(dynamic)
    #pragma omp parallel for default(none) shared(tmpFileNames, tmpFileNames2) schedule(dynamic)
    for (size_t i = 0; i < tmpFileNames.size(); ++i) {
        sortContacts(tmpFileNames[i], tmpFileNames2[i]);
        //progress2.increment();
        /*
        if (Progress::check_abort()) {
            return emptyDataFrame;
        }
        */
    }
    Rcout << "Merging buckets.\n";
    std::priority_queue < heap_data_t, std::vector < heap_data_t > > heap;
    std::vector < std::ifstream > files;
    IntegerVector ref1VectorR;
    IntegerVector ref2VectorR;
    IntegerVector pos1VectorR;
    IntegerVector pos2VectorR;
    IntegerVector countVectorR;
    files.reserve(tmpFileNames2.size());
    int chrId, chrId1, chrId2;
    int bin1, bin2;
    contact_t record;
    for (size_t fileId = 0; fileId < tmpFileNames2.size(); ++fileId) {
        files.emplace_back(tmpFileNames2[fileId].c_str(), std::fstream::binary);
        files[fileId].read(reinterpret_cast < char * >(&record), sizeof(record));
        //std::cout << "First (" << tmpFileNames2[fileId] << "): " << fileId << " -> " << record.chrId << ":" << record.bin1 << "-" << record.bin2 << "\n";
        heap.emplace(fileId, record.chrId, record.bin1, record.bin2);
    }
    unsigned int count = 0;
    chrId = -1;
    Progress progress3(maxRecords, true);
    heap_data_t heapElement;
    while (! heap.empty()) {
        heapElement = heap.top();
        heap.pop();
        //std::cout << "Third: " << heapElement.fileId << " -> " << heapElement.chrId << ":" << heapElement.bin1 << "-" << heapElement.bin2 << "\n";
        if (heapElement.fileId == 0) {
            progress3.increment();
        }
        if ((heapElement.chrId == chrId) && (heapElement.bin1 == bin1) && (heapElement.bin2 == bin2)) {
            ++count;
        }
        else {
            if (count > 0) {
            	std::tie(chrId1, chrId2) = chrIndices[chrId];
                ref1VectorR.push_back(chrId1 + 1);
                ref2VectorR.push_back(chrId2 + 1);
                pos1VectorR.push_back(bin1);
                pos2VectorR.push_back(bin2);
                countVectorR.push_back(count);
            }
            chrId = heapElement.chrId;
            bin1  = heapElement.bin1;
            bin2  = heapElement.bin2;
            count = 1;
        }
        if (files[heapElement.fileId].good()) {
            files[heapElement.fileId].read(reinterpret_cast < char * >(&record), sizeof(record));
            //std::cout << "Second: " << heapElement.fileId << " -> " << record.chrId << "-" << record.bin1 << "-" << record.bin2 << "\n";
            heap.emplace(heapElement.fileId, record.chrId, record.bin1, record.bin2);
        }
        if (Progress::check_abort()) {
            return emptyDataFrame;
        }
    }
    if (count > 0) {
    	std::tie(chrId1, chrId2) = chrIndices[chrId];
        ref1VectorR.push_back(chrId1 + 1);
        ref2VectorR.push_back(chrId2 + 1);
        pos1VectorR.push_back(bin1);
        pos2VectorR.push_back(bin2);
        countVectorR.push_back(count);
    }
    for (auto &tmpFileName: tmpFileNames2) {
        if (remove(tmpFileName.c_str()) != 0) {
            std::cerr << "Error, cannot remove temporary file '" << tmpFileName << "'\n";
        }
    }
    //Rcout << "Counting done...\n";
    ref1VectorR.attr("class") = "factor";
    ref1VectorR.attr("levels") = refs;
    ref2VectorR.attr("class") = "factor";
    ref2VectorR.attr("levels") = refs;
    DataFrame outputData = DataFrame::create(_["ref1"] = ref1VectorR, _["bin1"] = pos1VectorR, _["ref2"] = ref2VectorR, _["bin2"] = pos2VectorR, _["count"] = countVectorR);
    return outputData;
}
