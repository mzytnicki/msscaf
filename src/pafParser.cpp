#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]


using namespace Rcpp;

#include <fstream>
#include <sstream>
#include <string>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <progress.hpp>
#include <progress_bar.hpp>

struct interval_t {
  uint32_t chrId;
  uint32_t start;
  uint32_t end;
  interval_t (uint32_t c, uint32_t s, uint32_t e): chrId(c), start(s), end(e) {}
};


int median(std::vector<int> &v) {
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

class SparseMatrix {
public:
  size_t nRows;
  size_t nCols;
private:
  std::unordered_map<uint64_t, unsigned int> matrix;
public:
  void setSizes (size_t nr, size_t nc) {
    nRows = nr;
    nCols = nc;
  }
  SparseMatrix (size_t nr = 0, size_t nc = 0) {
    setSizes(nr, nc);
  }
  void addElement (uint32_t r, uint32_t c) {
    ++matrix[getKey(r, c)];
  }
  size_t getNElements () {
    return matrix.size();
  }
  void clear () {
    matrix.clear();
  }
  std::unordered_map<uint64_t, unsigned int>::const_iterator begin() {
    return matrix.cbegin();
  }
  std::unordered_map<uint64_t, unsigned int>::const_iterator end() {
    return matrix.cend();
  }
};


void updateMatrixCounts(std::unordered_map < uint64_t, int > &matrixCounts,
                       std::vector < interval_t >            &intervals) {
  for (size_t intervalId1 = 0; intervalId1 < intervals.size(); ++intervalId1) {
    interval_t &interval1 = intervals[intervalId1];
    uint32_t    size1     = interval1.end - interval1.start + 1;
    for (size_t intervalId2 = 0; intervalId2 <= intervalId1; ++intervalId2) {
      interval_t &interval2 = intervals[intervalId2];
      uint32_t    size2     = interval2.end - interval2.start + 1;
      if (interval2.chrId > interval1.chrId) {
        matrixCounts[getKey(interval2.chrId, interval1.chrId)] += size1 * size2;
      }
      else {
        matrixCounts[getKey(interval1.chrId, interval2.chrId)] += size1 * size2;
      }
    }
  }
}

void updateMatrices(std::unordered_map < uint64_t, SparseMatrix > &matrices,
                    std::unordered_set < uint64_t >               &validMatrices,
                    std::vector < interval_t >                    &intervals) {
  for (size_t intervalId1 = 0; intervalId1 < intervals.size(); ++intervalId1) {
    interval_t &interval1 = intervals[intervalId1];
    for (size_t intervalId2 = 0; intervalId2 <= intervalId1; ++intervalId2) {
      interval_t &interval2 = intervals[intervalId2];
      if (interval2.chrId > interval1.chrId) {
        if (validMatrices.find(getKey(interval2.chrId, interval1.chrId)) != validMatrices.end()) {
          for (uint32_t bin2 = interval2.start; bin2 <= interval2.end; ++bin2) {
            for (uint32_t bin1 = interval1.start; bin1 <= interval1.end; ++bin1) {
              if ((interval1.chrId == interval2.chrId) && (bin2 > bin1)) {
                matrices[getKey(interval2.chrId, interval1.chrId)].addElement(bin2, bin1);
              }
              else {
                matrices[getKey(interval2.chrId, interval1.chrId)].addElement(bin1, bin2);
              }
              //Rcerr << "Adding (1) " << bin2 << "-" << bin1 << " to matrix " << interval2.chrId << "-" << interval1.chrId << "\n";
            }
          }
        }
        //else { Rcerr << "1: " << interval2.chrId << "-" << interval1.chrId << " does not fit.\n"; }
      }
      else {
        if (validMatrices.find(getKey(interval1.chrId, interval2.chrId)) != validMatrices.end()) {
          for (uint32_t bin1 = interval1.start; bin1 <= interval1.end; ++bin1) {
            for (uint32_t bin2 = interval2.start; bin2 <= interval2.end; ++bin2) {
              if ((interval1.chrId == interval2.chrId) && (bin2 > bin1)) {
                matrices[getKey(interval1.chrId, interval2.chrId)].addElement(bin2, bin1);
              }
              else {
                matrices[getKey(interval1.chrId, interval2.chrId)].addElement(bin1, bin2);
              }
              //Rcerr << "Adding (2) " << bin1 << "-" << bin2 << " to matrix " << interval1.chrId << "-" << interval2.chrId << "\n";
            }
          }
        }
        //else { Rcerr << "2: " << interval1.chrId << "-" << interval2.chrId << " does not fit.\n"; }
      }
    }
  }
}

// [[Rcpp::export]]
List parsePafCpp(std::string &fname, int resolution, int minAlnLen, int minCount, int minNCells) {
  int nChrs = 0;
  std::unordered_map <std::string, int> chrIds;
  std::unordered_map < uint64_t, int > matrixCounts;
  std::unordered_set < uint64_t > validMatrices;
  std::unordered_map < uint64_t, SparseMatrix > matrices;
  IntegerVector chrSizes;
  std::vector < int > moleculeSizes;
  std::vector < int > stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts;
  IntegerVector chrs1, chrs2;
  IntegerVector bins1, bins2, counts;
  CharacterVector chrs;
  std::string currentReadName;
  std::vector < interval_t > intervals;
  std::ifstream infile(fname);
  std::string line;
  long nLines      = 0;
  long nLongEnough = 0;
  // First pass: get the size of the chromosomes, count the number of reads per matrix
  for (nLines = 0; std::getline(infile, line); ++nLines) {
    std::string chr, readName;
    int chrId;
    long start, end, size, dummyLong;
    char strand;
    std::istringstream iss(line);
    if (nLines % 10000000 == 0) {
      Rcerr << "\t" << nLines << " lines read.\n";
    }
    if (! (iss >> readName)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> dummyLong)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> dummyLong)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> dummyLong)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> strand)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> chr)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> dummyLong)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> start)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    if (! (iss >> end)) {
      Rcpp::stop("Cannot read line " + line + ".");
    }
    size  = end - start + 1;
    start = start / resolution;
    end   = end   / resolution;
    auto chrIt = chrIds.find(chr);
    if (chrIt == chrIds.end()) {
      chrId = nChrs;
      chrIds[chr] = chrId;
      ++nChrs;
      chrSizes.push_back(end, chr);
      //Rcerr << "new " << chr << " => " << chrId << "\n";
    }
    else if (end > chrSizes[chrIt->second]) {
      chrId = chrIt->second;
      chrSizes[chrId] = end;
    }
    //Rcerr << "reads: " << readName << " " << currentReadName;
    //for (auto &i: intervals) Rcerr << " " << i.chrId << ":" << i.start << "-" << i.end;
    //Rcerr << "\n";
    if (currentReadName != readName) {
      updateMatrixCounts(matrixCounts, intervals);
      intervals.clear();
      currentReadName = readName;
    }
    if (size >= minAlnLen) {
      //Rcerr << "adding " << chrId << "\n";
      ++nLongEnough;
      intervals.emplace_back(chrId, start, end);
    }
  }
  updateMatrixCounts(matrixCounts, intervals);
  Rcerr << "\t" << nLines << " lines read, " << nLongEnough << " passed the length threshold.\n";
  intervals.clear();
  Rcerr << "First pass over, " << chrSizes.size() << " refs found.\nStarting second pass.\n";
  for (auto &chrIdPair: matrixCounts) {
    if (chrIdPair.second >= minCount) {
      uint32_t chrId1 = getFirstKey(chrIdPair.first);
      uint32_t chrId2 = getSecondKey(chrIdPair.first);
      //Rcerr << "valid matrix: " << chrId1 << "-" << chrId2 << ": " << chrIdPair.second << "/" << minCount << "\n";
      validMatrices.insert(chrIdPair.first);
      matrices[chrIdPair.first] = SparseMatrix(chrSizes[chrId1]+1, chrSizes[chrId2]+1);
    }
  }
  //for (auto &c: chrIds) { Rcerr << c.second << ": " << c.first << "\n"; }
  Rcerr << matrixCounts.size() << " matrices in total, " << validMatrices.size() << " are valid.\n";
  currentReadName.clear();
  infile.clear();
  infile.seekg(0);
  long previousSize = 0;
  Progress progress(nLines, true);
  while(std::getline(infile, line)) {
    std::string chr, readName;
    int chrId;
    long start, end, size, dummyLong;
    char strand;
    std::istringstream iss(line);
    iss >> readName >> dummyLong >> dummyLong >> dummyLong >> strand >> chr >> dummyLong >> start >> end;
    chrId = chrIds[chr];
    size  = end - start + 1;
    start = start / resolution;
    end   = end   / resolution;
    if (currentReadName != readName) {
      if (intervals.size() == 1) {
        //Rcerr << "start: " << start << ", end: " << end << " resolution: " << resolution << ", size: " << (previousSize / resolution) << "\n";
        moleculeSizes.push_back(previousSize / resolution);
      }
      updateMatrices(matrices, validMatrices, intervals);
      intervals.clear();
      currentReadName = readName;
    }
    if (size >= minAlnLen) {
      intervals.emplace_back(chrId, start, end);
    }
    previousSize = size;
    progress.increment();
  }
  updateMatrices(matrices, validMatrices, intervals);
  chrs = CharacterVector(chrIds.size());
  for (auto chrId: chrIds) {
    chrs[chrId.second] = chrId.first;
  }
  int nMatrices = 0;
  for (auto &matrixId: matrices) {
    if (nMatrices % 100000 == 0) {
      Rcerr << "Exporting matrix # " << nMatrices << " / " << matrices.size() << "\n";
    }
    uint32_t chrId1 = getFirstKey(matrixId.first);
    uint32_t chrId2 = getSecondKey(matrixId.first);
    SparseMatrix &matrix = matrixId.second;
    if (static_cast<int>(matrix.getNElements()) >= minNCells) {
      std::vector < int > tmpChrs1(matrix.getNElements(), chrId1);
      std::vector < int > tmpChrs2(matrix.getNElements(), chrId2);
      std::vector < int > tmpBins1;
      std::vector < int > tmpBins2;
      std::vector < int > tmpCounts;
      tmpBins1.reserve(matrix.getNElements());
      tmpBins2.reserve(matrix.getNElements());
      tmpCounts.reserve(matrix.getNElements());
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
    //else { Rcerr << "Matrix has too few elements: " << matrix.getNElements() << "/" << minNCells << "\n"; }
    ++nMatrices;
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
  DataFrame outputDataFrame = DataFrame::create(_["ref1"]  = chrs1,
                                                _["bin1"]  = bins1,
                                                _["ref2"]  = chrs2,
                                                _["bin2"]  = bins2,
                                                _["count"] = counts);
  return List::create(_["data"] = outputDataFrame,
                      _["size"] = median(moleculeSizes),
                      _["sizes"] = chrSizes);
}
