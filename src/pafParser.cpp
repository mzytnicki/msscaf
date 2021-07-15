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
  interval_t (uint32_t c, uint32_t s, uint32_t e): chrId(c), start(s), end(e) {
    if (start > end) {
      Rcpp::stop("Error, end (" + std::to_string(end) + ") is before start (" + std::to_string(start) + ").");
    }
  }
  interval_t (): chrId(-1), start(-1), end(-1) {}
};

bool operator<(const interval_t &lhs, const interval_t &rhs) {
  if (lhs.chrId < rhs.chrId) return true;
  if (lhs.chrId > rhs.chrId) return false;
  if (lhs.start < rhs.start) return true;
  if (lhs.start > rhs.start) return false;
  return (lhs.end < rhs.end);
}

bool operator==(const interval_t &lhs, const interval_t &rhs) {
  return ((lhs.chrId == rhs.chrId) && (lhs.start == rhs.start) && (lhs.start == rhs.start));
}

std::ostream& operator<<(std::ostream& os, interval_t const &interval) {
  return os << interval.chrId << ":" << interval.start << "-" << interval.end;
}

int median(std::vector<int> &v) {
  if (v.size() == 0) return 0;
  size_t n = v.size() * 0.9;
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

using matrix_t = std::unordered_map < uint64_t, unsigned int >;

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

void updateMatrices(std::unordered_map < uint64_t, matrix_t > &matrices,
                    std::unordered_set < uint64_t >           &validMatrices,
                    std::vector < interval_t >                &intervals) {
  // Rcerr << "Entering UM with " << intervals.size() << " intervals\n";
  // for (auto &i: intervals) Rcerr << "\t" << i << "\n";
  std::sort(intervals.begin(), intervals.end());
  auto it = std::unique(intervals.begin(), intervals.end());
  intervals.resize(std::distance(intervals.begin(), it));
  for (size_t intervalId1 = 0; intervalId1 < intervals.size(); ++intervalId1) {
    interval_t &interval1 = intervals[intervalId1];
    if (interval1.start > interval1.end) Rcerr << "Problem in 'updateMatrices': [" << intervalId1 << "] " << interval1 << "\n";
    //Rcerr << "\tInt 1 ["  << intervalId1 << "]: " << interval1.chrId << ": " << interval1.start << "-" << interval1.end << "\n";
    for (size_t intervalId2 = 0; intervalId2 <= intervalId1; ++intervalId2) {
      interval_t &interval2 = intervals[intervalId2];
      //Rcerr << "\t\tInt 2 ["  << intervalId2 << "]: " << interval2.chrId << ": " << interval2.start << "-" << interval2.end << "\n";
      if (interval1.chrId < interval2.chrId) Rcerr << "Problem in 'updateMatrices': [" << intervalId1 << "] " << interval1 << ", [" << intervalId2 << "] " << interval2 << "\n";
      if (interval1.chrId > interval2.chrId) {
        if (validMatrices.find(getKey(interval1.chrId, interval2.chrId)) != validMatrices.end()) {
          for (uint32_t bin1 = interval1.start; bin1 <= interval1.end; ++bin1) {
            for (uint32_t bin2 = interval2.start; bin2 <= interval2.end; ++bin2) {
              ++matrices[getKey(interval1.chrId, interval2.chrId)][getKey(bin1, bin2)];
              //Rcerr << "Adding (2) " << bin1 << "-" << bin2 << " to matrix " << interval1.chrId << "-" << interval2.chrId << "\n";
            }
          }
        }
        //else { Rcerr << "2: " << interval1.chrId << "-" << interval2.chrId << " does not fit.\n"; }
      }
      else {
        if (validMatrices.find(getKey(interval1.chrId, interval2.chrId)) != validMatrices.end()) {
          for (uint32_t bin1 = interval1.start; bin1 <= interval1.end; ++bin1) {
            for (uint32_t bin2 = interval2.start; bin2 <= interval2.end; ++bin2) {
              ++matrices[getKey(interval1.chrId, interval2.chrId)][getKey(std::max<uint32_t>(bin1, bin2), std::min<uint32_t>(bin1, bin2))];
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
List parsePafCpp(std::string &fname, uint32_t resolution, int minAlnLen, int minCount, int minNCells) {
  uint32_t nChrs = 0;
  std::unordered_map <std::string, uint32_t> chrIds;
  std::unordered_map < uint64_t, int > matrixCounts;
  std::unordered_set < uint64_t > validMatrices;
  std::unordered_map < uint64_t, matrix_t > matrices;
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
    uint32_t chrId;
    uint32_t start, end, size, dummyLong;
    char strand;
    std::istringstream iss(line);
    if (nLines % 1000000 == 0) {
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
    if (end < start) {
      std::swap<uint32_t>(start, end);
    }
    size  = end - start + 1;
    start = start / resolution;
    end   = end   / resolution;
    if (start > end) {
      Rcpp::stop("Error!  End (" + std::to_string(end) + ") is before start (" + std::to_string(start) + ").");
    }
    auto chrIt = chrIds.find(chr);
    if (chrIt == chrIds.end()) {
      chrId = nChrs;
      chrIds[chr] = chrId;
      ++nChrs;
      chrSizes.push_back(end, chr);
      //Rcerr << "new " << chr << " => " << chrId << "\n";
    }
    else if (end > static_cast<uint32_t>(chrSizes[chrIt->second])) {
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
    if (size >= static_cast<uint32_t>(minAlnLen)) {
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
      //Rcerr << "valid matrix: " << chrId1 << "-" << chrId2 << ": " << chrIdPair.second << "/" << minCount << "\n";
      validMatrices.insert(chrIdPair.first);
    }
  }
  //for (auto &c: chrIds) { Rcerr << c.second << ": " << c.first << "\n"; }
  Rcerr << matrixCounts.size() << " matrices in total, " << validMatrices.size() << " are valid.\n";
  currentReadName.clear();
  infile.clear();
  infile.seekg(0);
  uint32_t previousSize = 0;
  Progress progress(nLines, true);
  while(std::getline(infile, line)) {
    std::string chr, readName;
    uint32_t chrId;
    uint32_t start, end, size, dummyLong;
    char strand;
    std::istringstream iss(line);
    iss >> readName >> dummyLong >> dummyLong >> dummyLong >> strand >> chr >> dummyLong >> start >> end;
    chrId = chrIds[chr];
    size  = end - start + 1;
    start = start / resolution;
    end   = end   / resolution;
    if (start > end) {
      Rcpp::stop("Error!  End (" + std::to_string(end) + ") is before start (" + std::to_string(start) + ").");
    }
    if (currentReadName != readName) {
      if (intervals.size() == 1) {
        //Rcerr << "start: " << start << ", end: " << end << " resolution: " << resolution << ", size: " << (previousSize / resolution) << "\n";
        moleculeSizes.push_back(previousSize / resolution);
      }
      updateMatrices(matrices, validMatrices, intervals);
      intervals.clear();
      currentReadName = readName;
    }
    if (size >= static_cast<uint32_t>(minAlnLen)) {
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
  Rcerr << "Merging matrices.\n";
  Progress progress2(matrices.size(), true);
  for (auto &matrixId: matrices) {
    uint32_t chrId1 = getFirstKey(matrixId.first);
    uint32_t chrId2 = getSecondKey(matrixId.first);
    matrix_t &matrix = matrixId.second;
    if (static_cast<int>(matrix.size()) >= minNCells) {
      std::vector < int > tmpChrs1(matrix.size(), chrId1);
      std::vector < int > tmpChrs2(matrix.size(), chrId2);
      std::vector < int > tmpBins1;
      std::vector < int > tmpBins2;
      std::vector < int > tmpCounts;
      tmpBins1.reserve(matrix.size());
      tmpBins2.reserve(matrix.size());
      tmpCounts.reserve(matrix.size());
      for (auto &it: matrix) {
        tmpBins1.push_back(getFirstKey(it.first));
        tmpBins2.push_back(getSecondKey(it.first));
        tmpCounts.push_back(it.second);
      }
      stdChrs1.insert(stdChrs1.end(), tmpChrs1.begin(), tmpChrs1.end());
      stdChrs2.insert(stdChrs2.end(), tmpChrs2.begin(), tmpChrs2.end());
      stdBins1.insert(stdBins1.end(), tmpBins1.begin(), tmpBins1.end());
      stdBins2.insert(stdBins2.end(), tmpBins2.begin(), tmpBins2.end());
      stdCounts.insert(stdCounts.end(), tmpCounts.begin(), tmpCounts.end());
    }
    //else { Rcerr << "Matrix has too few elements: " << matrix.size() << "/" << minNCells << "\n"; }
    progress2.increment();
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
