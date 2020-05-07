#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#include <fstream>
#include <sstream>
#include <string>
#include <cstdint>
#include <vector>
#include <unordered_map>


uint64_t getKey (const uint32_t i, const uint32_t j) {
  return (static_cast<uint64_t>(i) << 32) | j;
}

uint32_t getFirstKey (const uint64_t k) {
  return static_cast<uint32_t>(k >> 32);
}

uint32_t getSecondKey (const uint64_t k) {
  return static_cast<uint32_t>(k & 0xFFFFFFFF);
}

void getValues (const uint64_t k, uint32_t &i, uint32_t &j) {
  i = static_cast<uint32_t>(k >> 32);
  j = static_cast<uint32_t>(k);
}

void updateMatrixCounts(std::unordered_map < uint64_t, int > &matrixCounts,
                       std::vector < int >                   &currentChrIds) {
  if (currentChrIds.size() <= 1) {
    return;
  }
  int chrId1, chrId2;
  for (size_t i1 = 0; i1 < currentChrIds.size(); ++i1) {
    for (size_t i2 = 0; i2 < i1; ++i2) {
      chrId1 = currentChrIds[i1];
      chrId2 = currentChrIds[i2];
      if (chrId2 > chrId1) {
        std::swap(chrId1, chrId2);
      }
      if (chrId1 != chrId2) {
        ++matrixCounts[getKey(chrId1, chrId2)];
      }
    }
  }
}

void updateMatrices(std::unordered_map < uint64_t, arma::SpMat < int > > &matrices,
                    std::unordered_set < uint64_t >                      &validMatrices,
                    std::vector < int >                                  &currentChrIds,
                    std::vector < int >                                  &currentBins) {
  if (currentChrIds.size() <= 1) {
    return;
  }
  int chrId1, chrId2, bin1, bin2;
  for (size_t i1 = 0; i1 < currentChrIds.size(); ++i1) {
    for (size_t i2 = 0; i2 < i1; ++i2) {
      chrId1 = currentChrIds[i1];
      chrId2 = currentChrIds[i2];
      bin1   = currentBins[i1];
      bin2   = currentBins[i2];
      if (chrId1 != chrId2) {
        if (chrId2 > chrId1) {
          std::swap(chrId1, chrId2);
          std::swap(bin1, bin2);
        }
        if (validMatrices.find(getKey(chrId1, chrId2)) != validMatrices.end()) {
          ++matrices[getKey(chrId1, chrId2)](bin1, bin2);
        }
      }
    }
  }
}

// [[Rcpp::export]]
DataFrame parsePafCpp(std::string &fname, int resolution, int minAlnLen, int minCount, int minNCells) {
  int nChrs = 0;
  std::unordered_map <std::string, int> chrIds;
  std::unordered_map < uint64_t, int > matrixCounts;
  std::unordered_set < uint64_t > validMatrices;
  std::unordered_map < uint64_t, arma::SpMat < int > > matrices;
  std::vector < int > chrSizes;
  std::vector < int > stdChrs1, stdChrs2, stdBins1, stdBins2, stdCounts;
  IntegerVector chrs1, chrs2;
  IntegerVector bins1, bins2, counts;
  CharacterVector chrs;
  std::string currentReadName;
  std::vector < int > currentChrIds, currentBins;
  std::ifstream infile(fname);
  std::string line;
  for (long lineId = 0; std::getline(infile, line); ++lineId) {
    std::string chr, readName;
    int chrId;
    long start, end, size, dummyLong;
    char strand;
    std::istringstream iss(line);
    if (lineId % 1000000 == 0) {
      Rcerr << lineId << " lines read.\n";
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
    auto chrIt = chrIds.find(chr);
    if (chrIt == chrIds.end()) {
      chrId = nChrs;
      chrIds[chr] = chrId;
      ++nChrs;
    }
    else {
      chrId = chrIt->second;
    }
    size  = end - start + 1;
    start = start / resolution;
    end   = end   / resolution;
    if (chrId >= static_cast<int>(chrSizes.size())) {
      chrSizes.resize(chrId+1);
    }
    if (end > chrSizes[chrId]) {
      chrSizes[chrId] = end;
    }
    if (currentReadName != readName) {
      updateMatrixCounts(matrixCounts, currentChrIds);
      currentChrIds.clear();
      currentReadName = readName;
    }
    if (size >= minAlnLen) {
      currentChrIds.push_back(chrId);
    }
  }
  updateMatrixCounts(matrixCounts, currentChrIds);
  currentChrIds.clear();
  Rcerr << "First pass over, " << chrSizes.size() << " refs found.\nStarting second pass.\n";
  for (auto &chrIdPair: matrixCounts) {
    if (chrIdPair.second >= minCount) {
      uint32_t chrId1 = getFirstKey(chrIdPair.first);
      uint32_t chrId2 = getSecondKey(chrIdPair.first);
      validMatrices.insert(chrIdPair.first);
      matrices[chrIdPair.first] = arma::SpMat <int> (chrSizes[chrId1]+1, chrSizes[chrId2]+1);
    }
  }
  Rcerr << matrices.size() << " valid matrices.\n";
  currentReadName.clear();
  infile.clear();
  infile.seekg(0);
  for (long lineId = 0; std::getline(infile, line); ++lineId) {
    std::string chr, readName;
    int chrId;
    long start, end, bin, size, dummyLong;
    char strand;
    std::istringstream iss(line);
    if (lineId % 1000000 == 0) {
      Rcerr << lineId << " lines read.\n";
    }
    iss >> readName >> dummyLong >> dummyLong >> dummyLong >> strand >> chr >> dummyLong >> start >> end;
    chrId = chrIds[chr];
    bin   = ((start + end) / 2) / resolution;
    size  = end - start + 1;
    if (currentReadName != readName) {
      updateMatrices(matrices, validMatrices, currentChrIds, currentBins);
      currentChrIds.clear();
      currentBins.clear();
      currentReadName = readName;
    }
    if (size >= minAlnLen) {
      currentChrIds.push_back(chrId);
      currentBins.push_back(bin);
    }
  }
  updateMatrices(matrices, validMatrices, currentChrIds, currentBins);
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
    arma::SpMat < int > &matrix = matrixId.second;
    if (static_cast<int>(matrix.n_nonzero) >= minNCells) {
      std::vector < int > tmpChrs1(matrix.n_nonzero, chrId1);
      std::vector < int > tmpChrs2(matrix.n_nonzero, chrId2);
      std::vector < int > tmpBins1;
      std::vector < int > tmpBins2;
      std::vector < int > tmpCounts;
      tmpBins1.reserve(matrix.n_nonzero);
      tmpBins2.reserve(matrix.n_nonzero);
      tmpCounts.reserve(matrix.n_nonzero);
      for (arma::SpMat < int >::const_iterator it = matrix.begin(); it != matrix.end(); ++it) {
        tmpBins1.push_back(it.row());
        tmpBins2.push_back(it.col());
        tmpCounts.push_back(*it);
      }
      stdChrs1.insert(stdChrs1.end(), tmpChrs1.begin(), tmpChrs1.end());
      stdChrs2.insert(stdChrs2.end(), tmpChrs2.begin(), tmpChrs2.end());
      stdBins1.insert(stdBins1.end(), tmpBins1.begin(), tmpBins1.end());
      stdBins2.insert(stdBins2.end(), tmpBins2.begin(), tmpBins2.end());
      stdCounts.insert(stdCounts.end(), tmpCounts.begin(), tmpCounts.end());
    }
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
  return DataFrame::create(_["ref1"]  = chrs1,
                           _["bin1"]  = bins1,
                           _["ref2"]  = chrs2,
                           _["bin2"]  = bins2,
                           _["count"] = counts);
}