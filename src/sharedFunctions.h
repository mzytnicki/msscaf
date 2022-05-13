#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include <Rcpp.h>

using namespace Rcpp;                                                                                                                                                                                              

void computeRefsBins (IntegerVector &sizes, int nElements, IntegerVector &refs, IntegerVector &bins);
void computeRefsBinsMeta (IntegerVector &sizes, int nElements, int metaSize, bool metaBins, IntegerVector &refs, IntegerVector &bins);
long int computeOffsets (IntegerVector &sizes, std::vector < long int > &offsets);
long int computeOffsetsMeta (IntegerVector &sizes, int metaSize, std::vector < long int > &offsets);
bool caseInsensitiveEqual(std::string &s1, std::string &s2);
void setOutlierBinsToBool (IntegerVector &refs, IntegerVector &bins, std::vector < long int > &offsets, long int genomeSize, std::vector < bool > &outlierBins);
void setFullMatrixMeta (IntegerVector &refs1, IntegerVector &refs2, IntegerVector &bins1, IntegerVector &bins2, IntegerVector &counts, std::vector < long int > &offsets, long int metaGenomeSize, std::vector < bool > &outlierBinsBool, int distance, int metaSize, std::vector < std::vector < double > > &fullMatrix);
void setFullMatrix (IntegerVector &refs1, IntegerVector &refs2, IntegerVector &bins1, IntegerVector &bins2, IntegerVector &counts, std::vector < long int > &offsets, long int genomeSize, std::vector < bool > &outlierBinsBool, int distance, std::vector < std::vector < double > > &fullMatrix);
void computeDistanceCounts (IntegerVector &refs1, IntegerVector &refs2, IntegerVector &bins1, IntegerVector &bins2, IntegerVector &counts, std::vector < long int > &offsets, std::vector < bool > &outlierBinsBool, int distance, IntegerVector &sizes, std::vector < long int > &sums, std::vector < int > &nCounts);

#endif                                                                                                                                                                                                             
