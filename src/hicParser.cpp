// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

/*
 The MIT License (MIT)
 
 Copyright (c) 2011-2016 Broad Institute, Aiden Lab
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */
#include <cstring>
#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <streambuf>
#include "zlib.h"

/*
 Straw: fast C++ implementation of dump. Not as fully featured as the
 Java version. Reads the .hic file, finds the appropriate matrix and slice
 of data, and outputs as text in sparse upper triangular format.
 
 Currently only supporting matrices.
 */
// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf
{
    membuf(char* begin, char* end) {
        this->setg(begin, begin, end);
    }
};

// stores input information
struct hicInfo {
    long footerPos;
    std::vector <long> availableResolutions;
    int resolution;
    int resolutionIdSelected;
    int version;
    CharacterVector chrs;
    int nChrs;
    bool firstChromosomeAll;
};

// stores output information
struct outputStr {
    std::vector<int> chr1;
    std::vector<int> chr2;
    std::vector<int> bin1;
    std::vector<int> bin2;
    std::vector<int> count;
    /**
     * Transform chr id vectors to R factors
     */
    void append(int ch1, int b1, int ch2, int b2, int co) {
        // ch 0 should be ALL
        if (ch1 == 0) {
            if (ch2 != 0) {
                stop("The chromosome ALL should only match with an ALL chromosome");
            }
            return;
        }
        // factors start with 1 in R
        chr1.push_back(ch1);
        chr2.push_back(ch2);
        bin1.push_back(b1);
        bin2.push_back(b2);
        count.push_back(co);
    }
};


// returns whether or not this is valid HiC file
bool readMagicString(std::istream& fin) {
    std::string str;
    getline(fin, str, '\0');
    return ((str[0]=='H') && (str[1]=='I') && (str[2]=='C'));
}

char readCharFromFile(std::istream &fin) {
    char tempChar;
    fin.read(&tempChar, sizeof(char));
    return tempChar;
}

int16_t readInt16FromFile(std::istream &fin) {
    int16_t tempInt16;
    fin.read((char *) &tempInt16, sizeof(int16_t));
    return tempInt16;
}

int32_t readInt32FromFile(std::istream &fin) {
    int32_t tempInt32;
    fin.read((char *) &tempInt32, sizeof(int32_t));
    return tempInt32;
}

int64_t readInt64FromFile(std::istream &fin) {
    int64_t tempInt64;
    fin.read((char *) &tempInt64, sizeof(int64_t));
    return tempInt64;
}

float readFloatFromFile(std::istream &fin) {
    float tempFloat;
    fin.read((char *) &tempFloat, sizeof(float));
    return tempFloat;
}

double readDoubleFromFile(std::istream &fin) {
    double tempDouble;
    fin.read((char *) &tempDouble, sizeof(double));
    return tempDouble;
}

void readHeader(std::istream &fin, hicInfo &info) {
    info.resolutionIdSelected = -1;
    if (!readMagicString(fin)) {
        stop("Hi-C magic string is missing, does not appear to be a hic file.");
    }
    info.version = readInt32FromFile(fin);
    if (info.version < 6) {
        stop("Version " + std::to_string(info.version) + " no longer supported.");
    }
    std::string genome;
    info.footerPos = readInt64FromFile(fin);
    getline(fin, genome, '\0');
    if (info.version > 8) {
        readInt64FromFile(fin); // nviPosition
        readInt64FromFile(fin); // nviLength
    }
    int32_t totalAttributes = readInt32FromFile(fin);
    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < totalAttributes; i++) {
        std::string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
    }
    info.nChrs = readInt32FromFile(fin);
    // chromosome map for finding matrix
    for (int i = 0; i < info.nChrs; i++) {
        std::string name;
        int32_t length;
        getline(fin, name, '\0');
        if (info.version > 8) {
            length = readInt64FromFile(fin);
        } else {
            length = (int64_t) readInt32FromFile(fin);
        }
        info.chrs.push_back(name);
    }
    int32_t totalResolutions = readInt32FromFile(fin);
    for (int i = 0; i < totalResolutions; i++) {
        int32_t resolution = readInt32FromFile(fin);
        info.availableResolutions.push_back(resolution);
        if (resolution == info.resolution) {
            info.resolutionIdSelected = i;
        }
    }
    info.firstChromosomeAll = (
        info.chrs[0] == "ALL" || info.chrs[0] == "All"
    );
}



/*
void readHeader(std::istream& fin, hicInfo &info) {
    info.resolutionIdSelected = -1;
    if (! readMagicString(fin)) {
        stop("Hi-C magic string is missing, does not appear to be a hic file.");
    }
    fin.read((char*)& info.version, sizeof(int));
    Rcerr << "Hi-C format version: " << info.version << "\n";
    if (info.version < 6) {
        stop("Version " + std::to_string(info.version) + " no longer supported.");
    }
    std::string genome;
    int nattributes;
    fin.read((char*) &info.footerPos, sizeof(long));
    getline(fin, genome, '\0' );
    Rcerr << "Genome: " << genome << "\n";
    if (info.version >= 9) {
        long normVectorIndexPosition;
        fin.read((char*) &normVectorIndexPosition, sizeof(long));
        long normVectorIndexLength;
        fin.read((char*) &normVectorIndexLength, sizeof(long));
    }
    fin.read((char*)&nattributes, sizeof(int));
    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < nattributes; i++) {
        std::string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
    }
    fin.read((char*) &info.nChrs, sizeof(int));
    // chromosome map for finding matrix
    for (int i = 0; i < info.nChrs; i++) {
        std::string name;
        int length;
        getline(fin, name, '\0');
        fin.read((char*) &length, sizeof(int));
        info.chrs.push_back(name);
    }
    int nBpResolutions;
    fin.read((char*) &nBpResolutions, sizeof(int));
    for (int i = 0; i < nBpResolutions; i++) {
        if (info.version < 9) {
            int bpResolution;
            fin.read((char*) &bpResolution, sizeof(int));
            info.availableResolutions.push_back(bpResolution);
            if (bpResolution == info.resolution) {
                info.resolutionIdSelected = i;
            }
        }
        else {
            long bpResolution;
            fin.read((char*) &bpResolution, sizeof(long));
            info.availableResolutions.push_back(bpResolution);
            if (bpResolution == info.resolution) {
                info.resolutionIdSelected = i;
            }
        }
    }
    // The 'ALL' chromosome is useless.
    info.firstChromosomeAll = (info.chrs[0] == "ALL");
    // fragment resolutions: not used
    int nFragResolutions;
    fin.read((char*) &nFragResolutions, sizeof(int));
    Rcerr << "# frag. resolutions: " << nFragResolutions << "\n";
    for (int i = 0; i < nFragResolutions; i++) {
        int resFrag;
        fin.read((char*) &resFrag, sizeof(int));
    }
    // fragment site positions per chromosome: not used
    if (nFragResolutions > 0) {
        for (int i = 0; i < info.nChrs; i++) {
            int nSites;
            fin.read((char*) &nSites, sizeof(int));
            for (int i = 0; i < nSites; i++) {
                int sitePosition;
                fin.read((char*) &sitePosition, sizeof(int));
            }
        }
    }
}
*/


void readBlock(
        std::istream &fin,
        int64_t position,
        int32_t size,
        int32_t chrId1,
        int32_t chrId2,
        hicInfo &info,
        outputStr &output
) {
    
    if (size == 0) {
        return;
    }
    std::vector<int> chrIds1, chrIds2, bins1, bins2, counts;
    
    char* compressedBytes = new char[size];
    char* uncompressedBytes = new char[size*10]; // biggest seen so far is 3
    
    fin.seekg(position, std::ios::beg);
    fin.read(compressedBytes, size);
    
    // Decompress the block
    // zlib struct
    z_stream infstream;
    infstream.zalloc    = Z_NULL;
    infstream.zfree     = Z_NULL;
    infstream.opaque    = Z_NULL;
    infstream.avail_in  = (uInt)(size); // size of input
    infstream.next_in   = (Bytef *) compressedBytes; // input char array
    infstream.avail_out = (uInt)size*10; // size of output
    infstream.next_out  = (Bytef *)uncompressedBytes; // output char array
    // the actual decompression work
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    int uncompressedSize = infstream.total_out;
    
    // create stream from buffer for ease of use
    membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
    std::istream bufferIn(&sbuf);
    int32_t totalRecords = readInt32FromFile(bufferIn);
    bins1.reserve(totalRecords);
    bins2.reserve(totalRecords);
    counts.reserve(totalRecords);
    // different versions have different specific formats
    if (info.version < 7) {
        for (int i = 0; i < totalRecords; i++) {
            int32_t binX = readInt32FromFile(bufferIn);
            int32_t binY = readInt32FromFile(bufferIn);
            float  c     = readFloatFromFile(bufferIn);
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(c);
        }
    } else {
        int32_t binXOffset = readInt32FromFile(bufferIn);
        int32_t binYOffset = readInt32FromFile(bufferIn);
        bool    useShort   = readCharFromFile(bufferIn) == 0; // yes this is opposite of usual
        bool useShortBinX = true;
        bool useShortBinY = true;
        if (info.version > 8) {
            useShortBinX = readCharFromFile(bufferIn) == 0;
            useShortBinY = readCharFromFile(bufferIn) == 0;
        }
        char type = readCharFromFile(bufferIn);
        if (type == 1) {
            if (useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferIn);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferIn);
                    int16_t colCount = readInt16FromFile(bufferIn);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt16FromFile(bufferIn);
                        float c;
                        if (useShort) {
                            c = readInt16FromFile(bufferIn);
                        } else {
                            c = readFloatFromFile(bufferIn);
                        }
                        bins1.push_back(binX);
                        bins2.push_back(binY);
                        counts.push_back(c);
                    }
                }
            } else if (useShortBinX && !useShortBinY) {
                int32_t rowCount = readInt32FromFile(bufferIn);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferIn);
                    int16_t colCount = readInt16FromFile(bufferIn);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt16FromFile(bufferIn);
                        float c;
                        if (useShort) {
                            c = readInt16FromFile(bufferIn);
                        } else {
                            c = readFloatFromFile(bufferIn);
                        }
                        bins1.push_back(binX);
                        bins2.push_back(binY);
                        counts.push_back(c);
                    }
                }
            } else if (!useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferIn);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferIn);
                    int32_t colCount = readInt32FromFile(bufferIn);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt32FromFile(bufferIn);
                        float c;
                        if (useShort) {
                            c = readInt16FromFile(bufferIn);
                        } else {
                            c = readFloatFromFile(bufferIn);
                        }
                        bins1.push_back(binX);
                        bins2.push_back(binY);
                        counts.push_back(c);
                    }
                }
            } else {
                int32_t rowCount = readInt32FromFile(bufferIn);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferIn);
                    int32_t colCount = readInt32FromFile(bufferIn);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt32FromFile(bufferIn);
                        float c;
                        if (useShort) {
                            c = readInt16FromFile(bufferIn);
                        } else {
                            c = readFloatFromFile(bufferIn);
                        }
                        bins1.push_back(binX);
                        bins2.push_back(binY);
                        counts.push_back(c);
                    }
                }
            }
        } else if (type == 2) {
            int32_t nPts = readInt32FromFile(bufferIn);
            int16_t w = readInt16FromFile(bufferIn);
            for (int i = 0; i < nPts; i++) {
                int32_t row = i / w;
                int32_t col = i - row * w;
                int32_t bin1 = binXOffset + col;
                int32_t bin2 = binYOffset + row;
                if (useShort) {
                    int16_t c = readInt16FromFile(bufferIn);
                    if (c != -32768) {
                        bins1.push_back(bin1);
                        bins2.push_back(bin2);
                        counts.push_back(c);
                    }
                } else {
                    float c = readFloatFromFile(bufferIn);
                    if (!std::isnan(c)) {
                        bins1.push_back(bin1);
                        bins2.push_back(bin2);
                        counts.push_back(c);
                    }
                }
            }
        }
    }
    chrIds1 = std::vector<int>(bins1.size(), chrId1);
    chrIds2 = std::vector<int>(bins2.size(), chrId2);
    output.chr1.insert(output.chr1.end(),   chrIds1.begin(), chrIds1.end());
    output.chr2.insert(output.chr2.end(),   chrIds2.begin(), chrIds2.end());
    output.bin1.insert(output.bin1.end(),   bins1.begin(),   bins1.end());
    output.bin2.insert(output.bin2.end(),   bins2.begin(),   bins2.end());
    output.count.insert(output.count.end(), counts.begin(),  counts.end());
    delete[] compressedBytes;
    delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
}


/*
// this is the meat of reading the data.
// takes in the block number and returns the set of contact records corresponding to
// that block.
// the block data is compressed and must be decompressed using the zlib library functions
void readBlock(std::istream& fin, long position, int size, int chrId1, int chrId2, hicInfo &info, outputStr &output) {
    if (size == 0) {
        return;
    }
    std::vector<int> chrIds1, chrIds2, bins1, bins2, counts;
    
    char* compressedBytes = new char[size];
    char* uncompressedBytes = new char[size*10]; //biggest seen so far is 3
    
    fin.seekg(position, std::ios::beg);
    fin.read(compressedBytes, size);
    
    // Decompress the block
    // zlib struct
    z_stream infstream;
    infstream.zalloc    = Z_NULL;
    infstream.zfree     = Z_NULL;
    infstream.opaque    = Z_NULL;
    infstream.avail_in  = (uInt)(size); // size of input
    infstream.next_in   = (Bytef *) compressedBytes; // input char array
    infstream.avail_out = (uInt)size*10; // size of output
    infstream.next_out  = (Bytef *)uncompressedBytes; // output char array
    // the actual decompression work.
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    int uncompressedSize = infstream.total_out;
    
    // create stream from buffer for ease of use
    membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
    std::istream bufferin(&sbuf);
    int nRecords;
    bufferin.read((char*) &nRecords, sizeof(int));
    bins1.reserve(nRecords);
    bins2.reserve(nRecords);
    counts.reserve(nRecords);
    // different versions have different specific formats
    //Rcerr << "        Reading block with " << nRecords << " records\n";
    if (info.version < 7) {
        for (int i = 0; i < nRecords; i++) {
            int binX, binY;
            float count;
            bufferin.read((char*) &binX,   sizeof(int));
            bufferin.read((char*) &binY,   sizeof(int));
            bufferin.read((char*) &count, sizeof(float));
            //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (binX * info.resolution) << "\t" << (binY * info.resolution) << "\t" << counts << "\n";
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(count);
        }
    }
    else {
        int binXOffset, binYOffset;
        char useShort;
        char type;
        bufferin.read((char*) &binXOffset, sizeof(int));
        bufferin.read((char*) &binYOffset, sizeof(int));
        bufferin.read((char*) &useShort,   sizeof(char));
        bufferin.read((char*) &type,       sizeof(char));
        //Rcerr << "    Block " << binXOffset << "-" << binYOffset << " with format short " << static_cast<int>(useShort) << " and type " << static_cast<int>(type) << "\n";
        bool useShortBinX = true;
        bool useShortBinY = true;
        if (info.version > 8) {
            bufferin.read((char*) &useShortBinX, sizeof(char));
            bufferin.read((char*) &useShortBinY, sizeof(char));
        }
        if (type == 1) {
            // List-of-rows representation
            int rowCount;
            if (useShortBinY) {
                short tmp;
                bufferin.read((char*) &tmp, sizeof(short));
                rowCount = tmp;
            }
            else {
                bufferin.read((char*) &rowCount, sizeof(int));
            }
            // Rcerr << "      # rows: " << rowCount << "\n";
            for (int i = 0; i < rowCount; i++) {
                short y;
                int binY;
                short colCount;
                if (useShortBinY) {
                    short tmp;
                    bufferin.read((char*) &tmp, sizeof(short));
                    y = tmp;
                }
                else {
                    bufferin.read((char*) &y, sizeof(int));
                }
                binY = y + binYOffset;
                if (useShortBinX) {
                    short tmp;
                    bufferin.read((char*) &tmp, sizeof(short));
                    colCount = tmp;
                }
                else {
                    bufferin.read((char*) &colCount, sizeof(int));
                }
                // Rcerr << "        # cols: " << colCount << "\n";
                for (int j = 0; j < colCount; j++) {
                    short x;
                    int binX;
                    short c;
                    float count;
                    if (useShortBinX) {
                        short tmp;
                        bufferin.read((char*) &tmp, sizeof(short));
                        x = tmp;
                    }
                    else {
                        bufferin.read((char*) &x, sizeof(int));
                    }
                    binX = binXOffset + x;
                    if (useShort == 0) { // yes this is opposite of usual
                        bufferin.read((char*) &c, sizeof(short));
                        count = c;
                        bins1.push_back(binX);
                        bins2.push_back(binY);
                        counts.push_back(c);
                    }
                    else {
                        bufferin.read((char*) &count, sizeof(float));
                        bins1.push_back(binX);
                        bins2.push_back(binY);
                        counts.push_back(count);
                    }
                    //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (binX * info.resolution) << "\t" << (binY * info.resolution) << "\t" << counts << "\n";
                }
            }
        }
        else if (type == 2) { // have yet to find test file where this is true, possibly entirely deprecated
            int nPts;
            short w;
            bufferin.read((char*) &nPts, sizeof(int));
            bufferin.read((char*) &w, sizeof(short));
            
            for (int i = 0; i < nPts; i++) {
                int row = i / w;
                int col = i - row * w;
                int bin1 = binXOffset + col;
                int bin2 = binYOffset + row;
                float count;
                short c;
                if (useShort == 0) { // yes this is opposite of the usual
                    bufferin.read((char*) &c, sizeof(short));
                    if (c != -32768) {
                        //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (bin1 * info.resolution) << "\t" << (bin2 * info.resolution) << "\t" << c << "\n";
                        bins1.push_back(bin1);
                        bins2.push_back(bin2);
                        counts.push_back(c);
                    }
                }
                else {
                    bufferin.read((char*) &count, sizeof(float));
                    if (count != 0x7fc00000) { // not sure this works
                        //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (bin1 * info.resolution) << "\t" << (bin2 * info.resolution) << "\t" << counts << "\n";
                        bins1.push_back(bin1);
                        bins2.push_back(bin2);
                        counts.push_back(count);
                    }
                }
            }
        }
    }
    chrIds1 = std::vector<int>(bins1.size(), chrId1);
    chrIds2 = std::vector<int>(bins2.size(), chrId2);
    output.chr1.insert(output.chr1.end(),   chrIds1.begin(), chrIds1.end());
    output.chr2.insert(output.chr2.end(),   chrIds2.begin(), chrIds2.end());
    output.bin1.insert(output.bin1.end(),   bins1.begin(),   bins1.end());
    output.bin2.insert(output.bin2.end(),   bins2.begin(),   bins2.end());
    output.count.insert(output.count.end(), counts.begin(),  counts.end());
    
    //Rcerr << "  Done.\n";
    delete[] compressedBytes;
    delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
}
*/

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
void readMatrix(
        std::istream &fin,
        int64_t start,
        hicInfo &info,
        outputStr &output
) {
    
    std::streampos pos;
    if (start != -1) {
        fin.seekg(start, std::ios::beg);
        int32_t chromosomeId1 = readInt32FromFile(fin);
        int32_t chromosomeId2 = readInt32FromFile(fin);
        int32_t totalResolutions = readInt32FromFile(fin);
        if ((! info.firstChromosomeAll) || (chromosomeId1 != 0)) {
            for (int resolutionId = 0; resolutionId < totalResolutions; ++resolutionId) {
                std::string unit;
                getline(fin, unit, '\0');
                readInt32FromFile(fin); // resIdx
                readFloatFromFile(fin); // sumCounts
                readFloatFromFile(fin); // occupiedCellCount
                readFloatFromFile(fin); // stdDev
                readFloatFromFile(fin); // percent95
                readInt32FromFile(fin); // binSize
                readInt32FromFile(fin); // totalBlockBins
                readInt32FromFile(fin); // totalBlockColumns
                int32_t totalBlocks = readInt32FromFile(fin);
                for (int i = 0; i < totalBlocks; i++) {
                    readInt32FromFile(fin); // blockId
                    int64_t blockPosition = readInt64FromFile(fin);
                    int32_t blockSize     = readInt32FromFile(fin);
                    if (resolutionId == info.resolutionIdSelected) {
                        pos = fin.tellg();
                        readBlock(fin, blockPosition, blockSize, chromosomeId1, chromosomeId2, info, output);
                        fin.seekg(pos, std::ios::beg);
                    }
                }
            }
        }
    }
}


/*
void readMatrix(std::istream& fin, long start, int size, hicInfo &info, outputStr &output) {
    std::streampos pos;
    //Rcerr << "  Reading matrix start: " << start << ".\n";
    if (start != -1) {
        fin.seekg(start, std::ios::beg);
        int chrId1, chrId2, nResolutions;
        fin.read((char*) &chrId1, sizeof(int));
        fin.read((char*) &chrId2, sizeof(int));
        fin.read((char*) &nResolutions, sizeof(int));
        //Rcerr << "  Reading matrix " << info.chrs[chrId1] << "/" << info.chrs[chrId2] << ": " << nResolutions << " resolutions\n";
        for (int resolutionId = 0; resolutionId < nResolutions; ++resolutionId) {
            std::string unit;
            int resIdx;
            float tmp2;
            int binSize;
            int blockBinCount;
            int blockColumnCount;
            int blockCount;
            getline(fin, unit, '\0');
            fin.read((char*) &resIdx,           sizeof(int));
            //Rcerr << "    Resolution # " << resIdx << " in " << unit << "\n";
            fin.read((char*) &tmp2,             sizeof(float)); // sumCounts
            fin.read((char*) &tmp2,             sizeof(float)); // occupiedCellCount
            fin.read((char*) &tmp2,             sizeof(float)); // stdDev
            fin.read((char*) &tmp2,             sizeof(float)); // percent95
            fin.read((char*) &binSize,          sizeof(int));
            fin.read((char*) &blockBinCount,    sizeof(int));
            fin.read((char*) &blockColumnCount, sizeof(int));
            fin.read((char*) &blockCount,       sizeof(int));
            //Rcerr << "    # c blocks " << blockCount << "\n";
            for (int blockCountId = 0; blockCountId < blockCount; ++blockCountId) {
                int blockId, blockSize;
                long blockPosition;
                fin.read((char*) &blockId,       sizeof(int));
                fin.read((char*) &blockPosition, sizeof(long));
                fin.read((char*) &blockSize,     sizeof(int));
                //Rcerr << "      block # " << blockId << ": " << blockPosition << " / " << blockSize << "\n";
                if (resolutionId == info.resolutionIdSelected) {
                    //Rcerr << "        selected\n";
                    pos = fin.tellg();
                    readBlock(fin, blockPosition, blockSize, chrId1, chrId2, info, output);
                    fin.seekg(pos, std::ios::beg);
                }
            }
        }
    }
}
 */

// reads the footer from the footer pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
void readFooter(std::istream& fin, hicInfo &info, outputStr &output) {
    std::streampos pos;
    fin.seekg(info.footerPos, std::ios::beg);
    if (info.version > 8) {
        readInt64FromFile(fin); // totalBytes
    } else {
        readInt32FromFile(fin); // totalBytes
    }
    int32_t totalEntries = readInt32FromFile(fin);
    Progress progress(totalEntries, true);
    for (int i = 0; i < totalEntries; i++) {
        progress.increment();
        std::string str;
        getline(fin, str, '\0');
        int64_t fpos = readInt64FromFile(fin);
        readInt32FromFile(fin); // sizeInBytes
        pos = fin.tellg();
        readMatrix(fin, fpos, info, output);
        fin.seekg(pos, std::ios::beg);
    }
}

/*
void readFooter(std::istream& fin, hicInfo &info, outputStr &output) {
    std::streampos pos;
    fin.seekg(info.footerPos, std::ios::beg);
    if (info.version < 9) {
        int nBytes;
        fin.read((char*) &nBytes, sizeof(int));
    }
    else {
        long nBytes;
        fin.read((char*) &nBytes, sizeof(long));
    }
    int nEntries;
    fin.read((char*) &nEntries, sizeof(int));
    Progress progress(nEntries, true);
    Rcerr << "Found " << nEntries << " entries\n";
    for (int i = 0; i < nEntries; i++) {
        progress.increment();
        // if (i % 1000 == 0) {
        //   Rcerr << "Reading " << i << "/" << nEntries << " entries.\n";
        // }
        std::string str;
        getline(fin, str, '\0');
        long fpos;
        fin.read((char*)& fpos, sizeof(long));
        int sizeinbytes;
        fin.read((char*)& sizeinbytes, sizeof(int));
        pos = fin.tellg();
        Rcerr << "  Reading matrix: " << str << ", " << fpos << ", " << sizeinbytes << ", " << pos << "\n";
        readMatrix(fin, fpos, sizeinbytes, info, output);
        fin.seekg(pos, std::ios::beg);
    }
}
*/


// [[Rcpp::export]]
DataFrame parseHicCpp(std::string &fname, int resolution) {
    hicInfo info;
    outputStr output;
    std::ifstream fin;
    fin.open(fname, std::fstream::in);
    if (! fin) {
        stop("File " + fname + " cannot be opened for reading.");
    }
    info.resolution = resolution;
    readHeader(fin, info);
    Rcerr << "Available resolutions:\n";
    for (int resolution: info.availableResolutions) {
        Rcerr << "\t" << resolution << "\n";
    }
    Rcerr << info.nChrs << " chromosomes.\n";
    if (info.resolutionIdSelected == -1) {
        Rcerr << "Cannot find resolution " << resolution << ".\nTerminating here.\n";
        stop("Exiting.");
    }
    readFooter(fin, info, output);
    // Transform C++ vectors to R vectors and factors
    Rcerr << "Copying data..." << "\n";
    IntegerVector chrs1, chrs2, bins1, bins2, counts;
    chrs1  = wrap(output.chr1);
    chrs2  = wrap(output.chr2);
    bins1  = wrap(output.bin1);
    bins2  = wrap(output.bin2);
    counts = wrap(output.count);
    if (info.firstChromosomeAll) {
        // the first chr can be 'ALL'; remove it
        info.chrs.erase(0);
    }
    else {
        // factors start with in R
        chrs1 = chrs1 - 1;
        chrs2 = chrs2 - 1;
    }
    chrs1.attr("class") = "factor";
    chrs2.attr("class") = "factor";
    chrs1.attr("levels") = info.chrs;
    chrs2.attr("levels") = info.chrs;
    return DataFrame::create(_["ref1"]  = chrs1,
                             _["bin1"]  = bins1,
                             _["ref2"]  = chrs2,
                             _["bin2"]  = bins2,
                             _["count"] = counts);
}
