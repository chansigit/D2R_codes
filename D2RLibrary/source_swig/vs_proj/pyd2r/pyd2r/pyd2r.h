#ifndef D2R_H
#define D2R_H
#include <unordered_map>
//#include <algorithm>

using namespace std;
double calcWordProb(const char* w, double pt, double pc, double pg, double pa);


// Binary Overlapping Indicator
// We define the binary variable epsilon(w, u) for u in [1,k]
// to indicate that the word w can overlap itself with u letters
// i.d. its last u letters are the same as its first u letters
int isKOverlapped(const char* w, int u);

// Overlapping coefficient
double overlappingCoef(const char* w,
	double pt, double pc, double pg, double pa);

void countKmers(const char* seq, int k, unordered_map<std::string, int> & kmerCnt);
void printKmerCount(const unordered_map<std::string, int>& kmerCnt);

double D2R_varNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa);

double D2R_varNorm_nNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa);

double D2R_nNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa);

double D2R_nNorm_exactMean(const char* seq, int k,
	double pt, double pc, double pg, double pa);
#endif