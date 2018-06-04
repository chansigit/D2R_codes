#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "pyd2r.h"

using namespace std;

// Word Probability under the i.i.d. model
// letter occurrence probability is specified by passing the parameters
double calcWordProb(const char* w, double pt, double pc, double pg, double pa) {
	double likelihood = 1, pr;
	size_t n = strlen(w);
	for (size_t i = 0; i < n; ++i) {
		switch (w[i]){
		case 'T':
		case 't':
			pr = pt;
			break;
		case 'C':
		case 'c':
			pr = pc;
			break;
		case 'G':
		case 'g':
			pr = pg;
			break;
		case 'A':
		case 'a':
			pr = pa;
			break;
		default:
			break;
		}
		likelihood *= pr;
		//debug
		//cout << w[i] << "=" << pr<<" likelihood="<<likelihood<<endl;
	}
	return likelihood;
}

// Overlapping indicator of a word
int isKOverlapped(const char* w, int u){
	size_t n=strlen(w);
	for (size_t i=0; i<u; ++i){	
		if (w[i]!=w[n-u+i])
			return 0;
	}
	return 1;
}

// Overlapping coefficient
double overlappingCoef(const char* w, 
	double pt, double pc, double pg, double pa) {
	double result = 0;
	size_t k = strlen(w);
	for (size_t d = 1; d <= k - 1; ++d) {
		int overlapflag = (int)isKOverlapped(w, k - d);
		if (overlapflag != 0) {
			result += calcWordProb(w+(k-d)*sizeof(char),pt, pc, pg, pa);
		}
	}
	return result;
}


// k-mer counting with hash table
// You must guarantee kmerCnt is empty before calling this function
void countKmers(const char* seq, int k, unordered_map <std::string, int> & kmerCnt) {
	size_t seqLen = strlen(seq);

	// a buffer space holding a k-mer substring, MUST bigger than K
	const int MAXIMUM_KMER_SIZE = 100;
	char buff[MAXIMUM_KMER_SIZE] = { 0 };
	for (size_t pos = 0; pos <= seqLen - k; ++pos) { //enumerate k-mer starting position
													 // populate a k-mer 
		strncpy(buff, seq + pos, k);
		buff[k] = '\0';
		//cout << string(buff) << endl;
		std::string kmer = std::string(buff);
		// increate k-mer count, a zero will be automatically initialized if kmer unmet previously.
		kmerCnt[kmer]++;
	}
}


// print k-mer count
// for debug use
void printKmerCount(const unordered_map <std::string, int>& kmerCnt) {
	for (auto& elem : kmerCnt) {
		std::string kmer = elem.first;
		int cnt = elem.second;
		cout << kmer << " : " << cnt << endl;
	}
}


double D2R_varNorm(const char* seq, int k, 
				   double pt, double pc, double pg, double pa) {
	// Counting K-mers
	unordered_map<std::string, int> kmerCnt;
	countKmers(seq, k, kmerCnt);
	// Adjusting length
	double n = (double)strlen(seq) - k + 1; 
	// Computing statistic
	double result = 0.0;
	for (auto& elem : kmerCnt) {
		std::string kmer   =  elem.first;
		double kmerCnt     = (double)elem.second;
		double kmerProb    = calcWordProb( kmer.c_str(), pt, pc, pg, pa);
		double matchCnt    = kmerCnt * (kmerCnt - 1);
		double EMatchCnt   = n * n   *  kmerProb*kmerProb;
		double StdMatchCnt = n*kmerProb * sqrt(4*n*kmerProb+2);
		result += (matchCnt - EMatchCnt) / StdMatchCnt;
	}
	return result;
}


double D2R_varNorm_nNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa) {
	// Counting K-mers
	unordered_map<std::string, int> kmerCnt;
	countKmers(seq, k, kmerCnt);
	// Adjusting length
	double n = (double)strlen(seq) - k + 1;
	// Computing statistic
	double result = 0.0;
	for (auto& elem : kmerCnt) {
		std::string kmer = elem.first;
		double kmerCnt = (double)elem.second;
		double kmerProb = calcWordProb(kmer.c_str(), pt, pc, pg, pa);
		double matchCnt = kmerCnt * (kmerCnt - 1);
		double EMatchCnt = n * n   *  kmerProb*kmerProb;
		double StdMatchCnt = n * kmerProb * sqrt(4 * n*kmerProb + 2);
		result += (matchCnt - EMatchCnt) / StdMatchCnt;
	}
	result = result / ( n*(n - 1) );
	return result;
}


double D2R_nNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa) {
	// Counting K-mers
	unordered_map<std::string, int> kmerCnt;
	countKmers(seq, k, kmerCnt);
	// Adjusting length
	double n = (double)strlen(seq) - k + 1;
	// Computing statistic
	double result = 0.0;
	for (auto& elem : kmerCnt) {
		std::string kmer = elem.first;
		double kmerCnt = (double)elem.second;
		double kmerProb = calcWordProb(kmer.c_str(), pt, pc, pg, pa);
		double matchCnt = kmerCnt * (kmerCnt - 1);
		double EMatchCnt = n * n   *  kmerProb*kmerProb;
		//double StdMatchCnt = n * kmerProb * sqrt(4 * n*kmerProb + 2);
		result += (matchCnt - EMatchCnt) ;
	}
	result = result / (n*(n - 1));
	return result;
}


double D2R_nNorm_exactMean(const char* seq, int k,
	double pt, double pc, double pg, double pa) {
	// Counting K-mers
	unordered_map<std::string, int> kmerCnt;
	countKmers(seq, k, kmerCnt);
	// Adjusting length
	double n = (double)strlen(seq) - k + 1;
	// Computing statistic
	double result = 0.0;
	for (auto& elem : kmerCnt) {
		std::string kmer = elem.first;
		double kmerCnt = (double)elem.second;
		double kmerProb = calcWordProb(kmer.c_str(), pt, pc, pg, pa);
		double matchCnt = kmerCnt * (kmerCnt - 1);
		double EMatchCnt = (n*n-2*k*n+n)*kmerProb*kmerProb  
			              +2*n*kmerProb* overlappingCoef(kmer.c_str(), pt,pc,pg,pa) ;
		//double StdMatchCnt = n * kmerProb * sqrt(4 * n*kmerProb + 2);
		result += (matchCnt - EMatchCnt);
	}
	result = result / (n*(n - 1));
	return result;
}
