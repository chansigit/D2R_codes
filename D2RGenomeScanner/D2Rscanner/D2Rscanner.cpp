#include "stdafx.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <ctime>
#include <map>
using namespace std;

float pr_ch[256] ;
int winsize, k;
string genome_id;
ifstream genome_file;

float calcWordProb(string w) {
	float prob = 1.0;
	for (size_t pos = 0; pos != w.size(); ++pos) {
		char ch = toupper(w[pos]);
		if (ch == 'T' || ch == 'C' || ch == 'G' || ch == 'A') {
			prob *= pr_ch[size_t(ch)];
		}
		else {
			prob *= 0.25;
		}
	}
	return prob;
}


unordered_map <string, int> kmerCnt;
float D2R(string seq, int k) {
	
	//Count k-mer frequencies
	for (size_t pos = 0; pos <= seq.size() - k; ++pos) {
		string kmer=seq.substr(pos, k);
		kmerCnt[kmer]++;
	}

	float res = 0, n=(float)seq.size()-k+1;
	for (auto& elem : kmerCnt) {
		string kmer = elem.first;
		float cnt = (float)elem.second;
		//std::cout << kmer<< "  " << cnt<< '\n';

		float sqrte = n * calcWordProb(kmer);
		res += cnt * (cnt - 1) - sqrte * sqrte;
	}
	return res / (n*(n - 1));
}

float D2R_nosideeffect(string seq, int k) {
	unordered_map <string, int> kmerCnt;
	//Count k-mer frequencies
	for (size_t pos = 0; pos <= seq.size() - k; ++pos) {
		string kmer = seq.substr(pos, k);
		kmerCnt[kmer]++;
	}

	float res = 0, n = (float)seq.size() - k + 1;
	for (auto& elem : kmerCnt) {
		string kmer = elem.first;
		float cnt = (float)elem.second;
		//std::cout << kmer << "  " << cnt << '\n';

		float e = n * calcWordProb(kmer)*n * calcWordProb(kmer);
		res += cnt * (cnt - 1) - e;
	}
	return res / (n*(n - 1));
}
// Usage:
// D2Rscanner E:\code\D2R\CRISPR\genomes\NC_010162.fa  1000  7
int main(int argc, char* argv[]){
	assert(argc == 4);
	if (argc != 4) {
		cout << "Example Usage:\nD2Rscanner E:/code/D2R/CRISPR/genomes/NC_010162.fa  1000  7" << endl;
		return 0;
	}
	clock_t startT,endT;
	
	winsize = atoi(argv[2]);
	k       = atoi(argv[3]);
	assert(k < winsize);
	// ============================= loading genome file =============================
	startT = clock();
	
	genome_file.open(argv[1]);
	getline(genome_file, genome_id);
	std::cout << genome_id << endl;

	string line, genome;
	while (!genome_file.eof()) {
		genome_file >> line;
		genome += line;
	}
	genome_file.close();

	endT = clock();
	std::cout << "Loading Time : " <<fixed<<setprecision(4)
		 << (float)1000*(endT - startT) / CLOCKS_PER_SEC <<"ms"<< endl;
	std::cout << "Genome Length = "<<genome.size() << endl;
	
	// ============================= basic stats =============================
	startT = clock();
	unordered_map<char, unsigned int> frequency;
	
	for (auto& elem : genome)
	{
		if (isalpha(elem)){
			frequency[toupper(elem)]++;
		}
	}
	for (auto& elem : frequency){
		std::cout << elem.first << "  " << elem.second << '\n';
	}
	
	pr_ch[size_t('T')] = pr_ch[size_t('t')] = (float)frequency['T'] / genome.size();
	pr_ch[size_t('C')] = pr_ch[size_t('c')] = (float)frequency['C'] / genome.size();
	pr_ch[size_t('G')] = pr_ch[size_t('g')] = (float)frequency['G'] / genome.size();
	pr_ch[size_t('A')] = pr_ch[size_t('a')] = (float)frequency['A'] / genome.size();

	std::cout << "p(T)=" << pr_ch['T'] << endl;
	std::cout << "p(C)=" << pr_ch['C'] << endl;
	std::cout << "p(G)=" << pr_ch['G'] << endl;
	std::cout << "p(A)=" << pr_ch['A'] << endl;
	endT = clock();
	std::cout << "Letter Frequency Computing Time : " << fixed << setprecision(4)
		<< (float)1000 * (endT - startT) / CLOCKS_PER_SEC << "ms" << endl;

	// ============================= begin computation =============================
	startT = clock();
	
	vector<float> D2RHub;
	
	float currentD2R = D2R(genome.substr(0, winsize), k);
	D2RHub.push_back(currentD2R);
	std::cout << currentD2R << endl;
	
	for (int winpos = 1; winpos <= genome.size() - winsize; ++winpos) {
		string headKmer = genome.substr(winpos - 1, k);
		string tailKmer = genome.substr(winpos + winsize - k, k);
		if (headKmer == tailKmer)
			continue;

		float  headKmerPr = calcWordProb(headKmer);
		float  tailKmerPr = calcWordProb(tailKmer);
		
		float headChange=0, tailChange=0;
		if (kmerCnt[headKmer] != 1) {
			headChange = -2 * (kmerCnt[headKmer] - 1) / float((winsize - k + 1)*(winsize - k));
		}
		else {
			headChange = float(winsize - k + 1) / float(winsize - k) * headKmerPr*headKmerPr;		
		}

		if (kmerCnt[tailKmer] != 0) {
			tailChange = 2 * kmerCnt[tailKmer] / float(winsize - k + 1) / float(winsize - k);
		}
		else {
			tailChange = -float(winsize - k + 1) / float(winsize - k) * tailKmerPr*tailKmerPr;
		}
		
		currentD2R += headChange + tailChange;
		D2RHub.push_back(currentD2R);
		kmerCnt[headKmer] = kmerCnt[headKmer]- 1;
		kmerCnt[tailKmer] = kmerCnt[tailKmer]+ 1;
		
	}

	endT = clock();
	std::cout << "Stat Computing Time : " << fixed << setprecision(4)
		<< (float)1000 * (endT - startT) / CLOCKS_PER_SEC << "ms" << endl;

	// ============================= save results =============================
	startT = clock();
	ofstream txt("result.txt");
	for (int i = 0; i < D2RHub.size(); ++i)
		txt << D2RHub[i] << endl;
	txt.close();
	endT = clock();
	std::cout << "Storing Time : " << fixed << setprecision(4)
		<< (float)1000 * (endT - startT) / CLOCKS_PER_SEC << "ms" << endl;
    return 0;
}
