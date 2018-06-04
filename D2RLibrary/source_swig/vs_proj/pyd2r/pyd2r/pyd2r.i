%module pyd2r

%inline %{
#include "pyd2r.h"
%}
int    isKOverlapped(char* w, int u);

double calcWordProb(char* w, double pt, double pc, double pg, double pa);

int overlappingCoef(char* w, double pt, double pc, double pg, double pa);

double D2R_varNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa);

double D2R_varNorm_nNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa);

double D2R_nNorm(const char* seq, int k,
	double pt, double pc, double pg, double pa);

double D2R_nNorm_exactMean(const char* seq, int k,
	double pt, double pc, double pg, double pa);