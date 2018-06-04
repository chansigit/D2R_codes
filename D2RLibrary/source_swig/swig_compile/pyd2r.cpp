#include <iostream>
#include <cstring>
#include "pyd2r.h"
using namespace std;

int isKOverlapped(char* w, int u){
	int n=strlen(w);
	for (size_t i=0; i!=n; ++i){
		if (w[i]!=w[n-u+i])
			return false;
	}
	return true;
}