%module pyd2r

%inline %{
#include "pyd2r.h"
%}
int isKOverlapped(char* w, int u);