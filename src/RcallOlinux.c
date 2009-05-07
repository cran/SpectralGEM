#include<stdlib.h>
#include<stdio.h>

int mainl()
{
  FILE *fp;
  fp=fopen("./Spectral-GEM","r");
  if (!fp) {
    system("ifort ./GEM_sub.f ./SpectralGEM_v2.1.f90 -o ./Spectral-GEM");
  }
  system("./Spectral-GEM < tmp00001.txt");
  system("rm -f tmp00001.txt");
}
