#ifndef MOLEC_GEO_H
#define MOLEC_GEO_H

#include <iostream>
#include <fstream>
//#include <iomanip>

class Molecule{
public:
  /*constructor and destructor
  *since I wont use overloaded= in any case,
  *I may leave that one out from rule of three
  */
  Molecule(const char *filename, int charge);
  ~Molecule();
  int natom;
  /*varies functions*/
  /*printfn*/
  void printfn();

  /*bond distance calculation*/
  double bondLength(int i, int j);

  /*bond angle*/
  double bondAng(int i, int j, int k);

  /*bond out of plane angles*/
  double bondplanA(int i, int j,int k, int l);
  /*unit vector helper function, which I want to leave out initially*/
  double unit(int judge, int i, int j);

private:
  int _atomnum;
  int _charge;
  int *_zval;
  double *_x;
  double *_y;
  double *_z;
  //double **_distc; So we only store necessary data like coordinates and zvalue
};
#endif