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

  /*unit vector helper function, which I want to ignore initially*/
  double unit(int judge, int i, int j);

  /*torsion Angles of molecuels*/
  double torsionAng(int i, int j, int k, int l);

  /*to get the cent of mass of the molecule*/
  double centMass(int judge);

  /*Principle Moments of Inertia*/
  double *momntInert();
private:
  int _atomnum;
  int _charge;
  int *_zval;
  double *_x;//use coor[j][i] j={0,1,2}for x, y z would make code more concise 
  double *_y;
  double *_z;
  double _amass[9]= {0.0, 1.00782503223,  4.00260325413,  7.0160034366,  9.012183065, 11.00930536, 12.0000000, 14.00307400443, 15.99491461957};
  double *_momntInert;
  double *_xm;
  double *_ym;
  double *_zm;
  //double **_distc; So we only store necessary data like coordinates and zvalue，
  /*that depends actually. Usually for moment of inertia and cent of mass run for O(n)
  * we prefer to store them in class. Because sesearch them takes O(1)time.
  */
};
#endif
