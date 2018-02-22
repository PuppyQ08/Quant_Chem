/*@copyright 2018, Xiuyi Qin, all right reserved
*Metropolis Monte Carlo method for x^2*exp(-x^2) integrate
*From CHEM 550 UIUC class
*/
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "boost/random.hpp"
using namespace std;
//sampling function to filter random walk dots
double sampleFn(double x){
  double pi = 3.14159;
  return sqrt(1./pi)*exp(-x*x);
}
//generate random walk and use sampleFn to filter it
bool accept(double x, double xnew,boost::random::mt19937 &gen,  boost::random::uniform_int_distribution<> &dist){
  //double randn = (double)rand()/(RAND_MAX); Never never use rand()!
  double a = double(dist(gen))/ 1000000.0;
  double ratio = sampleFn(xnew)/sampleFn(x);
  if(ratio > a)
  return true;
  return false;
}
double integrand(double var){
  return var*var;
}
//to use filtered coordinate to generate integrate by num of walk times
void MetropolisMC_sampling(double &integOnetime, double &ratio,double numWalk, double stepsize,
  boost::random::mt19937 &gen,  boost::random::uniform_int_distribution<> &dist ){
    double x = 0.0, xnew = 0.0, sum_integrand = 0.0, pi = 3.14159;
    double count = 0.0;
    for(int i = 0; i < numWalk; i++){
      double r = double(dist(gen))/ 1000000.0;
      //this is random walk generation
      xnew = x + (2*r - 1) * stepsize;
      // this is using sampling fn to filter those coordinates
      if(accept(x, xnew, gen, dist)){
        sum_integrand += integrand(xnew);
        count++;
        x = xnew;
      }
      else{
        //this is really important, if rejected, still need to keep x and get integrand into sum as well
        sum_integrand += integrand(x);
      }
    }
     ratio = count/numWalk;//this is ratio of acceptance
     integOnetime = sqrt(pi) * sum_integrand/ numWalk;//remember to be divided by num of walk not count num
  }
  //to get avg result we need to measure more than one time
void MeasureIntegStd(double numMeasure, double &integOnetime, double &ratio,
    double numWalk, double stepsize, boost::random::mt19937 &gen,  boost::random::uniform_int_distribution<> &dist){
    double sum_I2avg = 0.0, sum_Iavg2 = 0.0, sumRatio = 0.0, IntegAvg = 0.0, ratioAvg = 0.0, stdv = 0.0;
    for(int i = 0; i < numMeasure; i++){
      MetropolisMC_sampling(integOnetime, ratio, numWalk, stepsize, gen, dist);
      sum_I2avg += integOnetime*integOnetime;
      sum_Iavg2 += integOnetime;
      sumRatio += ratio;
    }
    sum_I2avg =sum_I2avg/numMeasure;
    sum_Iavg2 = sum_Iavg2/numMeasure;
    IntegAvg = sum_Iavg2;
    sum_Iavg2 = sum_Iavg2*sum_Iavg2;
    stdv = sqrt(sum_I2avg - sum_Iavg2);
    ratioAvg = sumRatio/numMeasure;
    cout<< "num of measurement is " << numMeasure << "num of Walk is "<< numWalk <<endl;
    cout<< " step size is " << stepsize<<" avg ratio is " << ratioAvg <<endl;
    cout<< "avg integrate is "<< IntegAvg << " the standard dev is "<< stdv <<endl;
  }

int main(){
  double stepsize = 2.0;
  double numWalk = 100.0;
  double numMeasure = 100.0;
  double ratio = 0.0, integOnetime = 0.0, IntegAvg = 0.0, ratioAvg = 0.0;
  std::time_t now = std::time(0);
  boost::random::mt19937 gen(now);
  boost::random::uniform_int_distribution<> dist(1, 1000000);
  //increase num of walk by 10 time each time
  for(;numWalk<1000000000; numWalk = numWalk*10.0)
  MeasureIntegStd(numMeasure, integOnetime, ratio, numWalk, stepsize, gen, dist);

}
