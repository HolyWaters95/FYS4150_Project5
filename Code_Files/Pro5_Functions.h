#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

using namespace std;
using namespace arma;

//Functions
vec m_vector(double min, double max,double step_length);
void transaction(int i,int j,vec& M);
void Financial_analysis(int Ex, int Cycles, string file1, string file2, double lambda);
void transaction_savings(int i,int j, double lambda, vec& M);
void transaction_taxes(int i, int j, double t, vec& M);


//Runfunctions
void no_savings(int Ex, int Cycles);
void savings(int Ex, int Cycles);



//Tests
void test_stability();
void test_stability_savings();
void test_stability_taxes();
