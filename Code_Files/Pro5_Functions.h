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
void Financial_analysis(int Ex, int Cycles, int N, string file1, string file2, double lambda = 0, double alpha = 0, double gamma = 0);
void transaction(int i,int j, double lambda, vec& M);
void transaction_taxes(int i, int j, double t, vec& M);
vector<int> Sampling_Rule(vec M, mat& c, double alpha = 0, double gamma = 0);

//Runfunctions
void task_a(int Ex, int Cycles);
void task_c(int Ex, int Cycles);
void task_d(int Ex, int Cycles);
void task_e(int Ex, int Cycles);



//Tests
void test_sampling();
void test_stability();
void test_stability_savings();
void test_stability_taxes();
void test_P();
