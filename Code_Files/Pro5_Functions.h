#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

using namespace std;
using namespace arma;

vec m_vector(double min, double max,double step_length);
void transaction(int i,int j,vec& M);
