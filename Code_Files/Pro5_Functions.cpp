#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

using namespace std;
using namespace arma;


vec m_vector(double min, double max,double step_length){
    int num_steps = static_cast<int>((max-min)/step_length);
    vec m = vec(num_steps,fill::zeros);
    m(0) = min;
    for (uword i = 1;i<num_steps;i++){
        m(i) = m(i-1) + step_length;
    }
    return m;
}

void transaction(int i,int j, vec& M){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    double e = generate_canonical< double, 128 > (rng);
    double S = M(i)+M(j);
    M(i) = e*S;
    M(j) = (1-e)*S;
    return;
}
