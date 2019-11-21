#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

#include "Pro5_Functions.h"

using namespace std;
using namespace arma;

int main(){

    int Ex = 0;
    int Cycles = 0;
    cout << "Number of experiment runs? \n";
    cin >> Ex;
    cout << "Number of MC Cycles? \n";
    cin >> Cycles;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 5;
    int m0 = 100;
    vec M = m0*vec(N,fill::ones);
    vec m = m_vector(0,N*m0,0.05);
    vec counter = vec(m.n_elem,fill::zeros);
    // Start Experiment loop
    for (int i = 0;i<Ex;i++){
        // Start MC loop
        for (int j = 0; j<Cycles;j++){

            int i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N) );
            int i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N) );
            transaction(i1,i2,M);
            cout << M << endl;
            cout << sum(M) << endl;
        } //end MC loop
        vec temp = vec(M.n_elem,fill::zeros);
        for (int j = 0; j<N*m0;j++){
            temp.fill(1); temp *= m(j);
            temp = M-temp;
            for (uword i = 0;i<N;i++){
                if (abs(temp(i)) < 0.05){
                    counter(j) += 1;

                }
            }

        }
        //cout << counter << endl;
    } //end Experiment loop



    return 0;
} // end main
