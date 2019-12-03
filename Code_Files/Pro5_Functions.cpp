#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>



using namespace std;
using namespace arma;

#define EPS pow(10,-8)

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

void transaction_savings(int i,int j, double lambda, vec& M){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    double e = generate_canonical< double, 128 > (rng);
    double S = (1-lambda)*(M(i)+M(j));
    M(i) = lambda*M(i) + e*S;
    M(j) = lambda*M(j) + (1-e)*S;
    return;
}

void transaction_taxes(int i, int j, double t, vec& M){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    double e = generate_canonical< double, 128 > (rng);
    double S = M(i)+M(j);
    double mi = e*S;
    double mj = (1-e)*S;
    /*
    M(i) = mi ; M(j) = mj;
    double T = 0;
    for (uword i=0;i<M.n_elem;i++){
        T += t*M(i);
        M(i) = (1-t)*M(i);
    }
    M += T/M.n_elem;
    */

    double trans = mi - M(i);
    double Tax = abs(trans)/(1+1/t);
    M(i) = mi; M(j) = mj;

    if (trans > 0){
        M(i) -= Tax;    }
    else{M(j) -= Tax;}
    M += Tax/M.n_elem;

    return;
} // end of transaction_taxes

vector<int> Sampling_Rule(vec M, double alpha){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    vector<int> index = vector<int>{0, 0};
    int N = M.n_elem;

    int i1 = 0; int i2 = 0;
    int confirm = 0;

    double P = 0;
    double r = 0;

    while (i1 == i2 or confirm == 0){
    i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N) );
    i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N) );
    P = pow(abs(M(i1) - M(i2)), -alpha);
    r = generate_canonical< double, 128 > (rng);

    if (r < P or abs(M(i1)-M(i2)) < EPS){
        confirm = 1;
    } //end if

    } //end while
    index[0] = i1 ; index[1] = i2;
    return index;
} // end function Sampling_Rule

void Financial_analysis(int Ex, int Cycles, string file1, string file2, double lambda, double alpha){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 500;
    int m0 = 100;
    vec M;
    vec m = m_vector(0,N*m0,(N*m0)/500);

    int thread_id;
    int counter;
    string filename1;
    string filename2 = file2 + ".txt";

    #pragma omp parallel private(thread_id, counter,filename1, M)
    {
    thread_id = omp_get_thread_num();
    counter = 0;
    filename1 = file1 + "_" + to_string(thread_id) + ".txt";
    #pragma omp for
    // Start Experiment loop
    for (int i = 0;i<Ex;i++){
        M = m0*vec(N,fill::ones);

        // Start MC loop
        for (int j = 0; j<Cycles;j++){


            vector<int> index = Sampling_Rule(M,alpha);

            transaction(index[0],index[1],M);
            //transaction_taxes(i1,i2,lambda,M);
            if(thread_id==0 and counter==0){
                if(j == 0){
                vec M_temp = sort(M);
                double median = (M_temp(249) + M_temp(250))/2;
                ofstream output;
                output.open(filename2,ios::out);
                output << median << endl;
                output.close();
                }
                else if (j % (Cycles/100) == 0) {
                vec M_temp = sort(M);
                double median = (M_temp(249) + M_temp(250))/2;
                ofstream output;
                output.open(filename2,ios::app);
                output << median << endl;
                output.close();
                }

            }

        } //end MC loop

        if(thread_id == 0 and counter == 0){
            ofstream output;
            output.open(filename2,ios::app);
            output << endl;
            output.close();
        }


        //if(i % 100 == 0){
        if(counter == 0){
        ofstream output;
        output.open(filename1,ios::out);
        output << M << endl;
        output.close();
        }
        else{
        ofstream output;
        output.open(filename1,ios::app);
        output << M << endl;
        output.close();
        }
        //}

        counter += 1;



    } //end Experiment loop
    } //end pragma
} // end of function Financial_analysis

void no_savings(int Ex, int Cycles){

    time_t start, finish;
    start = clock();
    Financial_analysis(Ex,Cycles,"Money_distributions_no_savings","Median",0,0);
    finish = clock();
    cout << "time used by function Financial_analysis: " << (double) (finish-start)/CLOCKS_PER_SEC << " seconds" << endl;
    return;
} //end of no_savings

void savings(int Ex, int Cycles){

    vec lambdas = vec("0.25 0.5 0.9");
    for (uword i = 0;i<lambdas.n_elem;i++){
        double L = lambdas(i);
        time_t start, finish;
        start = clock();
        Financial_analysis(Ex,Cycles,"Money_distributions_savings_L_" + to_string(L),"Median_L_" + to_string(L),L,0);
        finish = clock();
        cout << "time used by function Financial_analysis: " << (double) (finish-start)/CLOCKS_PER_SEC << " seconds" << endl;
    }
    return;
} //end of no_savings
