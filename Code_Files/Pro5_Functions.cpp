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

void transaction(int i,int j, double lambda, vec& M){
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

vector<int> Sampling_Rule(vec M, mat& c, double alpha = 0, double gamma = 0){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    vector<int> index = vector<int>{0, 0};
    int N = M.n_elem;

    int i1 = 0; int i2 = 0;
    int confirm = 0;

    double P = 0;
    double r = 0;

    while (i1 == i2 or confirm == 0){
    confirm = 0;
    i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N) );
    i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N) );
    P = pow(abs(M(i1) - M(i2)), -alpha) * pow(c(i1,i2)+1,gamma);
    r = generate_canonical< double, 128 > (rng);

    if (r < P or abs(M(i1)-M(i2)) < EPS){
        confirm = 1;
    } //end if

    } //end while
    index[0] = i1 ; index[1] = i2;
    c(i1,i2) += 1;
    return index;
} // end function Sampling_Rule

void Financial_analysis(int Ex, int Cycles, int N, string file1, string file2, double lambda = 0, double alpha = 0, double gamma = 0){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int m0 = 1;
    vec M;

    vec m = m_vector(0,N*m0,(N*m0)/500);

    int thread_id;
    int counter;
    string filename1;
    string filename2 = "../Results/" + file2 + ".txt";

    #pragma omp parallel private(thread_id, counter,filename1, M)
    {
    thread_id = omp_get_thread_num();
    counter = 0;
    filename1 = "../Results/" + file1 + "_" + to_string(thread_id) + ".txt";
    #pragma omp for
    // Start Experiment loop
    for (int i = 0;i<Ex;i++){
        M = m0*vec(N,fill::ones);
        mat c = mat(M.n_elem,M.n_elem,fill::zeros);
        // Start MC loop
        for (int j = 0; j<Cycles;j++){


            vector<int> index = Sampling_Rule(M,c,alpha,gamma);

            transaction(index[0],index[1],lambda,M);

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


        counter += 1;



    } //end Experiment loop
    } //end pragma
    return;
} // end of function Financial_analysis

void task_a(int Ex, int Cycles){
    cout << "Task a) \n ------------------" << endl;
    int N = 500;
    int D = static_cast<int>(log10(Cycles));
    string numbers = "_D_" + to_string(D) + "_N_" + to_string(N);

    time_t start, finish;
    start = clock();
    Financial_analysis(Ex,Cycles,N,"Money_distributions"+numbers,"Median"+numbers);
    finish = clock();
    cout << "time used by function Financial_analysis: " << (double) (finish-start)/CLOCKS_PER_SEC << " seconds" << endl;
    return;
} //end of task a

void task_c(int Ex, int Cycles){
    cout << "Task c) \n ------------------" << endl;
    int D = static_cast<int>(log10(Cycles));
    vec lambdas = vec("0.25 0.5 0.9");
    int N = 500;
    for (uword i = 0;i<lambdas.n_elem;i++){
        double L = lambdas(i);
        cout << "Running Financial Analysis for lambda = " << L << endl;

        string numbers = "_D_" + to_string(D) + "_N_" + to_string(N) + "_L_" + to_string(L);
        time_t start, finish;
        start = clock();
        Financial_analysis(Ex,Cycles,N,"Money_distributions" + numbers,"Median" + numbers,L);
        finish = clock();
        cout << "time used by function Financial_analysis: " << (double) (finish-start)/CLOCKS_PER_SEC << " seconds" << endl << endl;
    }
    return;
} //end of task c

void task_d(int Ex, int Cycles){
    cout << "Task d) \n ------------------" << endl;
    int D = static_cast<int>(log10(Cycles));
    vector<int> Nvalues = vector<int>{500,1000};
    vec lambdas = vec("0 0.5");
    vec alphas = vec("0.5 1.0 1.5 2.0");
    for (uword i = 0;i<lambdas.n_elem;i++){
        double L = lambdas(i);
        for (int j = 0;j<Nvalues.size();j++){
            int N = Nvalues[j];
            for (uword k = 0;k<alphas.n_elem;k++){
                double alpha = alphas(k);

                string numbers = "_D_" + to_string(D) + "_N_" + to_string(N) + "_L_" + to_string(L) + "_a_" + to_string(alpha);
                cout << "Running Financial Analysis for" << endl
                     << "N = " << N << endl
                     << "L = " << L << endl
                     << "a = " << alpha << endl;
                time_t start, finish;
                start = clock();
                Financial_analysis(Ex,Cycles,N,"Money_distributions" + numbers,"Median"+numbers,L, alpha);
                finish = clock();
                cout << "time used by function Financial_analysis: " << (double) (finish-start)/CLOCKS_PER_SEC << " seconds" << endl << endl;
            }
        }
    }
    return;
} //end of task d

void task_e(int Ex, int Cycles){
    cout << "Task e) \n ------------------" << endl;
    int D = static_cast<int>(log10(Cycles));
    int N = 1000;
    vec lambdas = vec("0 0.5");
    vec alphas = vec("1.0 2.0");
    vec gammas = vec("1.0 2.0 3.0 4.0");

    for (uword i = 0;i<lambdas.n_elem;i++){
        double L = lambdas(i);
        for (uword j = 0;j<alphas.n_elem;j++){
            double alpha = alphas(j);
            for (uword k = 0;k<gammas.n_elem;k++){
                double gamma = gammas(k);

                string numbers = "_D_" + to_string(D) + "_N_" + to_string(N) + "_L_" + to_string(L) + "_a_" + to_string(alpha) + "_g_" + to_string(gamma);
                cout << "Running Financial Analysis for" << endl
                     << "L = " << L << endl
                     << "a = " << alpha << endl
                     << "g = " << gamma << endl;
                time_t start, finish;
                start = clock();
                Financial_analysis(Ex,Cycles,N,"Money_distributions" + numbers,"Median"+numbers,L, alpha,gamma);
                finish = clock();
                cout << "time used by function Financial_analysis: " << (double) (finish-start)/CLOCKS_PER_SEC << " seconds" << endl << endl;
            }
        }
    }
    return;
} //end of task e
