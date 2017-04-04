//This is class for DDMRG
#include <iostream>
//#include "dmrg_Hubbard_model.h" just a check
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include "tensor_type.h"
#include "functions.h"
#include "reading_input.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif
#ifndef _PARALLELIZE
#define _PARALLELIZE_AT_SITES_LEVEL false
#define _PARALLELIZE_AT_MATRICES_LEVEL false
#endif

using namespace std;

class DDMRG{


public:
    bool DDMRG_bool;
    int loop_0_DDMRG;
    Mat_1_doub Vec_A;
    Mat_1_doub Vec_X;
    double omega;
    double eta;
    bool Targetting_omega_space;
    int site_A;
    Mat_1_int Finite_loops, Finite_states;

    double weight_GS, weight_A, weight_X;

    void Calculate_X_vector(Mat_2_doub & Unitary_Eig_vecs, Mat_2_doub & Krylov_space_vecs,
                            double & Energy, Mat_1_real & Evals_Lanczos);
    void Initialize_parameters();
    void test1();
    void test2();

};

void DDMRG::Initialize_parameters(){

    loop_0_DDMRG=4;
    DDMRG_bool=false;
    Targetting_omega_space=false;
    Finite_loops.resize(3);
    Finite_loops[0]=10;Finite_loops[1]=-10;Finite_loops[2]=10;

}

void DDMRG::Calculate_X_vector(Mat_2_doub & Unitary_Eig_vecs, Mat_2_doub & Krylov_space_vecs,
                               double & Energy, Mat_1_real & Evals_Lanczos)
{
    type_double zero;
    type_double one;
#ifndef WITH_COMPLEX
    one=1.0;
    zero=0.0;
#endif
#ifdef WITH_COMPLEX
    one=(1.0,0.0);
    zero=(0.0,0.0);
#endif

    type_double value;
    double value_real, value_imag,value_mag2;
    Vec_X.clear();
    Vec_X.resize(Vec_A.size());

for(int i=0;i<Vec_A.size();i++){
    Vec_X[i]=zero;
    for(int m=0;m<Vec_A.size();m++){
        for(int j=0;j<Unitary_Eig_vecs.size();j++){
           for(int k=0;k<Unitary_Eig_vecs.size();k++){
                value_mag2=((Energy+omega - Evals_Lanczos[k])*(Energy+omega - Evals_Lanczos[k]))+
                            (eta*eta);
                value_real=(Energy+omega - Evals_Lanczos[k])/value_mag2;
                value_imag=-eta/value_mag2;
#ifndef WITH_COMPLEX
                value=value_real; //take care later
#endif
#ifdef WITH_COMPLEX
                value=(value_real,value_imag);
#endif

               for(int l=0;l<Unitary_Eig_vecs.size();l++){
                Vec_X[i] = Vec_X[i] + value*Krylov_space_vecs[j][i]*Unitary_Eig_vecs[k][j]*
                                      conjugate(Krylov_space_vecs[k][l])*
                                      conjugate(Unitary_Eig_vecs[l][m])*Vec_A[m];

               }
           }

        }
    }

}

}

void DDMRG::test2(){
    int j=0;
}
