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

#ifndef DDMRG_H
#define DDMRG_H
class DDMRG{


public:
    bool DDMRG_bool;
    int loop_0_DDMRG;
    Mat_1_doub Vec_A;
    Mat_1_Complex_doub Vec_X;
    double Norm_A, Norm_X_real, Norm_X_imag;
    double omega;
    double eta;
    bool Targetting_omega_space;
    int site_A;
    Mat_1_int Finite_loops, Finite_states;
    bool Operator_applied;

    double weight_GS, weight_A, weight_X ;

    void Calculate_X_vector(Mat_2_doub & Unitary_Eig_vecs, Mat_2_doub & Krylov_space_vecs,
                            double & Energy, Mat_1_real & Evals_Lanczos);
    void Calculate_Norms_of_vecX_and_A();
    void Initialize_parameters();
    void test1();
    void test2();

};

#endif
