#include <iostream>
#include "tensor_type.h"
#include "functions.h"
//#include "dmrg_Hubbard_model.h" just a check
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include "reading_input.h"
#include "DDMRG_keeper.h"
#include <fstream>
#include <limits> 
#include <iomanip>
#include <stdio.h>
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
#ifndef DMRG_H
#define DMRG_H
class DMRG{

public:

    DDMRG DDMRG_;

    string inp_filename,symmetry_used;
    int Target_L;
    bool Finite_algo_bool, Wavefuntion_transformation, Local_Sz, Spin_Correlations;
    double Energy;

    int Target_state;
    Mat_1_int le,ls                 ;
    int max_iter				    ;
    int m_infinite,m_infinite_max,m_finite_max;
    int nstatestouse;
    Mat_1_int Finite_loops, Finite_states;
    Mat_1_int Finite_states_max;
    Hamiltonian_1_COO H_LB, H_RB;   // Hamiltonian_1_COO[iteration no.]
    Matrix_COO H_SB;                   //  Hamiltonian_1_COO[iteration no.]
    Hamiltonian_2_COO OPSz_LB,OPSp_LB,OPSm_LB,OPSz_RB,OPSp_RB,OPSm_RB; // Hamiltonian_3_COO[iteration no.][site_no]
    Mat_1_doub Eig_vec, Eig_vec_transformed;
    // Unitary_Eig_vecs[vex_ind][component] and Krylov_space_vecs[vec_index][component]
    Mat_2_doub Unitary_Eig_vecs, Krylov_space_vecs;
    Mat_3_doub Red_den_mat_system, Red_den_mat_enviroment;
    Mat_1_real Evals_Lanczos;
    double Truncation_Error_S,Truncation_Error_E, Lanc_Error, Lanc_Error_out, Truncation_Error_Target;
    int max_lanczos_states;
    string ALGO;
    string LOOP_DIRECTION;
    type_double one,zero;
    Mat_1_doub Mag_field;
    string Mag_field_file;
    Mat_2_doub J_zz_long_range,J_pm_long_range,J_mp_long_range,J_pp_long_range,J_mm_long_range;
    string J_zz_long_range_file,J_pm_long_range_file,J_mp_long_range_file,J_pp_long_range_file,J_mm_long_range_file;

    bool NPoint_Corrs;

    // no. of sets of correlations
    int N_corr_sets;

    //N_point[set_no]=n; "n" for n point correlator
    Mat_1_int N_point;

    //at which "sweep no" calculation for set_no^th correlation starts
    Mat_1_int sweep_no;

    //"sweep no" for specific n point correlation ; n_point_sweep_no[set_no]
    Mat_1_int one_point_sweep_no, two_point_sweep_no, four_point_sweep_no, six_point_sweep_no;

    //N_point_oprnames[set_no][i]="A_i" ; "A_i" \belongs to {"Sz","Sp","Sm",..} is operator.
    Mat_2_string N_point_oprnames;

    //N_point_range[set_no]="filename to get sites for set_no^th correlation"
    Mat_1_string N_point_range;

    //n_point_oprnames[set_no][i]="Sz or Sp or Sm .." ; i\belongs to {0,1,2,..n-1}
    Mat_2_string one_point_oprnames, two_point_oprnames, four_point_oprnames, six_point_oprnames;

    //after reading from N_point_range ===> sites_corrs[set_no][operator_no][sites_index]
    Mat_3_int sites_corrs;

    //for specific "n" point correlations, n_point_sites[set_no][operator_no][sites_index]
    Mat_3_int one_point_sites,two_point_sites, four_point_sites, six_point_sites;

    //No. of sets with "n" point correlations
    int one_point_obs,two_point_corrs,four_point_corrs,six_point_corrs;

    //operators for correlations
    Hamiltonian_3_COO One_point_opr_LB, One_point_opr_RB;
    Hamiltonian_3_COO Two_point_opr_LB, Two_point_opr_RB; //Hamiltonian_4_COO[set_no][iteration no.][sites_index]
    Hamiltonian_3_COO Four_point_opr_LB, Four_point_opr_RB;
    Hamiltonian_3_COO Six_point_opr_LB, Six_point_opr_RB;

    //onsite operators for correlations
    Hamiltonian_2_COO One_point_opr_onsite,Two_point_opr_onsite,
    Four_point_opr_onsite,Six_point_opr_onsite;  //Two(4,6)_point_opr_onsite[set_no][opr_no]

    //observables
    Mat_1_doub Sz_vec;
    Mat_2_doub OnePoint_Observable;  //OnePoint_Observable[set_no][site_no]
    Mat_3_doub TwoPointCorrs; //TwoPointCorrs[set_no][site_no][site_no]
    Mat_2_doub FourPointCorrs; //FourPointCorrs[set_no][sites_index]
    Mat_7_doub SixPointCorrs; //SixPointCorrs[set_no][site_no][site_no][][][][]

    //Saving or restart
    bool _SAVING;
    bool _RESTART;
    string saving_filename;
    string restart_filename;
    bool _BINARY;




    void read_INPUT();
    void Create_Q_Eff();
    void Initialize_Hamiltonians();
    void Initialize_Oprts();
    void Grow_LB_and_update_Spin_oprts(int iter);
    int Inverse_LB(int itr,int q);
    void Grow_RB_and_update_Spin_oprts(int iter);

    int Inverse_RB(int itr,int q);
    void Perform_LANCZOS(int sys_iter, int env_iter); //Have to Parallelized later
    void Subtract(Mat_1_doub &temp1, type_double x, Mat_1_doub &temp2);
    void Operate_H_SB(Mat_1_doub &vec_in, int sys_itr, int env_itr,Mat_1_doub &vec_out);
    void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG);
    void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG,
                     Mat_2_doub & Unit_E_vecs, Mat_1_real & Evals_Lanczos);
    type_double Inner_Product(Mat_3_doub  &temp1,Mat_3_doub  &temp2);
    void Do_RENORMALIZATION_of_S_and_E(int sys_iter, int env_iter, int loop, int loop_i, int loop_no);
    void Choosing_m_states(Mat_1_real Eval, int & m_sts, int iter, string sysorenv);
    void Measure_observables(int sys_iter, int env_iter, int loop);
    void Operate_Interface_interactions(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                        type_double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB);
    void Operate_Interface_interactions_f1(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                           type_double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB);
    void Operate_Interface_interactions_f2(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                           type_double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB);
    void Operate_SB_operator(Mat_1_doub &vec_in, Mat_1_doub &vec_out, type_double coeff_value, Matrix_COO &OP_LB, Matrix_COO &OP_RB);

    bool Test_Hermiticity(int sys_iter, int env_iter);
    void Initialize_corr_operators(int loop, int loop_iter);
    void Update_n_point_corr_oprs(int loop,int loop_iter,int env_i,int sys_i);
    void Initialize_onsite_oprs_for_corrs();
    void Writing_data();
    void Reading_data();





};

#endif
