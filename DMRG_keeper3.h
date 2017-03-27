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
#define _PARALLELIZE_AT_MATRICES_LEVEL true
#endif

using namespace std;

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
    void Diagonalize(Mat_1_doub &X ,Mat_1_doub &Y2 , type_double & EG, Mat_1_doub & vecG);
    void Diagonalize(Mat_1_doub &X ,Mat_1_doub &Y2 , type_double & EG, Mat_1_doub & vecG,
                     Mat_2_doub & Unit_E_vecs, Mat_1_real & Evals_Lanczos);
    type_double Inner_Product(Mat_3_doub  &temp1,Mat_3_doub  &temp2);
    void Do_RENORMALIZATION_of_S_and_E(int sys_iter, int env_iter, int loop, int loop_i, int loop_no);
    void Renormalize(Matrix_COO A, type_double* UL,type_double* UR, Matrix_COO & B, int m_UL, int m_UR);
    Matrix_COO Renormalize(type_double* UL,type_double* UR, Matrix_COO A, int m_UL, int m_UR);
    void Choosing_m_states(Mat_1_doub Eval, int & m_sts, int iter, string sysorenv);
    void Measure_observables(int sys_iter, int env_iter, int loop);
    void Operate_Interface_interactions(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                        type_double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB);
    void Operate_Interface_interactions_f1(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                           type_double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB);
    void Operate_Interface_interactions_f2(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                           type_double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB);

    bool Test_Hermiticity(int sys_iter, int env_iter);
    void Initialize_corr_operators(int loop, int loop_iter);
    void Update_n_point_corr_oprs(int loop,int loop_iter,int env_i,int sys_i);
    void Initialize_onsite_oprs_for_corrs();
    void Writing_data();
    void Reading_data();





};


void DMRG::read_INPUT(){

    reading_input(inp_filename,  Mag_field_file, J_zz_long_range_file,J_pm_long_range_file,J_mp_long_range_file,J_pp_long_range_file,J_mm_long_range_file,
                  m_infinite,m_infinite_max, Target_L, Target_state, Truncation_Error_Target, Lanc_Error, Finite_loops,
                  Finite_states,Finite_states_max,Finite_algo_bool, Wavefuntion_transformation,
                  max_lanczos_states,N_corr_sets,
                  N_point, N_point_oprnames, N_point_range, NPoint_Corrs, sweep_no);

    cout<<"|_______________________________________________________________________________________________|"<<endl;
    cout<<"|           DMRG Ground state calculation is being done for following parameters                |"<<endl;
    cout<<"|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|"<<endl;
    cout<<"Target_Length = " <<Target_L<<"\t"<<"  DMRG_Truncation_Error = "<<Truncation_Error_Target<<endl;
    cout<<"Lanczos_Error = "<<Lanc_Error<<"\t"<<endl;
    cout<<"Minf ="<<m_infinite<<"\t"<<"Minf_max ="<<m_infinite_max<<endl;
    cout<<"Finite loops = ";
    for(int i=0;i<Finite_loops.size();i++){cout<<Finite_loops[i]<<"\t";}cout<<endl;
    cout<<"Finite states = ";
    for(int i=0;i<Finite_states.size();i++){cout<<Finite_states[i]<<"\t";}cout<<endl;
    cout<<"Finite states max = ";
    for(int i=0;i<Finite_states_max.size();i++){cout<<Finite_states_max[i]<<"\t";}cout<<endl;
    cout<<"Finite_algo_bool = "<<Finite_algo_bool<<endl;

    cout<<"Wavefuntion_transformation = "<<Wavefuntion_transformation<<endl;



    Truncation_Error_S=0;
    Truncation_Error_E=0;

    cout<<"|_______________________________________________________________________________________________|"<<endl;
    cout<<"|                               DMRG Ground state calculation started                           |"<<endl;
    cout<<"|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|"<<endl;


    sites_corrs.resize(N_corr_sets);
    for(int i =0;i<N_corr_sets;i++){
        sites_corrs[i].resize(N_point[i]);
    }

    reading_sites_for_correlations(N_point_range,sites_corrs);
    reading_long_range_connections(Mag_field_file,J_zz_long_range_file,J_pm_long_range_file,
                                   J_mp_long_range_file,J_pp_long_range_file,J_mm_long_range_file,
                                   Mag_field,J_zz_long_range,J_pm_long_range,
                                   J_mp_long_range,J_pp_long_range,J_mm_long_range,
                                   Target_L);
    reading_restart_or_saving(_RESTART,_SAVING, saving_filename, restart_filename, inp_filename);




}

void DMRG::Create_Q_Eff(){

    max_iter =  (Target_L/2) -2;

    ls.resize(2*max_iter+2);le.resize(2*max_iter+2);

    for(int i=0;i<=2*max_iter+1;i++){

        //-------------------------Creating ls,le----------------------------

        ls[i]= i+1;
        le[i]= Target_L -2 -i;


        //-------------------------Created ls,le-----------------------------

    }


    Red_den_mat_system.resize(2*max_iter +2);
    Red_den_mat_enviroment.resize(2*max_iter +2);


}


void DMRG::Initialize_Hamiltonians(){

#ifndef WITH_COMPLEX
    one=1.0;
    zero=0.0;
#endif
#ifdef WITH_COMPLEX
    one=(1.0,0.0);
    zero=(0.0,0.0);
#endif


    H_LB.resize(2*max_iter+2);H_RB.resize(2*max_iter+2);

    //basis = {down_spin, up_spin}

#ifndef WITH_COMPLEX
    H_LB[0].value.push_back(0.5);H_LB[0].value.push_back(-0.5);
#endif
#ifdef WITH_COMPLEX
    H_LB[0].value.push_back(complex<double>(0.5,0.0));H_LB[0].value.push_back(complex<double>(-0.5,0.0));
#endif
    H_LB[0].rows.push_back(0);H_LB[0].rows.push_back(1);
    H_LB[0].columns.push_back(0);H_LB[0].columns.push_back(1);
    H_LB[0].nrows=2;H_LB[0].ncols=2;
    H_RB[0]=H_LB[0];

    int N_p = omp_get_max_threads();
    cout<<N_p<<" processors are used parallely"<<endl;

}

void DMRG::Initialize_onsite_oprs_for_corrs(){


    one_point_obs=0;
    two_point_corrs=0;
    four_point_corrs=0;
    six_point_corrs=0;


    for(int set_no=0;set_no<N_corr_sets;set_no++){
        if(sites_corrs[set_no].size()==1){
            one_point_obs = one_point_obs+1;
            one_point_sites.push_back(sites_corrs[set_no]);
            one_point_sweep_no.push_back(sweep_no[set_no]);
            one_point_oprnames.push_back(N_point_oprnames[set_no]);
        }
        if(sites_corrs[set_no].size()==2){
            two_point_corrs = two_point_corrs+1;
            two_point_sites.push_back(sites_corrs[set_no]);
            two_point_sweep_no.push_back(sweep_no[set_no]);
            two_point_oprnames.push_back(N_point_oprnames[set_no]);
        }
        if(sites_corrs[set_no].size()==4){
            four_point_corrs = four_point_corrs+1;
            four_point_sites.push_back(sites_corrs[set_no]);
            four_point_sweep_no.push_back(sweep_no[set_no]);
            four_point_oprnames.push_back(N_point_oprnames[set_no]);
        }
        if(sites_corrs[set_no].size()==6){
            six_point_corrs = six_point_corrs+1;
            six_point_sites.push_back(sites_corrs[set_no]);
            six_point_sweep_no.push_back(sweep_no[set_no]);
            six_point_oprnames.push_back(N_point_oprnames[set_no]);
        }
    }


    One_point_opr_onsite.resize(one_point_obs);
    for(int set_no=0;set_no<one_point_obs;set_no++){
        One_point_opr_onsite[set_no].resize(one_point_oprnames[set_no].size());
        for(int opr_no=0;opr_no<one_point_oprnames[set_no].size();opr_no++){

            if(one_point_oprnames[set_no][opr_no]=="Sz"){
                One_point_opr_onsite[set_no][opr_no]=OPSz_LB[0][0];
            }
            if(one_point_oprnames[set_no][opr_no]=="Sp"){
                One_point_opr_onsite[set_no][opr_no]=OPSp_LB[0][0];
            }
            if(one_point_oprnames[set_no][opr_no]=="Sm"){
                One_point_opr_onsite[set_no][opr_no]=OPSm_LB[0][0];
            }
        }
    }





    Two_point_opr_onsite.resize(two_point_corrs);
    for(int set_no=0;set_no<two_point_corrs;set_no++){
        Two_point_opr_onsite[set_no].resize(two_point_oprnames[set_no].size());
        for(int opr_no=0;opr_no<two_point_oprnames[set_no].size();opr_no++){

            if(two_point_oprnames[set_no][opr_no]=="Sz"){
                Two_point_opr_onsite[set_no][opr_no]=OPSz_LB[0][0];
            }
            if(two_point_oprnames[set_no][opr_no]=="Sp"){
                Two_point_opr_onsite[set_no][opr_no]=OPSp_LB[0][0];
            }
            if(two_point_oprnames[set_no][opr_no]=="Sm"){
                Two_point_opr_onsite[set_no][opr_no]=OPSm_LB[0][0];
            }
        }
    }

    Four_point_opr_onsite.resize(four_point_corrs);
    for(int set_no=0;set_no<four_point_corrs;set_no++){
        Four_point_opr_onsite[set_no].resize(four_point_oprnames[set_no].size());
        for(int opr_no=0;opr_no<four_point_oprnames[set_no].size();opr_no++){

            if(four_point_oprnames[set_no][opr_no]=="Sz"){
                Four_point_opr_onsite[set_no][opr_no]=OPSz_LB[0][0];
            }
            if(four_point_oprnames[set_no][opr_no]=="Sp"){
                Four_point_opr_onsite[set_no][opr_no]=OPSp_LB[0][0];
            }
            if(four_point_oprnames[set_no][opr_no]=="Sm"){
                Four_point_opr_onsite[set_no][opr_no]=OPSm_LB[0][0];
            }
        }
    }

    Six_point_opr_onsite.resize(six_point_corrs);
    for(int set_no=0;set_no<six_point_corrs;set_no++){
        Six_point_opr_onsite[set_no].resize(six_point_oprnames[set_no].size());
        for(int opr_no=0;opr_no<six_point_oprnames[set_no].size();opr_no++){

            if(six_point_oprnames[set_no][opr_no]=="Sz"){
                Six_point_opr_onsite[set_no][opr_no]=OPSz_LB[0][0];
            }
            if(six_point_oprnames[set_no][opr_no]=="Sp"){
                Six_point_opr_onsite[set_no][opr_no]=OPSp_LB[0][0];
            }
            if(six_point_oprnames[set_no][opr_no]=="Sm"){
                Six_point_opr_onsite[set_no][opr_no]=OPSm_LB[0][0];
            }
        }
    }


}

void DMRG::Initialize_Oprts(){

    OPSz_LB.resize(2*max_iter+2);OPSp_LB.resize(2*max_iter+2);OPSm_LB.resize(2*max_iter+2);
    OPSz_RB.resize(2*max_iter+2);OPSp_RB.resize(2*max_iter+2);OPSm_RB.resize(2*max_iter+2);

    for(int i=0;i<2*max_iter+2;i++){
        OPSz_LB[i].resize(Target_L);OPSp_LB[i].resize(Target_L);OPSm_LB[i].resize(Target_L);
        OPSz_RB[i].resize(Target_L);OPSp_RB[i].resize(Target_L);OPSm_RB[i].resize(Target_L);

    }



#ifndef WITH_COMPLEX
    OPSz_LB[0][0].value.push_back(-0.5);OPSz_LB[0][0].value.push_back(0.5);
    OPSp_LB[0][0].value.push_back(1.0);
    OPSm_LB[0][0].value.push_back(1.0);
    OPSz_RB[0][Target_L-1].value.push_back(-0.5);OPSz_RB[0][Target_L-1].value.push_back(0.5);
    OPSp_RB[0][Target_L-1].value.push_back(1.0);
    OPSm_RB[0][Target_L-1].value.push_back(1.0);
#endif
#ifdef WITH_COMPLEX
    OPSz_LB[0][0].value.push_back(complex<double>(-0.5,0.0));OPSz_LB[0][0].value.push_back(complex<double>(0.5,0.0));
    OPSp_LB[0][0].value.push_back(complex<double>(1.0,0.0));
    OPSm_LB[0][0].value.push_back(complex<double>(1.0,0.0));
    OPSz_RB[0][Target_L-1].value.push_back(complex<double>(-0.5,0.0));OPSz_RB[0][Target_L-1].value.push_back(complex<double>(0.5,0.0));
    OPSp_RB[0][Target_L-1].value.push_back(complex<double>(1.0,0.0));
    OPSm_RB[0][Target_L-1].value.push_back(complex<double>(1.0,0.0));
#endif



    OPSz_LB[0][0].rows.push_back(0);OPSz_LB[0][0].rows.push_back(1);
    OPSz_LB[0][0].columns.push_back(0);OPSz_LB[0][0].columns.push_back(1);
    OPSz_LB[0][0].nrows=2;OPSz_LB[0][0].ncols=2;


    OPSp_LB[0][0].rows.push_back(0);
    OPSp_LB[0][0].columns.push_back(1);
    OPSp_LB[0][0].nrows=2;OPSp_LB[0][0].ncols=2;


    OPSm_LB[0][0].rows.push_back(1);
    OPSm_LB[0][0].columns.push_back(0);
    OPSm_LB[0][0].nrows=2;OPSm_LB[0][0].ncols=2;





    OPSz_RB[0][Target_L-1].rows.push_back(0);OPSz_RB[0][Target_L-1].rows.push_back(1);
    OPSz_RB[0][Target_L-1].columns.push_back(0);OPSz_RB[0][Target_L-1].columns.push_back(1);
    OPSz_RB[0][Target_L-1].nrows=2;OPSz_RB[0][Target_L-1].ncols=2;


    OPSp_RB[0][Target_L-1].rows.push_back(0);
    OPSp_RB[0][Target_L-1].columns.push_back(1);
    OPSp_RB[0][Target_L-1].nrows=2;OPSp_RB[0][Target_L-1].ncols=2;


    OPSm_RB[0][Target_L-1].rows.push_back(1);
    OPSm_RB[0][Target_L-1].columns.push_back(0);
    OPSm_RB[0][Target_L-1].nrows=2;OPSm_RB[0][Target_L-1].ncols=2;





}


void DMRG::Initialize_corr_operators(int loop,int loop_iter){




    //------------- 1 point operators updated when moving to TO RIGHT-----------------------//
    for(int set_no=0;set_no<one_point_obs;set_no++){
        if(loop==one_point_sweep_no[set_no]-1 && LOOP_DIRECTION=="TO_LEFT"){
            if(loop_iter==abs(Finite_loops[loop])){
                One_point_opr_LB.resize(one_point_obs);
                for(int set_no=0;set_no<one_point_obs;set_no++){
                    One_point_opr_LB[set_no].resize(2*max_iter +2);

                    for(int iter=0;iter<2*max_iter+2;iter++){
                        One_point_opr_LB[set_no][iter].resize(one_point_sites[set_no][0].size());

                        for(int sites_index=0;sites_index<one_point_sites[set_no][0].size();sites_index++){

                            if(one_point_sites[set_no][0][sites_index]==0){

                                One_point_opr_LB[set_no][0][sites_index]=One_point_opr_onsite[set_no][0];
                            }
                        }

                    }

                }
            }
        }
    }

    //------------- 1 point operators TO RIGHT --END-----------------------//




    //------------- 2 point operators updated when moving to TO RIGHT-----------------------//
    for(int set_no=0;set_no<two_point_corrs;set_no++){
        if(loop==two_point_sweep_no[set_no]-1 && LOOP_DIRECTION=="TO_LEFT"){
            if(loop_iter==abs(Finite_loops[loop])){
                Two_point_opr_LB.resize(two_point_corrs);
                for(int set_no=0;set_no<two_point_corrs;set_no++){
                    Two_point_opr_LB[set_no].resize(2*max_iter +2);

                    for(int iter=0;iter<2*max_iter+2;iter++){
                        Two_point_opr_LB[set_no][iter].resize(two_point_sites[set_no][0].size());

                        for(int sites_index=0;sites_index<two_point_sites[set_no][0].size();sites_index++){

                            if(two_point_sites[set_no][0][sites_index]==0){

                                Two_point_opr_LB[set_no][0][sites_index]=Two_point_opr_onsite[set_no][0];
                            }
                        }

                    }

                }
            }
        }
    }

    //------------- 2 point operators TO RIGHT --END-----------------------//

    /*
    //------------- 2 point operators TO LEFT-----------------------//
    for(int set_no=0;set_no<two_point_corrs;set_no++){
        if(loop=two_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_LEFT"){

            Two_point_opr_RB.resize(two_point_corrs);
            for(int set_no=0;set_no<two_point_corrs;set_no++){
                Two_point_opr_RB[set_no].resize(2*max_iter +2);

                for(int iter=0;iter<2*max_iter+2;iter++){
                    Two_point_opr_RB[set_no][iter].resize(two_point_sites[set_no][0].size());

                    for(int sites_index=0;sites_index<two_point_sites[set_no][0].size();sites_index++){


                        if(two_point_sites[set_no][1][sites_index]==Target_L-1){

                            Two_point_opr_RB[set_no][0][sites_index]=Two_point_opr_onsite[set_no][1];
                        }
                    }

                }

            }

        }
    }
    //------------- 2 point operators TO LEFT--END-----------------------//
*/

    //------------- 4 point operators TO RIGHT -----------------------//
    for(int set_no=0;set_no<four_point_corrs;set_no++){
        if(loop==four_point_sweep_no[set_no]-1 && LOOP_DIRECTION=="TO_LEFT"){
            if(loop_iter==abs(Finite_loops[loop])){

                Four_point_opr_LB.resize(four_point_corrs);
                for(int set_no=0;set_no<four_point_corrs;set_no++){
                    Four_point_opr_LB[set_no].resize(2*max_iter +2);

                    for(int iter=0;iter<2*max_iter+2;iter++){
                        Four_point_opr_LB[set_no][iter].resize(four_point_sites[set_no][0].size());

                        for(int sites_index=0;sites_index<four_point_sites[set_no][0].size();sites_index++){

                            if(four_point_sites[set_no][0][sites_index]==0){
                                Four_point_opr_LB[set_no][0][sites_index]=Four_point_opr_onsite[set_no][0];
                            }

                        }

                    }

                }

            }
        }
    }
    //------------- 4 point operators TO RIGHT --END-----------------------//


    //------------- 4 point operators TO LEFT -----------------------//
    for(int set_no=0;set_no<four_point_corrs;set_no++){
        if(loop=four_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_LEFT"){

            Four_point_opr_RB.resize(four_point_corrs);
            for(int set_no=0;set_no<four_point_corrs;set_no++){
                Four_point_opr_RB[set_no].resize(2*max_iter +2);

                for(int iter=0;iter<2*max_iter+2;iter++){
                    Four_point_opr_RB[set_no][iter].resize(four_point_sites[set_no][0].size());

                    for(int sites_index=0;sites_index<four_point_sites[set_no][0].size();sites_index++){

                        if(four_point_sites[set_no][3][sites_index]==0){
                            Four_point_opr_RB[set_no][0][sites_index]=Four_point_opr_onsite[set_no][3];
                        }

                    }

                }

            }

        }
    }
    //------------- 4 point operators TO LEFT --END---------------------//


    //------------- 6 point operators TO RIGHT -----------------------//
    for(int set_no=0;set_no<six_point_corrs;set_no++){
        if(loop=six_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT"){

            Six_point_opr_LB.resize(six_point_corrs);
            for(int set_no=0;set_no<six_point_corrs;set_no++){
                Six_point_opr_LB[set_no].resize(2*max_iter +2);

                for(int iter=0;iter<2*max_iter+2;iter++){
                    Six_point_opr_LB[set_no][iter].resize(six_point_sites[set_no][0].size());

                    for(int sites_index=0;sites_index<six_point_sites[set_no][0].size();sites_index++){

                        if(six_point_sites[set_no][0][sites_index]==0){
                            Six_point_opr_LB[set_no][0][sites_index]=Six_point_opr_onsite[set_no][0];
                        }

                    }

                }

            }


        }
    }
    //------------- 6 point operators TO RIGHT --END-----------------------//


    //------------- 6 point operators TO LEFT-----------------------//
    for(int set_no=0;set_no<six_point_corrs;set_no++){
        if(loop=six_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_LEFT"){

            Six_point_opr_RB.resize(six_point_corrs);
            for(int set_no=0;set_no<six_point_corrs;set_no++){
                Six_point_opr_RB[set_no].resize(2*max_iter +2);

                for(int iter=0;iter<2*max_iter+2;iter++){
                    Six_point_opr_RB[set_no][iter].resize(six_point_sites[set_no][0].size());

                    for(int sites_index=0;sites_index<six_point_sites[set_no][0].size();sites_index++){

                        if(six_point_sites[set_no][5][sites_index]==Target_L-1){
                            Six_point_opr_RB[set_no][0][sites_index]=Six_point_opr_onsite[set_no][5];
                        }

                    }

                }

            }


        }
    }
    //------------- 6 point operators TO LEFT --END-----------------------//




}


void DMRG::Grow_LB_and_update_Spin_oprts(int iter){

    bool Parallelize_Creation_OPSzpm_LB;
    Parallelize_Creation_OPSzpm_LB=_PARALLELIZE_AT_SITES_LEVEL;
    bool add_connection, creating_oprs;
    int tmp;
    Matrix_COO temp_COO;
    tmp =0;
    type_double h_temp;
    type_double one;
#ifndef WITH_COMPLEX
    one=1.0;
#endif
#ifdef WITH_COMPLEX
    one=(1.0,0.0);
#endif

    H_LB[iter+1].ncols=0;
    H_LB[iter+1].nrows=0;
    H_LB[iter+1].columns.clear();H_LB[iter+1].rows.clear();H_LB[iter+1].value.clear();

    //-------------------Just doing LBXI + IXsite---------------------//

    if(H_LB[iter].nrows!=0){
        Direct_Product(Identity(H_LB[iter].nrows),H_LB[0], H_LB[iter+1]);   //I_LB X H_site
        Direct_Product(H_LB[iter],Identity(H_LB[0].nrows),temp_COO);        //H_LB X I_site
        if(iter==0){
            h_temp=Mag_field[0];
        }
        else{
            h_temp=one;
        }
        Sum(H_LB[iter+1],temp_COO,H_LB[iter+1],Mag_field[ls[iter]],h_temp);
    }


    //------------------- LBXI + IXsite done---------------------//


    //---------------------ADDING connections-----------------------//

    for(int siteLB=0;siteLB<ls[iter];siteLB++){


        //S_z(LB)S_z(site)

        if(J_zz_long_range[siteLB][ls[iter]]!=zero)
        {
            Direct_Product(OPSz_LB[iter][siteLB],OPSz_LB[0][0],temp_COO);
            Sum(H_LB[iter+1],temp_COO,H_LB[iter+1],one,J_zz_long_range[siteLB][ls[iter]]);
        }

        //S_p(LB)S_m(site)
        if(J_pm_long_range[siteLB][ls[iter]]!=zero){
            Direct_Product(OPSp_LB[iter][siteLB],OPSm_LB[0][0],temp_COO);
            Sum(H_LB[iter+1],temp_COO,H_LB[iter+1],one,J_pm_long_range[siteLB][ls[iter]]);}

        //S_m(LB)S_p(site)
        if(J_mp_long_range[siteLB][ls[iter]]!=zero){
            Direct_Product(OPSm_LB[iter][siteLB],OPSp_LB[0][0],temp_COO);
            Sum(H_LB[iter+1],temp_COO,H_LB[iter+1],one,J_mp_long_range[siteLB][ls[iter]]);}

        //S_p(LB)S_p(site)
        if(J_pp_long_range[siteLB][ls[iter]]!=zero){
            Direct_Product(OPSp_LB[iter][siteLB],OPSp_LB[0][0],temp_COO);
            Sum(H_LB[iter+1],temp_COO,H_LB[iter+1],one,J_pp_long_range[siteLB][ls[iter]]);}

        //S_m(LB)S_m(site)
        if(J_mm_long_range[siteLB][ls[iter]]!=0){
            Direct_Product(OPSm_LB[iter][siteLB],OPSm_LB[0][0],temp_COO);
            Sum(H_LB[iter+1],temp_COO,H_LB[iter+1],one,J_mm_long_range[siteLB][ls[iter]]);}


    }

    //--------------------Connections DONE------------------------//



    temp_COO.value.clear();temp_COO.rows.clear();temp_COO.columns.clear();


    //__________________________________________________________________________________________________________________________//
    //-----------------------------------------CREATING OPSz,p,m  OPRS for LB ----------------------------------//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    if(!Parallelize_Creation_OPSzpm_LB){goto skiploop_8;}
#pragma omp parallel for default(shared)
skiploop_8:
    for(int site_i=0;site_i<=ls[iter];site_i++){

        if(site_i<ls[iter])
        {
#ifndef WITH_COMPLEX
            creating_oprs=(J_zz_long_range[site_i][ls[iter]]!=0);
#endif
#ifdef WITH_COMPLEX
            creating_oprs=(J_zz_long_range[site_i][ls[iter]].real()!=0 || J_zz_long_range[site_i][ls[iter]].imag()!=0 );
#endif
            if(creating_oprs){
                Direct_Product(OPSz_LB[iter][site_i],Identity(2),OPSz_LB[iter+1][site_i]);
            }
#ifndef WITH_COMPLEX
            creating_oprs=(J_pp_long_range[site_i][ls[iter]]!=0 || J_pm_long_range[site_i][ls[iter]]!=0 || J_mp_long_range[site_i][ls[iter]]!=0);
#endif
#ifdef WITH_COMPLEX
            creating_oprs=(J_pp_long_range[site_i][ls[iter]].real()!=0 || J_pm_long_range[site_i][ls[iter]].real()!=0 || J_mp_long_range[site_i][ls[iter]].real()!=0 ||
                           J_pp_long_range[site_i][ls[iter]].imag()!=0 || J_pm_long_range[site_i][ls[iter]].imag()!=0 || J_mp_long_range[site_i][ls[iter]].imag()!=0);
#endif

            if(creating_oprs){
                Direct_Product(OPSp_LB[iter][site_i],Identity(2),OPSp_LB[iter+1][site_i]);
                Direct_Product(OPSm_LB[iter][site_i],Identity(2),OPSm_LB[iter+1][site_i]);
            }

        }

        else{

            Direct_Product(Identity(H_LB[iter].nrows),OPSz_LB[0][0],OPSz_LB[iter+1][site_i]);
            Direct_Product(Identity(H_LB[iter].nrows),OPSp_LB[0][0],OPSp_LB[iter+1][site_i]);
            Direct_Product(Identity(H_LB[iter].nrows),OPSm_LB[0][0],OPSm_LB[iter+1][site_i]);

        }

    }






    //__________________________________________________________________________________________________________________________//
    //-----------------------------------------CREATED OPSz,p,m OPRS for LB----------------------------------//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


}


void DMRG::Grow_RB_and_update_Spin_oprts(int iter){

    bool Parallelize_Creation_OPSzpm_RB;
    Parallelize_Creation_OPSzpm_RB=_PARALLELIZE_AT_SITES_LEVEL;
    bool add_connection, creating_oprs;
    int tmp;
    Matrix_COO temp_COO;
    tmp =0;
    type_double h_temp;

    H_RB[iter+1].ncols=0;
    H_RB[iter+1].nrows=0;
    H_RB[iter+1].columns.clear();H_RB[iter+1].rows.clear();H_RB[iter+1].value.clear();

    //-------------------Just doing RBXI + IXsite---------------------//

    if(H_RB[iter].nrows!=0){
        Direct_Product(H_RB[0],Identity(H_RB[iter].nrows), H_RB[iter+1]);  // H_site X I(RB)
        Direct_Product(Identity(H_RB[0].nrows),H_RB[iter],temp_COO);       // I(site) X H_RB
        if(iter==0){
            h_temp=Mag_field[Target_L-1];
        }else{
            h_temp=one;
        }
        Sum(H_RB[iter+1],temp_COO,H_RB[iter+1],Mag_field[le[iter]],h_temp);  //  H_mag_site.(site X I(RB)) + (I(site) X H_RB)
    }


    //------------------- RBXI + IXsite done---------------------//


    //---------------------ADDING connections-----------------------//

    for(int siteRB=le[iter]+1;siteRB<Target_L;siteRB++){

        //S_z(RB)S_z(site)
        if(J_zz_long_range[le[iter]][siteRB]!=0){
            Direct_Product(OPSz_LB[0][0],OPSz_RB[iter][siteRB],temp_COO);
            Sum(H_RB[iter+1],temp_COO,H_RB[iter+1],1.0,J_zz_long_range[le[iter]][siteRB]);}

        //S_p(RB)S_m(site)
        if(J_pm_long_range[le[iter]][siteRB]!=0){
            Direct_Product(OPSp_LB[0][0],OPSm_RB[iter][siteRB],temp_COO);
            Sum(H_RB[iter+1],temp_COO,H_RB[iter+1],1.0,J_pm_long_range[le[iter]][siteRB]);}

        //S_m(RB)S_p(site)
        if(J_mp_long_range[le[iter]][siteRB]!=0){
            Direct_Product(OPSm_LB[0][0],OPSp_RB[iter][siteRB],temp_COO);
            Sum(H_RB[iter+1],temp_COO,H_RB[iter+1],1.0,J_mp_long_range[le[iter]][siteRB]);}

        //S_p(RB)S_p(site)
        if(J_pp_long_range[le[iter]][siteRB]!=0){
            Direct_Product(OPSp_LB[0][0],OPSp_RB[iter][siteRB],temp_COO);
            Sum(H_RB[iter+1],temp_COO,H_RB[iter+1],1.0,J_pp_long_range[le[iter]][siteRB]);}

        //S_m(RB)S_m(site)
        if(J_mm_long_range[le[iter]][siteRB]!=0){
            Direct_Product(OPSm_LB[0][0],OPSm_RB[iter][siteRB],temp_COO);
            Sum(H_RB[iter+1],temp_COO,H_RB[iter+1],1.0,J_mm_long_range[le[iter]][siteRB]);}

    }

    //--------------------Connections DONE------------------------//


    temp_COO.value.clear();temp_COO.rows.clear();temp_COO.columns.clear();



    //__________________________________________________________________________________________________________________________//
    //-----------------------------------------CREATING OPSz,p,m  OPRS for RB ----------------------------------//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    if(!Parallelize_Creation_OPSzpm_RB){goto skiploop_9;}
#pragma omp parallel for default(shared)
skiploop_9:
    for(int site_i=le[iter];site_i<Target_L;site_i++){


        if(site_i>le[iter]) {

            if(J_zz_long_range[site_i][le[iter]]!=0){
                Direct_Product(Identity(2),OPSz_RB[iter][site_i],OPSz_RB[iter+1][site_i]);
            }
            if(J_pp_long_range[site_i][le[iter]]!=0
                    || J_pm_long_range[site_i][le[iter]]!=0 ||
                    J_mp_long_range[site_i][le[iter]]!=0){
                Direct_Product(Identity(2),OPSp_RB[iter][site_i],OPSp_RB[iter+1][site_i]);
                Direct_Product(Identity(2),OPSm_RB[iter][site_i],OPSm_RB[iter+1][site_i]);
            }


        }

        else{

            Direct_Product(OPSz_LB[0][0],Identity(H_RB[iter].nrows),OPSz_RB[iter+1][site_i]);
            Direct_Product(OPSp_LB[0][0],Identity(H_RB[iter].nrows),OPSp_RB[iter+1][site_i]);
            Direct_Product(OPSm_LB[0][0],Identity(H_RB[iter].nrows),OPSm_RB[iter+1][site_i]);

        }

    }






    //__________________________________________________________________________________________________________________________//
    //-----------------------------------------CREATED OPSz,p,m OPRS for RB----------------------------------//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


}


void DMRG::Update_n_point_corr_oprs(int loop,int loop_iter,int env_i,int sys_i){


    bool Parallelize_update_n_p_corr;
    Parallelize_update_n_p_corr=_PARALLELIZE_AT_SITES_LEVEL;

    int site_0, site_1, site_2, site_3;



    //One_point
    for(int set_no=0;set_no<one_point_obs;set_no++){
        if( (loop==one_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT")
                || (loop==one_point_sweep_no[set_no]-1 && loop_iter==abs(Finite_loops[loop])) ){

            if(!Parallelize_update_n_p_corr){goto skiploop_6;}
#pragma omp parallel for default(shared) private(site_0,site_1)
skiploop_6:
            for(int sites_index=0;sites_index<One_point_opr_LB[set_no][sys_i].size();sites_index++){

                site_0=one_point_sites[set_no][0][sites_index];

                if(site_0<ls[sys_i]){
                    One_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(One_point_opr_LB[set_no][sys_i][sites_index], Identity(H_LB[0].nrows));

                }

                else if(site_0==ls[sys_i]){
                    One_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Identity(H_LB[sys_i].nrows), One_point_opr_onsite[set_no][0]);
                }


            }
            One_point_opr_LB[set_no][sys_i].clear();
        }

    }

    //Two_points
    for(int set_no=0;set_no<two_point_corrs;set_no++){
        if( (loop==two_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT")
                || (loop==two_point_sweep_no[set_no]-1 && loop_iter==abs(Finite_loops[loop])) ){

            if(!Parallelize_update_n_p_corr){goto skiploop_66;}
#pragma omp parallel for default(shared) private(site_0,site_1)
skiploop_66:
            for(int sites_index=0;sites_index<Two_point_opr_LB[set_no][sys_i].size();sites_index++){

                site_0=two_point_sites[set_no][0][sites_index];
                site_1=two_point_sites[set_no][1][sites_index];

                if(site_0<ls[sys_i] && site_1==ls[sys_i]){
                    Two_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Two_point_opr_LB[set_no][sys_i][sites_index], Two_point_opr_onsite[set_no][1]);

                }
                else if(site_0<ls[sys_i] && (site_1<ls[sys_i] || site_1>ls[sys_i]) ){
                    Two_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Two_point_opr_LB[set_no][sys_i][sites_index], Identity(H_LB[0].nrows));

                }
                else if(site_0==ls[sys_i] && site_1>ls[sys_i]){
                    Two_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Identity(H_LB[sys_i].nrows), Two_point_opr_onsite[set_no][0]);
                }


            }
            Two_point_opr_LB[set_no][sys_i].clear();
        }

    }



    //    if(loop=two_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_LEFT"){
    //        for(int set_no;set_no<two_point_corrs;set_no++){

    //            for(int sites_index=0;sites_index<Two_point_opr_RB[set_no][env_i].size();sites_index++){

    //                site_0=two_point_sites[set_no][0][sites_index];
    //                site_1=two_point_sites[set_no][1][sites_index];
    //                if(site_0<ls[sys_i] && site_1==ls[sys_i]){
    //                    Two_point_opr_RB[set_no][sys_i+1][sites_index]=
    //                            Direct_Product(Two_point_opr_RB[set_no][sys_i][sites_index], Two_point_opr_onsite[set_no][1]);
    //                }
    //                else if(site_0<ls[sys_i] && site_1<ls[sys_i]){
    //                    Two_point_opr_RB[set_no][sys_i+1][sites_index]=
    //                            Direct_Product(Two_point_opr_RB[set_no][sys_i][sites_index], Identity(H_RB[0].nrows));
    //                }
    //                else if(site_0==ls[sys_i] && site_1>ls[sys_i]){
    //                    Two_point_opr_RB[set_no][sys_i+1][sites_index]=
    //                            Direct_Product(Identity(H_RB[sys_i].nrows), Two_point_opr_onsite[set_no][0]);
    //                }


    //            }
    //        }
    //    }


    //Four_points
    for(int set_no=0;set_no<four_point_corrs;set_no++)
    {
        if( (loop==four_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT") ||
                (loop==four_point_sweep_no[set_no]-1 && loop_iter==abs(Finite_loops[loop])) )
        {

            if(!Parallelize_update_n_p_corr){goto skiploop_7;}
#pragma omp parallel for default(shared) private(site_0,site_1,site_2,site_3)
skiploop_7:
            for(int sites_index=0;sites_index<Four_point_opr_LB[set_no][sys_i].size();sites_index++)
            {

                site_0=four_point_sites[set_no][0][sites_index];
                site_1=four_point_sites[set_no][1][sites_index];
                site_2=four_point_sites[set_no][2][sites_index];
                site_3=four_point_sites[set_no][3][sites_index];


                //condition-0 {s0 [] s1 s2 s3} or {s0 s1 [] s2 s3} or  {s0 s1 s2 [] s3} or  {s0 s1 s2 s3 []}
                if( (site_0<ls[sys_i] && site_1>ls[sys_i] && site_2>ls[sys_i] && site_3>ls[sys_i]) ||
                        (site_0<ls[sys_i] && site_1<ls[sys_i] && site_2>ls[sys_i] && site_3>ls[sys_i]) ||
                        (site_0<ls[sys_i] && site_1<ls[sys_i] && site_2<ls[sys_i] && site_3>ls[sys_i]) ||
                        (site_0<ls[sys_i] && site_1<ls[sys_i] && site_2<ls[sys_i] && site_3<ls[sys_i])  )
                {
                    Four_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Four_point_opr_LB[set_no][sys_i][sites_index], Identity(H_LB[0].nrows));


                }

                //condition-1 {[s0] s1 s2 s3}
                if(site_0==ls[sys_i] && site_1>ls[sys_i] && site_2>ls[sys_i] && site_3>ls[sys_i])
                {
                    Four_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Identity(H_LB[sys_i].nrows), Four_point_opr_onsite[set_no][0]);
                }

                //condition-2 {s0 [s1] s2 s3}
                if((site_0<ls[sys_i] && site_1==ls[sys_i] && site_2>ls[sys_i] && site_3>ls[sys_i]))
                {
                    Four_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Four_point_opr_LB[set_no][sys_i][sites_index], Four_point_opr_onsite[set_no][1]);


                }

                //condition-3 {s0 s1 [s2] s3}
                if((site_0<ls[sys_i] && site_1<ls[sys_i] && site_2==ls[sys_i] && site_3>ls[sys_i]))
                {
                    Four_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Four_point_opr_LB[set_no][sys_i][sites_index], Four_point_opr_onsite[set_no][2]);

                }

                //condition-4 {s0 s1 s2 [s3]}
                if((site_0<ls[sys_i] && site_1<ls[sys_i] && site_2<ls[sys_i] && site_3==ls[sys_i]))
                {
                    Four_point_opr_LB[set_no][sys_i+1][sites_index]=
                            Direct_Product(Four_point_opr_LB[set_no][sys_i][sites_index], Four_point_opr_onsite[set_no][3]);

                }




            }
            Four_point_opr_LB[set_no][sys_i].clear();
        }

    }


}


void DMRG::Perform_LANCZOS(int sys_iter, int env_iter){


    int tmp_sz;
    int lanc_iter=0;
    double eps=Lanc_Error;
    double diff_E;
    double temp1, temp3, E0, E0_old;
    Mat_1_doub B2,A, red_eig_vec,Norms;
    Mat_1_doub Kvector_n,Kvector_nm1,Kvector_np1 ; //[element] ; element = i*(Dim(E)) + j(i~Sys, j~Env)


    srand(10);
    tmp_sz=H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;

    if(Wavefuntion_transformation==true){Kvector_n=Eig_vec_transformed;}
    else{
        for(int j=0;j<H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;j++){
            temp1=(rand()%RAND_MAX);
            Kvector_n.push_back(temp1/RAND_MAX);
        }
    }



    double tmpnrm=sqrt(dot_product(Kvector_n,Kvector_n));
    for(int j =0;j<H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;j++){
        Kvector_n[j] = (Kvector_n[j]/(tmpnrm));//*1.0e-10;

    }


    E0_old=0;
    diff_E=1.0;


    while(diff_E>eps && lanc_iter<max_lanczos_states){
        clock_t Lanc_time = clock();
        temp1 =dot_product(Kvector_n,Kvector_n);
        Norms.push_back(sqrt(temp1));
        if(lanc_iter==0){B2.push_back(0);}
        else{

            B2.push_back(tmpnrm*tmpnrm);
        }

        clock_t oprt_SB_time = clock();

        Operate_H_SB(Kvector_n,sys_iter,env_iter,Kvector_np1); // saved in K_vector_np1




        cout<<"Time to operate SB : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;

        temp3 = dot_product(Kvector_n, Kvector_np1);


        A.push_back(dot_product(Kvector_n, Kvector_np1)/dot_product(Kvector_n,Kvector_n));


        Subtract(Kvector_np1, A[lanc_iter], Kvector_n);	//

        if(lanc_iter!=0){Subtract(Kvector_np1, sqrt(B2[lanc_iter]), Kvector_nm1);	}


        //Normalizaton of Knp1 , added by myself, not included in std. Lanczos
        tmpnrm =sqrt(dot_product(Kvector_np1,Kvector_np1)); //new
        for(int i=0;i<Kvector_np1.size();i++){
            Kvector_np1[i]=Kvector_np1[i]/tmpnrm;

        }

        Diagonalize(A,B2,E0,red_eig_vec);


        diff_E=	fabs(E0-E0_old);
        cout<<"diff_E = "<<diff_E<<" E0 = "<<E0<<" E0_old = "<<E0_old<<endl;
        if(lanc_iter==0){if(dot_product(Kvector_np1,Kvector_np1)==0){diff_E = 0;}}
        //cout<<"Energy for lanc_iter("<<lanc_iter<<") is "<<E0<<endl;

        E0_old=E0;

        if(DDMRG_.DDMRG_bool==true){
            Krylov_space_vecs.push_back(Kvector_n);
        }

        //Kvector_nm1.clear();
        Kvector_nm1=Kvector_n;
        //Kvector_n.clear();
        Kvector_n=Kvector_np1;

        if(lanc_iter<=Target_state){diff_E=1.0;}//doing altleast 2 iterations of Lanczos

        lanc_iter=lanc_iter+1;


        if(lanc_iter==tmp_sz){diff_E=0;}
        cout<<"Time for 1 LAnczos iter : "<<double( clock() - Lanc_time ) / (double)CLOCKS_PER_SEC<<endl;
    }

    if(DDMRG_.DDMRG_bool==true){
        Diagonalize(A,B2,E0,red_eig_vec,Unitary_Eig_vecs,Evals_Lanczos);
    }

    cout<<"Perform_LANCZOS: "<<"NO. of itearations required to get convergence in LANCZOS(pass 1) = "<<lanc_iter<<endl;
    cout<<"Perform_LANCZOS: "<<"Energy(GS["<<Target_state<<"] of SB) for length = "<<(sys_iter+env_iter+4)<<" is "<<scientific<<setprecision(20)<< "Energy = "<<E0<<"   "<<"Lanczos_error = "<<diff_E<<"   "<<"E/L = "<<E0/(sys_iter+env_iter+4)<<endl;
    Energy=E0;
    Lanc_Error_out=diff_E;
    cout<<"Perform_LANCZOS: "<<"LANCZOS(pass 2) STARTING FOR SUPERBLOCK Eigenvector, "<<", Size of Matrix(SB) = "<<tmp_sz<<endl<<endl;


    Kvector_n.clear();	Kvector_nm1.clear(); Kvector_np1.clear();
    srand(10);
    tmp_sz=H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;

    if(Wavefuntion_transformation==true){Kvector_n=Eig_vec_transformed;}
    else{
        for(int j=0;j<H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;j++){
            temp1=(rand()%RAND_MAX);
            Kvector_n.push_back(temp1/RAND_MAX);
        }
    }



    tmpnrm =sqrt(dot_product(Kvector_n,Kvector_n));
    for(int j =0;j<H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;j++){
        Kvector_n[j] = (Kvector_n[j]/(tmpnrm));//*1.0e-10;

    }



    Eig_vec.clear();


    for(int j=0;j<H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;j++){
        Eig_vec.push_back(0);
    }


    for(int lanc_iter2=0;lanc_iter2<lanc_iter;lanc_iter2=lanc_iter2+1){



        Subtract(Eig_vec, (-1.0*(red_eig_vec[lanc_iter2])), Kvector_n);
        //cout<<"NOrm = "<<Norms[lanc_iter2]<<endl<<endl;




        Operate_H_SB(Kvector_n,sys_iter,env_iter,Kvector_np1);// saved in K_vector_np1




        Subtract(Kvector_np1, A[lanc_iter2], Kvector_n);	//
        if(lanc_iter2!=0){Subtract(Kvector_np1, sqrt(B2[lanc_iter2]), Kvector_nm1);	}

        //Normalizaton of Knp1 , not included in std. Lanczos
        tmpnrm =sqrt(dot_product(Kvector_np1,Kvector_np1)); //new
        for(int i=0;i<Kvector_np1.size();i++){
            Kvector_np1[i]=Kvector_np1[i]/tmpnrm;
        }


        Kvector_nm1=Kvector_n;
        Kvector_n=Kvector_np1;



    }

    double norm_ev=dot_product(Eig_vec, Eig_vec);
    for(int j=0;j<H_LB[sys_iter+1].nrows*H_RB[env_iter+1].nrows;j++){
        Eig_vec[j]= Eig_vec[j]/sqrt(norm_ev);

    }

    B2.clear();A.clear(); red_eig_vec.clear();Norms.clear();
    Kvector_n.clear();Kvector_nm1.clear();Kvector_np1.clear();

}


void DMRG::Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2){

    if(temp1.size()==temp2.size()){
        for(int k=0;k<temp1.size();k++){
            temp1[k]=temp1[k] - x*temp2[k];

        }
    }
    else{cout<<"Problem in DMRG::Inner_Subtract -- no. of Q_eff sectors not equal in S and E"<<endl;}

}


void DMRG::Diagonalize(Mat_1_doub &X ,Mat_1_doub &Y2 , double & EG, Mat_1_doub & vecG){

    int LDA=X.size(), info;
    /* Local arrays */


    double* eval = new double[X.size()];
    double* mat = new double[X.size()*X.size()];

    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                mat[i*(X.size())+j]=X[j];}
            else if(j==i-1){mat[i*(X.size())+j]=sqrt(Y2[i]);}
            else{mat[i*(X.size())+j]=0;}
        }
    }

    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L', X.size(), mat , LDA, eval );


    /* Check for convergence */
    if( info > 0 ) {

        cout<< "The LAPACKE_dsyev failed to diagonalize."<<endl;

    }


    int T_s;
    if(Target_state>X.size()-1){
        T_s=X.size()-1;
    }
    else{
        T_s=Target_state;
    }
    EG=eval[T_s];
    free(eval);
    vecG.clear();
    vecG.resize(X.size());
    for(int i=0;i<X.size();i++){vecG[i]=mat[i*X.size()+Target_state];}
    free(mat);


}

void DMRG::Diagonalize(Mat_1_doub &X ,Mat_1_doub &Y2 , double & EG, Mat_1_doub & vecG,
                       Mat_2_doub & Unit_E_vecs, Mat_1_real & Evals_Lanczos){

    int LDA=X.size(), info;
    /* Local arrays */

    Evals_Lanczos.clear();
    Evals_Lanczos.resize(X.size());

    double* eval = new double[X.size()];
    double* mat = new double[X.size()*X.size()];

    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                mat[i*(X.size())+j]=X[j];}
            else if(j==i-1){mat[i*(X.size())+j]=sqrt(Y2[i]);}
            else{mat[i*(X.size())+j]=0;}
        }
    }

    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L', X.size(), mat , LDA, eval );


    /* Check for convergence */
    if( info > 0 ) {

        cout<< "The LAPACKE_dsyev failed to diagonalize."<<endl;

    }


    int T_s;
    if(Target_state>X.size()-1){
        T_s=X.size()-1;
    }
    else{
        T_s=Target_state;
    }
    EG=eval[T_s];

    for(int i=0;i<X.size();i++){
    Evals_Lanczos[i]=eval[i];
    }

    free(eval);
    vecG.clear();
    vecG.resize(X.size());
    Unit_E_vecs.clear();
    Unit_E_vecs.resize(X.size());

    for(int i=0;i<X.size();i++){
        Unit_E_vecs[i].resize(X.size());}

    for(int i=0;i<X.size();i++){
        vecG[i]=mat[i*X.size()+Target_state];
        for(int j_state=0;j_state<X.size();j_state++){
            Unit_E_vecs[j_state][i]=mat[i*X.size()+j_state];
        }
    }

    free(mat);


}

void DMRG::Operate_H_SB(Mat_1_doub &vec_in, int sys_itr, int env_itr, Mat_1_doub &vec_out){



    int lp,l,mp,m;

    vec_out.clear();
    vec_out.resize(vec_in.size());
    for(int l=0;l<vec_in.size();l++){
        vec_out[l]=0;
    }

    //----OPEN MP ---------------//
    int th_id;
    int X1, X2;
    Mat_2_doub vec_saved;
    bool f1=false,f2=true;
    int row_A, col_A;
    row_A = omp_get_max_threads();
    col_A = vec_in.size();

    vec_saved.resize(row_A);
    for(int ra=0;ra<row_A;ra++){
        vec_saved[ra].resize(col_A);
#pragma omp parallel for default(shared)
        for(int ca=0;ca<col_A;ca++){
            vec_saved[ra][ca]=0;
        }
    }
    //----OPEN MP ---------------//



    //-------Acting HS X I + I X HE -----------------------------------------------------//

    //----OPEN MP ---------------//
    if(f1==false){goto skiploopU1F1;}
#pragma omp parallel for default(shared) private(th_id,mp,lp,l) if(f1)
skiploopU1F1:
    //----OPEN MP ---------------//

    for(mp=0;mp<H_RB[env_itr+1].nrows;mp++){

        //----OPEN MP ---------------//
        if(f2==false){goto skiploopU1F2;}
#pragma omp parallel for default(shared) private(th_id,lp,l) if(f2)
skiploopU1F2:
        //----OPEN MP ---------------//

        for(int is=0;is<H_LB[sys_itr+1].value.size();is++){

            //----OPEN MP ---------------//
            th_id=omp_get_thread_num();
            //----OPEN MP ---------------//

            lp=H_LB[sys_itr+1].rows[is];
            l=H_LB[sys_itr+1].columns[is];


            vec_saved[th_id][lp*(H_RB[env_itr+1].nrows)+mp] = vec_saved[th_id][lp*(H_RB[env_itr+1].nrows)+mp]  +
                    vec_in[l*(H_RB[env_itr+1].nrows)+mp]*H_LB[sys_itr+1].value[is];

        }
    }



    //----OPEN MP ---------------//
    if(f1==false){goto skiploopU2F1;}
#pragma omp parallel for default(shared) private(th_id,lp,mp,m) if(f1)
skiploopU2F1:
    //----OPEN MP ---------------//

    for(lp=0;lp<H_LB[sys_itr+1].nrows;lp++){

        //----OPEN MP ---------------//
        if(f2==false){goto skiploopU2F2;}
#pragma omp parallel for default(shared) private(th_id,mp,m) if(f2)
skiploopU2F2:
        //----OPEN MP ---------------//

        for(int ie=0;ie<H_RB[env_itr+1].value.size();ie++){

            //----OPEN MP ---------------//
            th_id=omp_get_thread_num();
            //----OPEN MP ---------------//

            mp=H_RB[env_itr+1].rows[ie];
            m=H_RB[env_itr+1].columns[ie];


            vec_saved[th_id][lp*(H_RB[env_itr+1].nrows)+mp] = vec_saved[th_id][lp*(H_RB[env_itr+1].nrows)+mp]  +
                    vec_in[lp*(H_RB[env_itr+1].nrows)+m]*H_RB[env_itr+1].value[ie];

        }
    }




    //-------HS X I + I X HE DONE-----------------------------------------------------//

    //-------------Connections-----------------------------------------------------//


    for(int site_sys=0;site_sys<=ls[sys_itr];site_sys++){
        for(int site_env=le[env_itr];site_env<Target_L;site_env++){

            int site_env_p;
            site_env_p=site_env - (le[env_itr]-1-ls[sys_itr]);

            if(J_zz_long_range[site_sys][site_env_p]!=0){

                X1=OPSz_LB[sys_itr+1][site_sys].value.size();
                X2=OPSz_RB[env_itr+1][site_env].value.size();
                f1=((X1>X2)&&(X1>1)); f2=((X2>=X1)&&(X2>1));
                if(f1==true){
                    Operate_Interface_interactions_f1(vec_in, vec_saved, J_zz_long_range[site_sys][site_env_p],
                                                      OPSz_LB[sys_itr+1][site_sys], OPSz_RB[env_itr+1][site_env]);
                }
                else if(f2==true){
                    Operate_Interface_interactions_f2(vec_in, vec_saved, J_zz_long_range[site_sys][site_env_p],
                                                      OPSz_LB[sys_itr+1][site_sys], OPSz_RB[env_itr+1][site_env]);
                }
                else{
                    Operate_Interface_interactions(vec_in, vec_saved, J_zz_long_range[site_sys][site_env_p],
                                                   OPSz_LB[sys_itr+1][site_sys], OPSz_RB[env_itr+1][site_env]);
                }
            }

            if(J_pm_long_range[site_sys][site_env_p]!=0){

                X1=OPSp_LB[sys_itr+1][site_sys].value.size();
                X2=OPSm_RB[env_itr+1][site_env].value.size();
                f1=((X1>X2)&&(X1>1)); f2=((X2>=X1)&&(X2>1));
                if(f1==true){
                    Operate_Interface_interactions_f1(vec_in, vec_saved, J_pm_long_range[site_sys][site_env_p],
                                                      OPSp_LB[sys_itr+1][site_sys], OPSm_RB[env_itr+1][site_env]);
                }
                else if(f2==true){
                    Operate_Interface_interactions_f2(vec_in, vec_saved, J_pm_long_range[site_sys][site_env_p],
                                                      OPSp_LB[sys_itr+1][site_sys], OPSm_RB[env_itr+1][site_env]);
                }
                else{
                    Operate_Interface_interactions(vec_in, vec_saved, J_pm_long_range[site_sys][site_env_p],
                                                   OPSp_LB[sys_itr+1][site_sys], OPSm_RB[env_itr+1][site_env]);
                }
            }

            if(J_mp_long_range[site_sys][site_env_p]!=0){

                X1=OPSm_LB[sys_itr+1][site_sys].value.size();
                X2=OPSp_RB[env_itr+1][site_env].value.size();
                f1=((X1>X2)&&(X1>1)); f2=((X2>=X1)&&(X2>1));
                if(f1==true){
                    Operate_Interface_interactions_f1(vec_in, vec_saved, J_mp_long_range[site_sys][site_env_p],
                                                      OPSm_LB[sys_itr+1][site_sys], OPSp_RB[env_itr+1][site_env]);
                }
                else if(f2==true){
                    Operate_Interface_interactions_f2(vec_in, vec_saved, J_mp_long_range[site_sys][site_env_p],
                                                      OPSm_LB[sys_itr+1][site_sys], OPSp_RB[env_itr+1][site_env]);
                }
                else{
                    Operate_Interface_interactions(vec_in, vec_saved, J_mp_long_range[site_sys][site_env_p],
                                                   OPSm_LB[sys_itr+1][site_sys], OPSp_RB[env_itr+1][site_env]);
                }
            }

            if(J_pp_long_range[site_sys][site_env_p]!=0){

                X1=OPSp_LB[sys_itr+1][site_sys].value.size();
                X2=OPSp_RB[env_itr+1][site_env].value.size();
                f1=((X1>X2)&&(X1>1)); f2=((X2>=X1)&&(X2>1));
                if(f1==true){
                    Operate_Interface_interactions_f1(vec_in, vec_saved, J_pp_long_range[site_sys][site_env_p],
                                                      OPSp_LB[sys_itr+1][site_sys], OPSp_RB[env_itr+1][site_env]);
                }
                else if(f2==true){
                    Operate_Interface_interactions_f2(vec_in, vec_saved, J_pp_long_range[site_sys][site_env_p],
                                                      OPSp_LB[sys_itr+1][site_sys], OPSp_RB[env_itr+1][site_env]);

                }
                else{
                    Operate_Interface_interactions(vec_in, vec_saved, J_pp_long_range[site_sys][site_env_p],
                                                   OPSp_LB[sys_itr+1][site_sys], OPSp_RB[env_itr+1][site_env]);

                }
            }

            if(J_mm_long_range[site_sys][site_env_p]!=0){

                X1=OPSm_LB[sys_itr+1][site_sys].value.size();
                X2=OPSm_RB[env_itr+1][site_env].value.size();
                f1=((X1>X2)&&(X1>1)); f2=((X2>=X1)&&(X2>1));
                if(f1==true){
                    Operate_Interface_interactions_f1(vec_in, vec_saved, J_mm_long_range[site_sys][site_env_p],
                                                      OPSm_LB[sys_itr+1][site_sys], OPSm_RB[env_itr+1][site_env]);
                }
                else if(f2==true){
                    Operate_Interface_interactions_f2(vec_in, vec_saved, J_mm_long_range[site_sys][site_env_p],
                                                      OPSm_LB[sys_itr+1][site_sys], OPSm_RB[env_itr+1][site_env]);
                }
                else{
                    Operate_Interface_interactions(vec_in, vec_saved, J_mm_long_range[site_sys][site_env_p],
                                                   OPSm_LB[sys_itr+1][site_sys], OPSm_RB[env_itr+1][site_env]);
                }
            }

        }
    }



    //-----------Connections Done-------------------------------------------------//

    for(int ra=0;ra<row_A;ra++){

#pragma omp parallel for default(shared)
        for(int ca=0;ca<col_A;ca++){
            vec_out[ca] +=  vec_saved[ra][ca];
        }
    }

    vec_saved.clear();
}

void DMRG::Operate_Interface_interactions(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                          double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB){

    int lp,l,mp,m;

    for(int is=0;is<OP_LB.value.size();is++){
        for(int ie=0;ie<OP_RB.value.size();ie++){
            lp=OP_LB.rows[is];
            l=OP_LB.columns[is];
            mp=OP_RB.rows[ie];
            m=OP_RB.columns[ie];

            vec_out[0][lp*(OP_RB.nrows)+mp] = vec_out[0][lp*(OP_RB.nrows)+mp] +
                    J_coeff*vec_in[l*(OP_RB.nrows)+m]*OP_LB.value[is]*
                    OP_RB.value[ie];
        }

    }

}

void DMRG::Operate_Interface_interactions_f1(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                             double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB){

    int lp,l,mp,m;
    int th_id;

#pragma omp parallel for default(shared) private(th_id,mp,lp,m,l)
    for(int is=0;is<OP_LB.value.size();is++){
        th_id=omp_get_thread_num();
        for(int ie=0;ie<OP_RB.value.size();ie++){
            lp=OP_LB.rows[is];
            l=OP_LB.columns[is];
            mp=OP_RB.rows[ie];
            m=OP_RB.columns[ie];

            vec_out[th_id][lp*(OP_RB.nrows)+mp] = vec_out[th_id][lp*(OP_RB.nrows)+mp] +
                    J_coeff*vec_in[l*(OP_RB.nrows)+m]*OP_LB.value[is]*
                    OP_RB.value[ie];
        }

    }

}

void DMRG::Operate_Interface_interactions_f2(Mat_1_doub &vec_in, Mat_2_doub &vec_out,
                                             double J_coeff, Matrix_COO &OP_LB, Matrix_COO &OP_RB){

    int lp,l,mp,m;
    int th_id;


    for(int is=0;is<OP_LB.value.size();is++){
#pragma omp parallel for default(shared) private(th_id,mp,lp,m,l)

        for(int ie=0;ie<OP_RB.value.size();ie++){
            th_id=omp_get_thread_num();
            lp=OP_LB.rows[is];
            l=OP_LB.columns[is];
            mp=OP_RB.rows[ie];
            m=OP_RB.columns[ie];

            vec_out[th_id][lp*(OP_RB.nrows)+mp] = vec_out[th_id][lp*(OP_RB.nrows)+mp] +
                    J_coeff*vec_in[l*(OP_RB.nrows)+m]*OP_LB.value[is]*
                    OP_RB.value[ie];
        }

    }

}

double DMRG::Inner_Product(Mat_3_doub  &temp1,Mat_3_doub  &temp2){
    //cout<<"fine 1"<<endl;
    bool bool_to_proceed;
    if(symmetry_used=="Total_Sz"){bool_to_proceed = false;}else {bool_to_proceed = true;}

    double temp, temp_k;
    temp=0;

    if(temp1.size()==temp2.size()){
        for(int k=0;k<temp1.size();k++){

            if(temp1[k].size()==temp2[k].size()){
                temp_k=0;
                for( int l=0;l<temp1[k].size();l++){
                    //  if((bool_to_proceed==true) || (k==l)){
                    if(temp1[k][l].size()==temp2[k][l].size()){

                        for(int i=0;i<temp1[k][l].size();i++){

                            temp_k=temp_k + temp1[k][l][i]*temp2[k][l][i];



                        }

                    }
                    else{cout<<"Problem in DMRG::Inner_Product -- size of vecs in ["<<k<<"]["<<l<<"]th Q_eff not equal"<<endl;}
                    // }
                }

                temp=temp+temp_k;
            }

            else{cout<<"Problem in DMRG::Inner_Product -- size of vecs in "<<k<<"th Q_eff not equal"<<endl;}

        }
    }
    else{cout<<"Problem in DMRG::Inner_Product -- no. of Q_eff sectors not equal in S and E"<<endl;}



    //cout<<"fine 2"<<endl;
    return temp;

}


void DMRG::Do_RENORMALIZATION_of_S_and_E(int sys_iter, int env_iter, int loop, int loop_i, int loop_no){

    bool Parallelize_Renormalization;
    Parallelize_Renormalization=_PARALLELIZE_AT_SITES_LEVEL;
    double* Red_den_mat_sys=(double *) calloc(H_LB[sys_iter+1].nrows*H_LB[sys_iter+1].ncols, sizeof(double));
    double* eval_sys=(double *) calloc(H_LB[sys_iter+1].nrows, sizeof(double));
    double* Red_den_mat_env=(double *) calloc(H_RB[env_iter+1].nrows*H_RB[env_iter+1].ncols, sizeof(double));
    double* eval_env=(double *) calloc(H_RB[env_iter+1].nrows, sizeof(double));
    Mat_1_doub Evl_sys,Evl_env;
    Evl_sys.resize(H_LB[sys_iter+1].nrows);Evl_env.resize(H_RB[env_iter+1].nrows);
    int LDA,info;
    int m_states_env,m_states_sys;




    Red_den_mat_system[sys_iter].clear();

    clock_t rnrm_inside_0 = clock();

    if(DDMRG_.DDMRG_bool==true){

        DDMRG_.Vec_A.clear();
        DDMRG_.Vec_A.assign (Eig_vec.size(),0.0);
        //right now OPR_A is Sz only
        if(DDMRG_.site_A<=ls[sys_iter]){
            Operate_Interface_interactions(Eig_vec, DDMRG_.Vec_A,
                                           1.0, OPSz_LB[sys_iter+1][DDMRG_.site_A],
                    Identity(H_RB[env_iter+1].nrows));
        }
        else if(DDMRG_.site_A>=le[env_iter]){
            Operate_Interface_interactions(Eig_vec, DDMRG_.Vec_A,
                                           1.0, Identity(H_LB[sys_iter+1].nrows),
                    OPSz_RB[env_iter+1][DDMRG_.site_A]);
        }

    DDMRG_.Calculate_X_vector(Unitary_Eig_vecs, Krylov_space_vecs,Energy,Evals_Lanczos);

    }

    LDA = H_LB[sys_iter+1].nrows;
    for(int m=0;m<H_LB[sys_iter+1].nrows;m++){
        for(int n=0;n<=m;n++){    //n<=m must  be done
            for(int l=0;l<H_RB[env_iter+1].nrows;l++){
                Red_den_mat_sys[m*H_LB[sys_iter+1].ncols + n]=
                        Red_den_mat_sys[m*H_LB[sys_iter+1].ncols + n] +
                        (Eig_vec[m*(H_RB[env_iter+1].nrows) + l])*
                        (Eig_vec[n*(H_RB[env_iter+1].nrows) + l]);

            }

        }
    }



    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L', H_LB[sys_iter+1].nrows, Red_den_mat_sys , LDA, eval_sys );


    /* Check for convergence */
    if( info > 0 ) {
        cout<< "The algorithm failed to diagonalize Red_den_mat_sys."<<endl;

    }
    Red_den_mat_system[sys_iter].resize(H_LB[sys_iter+1].nrows);
    for(int m=0;m<H_LB[sys_iter+1].nrows;m++){
        Red_den_mat_system[sys_iter][m].resize(H_LB[sys_iter+1].nrows);
        for(int mtil=0;mtil<H_LB[sys_iter+1].nrows;mtil++){
            Red_den_mat_system[sys_iter][m][mtil] = Red_den_mat_sys[m*H_LB[sys_iter+1].ncols +H_LB[sys_iter+1].ncols  -1 -mtil];
        }
    }








    double tmp_error1=0;

    for(int l=0; l<H_LB[sys_iter+1].nrows;l++){
        Evl_sys[l]=abs(eval_sys[l]);
    }


    Choosing_m_states(Evl_sys,m_states_sys,sys_iter,"SYSTEM");




    Red_den_mat_enviroment[env_iter].clear();

    LDA = H_RB[env_iter+1].nrows;

    for(int m=0;m<H_RB[env_iter+1].nrows;m++){
        for(int n=0;n<=m;n++){
            for(int l=0;l<H_LB[sys_iter+1].nrows;l++){
                Red_den_mat_env[m*H_RB[env_iter+1].ncols + n]=Red_den_mat_env[m*H_RB[env_iter+1].ncols + n] +
                        (Eig_vec[l*(H_RB[env_iter+1].nrows) + m])*
                        (Eig_vec[l*(H_RB[env_iter+1].nrows) + n]);
            }
        }
    }




    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L', H_RB[env_iter+1].nrows, Red_den_mat_env , LDA, eval_env);


    //       Check for convergence
    if( info > 0 ) {
        cout<< "The algorithm failed to diagonalize Red_den_mat_E."<<endl;

    }

    Red_den_mat_enviroment[env_iter].resize(H_RB[env_iter+1].nrows);
    for(int m=0;m<H_RB[env_iter+1].nrows;m++){
        Red_den_mat_enviroment[env_iter][m].resize(H_RB[env_iter+1].nrows);
        for(int mtil=0;mtil<H_RB[env_iter+1].nrows;mtil++){
            Red_den_mat_enviroment[env_iter][m][mtil] = Red_den_mat_env[m*H_RB[env_iter+1].ncols +H_RB[env_iter+1].ncols  -1 -mtil] ;
        }
    }




    for(int l=0; l<H_RB[env_iter+1].nrows;l++){
        Evl_env[l]=abs(eval_env[l]);
    }



    Choosing_m_states(Evl_env,m_states_env,env_iter,"ENVIROMENT");


    Truncation_Error_S=0;


    cout<<"Total_states_sys = "<<m_states_sys<<endl;




    tmp_error1=1.0;


    for(int l=H_LB[sys_iter+1].nrows -1; l >= H_LB[sys_iter+1].nrows - m_states_sys ;l--){

        tmp_error1=tmp_error1 - Evl_sys[l];

    }


    Truncation_Error_S = max(tmp_error1,0.0);


    cout<<"Truncation_Error_S : "<<Truncation_Error_S<<endl;








    cout<<"Total_states_env = "<<m_states_env<<endl;

    Truncation_Error_E=0;


    tmp_error1=1.0;


    for(int l=H_RB[env_iter+1].nrows -1; l >= H_RB[env_iter+1].nrows - m_states_env ;l--){

        tmp_error1=tmp_error1 - Evl_env[l];

    }






    Truncation_Error_E = max(tmp_error1,0.0);


    cout<<"Truncation_Error_E : "<<Truncation_Error_E<<endl;

    cout<<"Time for doing Renormalization 0 : "<<double( clock() - rnrm_inside_0 ) / (double)CLOCKS_PER_SEC<<endl;

    clock_t wft_time = clock();




    //------------Wavefunction Transformation------------------------//
    //```````````````````````````````````````````````````````````````//

    if(Wavefuntion_transformation==true){

        //System size decreases, Enviroment increases
        if((loop==0 && loop_i==0) || ((loop<0 && loop_i<abs(loop))) ||(loop>0 && loop_i==loop)  ){
            cout<<"DOING WAVEFUNCTION TRANSFORMATION"<<endl;

            int dim_env_p,dim_sys_p;
            dim_env_p=H_LB[0].nrows*m_states_env;
            dim_sys_p=H_LB[0].nrows*H_LB[sys_iter-1].ncols;


            Eig_vec_transformed.clear();
            Eig_vec_transformed.assign(dim_sys_p*dim_env_p, 0);

            int lp,m_til,m;
            for(int l_til=0;l_til<m_states_env;l_til++){

                for(int mp=0;mp<dim_sys_p;mp++){

                    for(int i =0;i<H_LB[0].nrows;i++){

                        lp = i*m_states_env + l_til;
                        for(int l=0;l<H_RB[env_iter+1].nrows;l++){
                            for(m_til=0;m_til<H_LB[sys_iter].nrows;m_til++){
                                m =  m_til*H_LB[0].nrows + i;
                                Eig_vec_transformed[mp*dim_env_p + lp] =  Eig_vec_transformed[mp*dim_env_p + lp] +
                                        Red_den_mat_enviroment[env_iter][l][l_til]*Eig_vec[m*H_RB[env_iter+1].nrows + l]*
                                        Red_den_mat_system[sys_iter-1][mp][m_til];


                            }
                        }
                    }
                }
            }
        }

        else{

            cout<<"DOING WAVEFUNCTION TRANSFORMATION"<<endl;

            int dim_env_p,dim_sys_p;
            dim_env_p=H_RB[env_iter-1].ncols*H_LB[0].nrows;
            dim_sys_p=m_states_sys*H_LB[0].nrows;

            Eig_vec_transformed.clear();
            Eig_vec_transformed.assign(dim_sys_p*dim_env_p, 0);

            for(int m_til=0;m_til<m_states_sys;m_til++){

                for(int lp=0;lp<dim_env_p;lp++){

                    for(int i =0;i<H_LB[0].nrows;i++){

                        int mp = m_til*H_LB[0].nrows + i;     //dim_before_ie,dim_before_isp have not calulated yet(1)

                        for(int m=0;m<H_LB[sys_iter+1].nrows;m++){
                            for(int l_til=0;l_til<H_RB[env_iter].nrows;l_til++){
                                int l = i*H_RB[env_iter].nrows + l_til;
                                Eig_vec_transformed[mp*dim_env_p + lp] =  Eig_vec_transformed[mp*dim_env_p + lp] +
                                        Red_den_mat_enviroment[env_iter-1][lp][l_til]*Eig_vec[m*H_RB[env_iter+1].nrows + l]*
                                        Red_den_mat_system[sys_iter][m][m_til];
                            }
                        }
                    }
                }
            }
        }

    }
    //------------Wavefunction Transformation Done-------------------//
    //```````````````````````````````````````````````````````````````//



    cout<<"Time for WFT : "<<double( clock() - wft_time ) / (double)CLOCKS_PER_SEC<<endl;




    clock_t rnrm_inside_1 = clock();

    bool FOR_ALL=false;
    Matrix_COO temp_coo;
    if(!Parallelize_Renormalization){goto skiploop_1;}
#pragma omp parallel for default(shared)
skiploop_1:
    for(int site_i=le[env_iter];site_i<Target_L;site_i++){

        if(FOR_ALL || site_i==le[env_iter] || J_zz_long_range[site_i][le[env_iter]]!=0){
            Renormalize(OPSz_RB[env_iter+1][site_i],Red_den_mat_env, Red_den_mat_env, OPSz_RB[env_iter+1][site_i],m_states_env, m_states_env);
        }

        if(FOR_ALL || site_i==le[env_iter] || J_pp_long_range[site_i][le[env_iter]]!=0
                || J_pm_long_range[site_i][le[env_iter]]!=0 ||
                J_mp_long_range[site_i][le[env_iter]]!=0){
            Renormalize(OPSp_RB[env_iter+1][site_i],Red_den_mat_env, Red_den_mat_env, OPSp_RB[env_iter+1][site_i],m_states_env, m_states_env);

            Renormalize(OPSm_RB[env_iter+1][site_i],Red_den_mat_env, Red_den_mat_env, OPSm_RB[env_iter+1][site_i],m_states_env, m_states_env);
        }
    }

    //temp_coo=
    Renormalize(H_RB[env_iter+1],Red_den_mat_env, Red_den_mat_env, H_RB[env_iter+1],m_states_env, m_states_env);
    //H_RB[env_iter+1] = temp_coo;

    if(!Parallelize_Renormalization){goto skiploop_2;}
#pragma omp parallel for default(shared)
skiploop_2:
    for(int site_i=0;site_i<=ls[sys_iter];site_i++){

        if(FOR_ALL || site_i==ls[sys_iter] || J_zz_long_range[site_i][ls[sys_iter]]!=0){
            Renormalize(OPSz_LB[sys_iter+1][site_i],Red_den_mat_sys, Red_den_mat_sys, OPSz_LB[sys_iter+1][site_i],m_states_sys, m_states_sys);
        }

        if(FOR_ALL || site_i==ls[sys_iter] || J_pp_long_range[site_i][ls[sys_iter]]!=0
                || J_pm_long_range[site_i][ls[sys_iter]]!=0 ||
                J_mp_long_range[site_i][ls[sys_iter]]!=0){
            Renormalize(OPSp_LB[sys_iter+1][site_i],Red_den_mat_sys, Red_den_mat_sys, OPSp_LB[sys_iter+1][site_i],m_states_sys, m_states_sys);

            Renormalize(OPSm_LB[sys_iter+1][site_i],Red_den_mat_sys, Red_den_mat_sys, OPSm_LB[sys_iter+1][site_i],m_states_sys, m_states_sys);
        }
    }


    //temp_coo=
    Renormalize(H_LB[sys_iter+1],Red_den_mat_sys, Red_den_mat_sys, H_LB[sys_iter+1],m_states_sys, m_states_sys);
    //H_LB[sys_iter+1]=temp_coo;


    //--------------RENORMALIZING CORRELATION OPERATORS---------------------//
    cout<<"Time for doing Renormalization inside 1: "<<double( clock() - rnrm_inside_1 ) / (double)CLOCKS_PER_SEC<<endl;

    clock_t rnrm_inside_2 = clock();
    int site_0;
    if(ALGO=="FINITE"){

        //One_point
        for(int set_no=0;set_no<one_point_obs;set_no++){
            if( (loop_no==one_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT")
                    || (loop_no==one_point_sweep_no[set_no] -1 && loop_i==abs(Finite_loops[loop_no])) ){
                if(!Parallelize_Renormalization){goto skiploop_3;}
#pragma omp parallel for default(shared) private(site_0)
skiploop_3:
                for(int sites_index=0;sites_index<one_point_sites[set_no][0].size();sites_index++){

                    site_0=one_point_sites[set_no][0][sites_index];

                    if(site_0<=ls[sys_iter]){
                        //temp_coo =
                        Renormalize( One_point_opr_LB[set_no][sys_iter+1][sites_index],Red_den_mat_sys, Red_den_mat_sys, One_point_opr_LB[set_no][sys_iter+1][sites_index],m_states_sys, m_states_sys);
                        //Two_point_opr_LB[set_no][sys_iter+1][sites_index]= temp_coo;
                    }
                }
            }
        }

        //Two_point
        for(int set_no=0;set_no<two_point_corrs;set_no++){
            if( (loop_no==two_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT")
                    || (loop_no==two_point_sweep_no[set_no] -1 && loop_i==abs(Finite_loops[loop_no])) ){
                if(!Parallelize_Renormalization){goto skiploop_36;}
#pragma omp parallel for default(shared) private(site_0)
skiploop_36:
                for(int sites_index=0;sites_index<two_point_sites[set_no][0].size();sites_index++){

                    site_0=two_point_sites[set_no][0][sites_index];

                    if(site_0<=ls[sys_iter]){
                        //temp_coo =
                        Renormalize( Two_point_opr_LB[set_no][sys_iter+1][sites_index],Red_den_mat_sys, Red_den_mat_sys, Two_point_opr_LB[set_no][sys_iter+1][sites_index],m_states_sys, m_states_sys);
                        //Two_point_opr_LB[set_no][sys_iter+1][sites_index]= temp_coo;
                    }
                }
            }
        }

        //Four_point
        for(int set_no=0;set_no<four_point_corrs;set_no++){
            if( (loop_no==four_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT")
                    || (loop_no==four_point_sweep_no[set_no] -1 && loop_i==abs(Finite_loops[loop_no])) ){

                if(!Parallelize_Renormalization){goto skiploop_4;}
#pragma omp parallel for default(shared) private(site_0)
skiploop_4:
                for(int sites_index=0;sites_index<four_point_sites[set_no][0].size();sites_index++){

                    site_0=four_point_sites[set_no][0][sites_index];

                    if(site_0<=ls[sys_iter]){
                        //temp_coo =
                        Renormalize(Four_point_opr_LB[set_no][sys_iter+1][sites_index],Red_den_mat_sys, Red_den_mat_sys, Four_point_opr_LB[set_no][sys_iter+1][sites_index],m_states_sys, m_states_sys);
                        //Four_point_opr_LB[set_no][sys_iter+1][sites_index]= temp_coo;
                    }
                }
            }
        }


    }

    //--------------RENORMALIZING CORRELATION OPERATORS DONE---------------------//

    cout<<"Time for doing Renormalization inside 2 : "<<double( clock() - rnrm_inside_2 ) / (double)CLOCKS_PER_SEC<<endl;
    temp_coo.value.clear();temp_coo.rows.clear();temp_coo.columns.clear();

    free(Red_den_mat_sys);free(eval_sys);free(Red_den_mat_env);free(eval_env);
    Evl_env.clear();Evl_sys.clear();

}

Matrix_COO DMRG::Renormalize(double* UL,double* UR, Matrix_COO A, int m_UL, int m_UR){

    bool Parallelize_Renormalize_function;
    Parallelize_Renormalize_function=_PARALLELIZE_AT_MATRICES_LEVEL;
    int m_inf_r = min(m_UL, A.nrows);
    int m_inf_c = min(m_UR, A.ncols);
    int tmp=0;
    Matrix_COO B;
    //    B.value.clear();
    //    B.rows.clear();
    //    B.columns.clear();
    B.nrows=m_inf_r;
    B.ncols=m_inf_c;



    double temp;


    for(int i=0;i<m_inf_r;i++){

        for(int l=0;l<m_inf_c;l++){

            temp=0;
            if(!Parallelize_Renormalize_function){goto skiploop_14;}
#pragma omp parallel for default(shared) reduction(+:temp)
skiploop_14:
            for( int p=0;p<A.value.size();p++){

                temp=temp+UL[A.rows[p]*A.nrows + A.nrows - 1 - i]*
                        A.value[p]*UR[A.columns[p]*A.ncols + A.ncols - 1 - l];

            }

            if(fabs(temp)>Lanc_Error){
                B.value.push_back(temp);
                B.rows.push_back(i);
                B.columns.push_back(l);
                tmp=tmp+1;
            }
        }
    }

    if(tmp==0){B.value.push_back(0);B.rows.push_back(0);B.columns.push_back(0);}
    //if(m_inf_r==0 || m_inf_c==0){cout<<"some problem"<<endl;}

    if(m_inf_r==0 || m_inf_c==0){
        B.value.clear();
        B.rows.clear();
        B.columns.clear();
        B.nrows=0;
        B.ncols=0;

    }
    return B;
}

void DMRG::Renormalize(Matrix_COO A, double* UL,double* UR, Matrix_COO & B, int m_UL, int m_UR){

    bool Parallelize_Renormalize_function;
    Parallelize_Renormalize_function=_PARALLELIZE_AT_MATRICES_LEVEL;
    int m_inf_r = min(m_UL, A.nrows);
    int m_inf_c = min(m_UR, A.ncols);
    int tmp=0;
    B.value.clear();
    B.rows.clear();
    B.columns.clear();
    B.nrows=m_inf_r;
    B.ncols=m_inf_c;
    double temp;
    for(int i=0;i<m_inf_r;i++){

        for(int l=0;l<m_inf_c;l++){

            temp=0;
            if(!Parallelize_Renormalize_function){goto skiploop_14;}
#pragma omp parallel for default(shared) reduction(+:temp)
skiploop_14:
            for( int p=0;p<A.value.size();p++){

                temp=temp+UL[A.rows[p]*A.nrows + A.nrows - 1 - i]*
                        A.value[p]*UR[A.columns[p]*A.ncols + A.ncols - 1 - l];

            }

            if(fabs(temp)>Lanc_Error){
                B.value.push_back(temp);
                B.rows.push_back(i);
                B.columns.push_back(l);
                tmp=tmp+1;
            }
        }
    }

    if(tmp==0){B.value.push_back(0);B.rows.push_back(0);B.columns.push_back(0);}
    //if(m_inf_r==0 || m_inf_c==0){cout<<"some problem"<<endl;}

    if(m_inf_r==0 || m_inf_c==0){
        B.value.clear();
        B.rows.clear();
        B.columns.clear();
        B.nrows=0;
        B.ncols=0;

    }
}

void DMRG::Choosing_m_states(Mat_1_doub  Eval, int & m_sts, int iter, string sysorenv){

    bool check_while;
    double eps=0;
    double error=1.0;
    if(sysorenv == "SYSTEM" || sysorenv == "ENVIROMENT"){


        m_sts=0;

        int m=0;
        double max;
        int m_limit=Eval.size();


        m_limit=min(m_limit,m_infinite_max);


        check_while=true;
        while(check_while==true){
            max=-1;

            for(int i=0;i<Eval.size();i++){
                if(Eval[i]>max){max=Eval[i];}
            }
            for(int i=0;i<Eval.size();i++){
                if(Eval[i]==max){
                    if(m_sts<Eval.size()){
                        m_sts=m_sts+1;m=m+1;
                        error = error - Eval[i];
                    }
                    Eval[i]=0;
                }
            }


            if(m_limit>=nstatestouse)
            {   if((m<nstatestouse)||((m<m_limit)&&(error>=Truncation_Error_Target))){
                    check_while=true;}else{check_while=false;}
            }
            else{
                if(m<m_limit){
                    check_while=true;}else{check_while=false;}
            }

        }

    }




}




void DMRG::Measure_observables(int sys_iter, int env_iter, int loop){

    //This have to be parallelized at _PARALLELIZE_AT_SITES_LEVEL
    // and take care of dot_product properly by making 2 dot product, one is always parallelized, which goes in Lanczos
    // and other depends on _PARALLELIZE_AT_MATRICES_LEVEL

    OnePoint_Observable.resize(one_point_obs);
    for(int set_no=0;set_no<one_point_obs;set_no++){
        OnePoint_Observable[set_no].resize(Target_L);
    }


    TwoPointCorrs.resize(two_point_corrs);
    for(int i=0;i<two_point_corrs;i++){
        TwoPointCorrs[i].resize(Target_L);
        for(int j=0;j<Target_L;j++){
            TwoPointCorrs[i][j].resize(Target_L);
        }
    }
    FourPointCorrs.resize(four_point_corrs);
    for(int set_no=0;set_no<four_point_corrs;set_no++){
        FourPointCorrs[set_no].resize(four_point_sites[set_no][0].size());
    }
    Mat_2_doub vec_out;
    vec_out.resize(1);
    vec_out[0].resize(Eig_vec.size(),0);
    Sz_vec.resize(Target_L);
    Matrix_COO temp_coo_sys, temp_coo_env;

    int site_0,site_1,site_2,site_3;




    //-------------------------------Measuring 4 point correlations-------------------------------//

    for(int set_no=0;set_no<four_point_corrs;set_no++){
        if(loop==four_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT"){

            for(int sites_index=0;sites_index<four_point_sites[set_no][0].size();sites_index++){
                site_0=four_point_sites[set_no][0][sites_index];
                site_1=four_point_sites[set_no][1][sites_index];
                site_2=four_point_sites[set_no][2][sites_index];
                site_3=four_point_sites[set_no][3][sites_index];


                //        SYSTEM=LB+site(L-3)   site(L-2)  RB(L-1)
                //      [------s0---s1--s2---]    []        [s3]
                if(site_3==Target_L-1 && site_2<Target_L-2)
                {
                    temp_coo_env=Direct_Product(Identity(H_LB[0].nrows), Four_point_opr_onsite[set_no][3]);
                    temp_coo_sys=Four_point_opr_LB[set_no][sys_iter+1][sites_index];
                }

                //        SYSTEM=LB+site(L-3)   site(L-2)  RB(L-1)
                //       [------s0---s1-------]    [s2]      [s3]

                if(site_3==Target_L-1 && site_2==Target_L-2)
                {
                    temp_coo_env=Direct_Product(Four_point_opr_onsite[set_no][2], Four_point_opr_onsite[set_no][3]);
                    temp_coo_sys=Four_point_opr_LB[set_no][sys_iter+1][sites_index];
                }

                //        SYSTEM=LB+site(L-3)   site(L-2)  RB(L-1)
                //       [------s0---s1---s2--]    [s3]      []

                if(site_3==Target_L-2)
                {
                    temp_coo_env=Direct_Product(Four_point_opr_onsite[set_no][3], Identity(H_LB[0].nrows));
                    temp_coo_sys=Four_point_opr_LB[set_no][sys_iter+1][sites_index];
                }

                //        SYSTEM=LB+site(L-3)   site(L-2)  RB(L-1)
                //       [--s0---s1---s2--s3--]    []         []

                if(site_3<Target_L-2)
                {
                    temp_coo_env=Identity(H_RB[env_iter+1].nrows);
                    temp_coo_sys=Four_point_opr_LB[set_no][sys_iter+1][sites_index];
                }

                vec_out[0].assign(Eig_vec.size(),0);
                Operate_Interface_interactions(Eig_vec, vec_out, 1.0, temp_coo_sys, temp_coo_env);
                FourPointCorrs[set_no][sites_index]=dot_product(Eig_vec,vec_out[0]);

            }
        }
    }

    //--------------------------Measured 4 point correlations-----------------------------------------------------//



    //-------------------------------Measuring 2 point correlations-------------------------------//

    for(int set_no=0;set_no<two_point_corrs;set_no++){
        if(loop==two_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT"){

            for(int sites_index=0;sites_index<two_point_sites[set_no][0].size();sites_index++){
                site_0=two_point_sites[set_no][0][sites_index];
                site_1=two_point_sites[set_no][1][sites_index];

                if(site_0<=ls[sys_iter] && site_1<=ls[sys_iter]){
                    temp_coo_sys=Two_point_opr_LB[set_no][sys_iter+1][sites_index];
                    temp_coo_env=Identity(H_RB[env_iter+1].nrows);
                }

                if(site_0<=ls[sys_iter] && site_1>=le[env_iter]){
                    temp_coo_sys=Two_point_opr_LB[set_no][sys_iter+1][sites_index];
                    if(site_1==le[env_iter]){
                        temp_coo_env=Direct_Product(Two_point_opr_onsite[set_no][1],Identity(H_LB[0].nrows));}
                    if(site_1==le[env_iter]+1){
                        temp_coo_env=Direct_Product(Identity(H_LB[0].nrows), Two_point_opr_onsite[set_no][1]);}

                }

                if(site_0==le[env_iter] && site_1==le[env_iter]+1){
                    temp_coo_sys=Identity(H_LB[sys_iter+1].nrows);
                    temp_coo_env=Direct_Product(Two_point_opr_onsite[set_no][0],Two_point_opr_onsite[set_no][1]);
                }


                vec_out[0].assign(Eig_vec.size(),0);
                Operate_Interface_interactions(Eig_vec, vec_out, 1.0, temp_coo_sys, temp_coo_env);
                TwoPointCorrs[set_no][site_0][site_1]=dot_product(Eig_vec,vec_out[0]);

            }
        }
    }

    //--------------------------Measured 2 point correlations-----------------------------------------------------//

    //--------------Measuring one point observable-------------------------------------------------------------//
    for(int set_no=0;set_no<one_point_obs;set_no++){
        if(loop==one_point_sweep_no[set_no] && LOOP_DIRECTION=="TO_RIGHT"){


            for(int sites_index=0;sites_index<one_point_sites[set_no][0].size();sites_index++){
                site_0=one_point_sites[set_no][0][sites_index];


                if(site_0<=ls[sys_iter]){
                    temp_coo_sys=One_point_opr_LB[set_no][sys_iter+1][sites_index];
                    temp_coo_env=Identity(H_RB[env_iter+1].nrows);
                }

                if(site_0>ls[sys_iter]){
                    temp_coo_sys=Identity(H_LB[sys_iter+1].nrows);
                    if(site_0==le[env_iter]){
                        temp_coo_env=Direct_Product(One_point_opr_onsite[set_no][0],Identity(H_LB[0].nrows));}
                    if(site_0==le[env_iter]+1){
                        temp_coo_env=Direct_Product(Identity(H_LB[0].nrows), One_point_opr_onsite[set_no][0]);}
                }

                vec_out[0].assign(Eig_vec.size(),0);
                Operate_Interface_interactions(Eig_vec, vec_out, 1.0, temp_coo_sys, temp_coo_env);
                OnePoint_Observable[set_no][site_0]=dot_product(Eig_vec,vec_out[0]);




            }
        }
    }

    //--------------Measured one point observable-------------------------------------------------------------//


    vec_out.clear();

    temp_coo_sys.value.clear();temp_coo_sys.rows.clear();temp_coo_sys.columns.clear();
    temp_coo_env.value.clear();temp_coo_env.rows.clear();temp_coo_env.columns.clear();

}

void DMRG::Writing_data(){

    string filepath = saving_filename;
    //ofstream outputfile(filepath.c_str(),std::ios::out|std::ios::binary);



    if(_BINARY==true){

        FILE * outputfile;
        outputfile = fopen(filepath.c_str(), "w");
        int _value_size;

        //SYSTEM
        for(int iter=0;iter<H_LB.size();iter++){

            fwrite (&iter, sizeof(int), 1,outputfile);

            //writing H_LB
            fwrite(& H_LB[iter].nrows,sizeof(int),1, outputfile);
            fwrite(& H_LB[iter].ncols,sizeof(int),1, outputfile);
            _value_size=H_LB[iter].value.size();
            fwrite(& _value_size,sizeof(int),1,outputfile);
            fwrite(& H_LB[iter].value[0],sizeof(double),
                    H_LB[iter].value.size(), outputfile);
            fwrite(& H_LB[iter].rows[0],sizeof(int),
                    H_LB[iter].rows.size(), outputfile);
            fwrite(& H_LB[iter].columns[0],sizeof(int),
                    H_LB[iter].columns.size(), outputfile);

            //writing Red_den_mat_system
            int temp=Red_den_mat_system[iter].size();
            fwrite(&temp ,sizeof(int),
                   1, outputfile);

            for(int m=0;m<Red_den_mat_system[iter].size();m++)
            {
                fwrite(& Red_den_mat_system[iter][m][0],sizeof(double),
                        Red_den_mat_system[iter].size(), outputfile);
            }


            //writing operators for all sites
            for(int site=0;site<OPSz_LB[iter].size();site++){

                fwrite (&site, sizeof(int), 1,outputfile);
                //             write(OPSz_LB[iter][site],outputfile);
                fwrite(& OPSz_LB[iter][site].nrows,sizeof(int),1, outputfile);
                fwrite(& OPSz_LB[iter][site].ncols,sizeof(int),1, outputfile);

                _value_size=OPSz_LB[iter][site].value.size();
                fwrite(& _value_size,sizeof(int),1,outputfile);

                fwrite(& OPSz_LB[iter][site].value[0],sizeof(double),
                        OPSz_LB[iter][site].value.size(), outputfile);
                fwrite(& OPSz_LB[iter][site].rows[0],sizeof(int),
                        OPSz_LB[iter][site].rows.size(), outputfile);
                fwrite(& OPSz_LB[iter][site].columns[0],sizeof(int),
                        OPSz_LB[iter][site].columns.size(), outputfile);

                //             write(OPSp_LB[iter][site],outputfile);
                fwrite(& OPSp_LB[iter][site].nrows,sizeof(int),1, outputfile);
                fwrite(& OPSp_LB[iter][site].ncols,sizeof(int),1, outputfile);

                _value_size=OPSp_LB[iter][site].value.size();
                fwrite(& _value_size,sizeof(int),1,outputfile);

                fwrite(& OPSp_LB[iter][site].value[0],sizeof(double),
                        OPSp_LB[iter][site].value.size(), outputfile);
                fwrite(& OPSp_LB[iter][site].rows[0],sizeof(int),
                        OPSp_LB[iter][site].rows.size(), outputfile);
                fwrite(& OPSp_LB[iter][site].columns[0],sizeof(int),
                        OPSp_LB[iter][site].columns.size(), outputfile);

                //             write(OPSm_LB[iter][site],outputfile);
                fwrite(& OPSm_LB[iter][site].nrows,sizeof(int),1, outputfile);
                fwrite(& OPSm_LB[iter][site].ncols,sizeof(int),1, outputfile);

                _value_size=OPSm_LB[iter][site].value.size();
                fwrite(& _value_size,sizeof(int),1,outputfile);

                fwrite(& OPSm_LB[iter][site].value[0],sizeof(double),
                        OPSm_LB[iter][site].value.size(), outputfile);
                fwrite(& OPSm_LB[iter][site].rows[0],sizeof(int),
                        OPSm_LB[iter][site].rows.size(), outputfile);
                fwrite(& OPSm_LB[iter][site].columns[0],sizeof(int),
                        OPSm_LB[iter][site].columns.size(), outputfile);


            }


        }



        //ENVIROMENT

        for(int iter=0;iter<H_RB.size();iter++){

            fwrite (&iter, sizeof(int), 1,outputfile);

            //writing H_RB
            fwrite(& H_RB[iter].nrows,sizeof(int),1, outputfile);
            fwrite(& H_RB[iter].ncols,sizeof(int),1, outputfile);
            _value_size=H_RB[iter].value.size();
            fwrite(& _value_size,sizeof(int),1,outputfile);
            fwrite(& H_RB[iter].value[0],sizeof(double),
                    H_RB[iter].value.size(), outputfile);
            fwrite(& H_RB[iter].rows[0],sizeof(int),
                    H_RB[iter].rows.size(), outputfile);
            fwrite(& H_RB[iter].columns[0],sizeof(int),
                    H_RB[iter].columns.size(), outputfile);

            //writing Red_den_mat_enviroment
            int temp=Red_den_mat_enviroment[iter].size();
            fwrite(&temp ,sizeof(int),
                   1, outputfile);

            for(int m=0;m<Red_den_mat_enviroment[iter].size();m++)
            {
                fwrite(& Red_den_mat_enviroment[iter][m][0],sizeof(double),
                        Red_den_mat_enviroment[iter].size(), outputfile);
            }


            //writing operators for all sites
            for(int site=0;site<OPSz_RB[iter].size();site++){

                fwrite (&site, sizeof(int), 1,outputfile);
                //             write(OPSz_RB[iter][site],outputfile);
                fwrite(& OPSz_RB[iter][site].nrows,sizeof(int),1, outputfile);
                fwrite(& OPSz_RB[iter][site].ncols,sizeof(int),1, outputfile);

                _value_size=OPSz_RB[iter][site].value.size();
                fwrite(& _value_size,sizeof(int),1,outputfile);

                fwrite(& OPSz_RB[iter][site].value[0],sizeof(double),
                        OPSz_RB[iter][site].value.size(), outputfile);
                fwrite(& OPSz_RB[iter][site].rows[0],sizeof(int),
                        OPSz_RB[iter][site].rows.size(), outputfile);
                fwrite(& OPSz_RB[iter][site].columns[0],sizeof(int),
                        OPSz_RB[iter][site].columns.size(), outputfile);

                //             write(OPSp_RB[iter][site],outputfile);
                fwrite(& OPSp_RB[iter][site].nrows,sizeof(int),1, outputfile);
                fwrite(& OPSp_RB[iter][site].ncols,sizeof(int),1, outputfile);

                _value_size=OPSp_RB[iter][site].value.size();
                fwrite(& _value_size,sizeof(int),1,outputfile);

                fwrite(& OPSp_RB[iter][site].value[0],sizeof(double),
                        OPSp_RB[iter][site].value.size(), outputfile);
                fwrite(& OPSp_RB[iter][site].rows[0],sizeof(int),
                        OPSp_RB[iter][site].rows.size(), outputfile);
                fwrite(& OPSp_RB[iter][site].columns[0],sizeof(int),
                        OPSp_RB[iter][site].columns.size(), outputfile);

                //             write(OPSm_RB[iter][site],outputfile);
                fwrite(& OPSm_RB[iter][site].nrows,sizeof(int),1, outputfile);
                fwrite(& OPSm_RB[iter][site].ncols,sizeof(int),1, outputfile);

                _value_size=OPSm_RB[iter][site].value.size();
                fwrite(& _value_size,sizeof(int),1,outputfile);

                fwrite(& OPSm_RB[iter][site].value[0],sizeof(double),
                        OPSm_RB[iter][site].value.size(), outputfile);
                fwrite(& OPSm_RB[iter][site].rows[0],sizeof(int),
                        OPSm_RB[iter][site].rows.size(), outputfile);
                fwrite(& OPSm_RB[iter][site].columns[0],sizeof(int),
                        OPSm_RB[iter][site].columns.size(), outputfile);

            }

        }


        //Transformed Eigenvector
        _value_size=Eig_vec_transformed.size();
        fwrite(& _value_size,sizeof(int),1,outputfile);
        fwrite(& Eig_vec_transformed[0] ,sizeof(double),_value_size,outputfile);


    }


    else{

        ofstream outputfile(filepath.c_str(),std::ios::out|std::ios::binary);
        outputfile<<"------------------------"<<endl;

        outputfile<<"SYSTEM INFORMATION : "<<endl;
        outputfile<<"------------------------"<<endl;

        for(int iter=0;iter<H_LB.size();iter++){
            outputfile<<"iter"<<endl;
            outputfile<<iter<<endl;

            //writing H_LB
            outputfile<<"H_LB"<<endl;
            write(H_LB[iter],outputfile);

            //writing Red_den_mat_system
            outputfile<<"Red_den_mat_system"<<endl;
            write(Red_den_mat_system[iter],outputfile);


            //writing operators for all sites
            for(int site=0;site<OPSz_LB[iter].size();site++){
                outputfile<<"site"<<endl;
                outputfile<<site<<endl;

                outputfile<<"OPSz_LB"<<endl;
                write(OPSz_LB[iter][site],outputfile);
                outputfile<<"OPSp_LB"<<endl;
                write(OPSp_LB[iter][site],outputfile);
                outputfile<<"OPSm_LB"<<endl;
                write(OPSm_LB[iter][site],outputfile);

            }


        }


        outputfile<<"------------------------"<<endl;
        outputfile<<"ENVIROMENT INFORMATION : "<<endl;
        outputfile<<"------------------------"<<endl;

        for(int iter=0;iter<H_RB.size();iter++){
            outputfile<<"iter"<<endl;
            outputfile<<iter<<endl;

            //writing H_RB
            outputfile<<"H_RB"<<endl;
            write(H_RB[iter],outputfile);

            //writing Red_den_mat_enviroment
            outputfile<<"Red_den_mat_enviroment"<<endl;
            write(Red_den_mat_enviroment[iter],outputfile);

            //writing operators for all sites
            for(int site=0;site<OPSz_RB[iter].size();site++){
                outputfile<<"site"<<endl;
                outputfile<<site<<endl;

                outputfile<<"OPSz_RB"<<endl;
                write(OPSz_RB[iter][site],outputfile);
                outputfile<<"OPSp_RB"<<endl;
                write(OPSp_RB[iter][site],outputfile);
                outputfile<<"OPSm_RB"<<endl;
                write(OPSm_RB[iter][site],outputfile);
            }


        }

        outputfile<<"------------------------"<<endl;
        outputfile<<"EIGENVECTOR TRANSFORMED SUPERBLOCK: "<<endl;
        outputfile<<"------------------------"<<endl;

        for(int i=0;i<Eig_vec_transformed.size();i++){
            outputfile<<Eig_vec_transformed[i]<<endl;
        }

    }




    //outputfile<<"All the data is written in this file"<<endl;


}

void DMRG::Reading_data(){

    cout<<"Data is being read"<<endl;
    string filepath = restart_filename;




    if(_BINARY==true){

        FILE * readingfile;
        readingfile = fopen(filepath.c_str(), "r");

        int _value_size,row_i,col_i;
        double _value_doub;
        int no_of_rows;

        //SYSTEM
        for(int iter=0;iter<H_LB.size();iter++){

            fread (&iter, sizeof(int), 1,readingfile);

            //reading H_LB
            fread(& H_LB[iter].nrows,sizeof(int),1, readingfile);
            fread(& H_LB[iter].ncols,sizeof(int),1, readingfile);

            fread(& _value_size,sizeof(int),1,readingfile);
            H_LB[iter].value.resize(_value_size);
            H_LB[iter].rows.resize(_value_size);
            H_LB[iter].columns.resize(_value_size);

            for(int val_i=0;val_i<_value_size;val_i++){
                fread(& _value_doub,sizeof(double),
                      1, readingfile);
                H_LB[iter].value[val_i]=_value_doub;
            }

            for(int val_i=0;val_i<_value_size;val_i++){
                fread(& row_i,sizeof(int),
                      1, readingfile);
                H_LB[iter].rows[val_i]=row_i;
            }

            for(int val_i=0;val_i<_value_size;val_i++){
                fread(& col_i,sizeof(int),
                      1, readingfile);
                H_LB[iter].columns[val_i]=col_i;
            }



            //reading Red_den_mat_system

            fread(&no_of_rows ,sizeof(int), 1, readingfile);
            Red_den_mat_system[iter].resize(no_of_rows);

            for(int i=0;i<no_of_rows;i++){
                Red_den_mat_system[iter][i].resize(no_of_rows);

                for(int j=0;j<no_of_rows;j++){
                    fread(&_value_doub,sizeof(double),
                          1, readingfile);
                    Red_den_mat_system[iter][i][j]=_value_doub;
                }

            }


            //reading operators for all sites
            for(int site=0;site<OPSz_LB[iter].size();site++){

                fread (&site, sizeof(int), 1,readingfile);

                //read(OPSz_LB[iter][site],outputfile);

                fread(& OPSz_LB[iter][site].nrows,sizeof(int),1, readingfile);
                fread(& OPSz_LB[iter][site].ncols,sizeof(int),1, readingfile);

                fread(& _value_size,sizeof(int),1,readingfile);
                OPSz_LB[iter][site].value.resize(_value_size);
                OPSz_LB[iter][site].rows.resize(_value_size);
                OPSz_LB[iter][site].columns.resize(_value_size);

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& _value_doub,sizeof(double),
                          1, readingfile);
                    OPSz_LB[iter][site].value[val_i]=_value_doub;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& row_i,sizeof(int),
                          1, readingfile);
                    OPSz_LB[iter][site].rows[val_i]=row_i;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& col_i,sizeof(int),
                          1, readingfile);
                    OPSz_LB[iter][site].columns[val_i]=col_i;
                }


                //read(OPSp_LB[iter][site],outputfile);

                fread(& OPSp_LB[iter][site].nrows,sizeof(int),1, readingfile);
                fread(& OPSp_LB[iter][site].ncols,sizeof(int),1, readingfile);

                fread(& _value_size,sizeof(int),1,readingfile);
                OPSp_LB[iter][site].value.resize(_value_size);
                OPSp_LB[iter][site].rows.resize(_value_size);
                OPSp_LB[iter][site].columns.resize(_value_size);

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& _value_doub,sizeof(double),
                          1, readingfile);
                    OPSp_LB[iter][site].value[val_i]=_value_doub;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& row_i,sizeof(int),
                          1, readingfile);
                    OPSp_LB[iter][site].rows[val_i]=row_i;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& col_i,sizeof(int),
                          1, readingfile);
                    OPSp_LB[iter][site].columns[val_i]=col_i;
                }


                //read(OPSm_LB[iter][site],outputfile);

                fread(& OPSm_LB[iter][site].nrows,sizeof(int),1, readingfile);
                fread(& OPSm_LB[iter][site].ncols,sizeof(int),1, readingfile);

                fread(& _value_size,sizeof(int),1,readingfile);
                OPSm_LB[iter][site].value.resize(_value_size);
                OPSm_LB[iter][site].rows.resize(_value_size);
                OPSm_LB[iter][site].columns.resize(_value_size);

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& _value_doub,sizeof(double),
                          1, readingfile);
                    OPSm_LB[iter][site].value[val_i]=_value_doub;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& row_i,sizeof(int),
                          1, readingfile);
                    OPSm_LB[iter][site].rows[val_i]=row_i;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& col_i,sizeof(int),
                          1, readingfile);
                    OPSm_LB[iter][site].columns[val_i]=col_i;
                }


            }//site


        }//iter



        //ENVIROMENT
        for(int iter=0;iter<H_RB.size();iter++){

            fread (&iter, sizeof(int), 1,readingfile);

            //reading H_RB
            fread(& H_RB[iter].nrows,sizeof(int),1, readingfile);
            fread(& H_RB[iter].ncols,sizeof(int),1, readingfile);

            fread(& _value_size,sizeof(int),1,readingfile);
            H_RB[iter].value.resize(_value_size);
            H_RB[iter].rows.resize(_value_size);
            H_RB[iter].columns.resize(_value_size);

            for(int val_i=0;val_i<_value_size;val_i++){
                fread(& _value_doub,sizeof(double),
                      1, readingfile);
                H_RB[iter].value[val_i]=_value_doub;
            }

            for(int val_i=0;val_i<_value_size;val_i++){
                fread(& row_i,sizeof(int),
                      1, readingfile);
                H_RB[iter].rows[val_i]=row_i;
            }

            for(int val_i=0;val_i<_value_size;val_i++){
                fread(& col_i,sizeof(int),
                      1, readingfile);
                H_RB[iter].columns[val_i]=col_i;
            }



            //reading Red_den_mat_enviroment

            fread(&no_of_rows ,sizeof(int), 1, readingfile);
            Red_den_mat_enviroment[iter].resize(no_of_rows);

            for(int i=0;i<no_of_rows;i++){
                Red_den_mat_enviroment[iter][i].resize(no_of_rows);

                for(int j=0;j<no_of_rows;j++){
                    fread(&_value_doub,sizeof(double),
                          1, readingfile);
                    Red_den_mat_enviroment[iter][i][j]=_value_doub;
                }

            }


            //reading operators for all sites
            for(int site=0;site<OPSz_RB[iter].size();site++){

                fread (&site, sizeof(int), 1,readingfile);

                //read(OPSz_RB[iter][site],outputfile);

                fread(& OPSz_RB[iter][site].nrows,sizeof(int),1, readingfile);
                fread(& OPSz_RB[iter][site].ncols,sizeof(int),1, readingfile);

                fread(& _value_size,sizeof(int),1,readingfile);
                OPSz_RB[iter][site].value.resize(_value_size);
                OPSz_RB[iter][site].rows.resize(_value_size);
                OPSz_RB[iter][site].columns.resize(_value_size);

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& _value_doub,sizeof(double),
                          1, readingfile);
                    OPSz_RB[iter][site].value[val_i]=_value_doub;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& row_i,sizeof(int),
                          1, readingfile);
                    OPSz_RB[iter][site].rows[val_i]=row_i;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& col_i,sizeof(int),
                          1, readingfile);
                    OPSz_RB[iter][site].columns[val_i]=col_i;
                }


                //read(OPSp_RB[iter][site],outputfile);

                fread(& OPSp_RB[iter][site].nrows,sizeof(int),1, readingfile);
                fread(& OPSp_RB[iter][site].ncols,sizeof(int),1, readingfile);

                fread(& _value_size,sizeof(int),1,readingfile);
                OPSp_RB[iter][site].value.resize(_value_size);
                OPSp_RB[iter][site].rows.resize(_value_size);
                OPSp_RB[iter][site].columns.resize(_value_size);

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& _value_doub,sizeof(double),
                          1, readingfile);
                    OPSp_RB[iter][site].value[val_i]=_value_doub;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& row_i,sizeof(int),
                          1, readingfile);
                    OPSp_RB[iter][site].rows[val_i]=row_i;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& col_i,sizeof(int),
                          1, readingfile);
                    OPSp_RB[iter][site].columns[val_i]=col_i;
                }


                //read(OPSm_RB[iter][site],outputfile);

                fread(& OPSm_RB[iter][site].nrows,sizeof(int),1, readingfile);
                fread(& OPSm_RB[iter][site].ncols,sizeof(int),1, readingfile);

                fread(& _value_size,sizeof(int),1,readingfile);
                OPSm_RB[iter][site].value.resize(_value_size);
                OPSm_RB[iter][site].rows.resize(_value_size);
                OPSm_RB[iter][site].columns.resize(_value_size);

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& _value_doub,sizeof(double),
                          1, readingfile);
                    OPSm_RB[iter][site].value[val_i]=_value_doub;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& row_i,sizeof(int),
                          1, readingfile);
                    OPSm_RB[iter][site].rows[val_i]=row_i;
                }

                for(int val_i=0;val_i<_value_size;val_i++){
                    fread(& col_i,sizeof(int),
                          1, readingfile);
                    OPSm_RB[iter][site].columns[val_i]=col_i;
                }

            }//site

        }//iter


        //Read Transformed Eigenvector

        fread(& _value_size,sizeof(int),1,readingfile);
        Eig_vec_transformed.resize(_value_size);

        for(int i=0;i<_value_size;i++){
            fread(& _value_doub ,sizeof(double),1,readingfile);
            Eig_vec_transformed[i]=_value_doub;
        }

    }//if binary==true




}

//bool DMRG::Test_Hermiticity(int sys_iter, int env_iter){

//    bool bool_to_proceed;
//    if(symmetry_used =="Total_Sz"){bool_to_proceed = false;}else {bool_to_proceed = true;}

//    Mat_3_doub vec_in, vec_out;
//    Mat_2_doub Hamil_Mat;

//    vec_in.resize(Q_Eff_System[sys_iter].size());	vec_out.resize(Q_Eff_System[sys_iter].size());
//    int tmp_sz=0;
//    for(int is=0;is<Q_Eff_System[sys_iter].size();is++){
//        vec_in[is].resize(Q_Eff_Enviroment[env_iter].size());
//        vec_out[is].resize(Q_Eff_Enviroment[env_iter].size());
//        for(int ie=0;ie<Q_Eff_Enviroment[env_iter].size();ie++){
//            if((bool_to_proceed==true) || (Q_Eff_System[sys_iter][is]+Q_Eff_Enviroment[env_iter][ie]==Q_Eff_Super_Block[max(sys_iter,env_iter)])){
//                tmp_sz=tmp_sz + H_LB[sys_iter+1][is].nrows*H_RB[env_iter+1][ie].nrows;
//                vec_in[is][ie].resize(H_LB[sys_iter+1][is].nrows*H_RB[env_iter+1][ie].nrows);
//                for(int n=0;n<H_LB[sys_iter+1][is].nrows*H_RB[env_iter+1][ie].nrows;n++){
//                    vec_in[is][ie][n]=0;

//                }
//            }
//        }

//    }

//    Hamil_Mat.resize(tmp_sz);
//    for(int i=0;i<tmp_sz;i++){
//        Hamil_Mat[i].resize(tmp_sz);
//    }

//    int row_no=0;
//    int col_no=0;
//    for(int i=0;i<Q_Eff_System[sys_iter].size();i++){
//        for(int l=0;l<Q_Eff_Enviroment[env_iter].size();l++){
//            if((bool_to_proceed==true) ||
//                    (Q_Eff_System[sys_iter][i]+Q_Eff_Enviroment[env_iter][l]
//                     ==Q_Eff_Super_Block[max(sys_iter,env_iter)])){

//                for(int n=0;n<H_LB[sys_iter+1][i].nrows*H_RB[env_iter+1][l].nrows;n++){

//                    vec_in[i][l][n]=1;

//                    Operate_H_SB(vec_in,sys_iter,env_iter,vec_out);
//                    vec_in[i][l][n]=0;

//                    row_no=0;

//                    for(int i2=0;i2<Q_Eff_System[sys_iter].size();i2++){
//                        for(int l2=0;l2<Q_Eff_Enviroment[env_iter].size();l2++){
//                            if((bool_to_proceed==true) ||
//                                    (Q_Eff_System[sys_iter][i2]+Q_Eff_Enviroment[env_iter][l2]
//                                     ==Q_Eff_Super_Block[max(sys_iter,env_iter)])){

//                                for(int n2=0;n2<H_LB[sys_iter+1][i2].nrows*H_RB[env_iter+1][l2].nrows;n2++){
//                                    Hamil_Mat[row_no][col_no]=vec_out[i2][l2][n2];

//                                    row_no=row_no+1;
//                                }

//                            }
//                        }
//                    }





//                    col_no= col_no+1;
//                }
//            }
//        }
//    }

//    cout<<"PRINTING HAMIL_SB ::"<<endl;
//    bool check=true;
//    bool check_tmp;
//    for(int i=0;i<tmp_sz;i++){
//        for(int j=0;j<tmp_sz;j++){
//            cout<<Hamil_Mat[i][j]<<"    ";
//            check_tmp=(Hamil_Mat[i][j]==Hamil_Mat[j][i]);
//            check=(check && check_tmp);


//        }
//        cout<<endl;
//    }
//    cout<<"PRINTED HAMIL_SB ::"<<endl;
//    return check;

//}
