#include <iostream>  //for cin and cout
#include "functions.h"
//#include "dmrg_Hubbard_model.h"
#include "DMRG_keeper3.h"
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include "iconprint.h"
using namespace std;  // instead of std::cin, now just cin is enough

// To compile "g++ reading_input.cpp functions.cpp  dmrg_solver.cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm"



int main(int argc, char** argv  ){


    iconprint();

    cout<<"Using input file = "<<argv[1]<<endl;



    DMRG DM_RG;

    DM_RG._BINARY=true;

    DM_RG.inp_filename = argv[1];
    DM_RG.read_INPUT();

    bool Wft_trns;
    Wft_trns = DM_RG.Wavefuntion_transformation;
    DM_RG.Wavefuntion_transformation=false;


    DM_RG.Create_Q_Eff();


    DM_RG.Initialize_Hamiltonians();


    DM_RG.Initialize_Oprts();

    DM_RG.Initialize_onsite_oprs_for_corrs();


    DM_RG.LOOP_DIRECTION="CENTER";

    string out_string;
    out_string=DM_RG.inp_filename;
    if (out_string.size () > 0)  out_string.resize (out_string.size () - 4);

    string E_string = "Energy_"+out_string+".txt";
    ofstream Energy_vs_length(E_string.c_str());
    Energy_vs_length.precision(20);

    string Sz_string = "Sz_val_"+out_string+".txt";
    ofstream Sz_val(Sz_string.c_str());
    Sz_val.precision(10);

    if(DM_RG._RESTART==true){
        DM_RG.Reading_data();
        DM_RG.Wavefuntion_transformation=Wft_trns;
    }
    else{
        cout<<"_________________________________________________________________________________________________"<<endl;
        cout<<"============================STARTING INFINITE DMRG ALGORITHM====================================="<<endl;
        cout<<"-------------------------------------------------------------------------------------------------"<<endl;
        DM_RG.nstatestouse=DM_RG.m_infinite;
        DM_RG.ALGO="INFINITE";

        for(int iteration=0;iteration<=DM_RG.max_iter;iteration++){


            cout<<endl;
            cout <<"Infinite iteration no. is "<<iteration<<endl;





            clock_t LB_growth = clock();
            DM_RG.Grow_LB_and_update_Spin_oprts(iteration);
            cout<<"Time for LB growth : "<<double( clock() - LB_growth ) / (double)CLOCKS_PER_SEC<<endl;


            clock_t RB_growth = clock();
            DM_RG.Grow_RB_and_update_Spin_oprts(iteration);
            cout<<"Time for RB growth : "<<double( clock() - RB_growth ) / (double)CLOCKS_PER_SEC<<endl;








            clock_t lanc = clock();
            DM_RG.Perform_LANCZOS(iteration,iteration);
            cout<<"Time for doing lanczos : "<<double( clock() - lanc ) / (double)CLOCKS_PER_SEC<<endl;


            if(iteration!=DM_RG.max_iter){
                clock_t rnrm = clock();
                DM_RG.Do_RENORMALIZATION_of_S_and_E(iteration,iteration,0,0,0);
                cout<<"Time for doing Renormalization : "<<double( clock() - rnrm ) / (double)CLOCKS_PER_SEC<<endl;
            }





            Energy_vs_length<<2*(iteration)+4<<scientific<<"\t"<<DM_RG.Energy<<"\t"<<DM_RG.Lanc_Error_out<<"\t"<<DM_RG.Truncation_Error_S<<"\t"<<DM_RG.Truncation_Error_E<<endl;


            //        cout<<"stop infinite iter"<<iteration<<endl;
            //        getchar();

        }
        cout<<"_________________________________________________________________________________________________"<<endl;
        cout<<"================================INFINITE DMRG ALGORITHM DONE====================================="<<endl;
        cout<<"-------------------------------------------------------------------------------------------------"<<endl<<endl;



        if(DM_RG.Finite_algo_bool==true){
        DM_RG.Wavefuntion_transformation=Wft_trns;
        DM_RG.Do_RENORMALIZATION_of_S_and_E(DM_RG.max_iter,DM_RG.max_iter,0,0,0);}


        } //_RESTART==false









    if(DM_RG.Finite_algo_bool==true){

        clock_t rnrm = clock();




        cout<<"Time for doing Renormalization : "<<double( clock() - rnrm ) / (double)CLOCKS_PER_SEC<<endl;

        cout<<"_________________________________________________________________________________________________"<<endl;
        cout<<"================================STARTING FINITE DMRG ALGORITHM==================================="<<endl;
        cout<<"-------------------------------------------------------------------------------------------------"<<endl;
        DM_RG.ALGO="FINITE";
        cout<<"No. of Finite loops to do = "<<DM_RG.Finite_loops.size()<<endl<<endl;

        for(int loop=0;loop<DM_RG.Finite_loops.size();loop++){

            if(DM_RG.Finite_loops[loop]<0){
                DM_RG.LOOP_DIRECTION="TO_LEFT";
            }
            else if(DM_RG.Finite_loops[loop]>0){
                DM_RG.LOOP_DIRECTION="TO_RIGHT";
            }



            DM_RG.nstatestouse = DM_RG.Finite_states[loop];
            DM_RG.m_infinite_max = DM_RG.Finite_states_max[loop];

            for(int loop_iter=1;loop_iter<=abs(DM_RG.Finite_loops[loop]);loop_iter++){

                int env_i, sys_i;

                DM_RG.Initialize_corr_operators(loop, loop_iter);
                //                cout<<"stop finite loop"<<loop<<" loop_iter"<<loop_iter<<
                //                      "Initialize_corr_operators(loop, loop_iter) "<<endl;
                //                getchar();


                if(DM_RG._RESTART==true){

                    if(DM_RG.Finite_loops[loop] > 0){
                        env_i = 2*DM_RG.max_iter - loop_iter;
                        sys_i = loop_iter;
                        cout<<"SUPERBLOCK = "<< loop_iter + 2 <<"+"<< DM_RG.Target_L - 2 -loop_iter <<endl;

                    }
                    else if(DM_RG.Finite_loops[loop] < 0){
                        env_i = loop_iter;
                        sys_i = 2*DM_RG.max_iter - loop_iter;
                        cout<<"SUPERBLOCK = "<<DM_RG.Target_L - 2 -loop_iter  <<"+"<<  loop_iter + 2<<endl;

                    }

                }
                else{
                if(loop==0 && DM_RG.Finite_loops[loop] < 0){
                    env_i = DM_RG.max_iter + loop_iter;
                    sys_i = DM_RG.max_iter - loop_iter;
                    cout<<"SUPERBLOCK = "<< (int)(0.5*(DM_RG.Target_L) -loop_iter) <<"+"<<(int)(0.5*(DM_RG.Target_L)+loop_iter)<<endl;

                }
                else if(DM_RG.Finite_loops[loop] > 0){
                    env_i = 2*DM_RG.max_iter - loop_iter;
                    sys_i = loop_iter;
                    cout<<"SUPERBLOCK = "<< loop_iter + 2 <<"+"<< DM_RG.Target_L - 2 -loop_iter <<endl;

                }
                else if(DM_RG.Finite_loops[loop] < 0){
                    env_i = loop_iter;
                    sys_i = 2*DM_RG.max_iter - loop_iter;
                    cout<<"SUPERBLOCK = "<<DM_RG.Target_L - 2 -loop_iter  <<"+"<<  loop_iter + 2<<endl;

                }

                }





                clock_t RB_growth = clock();
                DM_RG.Grow_RB_and_update_Spin_oprts(env_i);
                cout<<"Time for RB growth : "<<double( clock() - RB_growth ) / (double)CLOCKS_PER_SEC<<endl;
                //                cout<<"stop finite loop"<<loop<<" loop_iter"<<loop_iter<<
                //                      "Grow_RB_and_update_Spin_oprts "<<endl;
                //                getchar();
                clock_t LB_growth = clock();
                DM_RG.Grow_LB_and_update_Spin_oprts(sys_i);
                cout<<"Time for LB growth : "<<double( clock() - LB_growth ) / (double)CLOCKS_PER_SEC<<endl;
                //                cout<<"stop finite loop"<<loop<<" loop_iter"<<loop_iter<<
                //                      "Grow_LB_and_update_Spin_oprts "<<endl;
                //                getchar();


                DM_RG.Update_n_point_corr_oprs(loop,loop_iter,env_i,sys_i);

                cout<<"AFTER GROWING :"<<endl;



                clock_t lanc = clock();
                DM_RG.Perform_LANCZOS(sys_i,env_i);
                cout<<"Time for doing lanczos : "<<double( clock() - lanc ) / (double)CLOCKS_PER_SEC<<endl;
                //                cout<<"stop finite loop"<<loop<<" loop_iter"<<loop_iter<<
                //                      "Perform_LANCZOS "<<endl;
                //                getchar();
                if(!((loop_iter==abs(DM_RG.Finite_loops[loop]) && (loop==DM_RG.Finite_loops.size()-1) ) )) {
                    rnrm = clock();
                    DM_RG.Do_RENORMALIZATION_of_S_and_E(sys_i,env_i,DM_RG.Finite_loops[loop],loop_iter,loop);
                    cout<<"Time for doing Renormalization and WFT : "<<double( clock() - rnrm ) / (double)CLOCKS_PER_SEC<<endl;
                    cout<<"Renormalization done"<<endl<<endl;

                }
                else{



                    DM_RG.Measure_observables(sys_i, env_i, loop);
                    cout<<"Observables calculation done"<<endl;

                    double Total_one_point_obs=0;
                    string one_p_file="one_point_obs_set_no_";

                    for(int set_no=0;set_no<DM_RG.one_point_obs;set_no++){

                        string setstr;
                        ostringstream temp;
                        temp<<set_no;
                        setstr=temp.str();

                        one_p_file=one_p_file + setstr +"_"+out_string+".txt";
                        ofstream One_p_obs(one_p_file.c_str());
                        One_p_obs.precision(10);
                        One_p_obs<<"<GS|"<<DM_RG.one_point_oprnames[set_no][0]<<"|GS>"<<endl;
                        Total_one_point_obs=0;
                        for(int i=0;i<DM_RG.Target_L;i++){
                                One_p_obs<<i<<"    "<<DM_RG.OnePoint_Observable[set_no][i]<<endl;
                                Total_one_point_obs = Total_one_point_obs +
                                                    DM_RG.OnePoint_Observable[set_no][i];


                        }
                        one_p_file="one_point_obs_set_no_";

                    cout<<"GROUND STATE LIES IN SECTOR, Total_"<<DM_RG.one_point_oprnames[set_no][0]
                       <<"(or you gave this sector) = "<<Total_one_point_obs<<endl;



                    }



                    cout<<"Spin-Spin Correlation Calculation done"<<endl;


                    string two_p_file="two_point_corr_set_no_";
                    for(int set_no=0;set_no<DM_RG.two_point_corrs;set_no++){

                        string setstr;
                        ostringstream temp;
                        temp<<set_no;
                        setstr=temp.str();

                        two_p_file=two_p_file + setstr +"_"+out_string+".txt";
                        ofstream Two_p_corr(two_p_file.c_str());
                        Two_p_corr.precision(5);
                        Two_p_corr<<"<GS|"<<DM_RG.two_point_oprnames[set_no][0]<<""<<DM_RG.two_point_oprnames[set_no][1]<<"|GS>"<<endl;
                        for(int i=0;i<DM_RG.Target_L;i++){
                            for(int j=0;j<DM_RG.Target_L;j++){
                                Two_p_corr<<DM_RG.TwoPointCorrs[set_no][i][j]<<"       ";
                            }
                            Two_p_corr<<endl<<endl;
                        }
                        two_p_file="two_point_corr_set_no_";
                    }

                    int site_0, site_1, site_2, site_3;
                    string four_p_file="four_point_corr_set_no_";
                    for(int set_no=0;set_no<DM_RG.four_point_corrs;set_no++){

                        string setstr;
                        ostringstream temp;
                        temp<<set_no;
                        setstr=temp.str();

                        four_p_file=four_p_file + setstr +"_"+out_string+".txt";
                        ofstream Four_p_corr(four_p_file.c_str());
                        Four_p_corr.precision(5);
                        Four_p_corr<<"<GS|"<<DM_RG.four_point_oprnames[set_no][0]<<" "<<DM_RG.four_point_oprnames[set_no][1]<<" "<<
                                     DM_RG.four_point_oprnames[set_no][2]<<" "<<DM_RG.four_point_oprnames[set_no][3]<<"|GS>"<<endl<<endl;

                        for(int sites_index=0;sites_index<DM_RG.four_point_sites[set_no][0].size();sites_index++){
                            site_0=DM_RG.four_point_sites[set_no][0][sites_index];
                            site_1=DM_RG.four_point_sites[set_no][1][sites_index];
                            site_2=DM_RG.four_point_sites[set_no][2][sites_index];
                            site_3=DM_RG.four_point_sites[set_no][3][sites_index];

                            Four_p_corr<<site_0<<" "<<site_1<<" "<<site_2<<" "<<site_3<<"    "<<DM_RG.FourPointCorrs[set_no][sites_index]<<endl;
                        }



                        four_p_file="four_point_corr_set_no_";
                    }

                    if(DM_RG._SAVING==true){
                        rnrm = clock();
                        DM_RG.Do_RENORMALIZATION_of_S_and_E(sys_i,env_i,DM_RG.Finite_loops[loop],loop_iter,loop);
                        cout<<"Time for doing Renormalization and WFT : "<<double( clock() - rnrm ) / (double)CLOCKS_PER_SEC<<endl;
                        cout<<"Renormalization done because Saving is done"<<endl<<endl;

                        cout<<"Now writing data to file "<< DM_RG.saving_filename<<endl;
                        DM_RG.Writing_data();
                    }
                    else{
                        cout<<"Renormalization is not done for the last iteration as SAVING is not done"<<endl<<endl;
                    }

                }


                Energy_vs_length<<DM_RG.Target_L<<scientific<<"\t"<<DM_RG.Energy<<"\t"<<DM_RG.Lanc_Error_out<<"\t"<<DM_RG.Truncation_Error_S<<"\t"<<DM_RG.Truncation_Error_E<<endl;



                //                cout<<"stop finite loop"<<loop<<" loop_iter"<<loop_iter<<endl;
                //                getchar();

            }



        }


        cout<<"_________________________________________________________________________________________________"<<endl;
        cout<<"================================FINITE DMRG ALGORITHM DONE======================================="<<endl;
        cout<<"-------------------------------------------------------------------------------------------------"<<endl<<endl;


    }
    else{cout<<"No FINITE ALGORITHM REQUESTED"<<endl<<endl;}

    cout<<"_________________________________________________________________________________________________"<<endl;
    cout<<"================================STARTING CALCULATION OF OBSERVABLES==============================="<<endl;
    cout<<"-------------------------------------------------------------------------------------------------"<<endl<<endl;


    cout<<"_________________________________________________________________________________________________"<<endl;
    cout<<"================================CALCULATION OF OBSERVABLES DONE==============================="<<endl;
    cout<<"-------------------------------------------------------------------------------------------------"<<endl<<endl;
    //getchar();

    return 0;
}









