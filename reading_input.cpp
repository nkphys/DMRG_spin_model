//#include "reading_input.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "tensor_type.h"
using namespace std;

//void reading_input(string filename, double & J1_p, double & J1_m, double & J1_z, double & J2_p, double & J2_m, double & J_z, double & H_mag, double & T_Sz, int & m_infinite);


void reading_input(string filename, string & Mag_field_file, string & J_zz_long_range_file,
                   string & J_pm_long_range_file, string & J_mp_long_range_file,
                   string & J_pp_long_range_file, string & J_mm_long_range_file,
                   int & m_infinite, int & m_infinite_max, int & Target_L,int & Target_state,
                   double & Truncation_Error_Target, double & Lanc_Error,
                   Mat_1_int & Finite_loops, Mat_1_int & Finite_states,
                   Mat_1_int & Finite_states_max, bool & Finite_algo_bool,
                   bool & Wavefuntion_transformation, int & max_lanczos_states, int & N_corr_sets,
                   Mat_1_int & N_point, Mat_2_string & N_point_filenames,
                   Mat_1_string & N_point_range, bool & NPoint_Corrs,
                   Mat_1_int & sweep_no){


    string filepath = filename;
    string target_length,Target_Length = "Target Length = ";
    string target_state,Target_State = "Target state in standard DMRG = ";
    string max_LANCZOS_states,Max_LANCZOS_States = "Max LANCZOS states = ";
    string dmrg_truncation_error,DMRG_Truncation_Error = "DMRG Truncation Error = ";
    string lanczos_error,Lanczos_Error = "LANCZOS Error = ";
    string jzz_longrange,Jzz_LongRange = "Jzz_longrange = ";
    string jpm_longrange,Jpm_LongRange = "Jpm_longrange = ";
    string jmp_longrange,Jmp_LongRange = "Jmp_longrange = ";
    string jpp_longrange,Jpp_LongRange = "Jpp_longrange = ";
    string jmm_longrange,Jmm_LongRange = "Jmm_longrange = ";
    string h_ext,H_Ext ="H_ext = ";
    string ncorrsets, NCorrSets= "Number of sets of correlations = ";
    string minf, Minf= "No.of states infinite = ";
    string minf_max, Minf_max= "No.of maxstates infinite = ";
    string finiteloops,FiniteLoops= "DMRG FINITE LOOPS = ";
    string finitestates, FiniteStates= "No.of states finite = ";
    string finitestates_max, FiniteStates_max= "No.of maxstates finite = ";
    string finitealgobool, FiniteAlgoBool="FINITE_ALGORITHM_BOOL = ";
    string wavefunctiontransformation, WavefunctionTransformation = "Wavefunction Transformation bool = ";
    string npointcorrs, NPointCorrs = "Calculating demanded set of Correlations = ";
    string npointf="Operators for N_points Observable set_";
    string npointr="N_points Observable site range set_";
    string sweepno_str="Finite Sweep no for set_";
    Mat_1_string npointfilenames, NPointFilenames;
    Mat_1_string sweepno, SweepNo;
    Mat_1_string npointrange, NpointRange;
    int offset;
    string line;
    ifstream inputfile(filepath.c_str());
    

    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);



            if ((offset = line.find(Jzz_LongRange, 0)) != string::npos) {
                jzz_longrange = line.substr (offset+Jzz_LongRange.length());				}
            if ((offset = line.find(Jpm_LongRange, 0)) != string::npos) {
                jpm_longrange = line.substr (offset+Jpm_LongRange.length());				}
            if ((offset = line.find(Jmp_LongRange, 0)) != string::npos) {
                jmp_longrange = line.substr (offset+Jmp_LongRange.length());				}
            if ((offset = line.find(Jpp_LongRange, 0)) != string::npos) {
                jpp_longrange = line.substr (offset+Jpp_LongRange.length());				}
            if ((offset = line.find(Jmm_LongRange, 0)) != string::npos) {
                jmm_longrange = line.substr (offset+Jmm_LongRange.length());				}
            if ((offset = line.find(H_Ext, 0)) != string::npos) {
                h_ext = line.substr (offset+H_Ext.length());				}


            if ((offset = line.find(Target_Length, 0)) != string::npos) {
                target_length = line.substr (offset + Target_Length.length());		}

            if ((offset = line.find(Target_State, 0)) != string::npos) {
                target_state = line.substr (offset + Target_State.length());		}

            if ((offset = line.find(Max_LANCZOS_States, 0)) != string::npos) {
                max_LANCZOS_states = line.substr (offset + Max_LANCZOS_States.length());		}


            if ((offset = line.find(DMRG_Truncation_Error, 0)) != string::npos) {
                dmrg_truncation_error = line.substr (offset + DMRG_Truncation_Error.length());}


            if ((offset = line.find(Lanczos_Error, 0)) != string::npos) {
                lanczos_error = line.substr (offset+Lanczos_Error.length());				}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file 1st time."<<endl;}


    inputfile.open(filepath.c_str());
    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Minf, 0)) != string::npos) {
                minf = line.substr (offset+Minf.length());				}

            if ((offset = line.find(Minf_max, 0)) != string::npos) {
                minf_max = line.substr (offset+Minf_max.length());				}

            if ((offset = line.find(FiniteAlgoBool, 0)) != string::npos) {
                finitealgobool = line.substr (offset+FiniteAlgoBool.length());				}

            if ((offset = line.find(WavefunctionTransformation, 0)) != string::npos) {
                wavefunctiontransformation = line.substr (offset+WavefunctionTransformation.length());	}


            if ((offset = line.find(NPointCorrs, 0)) != string::npos) {
                npointcorrs = line.substr (offset+NPointCorrs.length());	}


            if ((offset = line.find(FiniteLoops, 0)) != string::npos) {
                finiteloops = line.substr (offset+FiniteLoops.length());				}

            stringstream finiteloops_stream(finiteloops);
            int num;

            finiteloops_stream >> num ;
            Finite_loops.resize(num);

            for(int n=0;n<num;n++){
                finiteloops_stream >> Finite_loops[n];
            }

            if ((offset = line.find(FiniteStates, 0)) != string::npos) {
                finitestates = line.substr (offset+FiniteStates.length());				}

            stringstream finitestates_stream(finitestates);

            finitestates_stream >> num ;
            Finite_states.resize(num);

            for(int n=0;n<num;n++){
                finitestates_stream >> Finite_states[n];
            }

            if ((offset = line.find(NCorrSets, 0)) != string::npos) {
                ncorrsets = line.substr (offset+NCorrSets.length());}

            stringstream ncorrsets_stream(ncorrsets);
            N_corr_sets=0;
            ncorrsets_stream >> N_corr_sets;
            NPointFilenames.resize(N_corr_sets);npointfilenames.resize(N_corr_sets);
            NpointRange.resize(N_corr_sets);npointrange.resize(N_corr_sets);
            N_point_range.resize(N_corr_sets);
            N_point_filenames.resize(N_corr_sets);N_point.resize(N_corr_sets);
            sweep_no.resize(N_corr_sets);SweepNo.resize(N_corr_sets);
            sweepno.resize(N_corr_sets);

            for(int set_ind=0;set_ind<N_corr_sets;set_ind++){
                string strind ;
                stringstream strind_ss;
                strind_ss << set_ind+1;
                //cout<<"strind_ss : "<<strind_ss<<endl;

                strind = strind_ss.str();
                //cout<<"strind : "<<strind<<endl;
                NPointFilenames[set_ind] = npointf + strind + " = ";
                NpointRange[set_ind] = npointr + strind + " = ";
                SweepNo[set_ind] = sweepno_str + strind + " = ";

                //cout<<"NPointFilenames[set_ind] : "<<NPointFilenames[set_ind]<<endl;
                // cout<<"NpointRange[set_ind] : "<<NpointRange[set_ind]<<endl;

                if ((offset = line.find(NPointFilenames[set_ind], 0)) != string::npos) {
                    npointfilenames[set_ind] = line.substr (offset+NPointFilenames[set_ind].length());
                    stringstream npointfilenames_stream(npointfilenames[set_ind]);
                    N_point[set_ind]=0;
                    npointfilenames_stream >> N_point[set_ind];
                    // cout<<"npointfilenames : "<<npointfilenames[set_ind]<<endl;
                    //cout<<"npointfilenames_stream : "<<npointfilenames_stream<<endl;
                    //cout<<"N_point[set_ind] : "<<N_point[set_ind]<<endl;
                    N_point_filenames[set_ind].resize(N_point[set_ind]);

                    for(int n=0;n<N_point[set_ind];n++){
                        npointfilenames_stream >> N_point_filenames[set_ind][n];
                        // cout<<n<<" file : "<<N_point_filenames[set_ind][n]<<endl;
                    }
                }

                if ((offset = line.find(NpointRange[set_ind], 0)) != string::npos) {
                    npointrange[set_ind] = line.substr (offset+NpointRange[set_ind].length());
                    stringstream npointrange_stream(npointrange[set_ind]);
                    npointrange_stream >> N_point_range[set_ind];
                }

                if ((offset = line.find(SweepNo[set_ind], 0)) != string::npos) {
                    sweepno[set_ind] = line.substr (offset+SweepNo[set_ind].length());
                    stringstream sweepno_stream(sweepno[set_ind]);
                    sweepno_stream >> sweep_no[set_ind];
                }
            }



            if ((offset = line.find(FiniteStates_max, 0)) != string::npos) {
                finitestates_max = line.substr (offset+FiniteStates_max.length());				}

            stringstream finitestates_max_stream(finitestates_max);

            finitestates_max_stream >> num ;
            Finite_states_max.resize(num);

            for(int n=0;n<num;n++){
                finitestates_max_stream >> Finite_states_max[n];
            }
        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open the input file 2nd time."<<endl;}



    m_infinite =atoi(minf.c_str());
    m_infinite_max =atoi(minf_max.c_str());
    Target_L =atoi(target_length.c_str());
    Target_state=atoi(target_state.c_str());
    max_lanczos_states=atoi(max_LANCZOS_states.c_str());
    Truncation_Error_Target = atof(dmrg_truncation_error.c_str());
    Lanc_Error = atof(lanczos_error.c_str());
    J_zz_long_range_file=jzz_longrange;
    J_pm_long_range_file=jpm_longrange;
    J_mp_long_range_file=jmp_longrange;
    J_pp_long_range_file=jpp_longrange;
    J_mm_long_range_file=jmm_longrange;
    Mag_field_file=h_ext;



    if(finitealgobool=="true"){Finite_algo_bool= true ;}else{Finite_algo_bool=false;}
    if(wavefunctiontransformation=="true"){Wavefuntion_transformation= true ;}else{Wavefuntion_transformation=false;}
    if(npointcorrs=="true"){NPoint_Corrs= true ;}else{NPoint_Corrs=false;}

}


void reading_restart_or_saving(bool & _RESTART, bool & _SAVING, string & saving_filename,
                               string & restart_filename, string filename){

    string filepath = filename;
    ifstream inputfile(filepath.c_str());
    string line;
    int offset;

    string saving_file, SAVING_FILE = "Saving in this file = ";
    string restart_file, RESTART_FILE = "Restart from this file = ";
    string saving_bool, SAVING_BOOL = "Saving_bool = ";
    string restart_bool, RESTART_BOOL = "Restart_bool = ";


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(SAVING_BOOL, 0)) != string::npos) {
                saving_bool = line.substr (offset+SAVING_BOOL.length());				}

            if ((offset = line.find(RESTART_BOOL, 0)) != string::npos) {
                restart_bool = line.substr (offset+RESTART_BOOL.length());				}

            if ((offset = line.find(SAVING_FILE, 0)) != string::npos) {
                saving_file = line.substr (offset+SAVING_FILE.length());				}

            if ((offset = line.find(RESTART_FILE, 0)) != string::npos) {
                restart_file = line.substr (offset+RESTART_FILE.length());				}

        }
        inputfile.close();
    }

    saving_filename = saving_file;
    restart_filename = restart_file;
    if(saving_bool=="true"){_SAVING=true;}else{_SAVING=false;}
    if(restart_bool=="true"){_RESTART=true;}else{_RESTART=false;}



}

void reading_sites_for_correlations(Mat_1_string N_point_range, Mat_3_int & sites_corrs){

    string line;
    int site_tmp;

    for(int set_ind=0;set_ind<N_point_range.size();set_ind++){
        string filepath = N_point_range[set_ind];

        ifstream inputfile(filepath.c_str());

        if(inputfile.is_open())
        {
            while(getline(inputfile,line)){
                stringstream line_stream(line);

                for(int i=0;i<sites_corrs[set_ind].size();i++){

                    //sites_corrs[set_ind][i].resize(no_values);
                    line_stream >> site_tmp;
                    sites_corrs[set_ind][i].push_back(site_tmp);
                }
            }

            inputfile.close();
        }

    }

}

void reading_long_range_connections(string Mag_field_file, string J_zz_long_range_file,
                                    string J_pm_long_range_file, string J_mp_long_range_file,
                                    string J_pp_long_range_file, string J_mm_long_range_file,
                                    Mat_1_doub & Mag_field, Mat_2_doub & J_zz_long_range,
                                    Mat_2_doub & J_pm_long_range, Mat_2_doub & J_mp_long_range,
                                    Mat_2_doub & J_pp_long_range, Mat_2_doub & J_mm_long_range, int Target_L){

    string filearray[5]={J_zz_long_range_file,J_pm_long_range_file, J_mp_long_range_file,
                         J_pp_long_range_file,J_mm_long_range_file};
    string line;
    type_double tmp;
    Mat_2_doub tmp_mat;
    tmp_mat.resize(Target_L);

    for(int fileno=0;fileno<5;fileno++){
        tmp_mat.clear();
        tmp_mat.resize(Target_L);
        ifstream inputfile(filearray[fileno].c_str());
        if(inputfile.is_open())
        {
            while(getline(inputfile,line)){
                stringstream line_stream(line);
                for(int i=0;i<Target_L;i++){
#ifndef WITH_COMPLEX
                    line_stream >> tmp;
#endif
#ifdef WITH_COMPLEX
            line_stream >> tmp.real();
            tmp.imag()=0.0;
#endif
                    tmp_mat[i].push_back(tmp);

                }
            }

            inputfile.close();
        }

        if(filearray[fileno]==J_zz_long_range_file){
            J_zz_long_range = tmp_mat;
        }
        if(filearray[fileno]==J_pm_long_range_file){
            J_pm_long_range = tmp_mat;
        }
        if(filearray[fileno]==J_mp_long_range_file){
            J_mp_long_range = tmp_mat;
        }
        if(filearray[fileno]==J_pp_long_range_file){
            J_pp_long_range = tmp_mat;
        }
        if(filearray[fileno]==J_mm_long_range_file){
            J_mm_long_range = tmp_mat;
        }

    }

    Mag_field.clear();
    ifstream inputfile(Mag_field_file.c_str());
    if(inputfile.is_open())
    {
        while(getline(inputfile,line)){
            stringstream line_stream(line);
#ifndef WITH_COMPLEX
            line_stream >> tmp;
#endif
#ifdef WITH_COMPLEX
            line_stream >> tmp.real();
            tmp.imag()=0.0;
#endif
            Mag_field.push_back(tmp);
        }

        inputfile.close();
    }


}
