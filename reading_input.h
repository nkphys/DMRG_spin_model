void reading_input(string filename, string & Mag_field_file, string & J_zz_long_range_file,
                   string & J_pm_long_range_file, string & J_mp_long_range_file,
                   string & J_pp_long_range_file, string & J_mm_long_range_file,
                   int & m_infinite, int & m_infinite_max,
                   int & Target_L,int & Target_state, double & Truncation_Error_Target,
                   double & Lanc_Error,  Mat_1_int & Finite_loops,
                   Mat_1_int & Finite_states, Mat_1_int & Finite_states_max,
                   bool & Finite_algo_bool, bool & Wavefuntion_transformation,
                   int & max_lanczos_states, int & N_corr_sets,
                   Mat_1_int & N_point, Mat_2_string & N_point_filenames,
                   Mat_1_string & N_point_range, bool & NPoint_Corrs,
                   Mat_1_int & sweep_no);
void reading_sites_for_correlations(Mat_1_string N_point_range, Mat_3_int & sites_corrs);
void reading_long_range_connections(string Mag_field_file, string J_zz_long_range_file,
                                    string J_pm_long_range_file, string J_mp_long_range_file,
                                    string J_pp_long_range_file, string J_mm_long_range_file,
                                    Mat_1_doub & Mag_field, Mat_2_doub & J_zz_long_range,
                                    Mat_2_doub & J_pm_long_range, Mat_2_doub & J_mp_long_range,
                                    Mat_2_doub & J_pp_long_range, Mat_2_doub & J_mm_long_range,
                                    int Target_L);

void reading_restart_or_saving(bool & _Restart, bool & _SAVING, string & saving_filename, string & restart_filename, string filename);
