

void DMRG::Subtract( Mat_2_doub temp1, double x, Mat_2_doub temp2, Mat_2_doub & temp3){
	
	
	
	if(temp1.size()==temp2.size()){
	for(int k=0;k<temp1.size();k++){
		
		if(temp1[k].size()==temp2[k].size()){
		
		
		for(int i=0;i<temp1[k].size();i++){
			
				temp3[k][i]=temp1[k][i] - x*temp2[k][i];
				
												
			
											}
									
									}
									
		else{cout<<"Problem in DMRG::Inner_Product -- size of vecs in "<<k<<"th Q_eff not equal"<<endl;}
		
		}
								}
	else{cout<<"Problem in DMRG::Inner_Product -- no. of Q_eff sectors not equal in S and E"<<endl;} 
}
	
	


double DMRG::Inner_Product(Mat_2_doub temp1,Mat_2_doub temp2){
	
	
	
	double temp, temp_k;
	temp=0;
	if(temp1.size()==temp2.size()){
	for(int k=0;k<temp1.size();k++){
		
		if(temp1[k].size()==temp2[k].size()){
		temp_k=0;
		
		for(int i=0;i<temp1[k].size();i++){
			
				temp_k=temp_k + temp1[k][i]*temp2[k][i];
				
												
			
											}
										
		temp=temp+temp_k;
									}
									
		else{cout<<"Problem in DMRG::Inner_Product -- size of vecs in "<<k<<"th Q_eff not equal"<<endl;}
		
		}
								}
	else{cout<<"Problem in DMRG::Inner_Product -- no. of Q_eff sectors not equal in S and E"<<endl;}
	
	return temp;
	}


void DMRG::Diagonalize(Mat_1_doub X ,Mat_1_doub Y2 , double & EG, Mat_1_doub & vecG){
	
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
                cout<< "The algorithm failed to compute eigenvalues."<<endl;
                
        }
        
       
       
       EG=eval[0];
       vecG.clear();
       vecG.resize(X.size());
       for(int i=0;i<X.size();i++){vecG[i]=mat[i*X.size()];} 
			
	}


void DMRG::Operate_H_SB(Mat_2_doub vec_in, Mat_2_doub & vec_out, int itr){
	
		//cout<<"here"<<endl;
		
	
	
	
	int lp,l,mp,m,mp_old,lp_old;
	int j;
	int sign;
	
	vec_out.clear();
	vec_out.resize(vec_in.size());
	for(int i=0;i<vec_in.size();i++){
		vec_out[i].assign(vec_in[i].size(),0);
		}
		
		clock_t oprt_SB = clock();
	
	  // ith sec ---> ith sec
	for(int i=0;i<vec_in.size();i++){
		if(H_S[itr][i].nrows!=0 && H_E[itr][i].nrows!=0){
		
	/*if(i==4)	{
		
		cout<<"vec_in[i] = "<<endl;
		for(int j1=0;j1<vec_in[i].size();j1++){cout<<vec_in[i][j1]<<","<<"\t";}cout<<endl;
		cout<<"Q_S = "<<Q_Eff_System[itr][i]<<", Q_E = "<<Q_Eff_Enviroment[itr][i]<<endl;
		cout<<"System :"<<endl;
		Printing_COO_Matrix(H_S[itr][i]);
		cout<<"Enviroment :"<<endl;
		Printing_COO_Matrix(H_E[itr][i]);
		}*/
		

		mp_old=-1;			
		for(mp=0;mp<H_E[itr][i].nrows;mp++){
				
				for(int is=0;is<H_S[itr][i].value.size();is++){
		                lp=H_S[itr][i].rows[is];
						l=H_S[itr][i].columns[is];	
				
				
			    vec_out[i][lp*(H_E[itr][i].nrows)+mp] = vec_out[i][lp*(H_E[itr][i].nrows)+mp] +  vec_in[i][l*(H_E[itr][i].nrows)+mp]*H_S[itr][i].value[is]; 
			    
			    
			}
		}
		
		lp_old=-1;
		for(lp=0;lp<H_S[itr][i].nrows;lp++){		
			
		for(int ie=0;ie<H_E[itr][i].value.size();ie++){
				
				mp=H_E[itr][i].rows[ie];
				m=H_E[itr][i].columns[ie];
				

		       vec_out[i][lp*(H_E[itr][i].nrows)+mp] = vec_out[i][lp*(H_E[itr][i].nrows)+mp] +  vec_in[i][lp*(H_E[itr][i].nrows)+m]*H_E[itr][i].value[ie];
				  
					}
		    }
			    
		    
	/*	for(mp=0;mp<H_E[itr][i].nrows;mp++){
				
				for(int is=0;is<H_S[itr][i].value.size();is++){
	
				
				if(H_S[itr][i].columns[is]>=H_S[itr][i].rows[is]){
					    lp=H_S[itr][i].rows[is];
						l=H_S[itr][i].columns[is];
				if(lp==l){
					vec_out[i][lp*(H_E[itr][i].nrows)+mp] = vec_out[i][lp*(H_E[itr][i].nrows)+mp] +  vec_in[i][l*(H_E[itr][i].nrows)+mp]*H_S[itr][i].value[is];
					}
				else{
			    vec_out[i][lp*(H_E[itr][i].nrows)+mp] = vec_out[i][lp*(H_E[itr][i].nrows)+mp] +  vec_in[i][l*(H_E[itr][i].nrows)+mp]*H_S[itr][i].value[is]; 
			    vec_out[i][l*(H_E[itr][i].nrows)+mp] = vec_out[i][l*(H_E[itr][i].nrows)+mp] +  vec_in[i][lp*(H_E[itr][i].nrows)+mp]*H_S[itr][i].value[is];}
			 }
			    
			}
		}

		
			

		    
		    
		    
		  for(lp=0;lp<H_S[itr][i].nrows;lp++){		
			
		for(int ie=0;ie<H_E[itr][i].value.size();ie++){
				if(H_E[itr][i].columns[ie]>=H_E[itr][i].rows[ie]){
				mp=H_E[itr][i].rows[ie];
				m=H_E[itr][i].columns[ie];
				
				if(mp==m){
		       vec_out[i][lp*(H_E[itr][i].nrows)+mp] = vec_out[i][lp*(H_E[itr][i].nrows)+mp] +  vec_in[i][lp*(H_E[itr][i].nrows)+m]*H_E[itr][i].value[ie];}
		       else{
				   vec_out[i][lp*(H_E[itr][i].nrows)+mp] = vec_out[i][lp*(H_E[itr][i].nrows)+mp] +  vec_in[i][lp*(H_E[itr][i].nrows)+m]*H_E[itr][i].value[ie];
				   vec_out[i][lp*(H_E[itr][i].nrows)+m] = vec_out[i][lp*(H_E[itr][i].nrows)+m] +  vec_in[i][lp*(H_E[itr][i].nrows)+mp]*H_E[itr][i].value[ie];
				   }
				  
					}}
		    }   */		  	
				
				/*if((i==4)){
							
		//cout<<H_S[itr][i].value[is]<<"  "<<" values  "<<	H_E[itr][i].value[ie]<<"   "<<endl;				
		cout<<"vec_out[i] = "<<endl;
		for(int j1=0;j1<vec_out[i].size();j1++){cout<<vec_out[i][j1]<<","<<"\t";}cout<<endl;	
		
		}	*/
							
				}} cout<<"inside operate_H_SB : "<<double( clock() - oprt_SB ) / (double)CLOCKS_PER_SEC<<endl;	
		//	
								
		//C_up_S C_up_dag_E	   ith sector  ---> jth sector						
				for(int i=0;i<vec_in.size();i++){			
		j = Inverse_LB(itr+1,(Q_Eff_System[itr][i]-1000));
	    sign = 1;//(-1)*pow(-1,div(Q_Eff_System[itr][i],1000).quot);
		if(j!=1234567){		if(H_S[itr][i].nrows!=0 && H_E[itr][i].ncols!=0 && H_E[itr][j].nrows!=0 && H_S[itr][j].ncols!=0){				
		for(int is=0;is<C_up_LB[itr+1][i].value.size();is++){
		for(int ie=0;ie<C_up_dag_RB[itr+1][i].value.size();ie++){	
				lp=C_up_LB[itr+1][i].rows[is];
				l=C_up_LB[itr+1][i].columns[is];
				mp=C_up_dag_RB[itr+1][i].rows[ie];
				m=C_up_dag_RB[itr+1][i].columns[ie];
				//cout<<H_E[itr][j].nrows<<endl;	
				vec_out[j][lp*(H_E[itr][j].nrows)+mp] = vec_out[j][lp*(H_E[itr][j].nrows)+mp] +  (-1)*t_hop*sign*vec_in[i][l*(H_E[itr][i].nrows)+m]*C_up_LB[itr+1][i].value[is]*C_up_dag_RB[itr+1][i].value[ie];  
				}
			
		}							
					}}			
							}
									
									
		//C_dn_S C_dn_dag_E	  ith sec  ---> jth sec						
				for(int i=0;i<vec_in.size();i++){					
		j = Inverse_LB(itr+1,(Q_Eff_System[itr][i]-1));
	    sign = 1;//(-1)*pow(-1,(div(Q_Eff_Enviroment[itr][i],1000).quot + div(Q_Eff_System[itr][i],1000).quot + div(Q_Eff_System[itr][i],1000).rem));
		if(j!=1234567){			if(H_S[itr][i].nrows!=0 && H_E[itr][i].ncols!=0 && H_E[itr][j].nrows!=0 && H_S[itr][j].ncols!=0){			
		for(int is=0;is<C_dn_LB[itr+1][i].value.size();is++){
		for(int ie=0;ie<C_dn_dag_RB[itr+1][i].value.size();ie++){
				lp=C_dn_LB[itr+1][i].rows[is];
				l=C_dn_LB[itr+1][i].columns[is];
				mp=C_dn_dag_RB[itr+1][i].rows[ie];
				m=C_dn_dag_RB[itr+1][i].columns[ie];
				vec_out[j][lp*(H_E[itr][j].nrows)+mp] = vec_out[j][lp*(H_E[itr][j].nrows)+mp] +  (-1)*t_hop*sign*vec_in[i][l*(H_E[itr][i].nrows)+m]*C_dn_LB[itr+1][i].value[is]*C_dn_dag_RB[itr+1][i].value[ie];  
				}
			
		}							
					}}
					
				}
					
		//C_up_dag_S C_up_E	   ith sec ---> jth sec						
			for(int i=0;i<vec_in.size();i++){						
		j = Inverse_LB(itr+1,(Q_Eff_System[itr][i]+1000));
	    sign = 1;//pow(-1,div(Q_Eff_System[itr][i],1000).quot);
		if(j!=1234567){		if(H_S[itr][i].nrows!=0 && H_E[itr][i].ncols!=0 && H_E[itr][j].nrows!=0 && H_S[itr][j].ncols!=0){				
		for(int is=0;is<C_up_dag_LB[itr+1][i].value.size();is++){
		for(int ie=0;ie<C_up_RB[itr+1][i].value.size();ie++){
				lp=C_up_dag_LB[itr+1][i].rows[is];
				l=C_up_dag_LB[itr+1][i].columns[is];
				mp=C_up_RB[itr+1][i].rows[ie];
				m=C_up_RB[itr+1][i].columns[ie];
				vec_out[j][lp*(H_E[itr][j].nrows)+mp] = vec_out[j][lp*(H_E[itr][j].nrows)+mp] +  (-1)*t_hop*sign*vec_in[i][l*(H_E[itr][i].nrows)+m]*C_up_dag_LB[itr+1][i].value[is]*C_up_RB[itr+1][i].value[ie];  
				}
			
		}							
				}	}
				}
					
					
					
		//C_dn_dag_S C_dn_E	  ith sec  ---> jth sec						
			for(int i=0;i<vec_in.size();i++){						
		j = Inverse_LB(itr+1,(Q_Eff_System[itr][i]+1));
	    sign = 1;//pow(-1,(div(Q_Eff_Enviroment[itr][i],1000).quot + div(Q_Eff_System[itr][i],1000).quot + div(Q_Eff_System[itr][i],1000).rem));
		if(j!=1234567){		if(H_S[itr][i].nrows!=0 && H_E[itr][i].ncols!=0 && H_E[itr][j].nrows!=0 && H_S[itr][j].ncols!=0){
		for(int is=0;is<C_dn_dag_LB[itr+1][i].value.size();is++){
		for(int ie=0;ie<C_dn_RB[itr+1][i].value.size();ie++){
				lp=C_dn_dag_LB[itr+1][i].rows[is];
				l=C_dn_dag_LB[itr+1][i].columns[is];
				mp=C_dn_RB[itr+1][i].rows[ie];
				m=C_dn_RB[itr+1][i].columns[ie];
				vec_out[j][lp*(H_E[itr][j].nrows)+mp] = vec_out[j][lp*(H_E[itr][j].nrows)+mp] +  (-1)*t_hop*sign*vec_in[i][l*(H_E[itr][i].nrows)+m]*C_dn_dag_LB[itr+1][i].value[is]*C_dn_RB[itr+1][i].value[ie];  
				}
			
		}							
					}}											
									
									
									
									}
		
		
		//
		
		
		
		
	}


void DMRG::Perform_LANCZOS(int iter){
	
	int tmp_sz;
	int lanc_iter=0;
	double eps=0.000000000001,diff_E;
	double temp1, temp2, temp3, E0, E0_old;
	Mat_1_doub B2,A, red_eig_vec,Norms;
	Mat_2_doub Kvector_n,Kvector_nm1,Kvector_np1 ; //[q_S_index][element] element = i*(Dim(E)) + j(i~Sys, j~Env) 
	
	//Mat_2_doub temp_Wvec_1, temp_Wvec_2; //[q_S_index][element]
	
	
	Kvector_n.clear();	Kvector_nm1.clear(); Kvector_np1.clear();
	Kvector_n.resize(Q_Eff_System[iter].size());	Kvector_nm1.resize(Q_Eff_System[iter].size());	Kvector_np1.resize(Q_Eff_System[iter].size());
    
	
	//temp_Wvec_S.resize(Q_Eff_System[iter].size());
	//temp_Wvec_E.resize(Q_Eff_Enviroment[iter].size());
	    srand(1);
	    tmp_sz=0;
	for(int i=0;i<Q_Eff_System[iter].size();i++){
		//temp_Wvec_S[i].resize(H_LB[iter+1][i]);
		//cout<<H_S[iter][i].nrows<<"   "<<H_E[iter][i].nrows<<endl;
			tmp_sz=tmp_sz + H_S[iter][i].nrows*H_E[iter][i].nrows;
		for(int j=0;j<H_S[iter][i].nrows*H_E[iter][i].nrows;j++){
			temp1=(rand()%RAND_MAX);
			Kvector_n[i].push_back(temp1/RAND_MAX);

			}}
	
   
	
	E0_old=0;
	diff_E=1.0;
	cout<<"NO. of Q_EFF sectors required to create SB  = "<<Q_Eff_System[iter].size()<<endl<<endl;
	cout<<"LANCZOS(pass 1) STARTING FOR SUPERBLOCK for ITERATION = "<<iter<<", Size of Matrix(SB) = "<<tmp_sz<<endl;
	
	
	while(diff_E>eps){
	clock_t Lanc_time = clock();
	temp1 =Inner_Product(Kvector_n,Kvector_n);	
	Norms.push_back(sqrt(temp1));
	if(lanc_iter==0){B2.push_back(0);}
	else{
	temp2 = Inner_Product(Kvector_nm1,Kvector_nm1);
	B2.push_back(temp1/temp2);}
		
	/*	if(lanc_iter==0){
		for(int i=0;i<Q_Eff_System[iter].size();i++){

		for(int j=0;j<H_S[iter][i].nrows*H_E[iter][i].nrows;j++){
		
			cout<<Kvector_n[i][j]<<"  ";

			}cout<<endl;}}*/
	clock_t oprt_SB_time = clock();	//cout<<"here"<<endl;
		
	Operate_H_SB(Kvector_n,Kvector_np1,iter);// saved in K_vector_np1
	cout<<"Time to operate SB : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;
	/*		if(lanc_iter==0){
		for(int i=0;i<Q_Eff_System[iter].size();i++){

		for(int j=0;j<H_S[iter][i].nrows*H_E[iter][i].nrows;j++){
		
			cout<<Kvector_np1[i][j]<<"  ";

			}cout<<endl;}}*/
	
	temp3 =Inner_Product(Kvector_n, Kvector_np1);
	
	
	A.push_back(temp3/temp1);
	//cout<<temp3<<"  "<<temp1<<"  "<<A[lanc_iter]<<endl;
	
	Subtract(Kvector_np1, A[lanc_iter], Kvector_n, Kvector_np1);	//
	if(lanc_iter!=0){Subtract(Kvector_np1, B2[lanc_iter], Kvector_nm1, Kvector_np1);	}	
		

		
	Diagonalize(A,B2,E0,red_eig_vec);//
	
	diff_E=	fabs(E0-E0_old);
	//cout<<"Energy for lanc_iter("<<lanc_iter<<") is "<<E0<<endl;
	
	E0_old=E0;
	
	Kvector_nm1=Kvector_n;
	Kvector_n=Kvector_np1;	
	lanc_iter=lanc_iter+1;
	
	if(A.size()==tmp_sz){diff_E=0;}
	cout<<"Time for 1 LAnczos iter : "<<double( clock() - Lanc_time ) / (double)CLOCKS_PER_SEC<<endl;	
		}
		
	cout<<"NO. of itearations required to get convergence in LANCZOS(pass 1) = "<<lanc_iter<<endl;
	cout<<"Energy(GS of SB) = "<<scientific<< E0<<"   "<<eps<<"  "<<diff_E<<endl;	
	Energy=E0;Lanc_Error=diff_E;
	
	cout<<"LANCZOS(pass 2) STARTING FOR SUPERBLOCK Eigenvector = "<<iter<<", Size of Matrix(SB) = "<<tmp_sz<<endl;
	
	
	Kvector_n.clear();	Kvector_nm1.clear(); Kvector_np1.clear();
	Kvector_n.resize(Q_Eff_System[iter].size());	Kvector_nm1.resize(Q_Eff_System[iter].size());	Kvector_np1.resize(Q_Eff_System[iter].size());
    
	    srand(1);
	    
	for(int i=0;i<Q_Eff_System[iter].size();i++){
			
		for(int j=0;j<H_S[iter][i].nrows*H_E[iter][i].nrows;j++){
			temp1=(rand()%RAND_MAX);
			Kvector_n[i].push_back(temp1/RAND_MAX);

			}}
			
						
	
      Eig_vec.clear();
      Eig_vec.resize((Q_Eff_System[iter].size()));
      
      for(int i=0;i<Q_Eff_System[iter].size();i++){
			
		for(int j=0;j<H_S[iter][i].nrows*H_E[iter][i].nrows;j++){
			temp1=0;
		    Eig_vec[i].push_back(temp1);

			}}
	
		for(int lanc_iter2=0;lanc_iter2<lanc_iter;lanc_iter2=lanc_iter2+1){

			
	
	Subtract(Eig_vec, (-1.0*(red_eig_vec[lanc_iter2]))/Norms[lanc_iter2], Kvector_n, Eig_vec);
	
			

			
	Operate_H_SB(Kvector_n,Kvector_np1,iter);// saved in K_vector_np1
   	
	
	
	
	Subtract(Kvector_np1, A[lanc_iter2], Kvector_n, Kvector_np1);	//
	if(lanc_iter2!=0){Subtract(Kvector_np1, B2[lanc_iter2], Kvector_nm1, Kvector_np1);	}	
		
   
		
	
	Kvector_nm1=Kvector_n;
	Kvector_n=Kvector_np1;	
	
	
	
		}
				

			
			
			
			
			
			
		
	}


void DMRG::Choosing_m_states(Mat_2_doub  Eval, Mat_1_int & m_sts, int iter){
	
	m_sts.clear();
	m_sts.resize(Q_Eff_System[iter].size());
	for(int i=0;i<Q_Eff_System[iter].size();i++){m_sts[i]=0;}
	
	int m=0;
	double max;
	while(m<m_infinite){
	max=-1;
	
	for(int i=0;i<Eval.size();i++){if(H_S[iter][i].nrows!=0){
		for(int l=0;l<Eval[i].size();l++){
			
			if(Eval[i][l]>max){max=Eval[i][l];}
			
			}}}
			
	for(int i=0;i<Eval.size();i++){if(H_S[iter][i].nrows!=0){
	for(int l=0;l<Eval[i].size();l++){
			
		if(Eval[i][l]==max){Eval[i][l]=0;m_sts[i]=m_sts[i]+1;m=m+1;}
			
			}}}		
				
	 
	
    
	
	}
	
	
	
	}
	
	
void DMRG::Do_RENORMALIZATION_of_S_and_E(int iter){
	
	   double** Red_den_mat=(double **) malloc(sizeof(double *)*Q_Eff_System[iter].size());
	   
	   double** eval=(double **) malloc(sizeof(double *)*Q_Eff_System[iter].size());
	   Mat_2_doub Evl;
	   Evl.resize(Q_Eff_System[iter].size());
	   int j;
	   int LDA,info;
	   Mat_1_int m_states;
	
	for(int i=0;i<Q_Eff_System[iter].size();i++){ if(H_S[iter][i].nrows!=0){
		
		LDA = H_S[iter][i].nrows;
		Red_den_mat[i]= (double*) calloc (H_S[iter][i].nrows*H_S[iter][i].ncols ,sizeof(double));
	    eval[i] =  (double*) calloc(H_S[iter][i].nrows, sizeof(double));
	    Evl[i].resize(H_S[iter][i].nrows);
		
		
		
		for(int m=0;m<H_S[iter][i].nrows;m++){
			for(int n=0;n<=m;n++){
                  for(int l=0;l<H_E[iter][i].nrows;l++){
					  Red_den_mat[i][m*H_S[iter][i].ncols + n]=Red_den_mat[i][m*H_S[iter][i].ncols + n] + (Eig_vec[i][m*(H_E[iter][i].nrows) + l])*(Eig_vec[i][n*(H_E[iter][i].nrows) + l]);
					  }				
				}
			}
		
		
		
		
		info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L', H_S[iter][i].nrows, Red_den_mat[i] , LDA, eval[i] );
       
     
        /* Check for convergence */
        if( info > 0 ) {
                cout<< "The algorithm failed to diagonalize Red_den_mat_S."<<endl;
                
        }
		
		 }}   
		 
		 for(int i=0;i<Q_Eff_System[iter].size();i++){if(H_S[iter][i].nrows!=0){
			 for(int l=0; l<H_S[iter][i].nrows;l++){
				 Evl[i][l]=eval[i][l];
				 }
			 }} 
		 
		 //cout<<"fine 1"<<endl;
		 Choosing_m_states(Evl,m_states,iter);//cout<<"fine till here 2"<<endl;
		 Truncation_Error_S=1;
		 for(int i=0;i<Q_Eff_System[iter].size();i++){if(H_S[iter][i].nrows!=0){
			 for(int l=H_S[iter][i].nrows-1; l>H_S[iter][i].nrows-1-m_states[i];l--){
				 Truncation_Error_S=Truncation_Error_S - Evl[i][l];
				 }
			 }} 
		 //cout<<Q_Eff_System[iter].size()<<"   "<<m_states.size()<<endl;
		 for(int i=0;i<Q_Eff_System[iter].size();i++){if(H_S[iter][i].nrows!=0){
		
			 
		Renormalize(H_S[iter][i], Red_den_mat[i], Red_den_mat[i], H_LB[iter+1][i], m_states[i], m_states[i]);
		j=Inverse_LB(iter+1,Q_Eff_System[iter][i]-1000);
		if(j!=1234567){if(H_S[iter][j].nrows!=0){
		Renormalize(C_up_LB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_up_LB[iter+1][i],m_states[j], m_states[i] );}}
		j=Inverse_LB(iter+1,Q_Eff_System[iter][i]-1);
		if(j!=1234567){if(H_S[iter][j].nrows!=0){
		Renormalize(C_dn_LB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_dn_LB[iter+1][i], m_states[j], m_states[i]);}}
		j=Inverse_LB(iter+1,Q_Eff_System[iter][i]+1000);
		if(j!=1234567){if(H_S[iter][j].nrows!=0){
		Renormalize(C_up_dag_LB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_up_dag_LB[iter+1][i] , m_states[j], m_states[i]);}}
	    j=Inverse_LB(iter+1,Q_Eff_System[iter][i]+1);
		if(j!=1234567){if(H_S[iter][j].nrows!=0){
		Renormalize(C_dn_dag_LB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_dn_dag_LB[iter+1][i] , m_states[j], m_states[i]);}}
			 
			 
			} }
		 
		 for(int i=0;i<m_states.size();i++){if(m_states[i]==0){Q_Eff_Left_Block[iter+1][i]=Q_Eff_Super_Block[max_iter]+100000;}}     
		
		
		
		
		for(int i=0;i<Q_Eff_System[iter].size();i++){
		
		LDA = H_E[iter][i].nrows;
		Red_den_mat[i]= (double*) calloc(H_E[iter][i].nrows*H_E[iter][i].ncols, sizeof(double));
	    eval[i] = (double*) calloc(H_E[iter][i].nrows, sizeof(double));
		 Evl[i].resize(H_E[iter][i].nrows);
		
		
		for(int m=0;m<H_E[iter][i].nrows;m++){
			for(int n=0;n<=m;n++){
                  for(int l=0;l<H_S[iter][i].nrows;l++){
					  Red_den_mat[i][m*H_E[iter][i].ncols + n]=Red_den_mat[i][m*H_E[iter][i].ncols + n] + (Eig_vec[i][l*(H_E[iter][i].nrows) + m])*(Eig_vec[i][l*(H_E[iter][i].nrows) + n]);
					  }				
				}
			}
		
		
		
		
		info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L', H_E[iter][i].nrows, Red_den_mat[i] , LDA, eval[i] );
       
     
        /* Check for convergence */
        if( info > 0 ) {
                cout<< "The algorithm failed to diagonalize Red_den_mat_E."<<endl;
                
        }
		
		
				}
		
				for(int i=0;i<Q_Eff_Enviroment[iter].size();i++){
			 for(int l=0; l<H_E[iter][i].nrows;l++){
				 Evl[i][l]=eval[i][l];
				 }
			 }
		 
		 
		 Choosing_m_states(Evl,m_states,iter);
		 
		 
		 		 Truncation_Error_E=1;
		 		 
		 		for(int i=0;i<Q_Eff_Enviroment[iter].size();i++){if(H_E[iter][i].nrows!=0){
			 for(int l=H_E[iter][i].nrows-1; l>H_E[iter][i].nrows-1-m_states[i];l--){
				 Truncation_Error_E=Truncation_Error_E - Evl[i][l];
				 }
			} }
		 		 
		
			 for(int i=0;i<Q_Eff_System[iter].size();i++){
			 
		Renormalize(H_E[iter][i], Red_den_mat[i], Red_den_mat[i], H_RB[iter+1][i],m_states[i], m_states[i]);
		j=Inverse_RB(iter+1,Q_Eff_Enviroment[iter][i]-1000);
		if(j!=1234567){
		Renormalize(C_up_RB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_up_RB[iter+1][i],m_states[j], m_states[i]);}
		j=Inverse_RB(iter+1,Q_Eff_Enviroment[iter][i]-1);
		if(j!=1234567){
		Renormalize(C_dn_RB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_dn_RB[iter+1][i],m_states[j], m_states[i]);}
		j=Inverse_RB(iter+1,Q_Eff_Enviroment[iter][i]+1000);
		if(j!=1234567){
		Renormalize(C_up_dag_RB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_up_dag_RB[iter+1][i],m_states[j], m_states[i]);}
	    j=Inverse_RB(iter+1,Q_Eff_Enviroment[iter][i]+1);
		if(j!=1234567){
		Renormalize(C_dn_dag_RB[iter+1][i], Red_den_mat[j], Red_den_mat[i], C_dn_dag_RB[iter+1][i],m_states[j], m_states[i]);}
			 
			 
			 }
		
		for(int i=0;i<m_states.size();i++){if(m_states[i]==0){Q_Eff_Right_Block[iter+1][i]=Q_Eff_Super_Block[max_iter]+100000;}}
		
		
		
	
	
	
	}


void DMRG::Renormalize(Matrix_COO A, double* UL,double* UR, Matrix_COO & B, int m_UL, int m_UR){
	
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
			for( int p=0;p<A.value.size();p++){
				temp=temp+UL[A.rows[p]*A.nrows + A.nrows - 1 - i]*A.value[p]*UR[A.columns[p]*A.ncols + A.ncols - 1 - l];
				
				}
			if(fabs(temp)>0.0000000001){B.value.push_back(temp);B.rows.push_back(i);B.columns.push_back(l);tmp=tmp+1;}	
			}
		}
	
    if(tmp==0){B.value.push_back(0);B.rows.push_back(0);B.columns.push_back(0);}
	//if(m_inf_r==0 || m_inf_c==0){cout<<"some problem"<<endl;}
	}	
	
	
void DMRG::Grow_LB_and_update_C_LB_oprts(int iter){
	
	
int i_pr, k_prime;
int Q_i_prime, Q_j_prime, Q_k_prime;
int sign;
int tmp;
Matrix_COO temp_COO;
Mat_2_int cols_before;
Mat_2_int rows_before;
	
	
	tmp =0;
	cols_before.clear();
	rows_before.clear();
	cols_before.resize(Q_Eff_System[iter].size());
	rows_before.resize(Q_Eff_System[iter].size());
	
		//for(int k=0;k<Q_Eff_System[iter].size();k++){H_S[iter][k].nrows=1;H_S[iter][k].ncols=1;H_S[iter][k].value.assign(1,0);H_S[iter][k].rows.assign(1,0);H_S[iter][k].columns.assign(1,0);}
	
	for(int k=0;k<Q_Eff_System[iter].size();k++){

	cols_before[k].clear();
	rows_before[k].clear();	
	cols_before[k].resize(Q_Eff_Left_Block[iter].size());
	rows_before[k].resize(Q_Eff_Left_Block[iter].size());
		
	for(int i=0;i<Q_Eff_Left_Block[iter].size();i++){
		for(int j=0;j<Q_Eff_Site.size();j++){
								if(Q_Eff_Left_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_System[iter][k]){
									//cout<<Q_Eff_Left_Block[iter][i]<<"\t"<<Q_Eff_Site[j] <<"\t"<< Q_Eff_System[iter][k]<<endl;
									
										cols_before[k][i]=H_S[iter][k].ncols ;
										rows_before[k][i]=H_S[iter][k].nrows ;
										 
									if(tmp==0){
										 
										Direct_Product(Identity(H_LB[iter][i].nrows),H_LB[0][j], H_S[iter][k]);
										
										Sum(H_LB[iter][i],H_S[iter][k],H_S[iter][k]);
										
										tmp=tmp+1;
										}
										else{
										Sum(  H_LB[iter][i]     ,   Direct_Product(Identity(H_LB[iter][i].nrows),H_LB[0][j]), temp_COO );
										Direct_Sum(H_S[iter][k],temp_COO,H_S[iter][k]);}
																										
															
																									
																										}
										     }
													 }
													
										
							//---------------------ADDING FC-----------------------//
							
			for(int i=0;i<Q_Eff_Left_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){
							if(Q_Eff_Left_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_System[iter][k]){
							
								//C^_up(LB)C_up(site)
								Direct_Product(C_up_dag_LB[iter][i],C_up_LB[0][j],temp_COO);

								Q_i_prime=Q_Eff_Left_Block[iter][i] + 1000;  
								Q_j_prime=Q_Eff_Site[j]-1000;
								
							    i_pr = Inverse_LB(iter,Q_i_prime);
							    sign = 1;//pow(-1,div(Q_Eff_Left_Block[iter][i],1000).quot);
							    if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_S[iter][k], sign, t_hop);}
																																	
								
								
								//C^_dn(LB)C_dn(site)
								Direct_Product(C_dn_dag_LB[iter][i],C_dn_LB[0][j],temp_COO);
								
								Q_i_prime=Q_Eff_Left_Block[iter][i] + 1;  
								Q_j_prime=Q_Eff_Site[j]-1;
							
								i_pr = Inverse_LB(iter,Q_i_prime);
								sign = 1;//pow(-1,(div(Q_Eff_Left_Block[iter][i],1000).quot  + div(Q_Eff_Left_Block[iter][i],1000).rem  + div(Q_Eff_Site[j],1000).quot));
								//Printing_COO_Matrix(temp_COO);
								if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_S[iter][k], sign, t_hop);}
																																	
								
								//C^_dn(site)C_dn(LB)
								Direct_Product(C_dn_LB[iter][i],C_dn_dag_LB[0][j],temp_COO);
								
								Q_i_prime=Q_Eff_Left_Block[iter][i] - 1;  
								Q_j_prime=Q_Eff_Site[j]+1;
							
								i_pr = Inverse_LB(iter,Q_i_prime);
								sign = 1;// -1*pow(-1,(div(Q_Eff_Left_Block[iter][i],1000).quot  + div(Q_Eff_Left_Block[iter][i],1000).rem  + div(Q_Eff_Site[j],1000).quot));
								//Printing_COO_Matrix(temp_COO);
								if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_S[iter][k], sign, t_hop);}
								
								
								//C^_up(site)C_up(LB)
								Direct_Product(C_up_LB[iter][i],C_up_dag_LB[0][j],temp_COO);
								
								Q_i_prime=Q_Eff_Left_Block[iter][i] - 1000;  
								Q_j_prime=Q_Eff_Site[j]+1000;
							
								i_pr = Inverse_LB(iter,Q_i_prime);
								sign = 1;//-1*pow(-1,div(Q_Eff_Left_Block[iter][i],1000).quot);
								//Printing_COO_Matrix(temp_COO);
								if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_S[iter][k], sign, t_hop);} //
								
								
								}
								}
								}
							
							
							
							
							
							
							//--------------------FC DONE------------------------//			
										
								tmp=0;			
													}
	
	
	
	
	 //__________________________________________________________________________________________________________________________//	
	//-----------------------------------------CREATING ANHILATION/CREATION OPRS for LB/System----------------------------------//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
		
	C_up_LB[iter+1].resize(Q_Eff_System[iter].size());
	C_dn_LB[iter+1].resize(Q_Eff_System[iter].size());
	C_up_dag_LB[iter+1].resize(Q_Eff_System[iter].size());
	C_dn_dag_LB[iter+1].resize(Q_Eff_System[iter].size());
	for(int k=0;k<Q_Eff_System[iter].size();k++){
		
		//C_up_LB
		Q_k_prime = Q_Eff_System[iter][k]-1000;
		k_prime = Inverse_LB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_up_LB[iter+1][k].ncols=H_S[iter][k].ncols;
		C_up_LB[iter+1][k].nrows=H_S[iter][k_prime].ncols; 
		C_up_LB[iter+1][k].value.push_back(0);C_up_LB[iter+1][k].rows.push_back(0);C_up_LB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Left_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){							
							if(Q_Eff_Left_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_System[iter][k]){
							if((Q_Eff_Site[j]==1000) || (Q_Eff_Site[j]==1001))	{
							sign = 1;//pow(-1,div(Q_Eff_Left_Block[iter][i],1000).quot);
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_LB[iter][i].nrows), C_up_LB[iter+1][k], sign*C_up_LB[0][j].value[0], (-1)*C_up_LB[0][j].value[0]);
						}
							}}} 
		
							}
							
		//C_dn_LB
		Q_k_prime = Q_Eff_System[iter][k]-1;
		k_prime = Inverse_LB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_dn_LB[iter+1][k].ncols=H_S[iter][k].ncols;
		C_dn_LB[iter+1][k].nrows=H_S[iter][k_prime].ncols; 
		C_dn_LB[iter+1][k].value.push_back(0);C_dn_LB[iter+1][k].rows.push_back(0);C_dn_LB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Left_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){
							if(Q_Eff_Left_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_System[iter][k]){
							if((Q_Eff_Site[j]==1) || (Q_Eff_Site[j]==1001))	{
							sign = 1;//pow(-1,(div(Q_Eff_Left_Block[iter][i],1000).quot + div(Q_Eff_Left_Block[iter][i],1000).rem));
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_LB[iter][i].nrows), C_dn_LB[iter+1][k], sign*C_dn_LB[0][j].value[0], (-1)*C_dn_LB[0][j].value[0]);
						}
							}}} 
		
							}
							
							
							
		//C_dn_dag_LB
		Q_k_prime = Q_Eff_System[iter][k]+1;
		k_prime = Inverse_LB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_dn_dag_LB[iter+1][k].ncols=H_S[iter][k].ncols;
		C_dn_dag_LB[iter+1][k].nrows=H_S[iter][k_prime].ncols; 
		C_dn_dag_LB[iter+1][k].value.push_back(0);C_dn_dag_LB[iter+1][k].rows.push_back(0);C_dn_dag_LB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Left_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){
							if(Q_Eff_Left_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_System[iter][k]){
							if((Q_Eff_Site[j]==0) || (Q_Eff_Site[j]==1000))	{
							sign = 1;//pow(-1,(div(Q_Eff_Left_Block[iter][i],1000).quot + div(Q_Eff_Left_Block[iter][i],1000).rem ));
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_LB[iter][i].nrows), C_dn_dag_LB[iter+1][k], sign*C_dn_dag_LB[0][j].value[0], (-1)*C_dn_dag_LB[0][j].value[0]);
						}
							}}} 
		
							}
							
		//C_up_dag_LB
		Q_k_prime = Q_Eff_System[iter][k]+1000;
		k_prime = Inverse_LB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_up_dag_LB[iter+1][k].ncols=H_S[iter][k].ncols;
		C_up_dag_LB[iter+1][k].nrows=H_S[iter][k_prime].ncols; 
		C_up_dag_LB[iter+1][k].value.push_back(0);C_up_dag_LB[iter+1][k].rows.push_back(0);C_up_dag_LB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Left_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){							
							if(Q_Eff_Left_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_System[iter][k]){
							if((Q_Eff_Site[j]==0) || (Q_Eff_Site[j]==1))	{
							sign = 1;//pow(-1,div(Q_Eff_Left_Block[iter][i],1000).quot);
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_LB[iter][i].nrows), C_up_dag_LB[iter+1][k], sign*C_up_dag_LB[0][j].value[0], (-1)*C_up_dag_LB[0][j].value[0]);
						}
							}}} 
		
							}
		
		}
	
	
	
	 //__________________________________________________________________________________________________________________________//	
	//-----------------------------------------CREATED ANHILATION/CREATION OPRS for LB/System----------------------------------//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
	
	
	
	}


void DMRG::Grow_RB_and_update_C_RB_oprts(int iter){
	int i_pr, k_prime;
int Q_i_prime, Q_j_prime, Q_k_prime;
int sign;
int tmp;
Matrix_COO temp_COO;
Mat_2_int cols_before;
Mat_2_int rows_before;
	
	
		tmp =0;
	cols_before.clear();
	rows_before.clear();
	cols_before.resize(Q_Eff_Enviroment[iter].size());
	rows_before.resize(Q_Eff_Enviroment[iter].size());
	//for(int k=0;k<Q_Eff_Enviroment[iter].size();k++){H_E[iter][k].nrows=1;H_E[iter][k].ncols=1;H_E[iter][k].value.assign(1,0);H_E[iter][k].rows.assign(1,0);H_E[iter][k].columns.assign(1,0);}
	
	for(int k=0;k<Q_Eff_Enviroment[iter].size();k++){

	cols_before[k].clear();
	rows_before[k].clear();	
	cols_before[k].resize(Q_Eff_Right_Block[iter].size());
	rows_before[k].resize(Q_Eff_Right_Block[iter].size());
		
	for(int i=0;i<Q_Eff_Right_Block[iter].size();i++){
		for(int j=0;j<Q_Eff_Site.size();j++){
								if(Q_Eff_Right_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_Enviroment[iter][k]){
									//cout<<Q_Eff_Left_Block[iter][i]<<"\t"<<Q_Eff_Site[j] <<"\t"<< Q_Eff_System[iter][k]<<endl;
										cols_before[k][i]=H_E[iter][k].ncols ;
										rows_before[k][i]=H_E[iter][k].nrows ;
									
									if(tmp==0){
										 
										Direct_Product(H_LB[0][j],Identity(H_RB[iter][i].nrows), H_E[iter][k]);
										
										Sum(H_RB[iter][i],H_E[iter][k],H_E[iter][k]);
										
										tmp=tmp+1;
										}
										else{
										Sum(  H_RB[iter][i]     ,   Direct_Product(H_LB[0][j],Identity(H_RB[iter][i].nrows)), temp_COO );
										Direct_Sum(H_E[iter][k],temp_COO,H_E[iter][k]);}
																										
															
																									
																										}
										     }
													 }
													
										
							//---------------------ADDING FC-----------------------//
							
			for(int i=0;i<Q_Eff_Right_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){
							if(Q_Eff_Right_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_Enviroment[iter][k]){
							
								//C^_up(site)C_up(RB)
								Direct_Product(C_up_dag_LB[0][j],C_up_RB[iter][i],temp_COO);

								Q_j_prime=Q_Eff_Site[j] + 1000;  
								Q_i_prime=Q_Eff_Right_Block[iter][i]-1000;
								
							    i_pr = Inverse_RB(iter,Q_i_prime);
							    sign = 1;//pow(-1,div(Q_Eff_Site[j],1000).quot);
							    if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_E[iter][k], sign, t_hop);}
																																	
								
								
								//C^_dn(site)C_dn(RB)
								Direct_Product(C_dn_dag_LB[0][j],C_dn_RB[iter][i],temp_COO);
								
								Q_i_prime=Q_Eff_Right_Block[iter][i] - 1;  
								Q_j_prime=Q_Eff_Site[j]+1;
							
								i_pr = Inverse_RB(iter,Q_i_prime);
								sign = 1;//pow(-1,(div(Q_Eff_Site[j],1000).quot  + div(Q_Eff_Site[j],1000).rem  + div(Q_Eff_Right_Block[iter][i],1000).quot));
								//Printing_COO_Matrix(temp_COO);
								if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_E[iter][k], sign, t_hop);}
																																	
								
								//C^_dn(RB)C_dn(site)
								Direct_Product(C_dn_LB[0][j],C_dn_dag_RB[iter][i],temp_COO);
								
								Q_i_prime=Q_Eff_Right_Block[iter][i] + 1;  
								Q_j_prime=Q_Eff_Site[j]-1;
							
								i_pr = Inverse_RB(iter,Q_i_prime);
								sign = 1;//-1*pow(-1,(div(Q_Eff_Site[j],1000).quot  + div(Q_Eff_Site[j],1000).rem  + div(Q_Eff_Right_Block[iter][i],1000).quot));
								//Printing_COO_Matrix(temp_COO);
								if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_E[iter][k], sign, t_hop);}
								
								
								//C^_up(RB)C_up(site)
								Direct_Product(C_up_LB[0][j],C_up_dag_RB[iter][i],temp_COO);
								
								Q_i_prime=Q_Eff_Right_Block[iter][i] + 1000;  
								Q_j_prime=Q_Eff_Site[j]-1000;
							
								i_pr = Inverse_RB(iter,Q_i_prime);
								sign = 1;//-1*pow(-1,div(Q_Eff_Site[j],1000).quot);
								//Printing_COO_Matrix(temp_COO);
								if(i_pr!=1234567){
								Add_FC(cols_before[k][i], rows_before[k][i_pr], temp_COO, H_E[iter][k], sign, t_hop);} //
								
								
								
								
								
								}
								}
								}
							
							
							
							
							
							
							//--------------------FC DONE------------------------//			
										
								tmp=0;			
													}
															
    //__________________________________________________________________________________________________________________________//	
	//-----------------------------------------CREATING ANHILATION/CREATION OPRS for RB/Enviroment------------------------------//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
		
	C_up_RB[iter+1].resize(Q_Eff_Enviroment[iter].size());
	C_dn_RB[iter+1].resize(Q_Eff_Enviroment[iter].size());
	C_up_dag_RB[iter+1].resize(Q_Eff_Enviroment[iter].size());
	C_dn_dag_RB[iter+1].resize(Q_Eff_Enviroment[iter].size());
	for(int k=0;k<Q_Eff_Enviroment[iter].size();k++){
		
		//C_up_RB
		Q_k_prime = Q_Eff_Enviroment[iter][k]-1000;
		k_prime = Inverse_RB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_up_RB[iter+1][k].ncols=H_E[iter][k].ncols;
		C_up_RB[iter+1][k].nrows=H_E[iter][k_prime].ncols; 
		C_up_RB[iter+1][k].value.push_back(0);C_up_RB[iter+1][k].rows.push_back(0);C_up_RB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Right_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){							
							if(Q_Eff_Right_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_Enviroment[iter][k]){
							if((Q_Eff_Site[j]==1000) || (Q_Eff_Site[j]==1001))	{
							sign = 1.0;
							//cout<<cols_before[k][i]<<"   "<<rows_before[k_prime][i]<<endl;
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_RB[iter][i].nrows), C_up_RB[iter+1][k], C_up_LB[0][j].value[0], (-1)*C_up_LB[0][j].value[0]);
						}
							}}} 
		
							}
							
		//C_dn_RB
		Q_k_prime = Q_Eff_Enviroment[iter][k]-1;
		k_prime = Inverse_RB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_dn_RB[iter+1][k].ncols=H_E[iter][k].ncols;
		C_dn_RB[iter+1][k].nrows=H_E[iter][k_prime].ncols; 
		C_dn_RB[iter+1][k].value.push_back(0);C_dn_RB[iter+1][k].rows.push_back(0);C_dn_RB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Right_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){
							if(Q_Eff_Right_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_Enviroment[iter][k]){
							if((Q_Eff_Site[j]==1) || (Q_Eff_Site[j]==1001))	{
							sign = 1;//pow(-1,(div(Q_Eff_Right_Block[iter][i],1000).quot));
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_RB[iter][i].nrows), C_dn_RB[iter+1][k], sign*C_dn_LB[0][j].value[0], (-1)*C_dn_LB[0][j].value[0]);
						}
							}}} 
		
							}
							
							
							
		//C_dn_dag_RB
		Q_k_prime = Q_Eff_Enviroment[iter][k]+1;
		k_prime = Inverse_RB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_dn_dag_RB[iter+1][k].ncols=H_E[iter][k].ncols;
		C_dn_dag_RB[iter+1][k].nrows=H_E[iter][k_prime].ncols; 
		C_dn_dag_RB[iter+1][k].value.push_back(0);C_dn_dag_RB[iter+1][k].rows.push_back(0);C_dn_dag_RB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Right_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){
							if(Q_Eff_Right_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_Enviroment[iter][k]){
							if((Q_Eff_Site[j]==0) || (Q_Eff_Site[j]==1000))	{
							sign = 1;//pow(-1,(div(Q_Eff_Right_Block[iter][i],1000).quot));
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_RB[iter][i].nrows), C_dn_dag_RB[iter+1][k], sign*C_dn_dag_LB[0][j].value[0], (-1)*C_dn_dag_LB[0][j].value[0]);
						}
							}}} 
		
							}
							
		//C_up_dag_RB
		Q_k_prime = Q_Eff_Enviroment[iter][k]+1000;
		k_prime = Inverse_RB(iter+1,Q_k_prime);
		if(k_prime!=1234567){	
		C_up_dag_RB[iter+1][k].ncols=H_E[iter][k].ncols;
		C_up_dag_RB[iter+1][k].nrows=H_E[iter][k_prime].ncols; 
		C_up_dag_RB[iter+1][k].value.push_back(0);C_up_dag_RB[iter+1][k].rows.push_back(0);C_up_dag_RB[iter+1][k].columns.push_back(0);
		 
		for(int i=0;i<Q_Eff_Right_Block[iter].size();i++){
						for(int j=0;j<Q_Eff_Site.size();j++){							
							if(Q_Eff_Right_Block[iter][i] + Q_Eff_Site[j] == Q_Eff_Enviroment[iter][k]){
							if((Q_Eff_Site[j]==0) || (Q_Eff_Site[j]==1))	{
							sign = 1;
							
							Add_FC(cols_before[k][i], rows_before[k_prime][i], Identity(H_RB[iter][i].nrows), C_up_dag_RB[iter+1][k], C_up_dag_LB[0][j].value[0], (-1)*C_up_dag_LB[0][j].value[0]);
						}
							}}} 
		
							}
		
		}
	
	
	
	//__________________________________________________________________________________________________________________________//	
	//-----------------------------------------CREATED ANHILATION/CREATION OPRS for RB/Enviroment-------------------------------//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//		
	
	
	}


void DMRG::Initialize_Hamiltonians(){
	
	//HERE

H_LB.resize(max_iter+1);H_RB.resize(max_iter+1);
H_S.resize(max_iter+1);H_E.resize(max_iter+1);


for(int i=0;i<=max_iter;i++){
	H_LB[i].resize(Q_Eff_Left_Block[i].size());
	H_RB[i].resize(Q_Eff_Right_Block[i].size());
	H_S[i].resize(Q_Eff_System[i].size());
	H_E[i].resize(Q_Eff_Enviroment[i].size());
}
H_LB[0][0].value.push_back(0);H_LB[0][0].rows.push_back(0);H_LB[0][0].columns.push_back(0);H_LB[0][0].nrows=1;H_LB[0][0].ncols=1;
H_LB[0][1].value.push_back(0);H_LB[0][1].rows.push_back(0);H_LB[0][1].columns.push_back(0);H_LB[0][1].nrows=1;H_LB[0][1].ncols=1;
H_LB[0][2].value.push_back(0);H_LB[0][2].rows.push_back(0);H_LB[0][2].columns.push_back(0);H_LB[0][2].nrows=1;H_LB[0][2].ncols=1;
H_LB[0][3].value.push_back(U);H_LB[0][3].rows.push_back(0);H_LB[0][3].columns.push_back(0);H_LB[0][3].nrows=1;H_LB[0][3].ncols=1;

//H_site = H_LB[0];

H_RB[0][0].value.push_back(U);H_RB[0][0].rows.push_back(0);H_RB[0][0].columns.push_back(0);H_RB[0][0].nrows=1;H_RB[0][0].ncols=1;
H_RB[0][1].value.push_back(0);H_RB[0][1].rows.push_back(0);H_RB[0][1].columns.push_back(0);H_RB[0][1].nrows=1;H_RB[0][1].ncols=1;
H_RB[0][2].value.push_back(0);H_RB[0][2].rows.push_back(0);H_RB[0][2].columns.push_back(0);H_RB[0][2].nrows=1;H_RB[0][2].ncols=1;
H_RB[0][3].value.push_back(0);H_RB[0][3].rows.push_back(0);H_RB[0][3].columns.push_back(0);H_RB[0][3].nrows=1;H_RB[0][3].ncols=1;
	
	
	
	
	
	
	}
	

void DMRG::Initialize_Creation_Anhihilation_Oprts(){
	
	C_up_LB.resize(max_iter+1);C_up_RB.resize(max_iter+1);
C_dn_LB.resize(max_iter+1);C_dn_RB.resize(max_iter+1);
C_up_dag_LB.resize(max_iter+1);C_up_dag_RB.resize(max_iter+1);
C_dn_dag_LB.resize(max_iter+1);C_dn_dag_RB.resize(max_iter+1);
C_up_LB[0].resize(4);C_up_dag_LB[0].resize(4);C_dn_LB[0].resize(4);C_dn_dag_LB[0].resize(4);
C_dn_RB[0].resize(4);C_dn_dag_RB[0].resize(4);C_up_RB[0].resize(4);


C_up_LB[0][0].value.push_back(0);C_up_LB[0][0].rows.push_back(0);C_up_LB[0][0].columns.push_back(0);C_up_LB[0][0].nrows=1;C_up_LB[0][0].ncols=1;
C_up_LB[0][1].value.push_back(0);C_up_LB[0][1].rows.push_back(0);C_up_LB[0][1].columns.push_back(0);C_up_LB[0][1].nrows=1;C_up_LB[0][1].ncols=1;
C_up_LB[0][2].value.push_back(1);C_up_LB[0][2].rows.push_back(0);C_up_LB[0][2].columns.push_back(0);C_up_LB[0][2].nrows=1;C_up_LB[0][2].ncols=1;
C_up_LB[0][3].value.push_back(1);C_up_LB[0][3].rows.push_back(0);C_up_LB[0][3].columns.push_back(0);C_up_LB[0][3].nrows=1;C_up_LB[0][3].ncols=1;


C_up_dag_LB[0][0].value.push_back(1);C_up_dag_LB[0][0].rows.push_back(0);C_up_dag_LB[0][0].columns.push_back(0);C_up_dag_LB[0][0].nrows=1;C_up_dag_LB[0][0].ncols=1;
C_up_dag_LB[0][1].value.push_back(1);C_up_dag_LB[0][1].rows.push_back(0);C_up_dag_LB[0][1].columns.push_back(0);C_up_dag_LB[0][1].nrows=1;C_up_dag_LB[0][1].ncols=1;
C_up_dag_LB[0][2].value.push_back(0);C_up_dag_LB[0][2].rows.push_back(0);C_up_dag_LB[0][2].columns.push_back(0);C_up_dag_LB[0][2].nrows=1;C_up_dag_LB[0][2].ncols=1;
C_up_dag_LB[0][3].value.push_back(0);C_up_dag_LB[0][3].rows.push_back(0);C_up_dag_LB[0][3].columns.push_back(0);C_up_dag_LB[0][3].nrows=1;C_up_dag_LB[0][3].ncols=1;

C_dn_LB[0][0].value.push_back(0);C_dn_LB[0][0].rows.push_back(0);C_dn_LB[0][0].columns.push_back(0);C_dn_LB[0][0].nrows=1;C_dn_LB[0][0].ncols=1;
C_dn_LB[0][1].value.push_back(1);C_dn_LB[0][1].rows.push_back(0);C_dn_LB[0][1].columns.push_back(0);C_dn_LB[0][1].nrows=1;C_dn_LB[0][1].ncols=1;
C_dn_LB[0][2].value.push_back(0);C_dn_LB[0][2].rows.push_back(0);C_dn_LB[0][2].columns.push_back(0);C_dn_LB[0][2].nrows=1;C_dn_LB[0][2].ncols=1;
C_dn_LB[0][3].value.push_back(1);C_dn_LB[0][3].rows.push_back(0);C_dn_LB[0][3].columns.push_back(0);C_dn_LB[0][3].nrows=1;C_dn_LB[0][3].ncols=1;
//C_dn_LB[0][3].value.push_back(1);C_dn_LB[0][3].rows.push_back(0);C_dn_LB[0][3].columns.push_back(0);C_dn_LB[0][3].nrows=1;C_dn_LB[0][3].ncols=1;

C_dn_dag_LB[0][0].value.push_back(1);C_dn_dag_LB[0][0].rows.push_back(0);C_dn_dag_LB[0][0].columns.push_back(0);C_dn_dag_LB[0][0].nrows=1;C_dn_dag_LB[0][0].ncols=1;
C_dn_dag_LB[0][1].value.push_back(0);C_dn_dag_LB[0][1].rows.push_back(0);C_dn_dag_LB[0][1].columns.push_back(0);C_dn_dag_LB[0][1].nrows=1;C_dn_dag_LB[0][1].ncols=1;
C_dn_dag_LB[0][2].value.push_back(1);C_dn_dag_LB[0][2].rows.push_back(0);C_dn_dag_LB[0][2].columns.push_back(0);C_dn_dag_LB[0][2].nrows=1;C_dn_dag_LB[0][2].ncols=1;
//C_dn_dag_LB[0][2].value.push_back(1);C_dn_dag_LB[0][2].rows.push_back(0);C_dn_dag_LB[0][2].columns.push_back(0);C_dn_dag_LB[0][2].nrows=1;C_dn_dag_LB[0][2].ncols=1;
C_dn_dag_LB[0][3].value.push_back(0);C_dn_dag_LB[0][3].rows.push_back(0);C_dn_dag_LB[0][3].columns.push_back(0);C_dn_dag_LB[0][3].nrows=1;C_dn_dag_LB[0][3].ncols=1;

C_up_RB[0][0].value.push_back(1);C_up_RB[0][0].rows.push_back(0);C_up_RB[0][0].columns.push_back(0);C_up_RB[0][0].nrows=1;C_up_RB[0][0].ncols=1;
C_up_RB[0][1].value.push_back(1);C_up_RB[0][1].rows.push_back(0);C_up_RB[0][1].columns.push_back(0);C_up_RB[0][1].nrows=1;C_up_RB[0][1].ncols=1;
C_up_RB[0][2].value.push_back(0);C_up_RB[0][2].rows.push_back(0);C_up_RB[0][2].columns.push_back(0);C_up_RB[0][2].nrows=1;C_up_RB[0][2].ncols=1;
C_up_RB[0][3].value.push_back(0);C_up_RB[0][3].rows.push_back(0);C_up_RB[0][3].columns.push_back(0);C_up_RB[0][3].nrows=1;C_up_RB[0][3].ncols=1;

C_up_dag_RB[0]=C_up_LB[0];

C_dn_RB[0][0].value.push_back(1);C_dn_RB[0][0].rows.push_back(0);C_dn_RB[0][0].columns.push_back(0);C_dn_RB[0][0].nrows=1;C_dn_RB[0][0].ncols=1;
//C_dn_RB[0][0].value.push_back(1);C_dn_RB[0][0].rows.push_back(0);C_dn_RB[0][0].columns.push_back(0);C_dn_RB[0][0].nrows=1;C_dn_RB[0][0].ncols=1;
C_dn_RB[0][1].value.push_back(0);C_dn_RB[0][1].rows.push_back(0);C_dn_RB[0][1].columns.push_back(0);C_dn_RB[0][1].nrows=1;C_dn_RB[0][1].ncols=1;
C_dn_RB[0][2].value.push_back(1);C_dn_RB[0][2].rows.push_back(0);C_dn_RB[0][2].columns.push_back(0);C_dn_RB[0][2].nrows=1;C_dn_RB[0][2].ncols=1;
C_dn_RB[0][3].value.push_back(0);C_dn_RB[0][3].rows.push_back(0);C_dn_RB[0][3].columns.push_back(0);C_dn_RB[0][3].nrows=1;C_dn_RB[0][3].ncols=1;

C_dn_dag_RB[0][0].value.push_back(0);C_dn_dag_RB[0][0].rows.push_back(0);C_dn_dag_RB[0][0].columns.push_back(0);C_dn_dag_RB[0][0].nrows=1;C_dn_dag_RB[0][0].ncols=1;
C_dn_dag_RB[0][1].value.push_back(1);C_dn_dag_RB[0][1].rows.push_back(0);C_dn_dag_RB[0][1].columns.push_back(0);C_dn_dag_RB[0][1].nrows=1;C_dn_dag_RB[0][1].ncols=1;
//C_dn_dag_RB[0][1].value.push_back(1);C_dn_dag_RB[0][1].rows.push_back(0);C_dn_dag_RB[0][1].columns.push_back(0);C_dn_dag_RB[0][1].nrows=1;C_dn_dag_RB[0][1].ncols=1;
C_dn_dag_RB[0][2].value.push_back(0);C_dn_dag_RB[0][2].rows.push_back(0);C_dn_dag_RB[0][2].columns.push_back(0);C_dn_dag_RB[0][2].nrows=1;C_dn_dag_RB[0][2].ncols=1;
C_dn_dag_RB[0][3].value.push_back(1);C_dn_dag_RB[0][3].rows.push_back(0);C_dn_dag_RB[0][3].columns.push_back(0);C_dn_dag_RB[0][3].nrows=1;C_dn_dag_RB[0][3].ncols=1;

	
	
	
	
	}	
	
	
int DMRG::Inverse_LB(int itr,int q){
	int tmp=1234567;
	for(int i=0;i<Q_Eff_Left_Block[itr].size();i++){
		if(Q_Eff_Left_Block[itr][i]==q){tmp=i;}
		}
	return tmp;
	
	}
	
	
int DMRG::Inverse_RB(int itr,int q){
	int tmp=1234567;
	for(int i=0;i<Q_Eff_Right_Block[itr].size();i++){
		if(Q_Eff_Right_Block[itr][i]==q){tmp=i;}
		}
	return tmp;
	
	}


void DMRG::Create_Q_Eff(){
	
	
	
	    double temp_up, temp_dn, temp, temp2; 
				
		Q_Eff_Site.push_back(0);Q_Eff_Site.push_back(1);Q_Eff_Site.push_back(1000);Q_Eff_Site.push_back(1001);
		Q_Eff_System.resize(max_iter+1);
		Q_Eff_Enviroment.resize(max_iter+1);
		Q_Eff_Left_Block.resize(max_iter+1);
		Q_Eff_Right_Block.resize(max_iter+1);
		Q_Eff_Super_Block.resize(max_iter+1);
		
		for(int i=0;i<=max_iter;i++){
			
			
			// --------------------Creating Q_SB---------------------------
			temp_up = (Nup/(2+2*(max_iter+1)))*(2 + 2*(i+1)) ;
			
			if(temp_up+0.5 >= ceil(temp_up)){temp_up = ceil(temp_up);}
			else{temp_up = floor(temp_up);}
			
			temp_dn = (Ndn/(2+2*(max_iter+1)))*(2 + 2*(i+1)) ;
			
			if(temp_dn+0.5 >= ceil(temp_dn)){temp_dn = ceil(temp_dn);}
			else{temp_dn = floor(temp_dn);}
			
			Q_Eff_Super_Block[i]= 1000*temp_up + temp_dn;
			// --------------------Created Q_SB----------------------------
			
			// -------------------Creating Q_E, Q_S----------------------------
			
			temp_up=min( (div(Q_Eff_Super_Block[i],1000).quot), i+2  );
			temp_dn=min( (div(Q_Eff_Super_Block[i],1000).rem), i+2  );
		    temp = temp_up*1000 + temp_dn;
			
			for(int x=0;x<=temp_up;x++){
				for(int y=0;y<=temp_dn;y++){
					
                    if( (div(  (Q_Eff_Super_Block[i]-(temp - x*(1000)- y)),1000  ).quot<=(i+2)) && (div(  (Q_Eff_Super_Block[i]-(temp - x*(1000)- y)),1000  ).rem<=(i+2))  ){
					Q_Eff_Enviroment[i].push_back(temp - x*(1000)- y);
					Q_Eff_System[i].push_back(Q_Eff_Super_Block[i]-(temp - x*(1000)- y));}
					
					}
				}
			
			// --------------------Created Q_E, Q_S----------------------------
			
			
			
			/*if(i==34){
			for(int j=0;j<Enviroment[i].size();j++){
				cout<<System[i][j]<<"\t"<<Enviroment[i][j]<<"\t"<<Enviroment[i][j]+System[i][j]<<"\t"<<Super_Block[i]<<endl;
												}
					cout<<"size   =  "<<Enviroment[i].size()<<endl;							
			
					}*/
			
			
			
			// -----------------Creating Q_LB, Q_RB ----------------------------
			if(i==0){Q_Eff_Left_Block[i].push_back(0);Q_Eff_Left_Block[i].push_back(1);Q_Eff_Left_Block[i].push_back(1000);Q_Eff_Left_Block[i].push_back(1001);
				   Q_Eff_Right_Block[i].push_back(1001);Q_Eff_Right_Block[i].push_back(1000);Q_Eff_Right_Block[i].push_back(1);Q_Eff_Right_Block[i].push_back(0);}
			else{
				Q_Eff_Left_Block[i].resize(Q_Eff_Enviroment[i-1].size());
				Q_Eff_Right_Block[i].resize(Q_Eff_Enviroment[i-1].size());
				
				for(int j=0;j<Q_Eff_Enviroment[i-1].size();j++){
				Q_Eff_Left_Block[i][j]=Q_Eff_System[i-1][j];
				Q_Eff_Right_Block[i][j]=Q_Eff_Enviroment[i-1][j];								}
												
			
					}
				
					   
			
			// -----------------Created Q_LB, Q_RB -----------------------------
			
			}
		
			
		
		
		
}
