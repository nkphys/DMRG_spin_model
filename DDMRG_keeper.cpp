//This is class for DDMRG
#include "DDMRG_keeper.h"
void DDMRG::Initialize_parameters(){

    loop_0_DDMRG=4;
    DDMRG_bool=true;
    Targetting_omega_space=false;
    Finite_loops.resize(3);
    Finite_loops[0]=-8;Finite_loops[1]=8;Finite_loops[2]=-8;
    site_A=1;
    omega=0.1;
    eta=0.01;
    weight_GS=0.4;weight_A=0.3; weight_X=0.3;
    Operator_applied=false;


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
    one.real(1.0);
    one.imag(0.0);
    zero.real(0.0);
    zero.imag(0.0);
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
                    value.real(value_real);
                    value.imag(value_imag);
#endif

                    for(int l=0;l<Unitary_Eig_vecs.size();l++){
                        Vec_X[i] = Vec_X[i] + value*Krylov_space_vecs[j][i]*Unitary_Eig_vecs[k][j]*
                                conjugate(Krylov_space_vecs[l][m])*
                                conjugate(Unitary_Eig_vecs[k][l])*Vec_A[m];//may be,problem is here

                    }
                }

            }
        }

    }

}

void DDMRG::Calculate_Norms_of_vecX_and_A(){
    type_double tmpnrm_type_double; //new
    double tmpnrm;

    tmpnrm_type_double =(dot_product(Vec_X,Vec_X));
#ifndef WITH_COMPLEX
    tmpnrm=sqrt(tmpnrm_type_double);
#endif
#ifdef WITH_COMPLEX
    tmpnrm=sqrt(tmpnrm_type_double.real());
#endif
    Norm_X=tmpnrm;

    tmpnrm_type_double =(dot_product(Vec_A,Vec_A));
#ifndef WITH_COMPLEX
    tmpnrm=sqrt(tmpnrm_type_double);
#endif
#ifdef WITH_COMPLEX
    tmpnrm=sqrt(tmpnrm_type_double.real());
#endif
    Norm_A=tmpnrm;

}
void DDMRG::test2(){
    int j=0;
}
