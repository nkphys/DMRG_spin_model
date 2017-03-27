#include "functions.h"
#include <iostream>
#include <fstream>
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


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx writing Mat_2_doub in file xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void write(Mat_2_doub A, ofstream & outputfile){

    //ofstream outputfile(filepath.c_str());

    outputfile<<"nrows"<<endl;
    outputfile<<A.size()<<endl;
    outputfile<<"ncols"<<endl;
    outputfile<<A.size()<<endl;

    outputfile<<"rows columns values"<<endl;
    for(int row=0;row<A.size();row++){
        for(int col=0;col<A.size();col++){
#ifndef WITH_COMPLEX
            outputfile<<row<<"  "<<col<<"  "<<A[row][col]<<endl;
#endif
#ifdef WITH_COMPLEX
            outputfile<<row<<"  "<<col<<"  "<<A[row][col].real()<<"   "<<A[row][col].imag()<<endl;
#endif
        }
    }


}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx writing Matrix_COO in file xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void write(Matrix_COO A, ofstream & outputfile){


    //ofstream outputfile(filepath.c_str());

    outputfile<<"nrows"<<endl;
    outputfile<<A.nrows<<endl;
    outputfile<<"ncols"<<endl;
    outputfile<<A.ncols<<endl;

    outputfile<<"rows columns values"<<endl;
    for(int val_size=0;val_size<A.value.size();val_size++){
        outputfile<< A.rows[val_size]<<"  "<<A.columns[val_size]<<"  ";
#ifndef WITH_COMPLEX
        outputfile<<A.value[val_size]<<endl;
#endif
#ifdef WITH_COMPLEX
        outputfile<<A.value[val_size].real()<<"   "<<A.value[val_size].imag()<<endl;
#endif

    }




}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxx--------------------------Dot Product--------------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
type_double dot_product(Mat_1_doub vec1, Mat_1_doub vec2){
    //This dot_product is parallelized always, create another one with _PARALLELIZE_AT_MATRICES_LEVEL
    type_double temp;
    type_double temp1=0;
    bool Parallelize_dot_product;
    Parallelize_dot_product=true;

    if(!Parallelize_dot_product){goto skiploop_143;}
#pragma omp parallel for default(shared) reduction(+:temp1)
skiploop_143:
    for(int i=0;i<vec1.size();i++){
#ifndef WITH_COMPLEX
        temp1 = temp1 + (vec1[i])*(vec2[i]);
#endif
#ifdef WITH_COMPLEX
        temp1.real() = temp1.real() + (vec1[i].real())*(vec2[i].real()) - (vec1[i].imag())*(vec2[i].imag());
        temp1.imag() = temp1.imag() + (vec1[i].real())*(vec2[i].imag()) + (vec1[i].imag())*(vec2[i].real());
#endif

    }

    (temp) = temp1;

    return temp;


}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//

DOUBLE_MKL_CSR_MAT  CSR_MAT_TO_MKL_SPARSE(DOUBLE_CSR_MAT mat1){

    int i;
    DOUBLE_MKL_CSR_MAT mat2;
    mat2.nrows=mat1.row_ind.size() - 1;
    mat2.ncols=mat2.nrows;


    mat2.value = new double [mat1.value.size()];
    mat2.columns = new MKL_INT[mat1.columns.size()];
    mat2.row_ind = new MKL_INT[mat1.row_ind.size()];

    for(i=0;i<mat1.value.size();i++){
        mat2.value[i]=mat1.value[i];
    }
    for(i=0;i<mat1.columns.size();i++){
        mat2.columns[i]=mat1.columns[i];
    }
    for(i=0;i<mat1.row_ind.size();i++){
        mat2.row_ind[i]=mat1.row_ind[i];
    }

    return mat2;
}

//-----------------------------------------------------------------------------------------------------------------------------//
void Direct_Sum(Matrix_COO A, Matrix_COO B, Matrix_COO &temp){

    int tmp=0;
    bool add_it;
    temp.nrows = A.nrows + B.nrows;
    temp.ncols = A.ncols + B.ncols;

    temp.value.clear();
    temp.rows.clear();
    temp.columns.clear();

    for(int i=0;i<A.value.size();i++){

#ifndef WITH_COMPLEX
        add_it = (A.value[i]!=0);
#endif
#ifdef WITH_COMPLEX
        add_it = (A.value[i]!=(0.0,0.0));
#endif

        if(add_it){
            temp.value.push_back(A.value[i]);tmp=tmp+1;
            temp.rows.push_back(A.rows[i]);
            temp.columns.push_back(A.columns[i]);}
    }

    for(int i=0;i<B.value.size();i++){

#ifndef WITH_COMPLEX
        add_it = (B.value[i]!=0);
#endif
#ifdef WITH_COMPLEX
        add_it = (B.value[i].real()!=(0.0,0.0));
#endif

        if(add_it){
            temp.value.push_back(B.value[i]);tmp=tmp+1;
            temp.rows.push_back(B.rows[i]+A.nrows);
            temp.columns.push_back(B.columns[i]+A.ncols);}
    }

    if (tmp==0){temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();
        temp.value.push_back(0);
        temp.rows.push_back(0);
        temp.columns.push_back(0);}



}



//-----------------------------------------------------------------------------------------------------------------------------//
Matrix_COO Direct_Sum(Matrix_COO A, Matrix_COO B){

    Matrix_COO temp;
    int tmp=0;
    temp.nrows = A.nrows + B.nrows;
    temp.ncols = A.ncols + B.ncols;

    temp.value.clear();
    temp.rows.clear();
    temp.columns.clear();

    for(int i=0;i<A.value.size();i++){

#ifndef WITH_COMPLEX
        add_it = (A.value[i]!=0);
#endif
#ifdef WITH_COMPLEX
        add_it = (A.value[i].real()!=(0.0,0.0));
#endif

        if(add_it){
            temp.value.push_back(A.value[i]);tmp=tmp+1;
            temp.rows.push_back(A.rows[i]);
            temp.columns.push_back(A.columns[i]);}
    }

    for(int i=0;i<B.value.size();i++){

#ifndef WITH_COMPLEX
        add_it = (B.value[i]!=0);
#endif
#ifdef WITH_COMPLEX
        add_it = (B.value[i].real()!=(0.0,0.0));
#endif

        if(add_it){
            temp.value.push_back(B.value[i]);tmp=tmp+1;
            temp.rows.push_back(B.rows[i]+A.nrows);
            temp.columns.push_back(B.columns[i]+A.ncols);}
    }

    if (tmp==0){temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();
        temp.value.push_back(0);
        temp.rows.push_back(0);
        temp.columns.push_back(0);}

    return temp;


}

//-------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------//
void Direct_Product(Matrix_COO A, Matrix_COO B, Matrix_COO & temp){
    //cout<<"fine"<<endl;
    bool Parallelize_this;
    Parallelize_this=_PARALLELIZE_AT_MATRICES_LEVEL;
    temp.nrows = A.nrows*B.nrows;
    temp.ncols = A.ncols*B.ncols;

    temp.value.clear();
    temp.value.resize(A.value.size()*B.value.size());
    temp.rows.clear();
    temp.rows.resize(A.value.size()*B.value.size());
    temp.columns.clear();
    temp.columns.resize(A.value.size()*B.value.size());


    int i_A;
    int i_B;
    Mat_1_int A_R_offset, A_L_offset;
    Mat_1_int B_R_offset, B_L_offset;
    A_L_offset.resize(A.value.size());
    A_R_offset.resize(A.value.size());
    B_L_offset.resize(B.value.size());
    B_R_offset.resize(B.value.size());


    int i_A_l=0;
    int i_B_l=0;
    int i_A_r=A.value.size()-1;
    int i_B_r=B.value.size()-1;


    //Left_offset
    for (i_A=0;i_A<A.value.size();i_A++){
        if(A.rows[i_A]!=A.rows[i_A_l]){
            i_A_l=i_A;
        }
        A_L_offset[i_A]=i_A_l;

    }

    //Right_offset
    for(i_A=A.value.size()-1;i_A>=0;i_A--){
        if(A.rows[i_A]!=A.rows[i_A_r]){
            i_A_r=i_A;
        }
        A_R_offset[i_A]=i_A_r+1;
    }


    //Left_offset
    for (i_B=0;i_B<B.value.size();i_B++){
        if(B.rows[i_B]!=B.rows[i_B_l]){
            i_B_l=i_B;
        }
        B_L_offset[i_B]=i_B_l;

    }

    //Right_offset
    for(i_B=B.value.size()-1;i_B>=0;i_B--){
        if(B.rows[i_B]!=B.rows[i_B_r]){
            i_B_r=i_B;
        }
        B_R_offset[i_B]=i_B_r+1;
    }

    bool f1, f2;
    f1 = (A.nrows>B.nrows && A.nrows>4);
    f2 = (B.nrows>=A.nrows && B.nrows>4);
    int temp_i;

    if(f1==false){goto skiploop101;}
    if(!Parallelize_this){goto skiploop101;}
#pragma omp parallel for default(shared) private(i_A,temp_i) if(f1)
skiploop101:
    for (i_A=0;i_A<A.value.size();i_A++){
        if(f2==false){goto skiploop102;}
        if(!Parallelize_this){goto skiploop102;}
#pragma omp parallel for default(shared) private(i_B,temp_i) if(f2)
skiploop102:
        for (i_B=0;i_B<B.value.size();i_B++){

            temp_i = A_L_offset[i_A]*B.value.size() + (A_R_offset[i_A] - A_L_offset[i_A])*B_L_offset[i_B] +
                    (i_A - A_L_offset[i_A])*(B_R_offset[i_B] - B_L_offset[i_B]) +
                    (i_B - B_L_offset[i_B]);
#ifndef WITH_COMPLEX
            temp.value[temp_i]=A.value[i_A]*B.value[i_B];
#endif
#ifdef WITH_COMPLEX
            temp.value[temp_i].real()=A.value[i_A].real()*B.value[i_B].real() - A.value[i_A].imag()*B.value[i_B].imag();
            temp.value[temp_i].real()=A.value[i_A].real()*B.value[i_B].imag() + A.value[i_A].imag()*B.value[i_B].real();
#endif
            temp.rows[temp_i]=(B.nrows*A.rows[i_A] + B.rows[i_B]);
            temp.columns[temp_i]=(B.ncols*A.columns[i_A] + B.columns[i_B]);

        }

    }




    if((A.value.size()==0 || B.value.size()==0)){
        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();
        temp.value.push_back(0);
        temp.rows.push_back(0);
        temp.columns.push_back(0);}
}
//---------------------------------------------------------------------------------------------------------------------------------//
Matrix_COO Direct_Product(Matrix_COO A, Matrix_COO B){

    bool Parallelize_this;
    Parallelize_this=_PARALLELIZE_AT_MATRICES_LEVEL;
    Matrix_COO temp;
    temp.nrows = A.nrows*B.nrows;
    temp.ncols = A.ncols*B.ncols;

    temp.value.clear();
    temp.value.resize(A.value.size()*B.value.size());
    temp.rows.clear();
    temp.rows.resize(A.value.size()*B.value.size());
    temp.columns.clear();
    temp.columns.resize(A.value.size()*B.value.size());


    int i_A;
    int i_B;
    Mat_1_int A_R_offset, A_L_offset;
    Mat_1_int B_R_offset, B_L_offset;
    A_L_offset.resize(A.value.size());
    A_R_offset.resize(A.value.size());
    B_L_offset.resize(B.value.size());
    B_R_offset.resize(B.value.size());


    int i_A_l=0;
    int i_B_l=0;
    int i_A_r=A.value.size()-1;
    int i_B_r=B.value.size()-1;


    //Left_offset
    for (i_A=0;i_A<A.value.size();i_A++){
        if(A.rows[i_A]!=A.rows[i_A_l]){
            i_A_l=i_A;
        }
        A_L_offset[i_A]=i_A_l;

    }

    //Right_offset
    for(i_A=A.value.size()-1;i_A>=0;i_A--){
        if(A.rows[i_A]!=A.rows[i_A_r]){
            i_A_r=i_A;
        }
        A_R_offset[i_A]=i_A_r+1;
    }


    //Left_offset
    for (i_B=0;i_B<B.value.size();i_B++){
        if(B.rows[i_B]!=B.rows[i_B_l]){
            i_B_l=i_B;
        }
        B_L_offset[i_B]=i_B_l;

    }

    //Right_offset
    for(i_B=B.value.size()-1;i_B>=0;i_B--){
        if(B.rows[i_B]!=B.rows[i_B_r]){
            i_B_r=i_B;
        }
        B_R_offset[i_B]=i_B_r+1;
    }

    bool f1, f2;
    f1 = (A.nrows>B.nrows && A.nrows>4);
    f2 = (B.nrows>=A.nrows && B.nrows>4);
    int temp_i;

    if(f1==false){goto skiploop111;}
    if(!Parallelize_this){goto skiploop111;}
#pragma omp parallel for default(shared) private(temp_i,i_A) if(f1)
skiploop111:
    for (i_A=0;i_A<A.value.size();i_A++){
        if(f2==false){goto skiploop112;}
        if(!Parallelize_this){goto skiploop112;}
#pragma omp parallel for default(shared) private(temp_i,i_B) if(f2)
skiploop112:
        for (i_B=0;i_B<B.value.size();i_B++){

            temp_i = A_L_offset[i_A]*B.value.size() + (A_R_offset[i_A] - A_L_offset[i_A])*B_L_offset[i_B] +
                    (i_A - A_L_offset[i_A])*(B_R_offset[i_B] - B_L_offset[i_B]) +
                    (i_B - B_L_offset[i_B]);

#ifndef WITH_COMPLEX
            temp.value[temp_i]=A.value[i_A]*B.value[i_B];
#endif
#ifdef WITH_COMPLEX
            temp.value[temp_i].real()=A.value[i_A].real()*B.value[i_B].real() - A.value[i_A].imag()*B.value[i_B].imag();
            temp.value[temp_i].real()=A.value[i_A].real()*B.value[i_B].imag() + A.value[i_A].imag()*B.value[i_B].real();
#endif
            temp.rows[temp_i]=(B.nrows*A.rows[i_A] + B.rows[i_B]);
            temp.columns[temp_i]=(B.ncols*A.columns[i_A] + B.columns[i_B]);

        }

    }




    if((A.value.size()==0 || B.value.size()==0)){
        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();
        temp.value.push_back(0);
        temp.rows.push_back(0);
        temp.columns.push_back(0);}

    return temp;
}


//---------------------------------------------------------------------------------------//
Matrix_COO Direct_Product_prll(Matrix_COO A, Matrix_COO B){
    Matrix_COO temp;
    //cout<<"fine"<<endl;
    temp.nrows = A.nrows*B.nrows;
    temp.ncols = A.ncols*B.ncols;

    temp.value.clear();
    temp.value.resize(A.value.size()*B.value.size());
    temp.rows.clear();
    temp.rows.resize(A.rows.size()*B.rows.size());
    temp.columns.clear();
    temp.columns.resize(A.columns.size()*B.columns.size());
    bool f1, f2;
    int counter=0;

    f1 = (A.nrows>B.nrows && A.nrows>4);
    f2 = (B.nrows>=A.nrows && B.nrows>4);

    if(f1==false){goto skiploop1;}
#pragma omp parallel for default(shared) private(counter) if(f1)
skiploop1:
    for(int i=0;i<A.value.size();i++){
        counter=i*B.value.size();
        if(f2==false){goto skiploop2;}
#pragma omp parallel for default(shared) if(f2)
skiploop2:
        for(int j=0;j<B.value.size();j++){


#ifndef WITH_COMPLEX
            temp.value[counter+j]=A.value[i]*B.value[j];
#endif
#ifdef WITH_COMPLEX
            temp.value[counter+j].real()=A.value[i].real()*B.value[j].real() - A.value[i].imag()*B.value[j].imag();
            temp.value[counter+j].real()=A.value[i].real()*B.value[j].imag() + A.value[i].imag()*B.value[j].real();
#endif

            temp.rows[counter+j]=(B.nrows*A.rows[i] + B.rows[j]);
            temp.columns[counter+j]=(B.ncols*A.columns[i] + B.columns[j]);

        }

    }

    if((A.value.size()==0 || B.value.size()==0)){
        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();
        temp.value.push_back(0);
        temp.rows.push_back(0);
        temp.columns.push_back(0);}

    return temp;

}
//---------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------//
//Matrix_COO Direct_Product_push_back(Matrix_COO A, Matrix_COO B){
//    Matrix_COO temp;
//    temp.nrows = A.nrows*B.nrows;
//    temp.ncols = A.ncols*B.ncols;
//    temp.value.clear();
//    temp.rows.clear();
//    temp.columns.clear();
//    bool f1, f2;

//    int tmp=0;


//    for(int i=0;i<A.value.size();i++){
//        for(int j=0;j<B.value.size();j++){
//            if(A.value[i]*B.value[j]!=0){
//                temp.value.push_back(A.value[i]*B.value[j]);
//                temp.rows.push_back(B.nrows*A.rows[i] + B.rows[j]);
//                temp.columns.push_back(B.ncols*A.columns[i] + B.columns[j]);
//                tmp=tmp+1;}

//        }
//    }

//    if(tmp==0){
//        temp.value.clear();
//        temp.rows.clear();
//        temp.columns.clear();
//        temp.value.push_back(0);
//        temp.rows.push_back(0);
//        temp.columns.push_back(0);}


//    return temp;

//}
//---------------------------------------------------------------------------------------------------------------------------------//
void Direct_Product_prll(Matrix_COO A, Matrix_COO B, Matrix_COO & temp){
    //cout<<"fine"<<endl;
    temp.nrows = A.nrows*B.nrows;
    temp.ncols = A.ncols*B.ncols;

    temp.value.clear();
    temp.value.resize(A.value.size()*B.value.size());
    temp.rows.clear();
    temp.rows.resize(A.rows.size()*B.rows.size());
    temp.columns.clear();
    temp.columns.resize(A.columns.size()*B.columns.size());
    bool f1, f2;
    int counter=0;

    f1 = (A.nrows>B.nrows && A.nrows>4);
    f2 = (B.nrows>=A.nrows && B.nrows>4);

    if(f1==false){goto skiploop1;}
#pragma omp parallel for default(shared) private(counter) if(f1)
skiploop1:
    for(int i=0;i<A.value.size();i++){
        counter=i*B.value.size();
        if(f2==false){goto skiploop2;}
#pragma omp parallel for default(shared) if(f2)
skiploop2:
        for(int j=0;j<B.value.size();j++){

#ifndef WITH_COMPLEX
                temp.value[counter+j]=A.value[i]*B.value[j];
#endif
#ifdef WITH_COMPLEX
                temp.value[counter+j].real()=A.value[i].real()*B.value[j].real() - A.value[i].imag()*B.value[j].imag();
                temp.value[counter+j].real()=A.value[i].real()*B.value[j].imag() + A.value[i].imag()*B.value[j].real();
#endif
            temp.rows[counter+j]=(B.nrows*A.rows[i] + B.rows[j]);
            temp.columns[counter+j]=(B.ncols*A.columns[i] + B.columns[j]);



        }

    }

    if((A.value.size()==0 || B.value.size()==0)){
        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();
        temp.value.push_back(0);
        temp.rows.push_back(0);
        temp.columns.push_back(0);}
}
//---------------------------------------------------------------------------------------------------------------------------------//
//*/---------------------------------------------------------------------------------------------------------------------------------//
//void Direct_Product_pr(Matrix_COO A, Matrix_COO B, Matrix_COO & temp){

//    //cout<<"fine"<<endl;
//    temp.nrows = A.nrows*B.nrows;
//    temp.ncols = A.ncols*B.ncols;

//    temp.value.clear();
//    temp.rows.clear();
//    temp.columns.clear();

//    int tmp=0;
//    for(int i=0;i<A.value.size();i++){
//        for(int j=0;j<B.value.size();j++){
//            if(/*A.value[i]*B.value[j]!=0*/ true){
//                temp.value.push_back(A.value[i]*B.value[j]);
//                temp.rows.push_back(B.nrows*A.rows[i] + B.rows[j]);
//                temp.columns.push_back(B.ncols*A.columns[i] + B.columns[j]);
//                tmp=tmp+1;}

//        }
//    }

//    if(tmp==0){
//        temp.value.clear();
//        temp.rows.clear();
//        temp.columns.clear();
//        temp.value.push_back(0);
//        temp.rows.push_back(0);
//        temp.columns.push_back(0);}


//}
//----------------------------------------------------------------------------------------------------------*/
//------------------------------------------------------------------------------------------------------------------------------//
//void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & temp, double value1, double value2){

//    if((A.nrows == B.nrows) && (A.ncols == B.ncols) ){


//        int x=0;
//        temp.nrows = A.nrows;
//        temp.ncols = A.ncols;

//        temp.value.clear();
//        temp.rows.clear();
//        temp.columns.clear();

//        int N= A.nrows + A.ncols;
//        int i=0;
//        int j=0;
//        int a_i, b_j;
//        while( i < A.value.size()){
//            if(j<B.value.size()){
//                while( j<B.value.size()){

//                    if(i<A.value.size()){
//                        a_i=N*(A.rows[i])+ (A.columns[i]);
//                        b_j=N*(B.rows[j])+ (B.columns[j]);



//                        //cout<<"i,j = "<<i<<", "<<j<<endl;
//                        if(a_i==b_j){temp.value.push_back(value1*A.value[i]+value2*B.value[j]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;j=j+1;}
//                        if(a_i>b_j){temp.value.push_back(value2*B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                            j=j+1;}
//                        if(a_i<b_j){temp.value.push_back(value1*A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;}
//                    }
//                    else{temp.value.push_back(value2*B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                        j=j+1;}


//                }}

//            else{//cout<<"i,j = "<<i<<", "<<j<<endl;
//                temp.value.push_back(value1*A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                i=i+1;}

//        }


//        for(int i=1;i<temp.value.size();i++){
//            if(temp.value[i]==0){temp.value.erase (temp.value.begin()+i);temp.rows.erase (temp.rows.begin()+i);temp.columns.erase (temp.columns.begin()+i);}
//        }
//        if((temp.value.size()>1)&&(temp.value[0]==0)){temp.value.erase (temp.value.begin());temp.rows.erase (temp.rows.begin());temp.columns.erase (temp.columns.begin());}

//    }

//    else{cout<<"Error in doing Sum"<<endl;}

//}
//----------------------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------------------//
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & temp, type_double value1, type_double value2){

    int a_i=0;
    int b_j=0;
    int row_a_i,col_a_i, row_b_j,col_b_j;
    type_double temp1;

    if((A.nrows == B.nrows) && (A.ncols == B.ncols)){

        temp.nrows = A.nrows;
        temp.ncols = A.ncols;

        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();

        while(a_i<A.value.size() || b_j<B.value.size()){

            if(a_i<A.value.size()){
                row_a_i=A.rows[a_i];
                col_a_i=A.columns[a_i];
            }
            else{
                row_a_i=A.nrows+1;
                col_a_i=A.ncols+1;
            }

            if(b_j<B.value.size()){
                row_b_j=B.rows[b_j];
                col_b_j=B.columns[b_j];
            }
            else{
                row_b_j=B.nrows+1;
                col_b_j=B.ncols+1;
            }

            //if element of B comes before element of A
            if( (row_b_j < row_a_i)  ||
                    (row_b_j == row_a_i && col_b_j < col_a_i)
                    )
            {
#ifndef WITH_COMPLEX
                    temp.value.push_back(value2*B.value[b_j]);
#endif
#ifdef WITH_COMPLEX
                    temp1.real()=value2.real()*B.value[b_j].real() - value2.imag()*B.value[b_j].imag();
                    temp1.imag()=value2.real()*B.value[b_j].imag() + value2.imag()*B.value[b_j].real();
                    temp.value.push_back(temp1);
#endif
                temp.rows.push_back(row_b_j);
                temp.columns.push_back(col_b_j);
                b_j=b_j+1;

            }
            //if element of A comes before element of B
            if( (row_b_j > row_a_i) ||
                    (row_b_j == row_a_i && col_b_j > col_a_i)
                    )

            {
                temp.value.push_back(value1*A.value[a_i]);
                temp.rows.push_back(row_a_i);
                temp.columns.push_back(col_a_i);
                a_i=a_i+1;

            }
            //if elements of B, A comes at same place
            if(row_b_j==row_a_i && col_b_j==col_a_i){
                temp.value.push_back(value1*A.value[a_i] + value2*B.value[b_j]);
                temp.rows.push_back(row_a_i);
                temp.columns.push_back(col_a_i);
                a_i=a_i+1;
                b_j=b_j+1;
            }
        }

    }
    else{cout<<"Error in doing Sum"<<endl;}
}
//----------------------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------------------//
//void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & temp){

//    if((A.nrows == B.nrows) && (A.ncols == B.ncols) ){


//        int x=0;
//        temp.nrows = A.nrows;
//        temp.ncols = A.ncols;

//        temp.value.clear();
//        temp.rows.clear();
//        temp.columns.clear();

//        int N= A.nrows + A.ncols;
//        int i=0;
//        int j=0;
//        int a_i, b_j;
//        while( i < A.value.size()){
//            if(j<B.value.size()){
//                while( j<B.value.size()){

//                    if(i<A.value.size()){
//                        a_i=N*(A.rows[i])+ (A.columns[i]);
//                        b_j=N*(B.rows[j])+ (B.columns[j]);



//                        //cout<<"i,j = "<<i<<", "<<j<<endl;
//                        if(a_i==b_j){temp.value.push_back(A.value[i]+B.value[j]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;j=j+1;}
//                        if(a_i>b_j){temp.value.push_back(B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                            j=j+1;}
//                        if(a_i<b_j){temp.value.push_back(A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;}
//                    }
//                    else{temp.value.push_back(B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                        j=j+1;}


//                }}

//            else{//cout<<"i,j = "<<i<<", "<<j<<endl;
//                temp.value.push_back(A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                i=i+1;}

//        }


//        for(int i=1;i<temp.value.size();i++){
//            if(temp.value[i]==0){temp.value.erase (temp.value.begin()+i);temp.rows.erase (temp.rows.begin()+i);temp.columns.erase (temp.columns.begin()+i);}
//        }
//        if((temp.value.size()>1)&&(temp.value[0]==0)){temp.value.erase (temp.value.begin());temp.rows.erase (temp.rows.begin());temp.columns.erase (temp.columns.begin());}

//    }

//    else{cout<<"Error in doing Sum"<<endl;}

//}
//----------------------------------------------------------------------------------------------------------------------------//
//Matrix_COO Sum(Matrix_COO A, Matrix_COO B){

//    if((A.nrows == B.nrows) && (A.ncols == B.ncols) ){


//        Matrix_COO temp;
//        temp.nrows = A.nrows;
//        temp.ncols = A.ncols;

//        temp.value.clear();
//        temp.rows.clear();
//        temp.columns.clear();

//        int N= A.nrows + A.ncols;
//        int i=0;
//        int j=0;
//        int a_i, b_j;
//        while( i < A.value.size()){
//            if(j<B.value.size()){
//                while( j<B.value.size()){

//                    if(i<A.value.size()){
//                        a_i=N*(A.rows[i])+ (A.columns[i]);
//                        b_j=N*(B.rows[j])+ (B.columns[j]);



//                        //cout<<"i,j = "<<i<<", "<<j<<endl;
//                        if(a_i==b_j){temp.value.push_back(A.value[i]+B.value[j]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;j=j+1;}
//                        if(a_i>b_j){temp.value.push_back(B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                            j=j+1;}
//                        if(a_i<b_j){temp.value.push_back(A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;}
//                    }
//                    else{temp.value.push_back(B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                        j=j+1;}


//                }}

//            else{//cout<<"i,j = "<<i<<", "<<j<<endl;
//                temp.value.push_back(A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                i=i+1;}

//        }


//        for(int i=1;i<temp.value.size();i++){
//            if(temp.value[i]==0){temp.value.erase (temp.value.begin()+i);temp.rows.erase (temp.rows.begin()+i);temp.columns.erase (temp.columns.begin()+i);}
//        }
//        if((temp.value.size()>1)&&(temp.value[0]==0)){temp.value.erase (temp.value.begin());temp.rows.erase (temp.rows.begin());temp.columns.erase (temp.columns.begin());}

//        return temp;
//    }

//    else{cout<<"Error in doing Sum"<<endl;}

//}
//----------------------------------------------------------------------------------------------------------------------------//
//Matrix_COO Sum(Matrix_COO A, Matrix_COO B, double value1, double value2){

//    if((A.nrows == B.nrows) && (A.ncols == B.ncols) ){


//        Matrix_COO temp;
//        temp.nrows = A.nrows;
//        temp.ncols = A.ncols;

//        temp.value.clear();
//        temp.rows.clear();
//        temp.columns.clear();

//        int N= A.nrows + A.ncols;
//        int i=0;
//        int j=0;
//        int a_i, b_j;
//        while( i < A.value.size()){
//            if(j<B.value.size()){
//                while( j<B.value.size()){

//                    if(i<A.value.size()){
//                        a_i=N*(A.rows[i])+ (A.columns[i]);
//                        b_j=N*(B.rows[j])+ (B.columns[j]);



//                        //cout<<"i,j = "<<i<<", "<<j<<endl;
//                        if(a_i==b_j){temp.value.push_back(value1*A.value[i]+value2*B.value[j]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;j=j+1;}
//                        if(a_i>b_j){temp.value.push_back(value2*B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                            j=j+1;}
//                        if(a_i<b_j){temp.value.push_back(value1*A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                            i=i+1;}
//                    }
//                    else{temp.value.push_back(value2*B.value[j]); temp.rows.push_back(B.rows[j]); temp.columns.push_back(B.columns[j]);
//                        j=j+1;}


//                }}

//            else{//cout<<"i,j = "<<i<<", "<<j<<endl;
//                temp.value.push_back(value1*A.value[i]); temp.rows.push_back(A.rows[i]); temp.columns.push_back(A.columns[i]);
//                i=i+1;}

//        }


//        for(int i=1;i<temp.value.size();i++){
//            if(temp.value[i]==0){temp.value.erase (temp.value.begin()+i);temp.rows.erase (temp.rows.begin()+i);temp.columns.erase (temp.columns.begin()+i);}
//        }
//        if((temp.value.size()>1)&&(temp.value[0]==0)){temp.value.erase (temp.value.begin());temp.rows.erase (temp.rows.begin());temp.columns.erase (temp.columns.begin());}

//        return temp;
//    }

//    else{cout<<"Error in doing Sum"<<endl;}

//}
//----------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------//
Matrix_COO Identity(int N){

    Matrix_COO Iden;
    type_double temp1;
    for(int i=0;i<N;i++){
#ifndef WITH_COMPLEX
        Iden.value.push_back(1.0);
#endif
#ifdef WITH_COMPLEX
        temp1=(1.0,0.0);
        Iden.value.push_back(temp1);
#endif
        Iden.rows.push_back(i);
        Iden.columns.push_back(i);
    }
    Iden.nrows = N;
    Iden.ncols = N;


    return Iden;


}
//----------------------------------------------------------------------------------------------------------------------------//
//void Add_FC(int c, int r, Matrix_COO temp, Matrix_COO & mat, double sgn, double thop){

//    if(temp.value.size()==0){}
//    else{
//        Matrix_COO A;
//        A=temp;
//        for(int i=0;i<A.rows.size();i++){
//            A.rows[i]=A.rows[i]+r;
//            A.columns[i]=A.columns[i]+c;
//            A.value[i]=sgn*thop*A.value[i];}
//        A.nrows=mat.nrows;
//        A.ncols=mat.ncols;


//        Sum(A, mat, mat ,1, 1);
//        A.value.clear();A.rows.clear();A.columns.clear();}


//}
//-----------------------------------------------------------------------------------------------------------------------------//
void Printing_COO_Matrix(Matrix_COO A){

    cout<< "no. of rows = "<<A.nrows<<", "<< "no. of columns = "<<A.ncols<<endl;

    int vs,rs,cs;
    vs=A.value.size();
    rs=A.rows.size();
    cs=A.columns.size();
    cout<<"value = "<<"\t"<<"{";
    for(int i=0;i<vs-1;i++){cout<<A.value[i]<<","<<"\t";}
    cout<<A.value[vs-1]<<"}"<<endl;

    cout<<"rows = "<<"\t"<<"\t"<<"{";
    for(int i=0;i<A.rows.size()-1;i++){cout<<A.rows[i]<<","<<"\t";}
    cout<<A.rows[rs-1]<<"}"<<endl;

    cout<<"columns = "<<"\t"<<"{";
    for(int i=0;i<A.columns.size()-1;i++){cout<<A.columns[i]<<","<<"\t";}
    cout<<A.columns[cs-1]<<"}"<<endl;

}
//----------------------------------------------------------------------------------------------------------------------------//
