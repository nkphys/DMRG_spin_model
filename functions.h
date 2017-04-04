#include "tensor_type.h"

type_double dot_product(Mat_1_doub vec1, Mat_1_doub vec2);
DOUBLE_MKL_CSR_MAT CSR_MAT_TO_MKL_SPARSE(DOUBLE_CSR_MAT mat1 );
//void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & temp);
//Matrix_COO Sum(Matrix_COO A, Matrix_COO B);
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & temp, type_double  value1, type_double value2);
//Matrix_COO Sum(Matrix_COO A, Matrix_COO B, double value1, double value2);
Matrix_COO Direct_Sum(Matrix_COO A, Matrix_COO B);
void Direct_Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & temp);
Matrix_COO Direct_Product(Matrix_COO A, Matrix_COO B);
Matrix_COO Direct_Product_push_back(Matrix_COO A, Matrix_COO B);
void Direct_Product(Matrix_COO A, Matrix_COO B, Matrix_COO & temp);
void Printing_COO_Matrix(Matrix_COO A);
Matrix_COO Identity(int N);
void Add_FC(int c, int r, Matrix_COO temp, Matrix_COO & mat, double sgn, double thop);
void write(Matrix_COO A, ofstream & outputfile);
void write(Mat_2_doub A, ofstream & outputfile);
complex<double> conjugate(complex<double> temp);
double conjugate(double temp);
Matrix_COO Renormalize(double* UL,double* UR, Matrix_COO A, int m_UL, int m_UR, double Lanc_Error);
void Renormalize(Matrix_COO A, double* UL,double* UR, Matrix_COO & B, int m_UL, int m_UR, double Lanc_Error);
Matrix_COO Renormalize(lapack_complex_double* UL,lapack_complex_double* UR, Matrix_COO A, int m_UL, int m_UR, double Lanc_Error);
void Renormalize(Matrix_COO A, lapack_complex_double* UL,lapack_complex_double* UR, Matrix_COO & B, int m_UL, int m_UR, double Lanc_Error);
