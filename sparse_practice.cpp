#include <iostream> 
//#include "mkl.h"
//#include "mkl_types.h"
//#include "mkl_spblas.h" //for sparse matrices subroutines.
//#include <vector>
#include "functions.h"
#include "dmrg_Hubbard_model.h"


//*********************************************************//

/*		TO COMPILE THIS PROGRAM USE FOLLOWING
		
		g++ sparse_practice.cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core               */
		
	 
//*********************************************************//


using namespace std;


int main() {
	
int i;

double input_vec[4] = {11.0,1.0,9.0,1.0};
double output_vec[4]; 





Single_Site site;

site.Site_construction();


DOUBLE_MKL_CSR_MAT mat2 = CSR_MAT_TO_MKL_SPARSE(site.C_up_dag);
DOUBLE_MKL_CSR_MAT mat3 = CSR_MAT_TO_MKL_SPARSE(site.C_up);





char transa = 'n';
MKL_INT m=4;
cout<<"fine till here"<<endl;
mkl_cspblas_dcsrgemv(&transa, &m, mat2.value, mat2.row_ind, mat2.columns, input_vec, output_vec);

for (i=0;i<4;i++){
	cout<<output_vec[i]<<endl;
	
	}

return 0;


}
