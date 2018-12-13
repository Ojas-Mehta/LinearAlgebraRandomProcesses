#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <math.h>

using namespace std;
int m=0,n=0, mtrain, mval, rowind=0, colind=0, mtest=0;
double **test_long, **train_long,**val_long,**train1, **train2, **val1, **val2, **sigma1, **sigma2,**sigma1inv, **sigma2inv, *mu1, *mu2, accuracy1, accuracy2, *mu_long, **sigma_long;
double det1, det2;
double **test;
ofstream op ("output.txt");
void printArr(double **M, int mm, int nn){
	for(int i =0; i<mm; i++){
		for(int j=0; j<nn; j++){
			cout<<M[i][j]<<" ";
		}
		cout<<'\n';
	}
}

void calMean(double **arr, int row, double *tempmean){

	for(int j=0; j<n; j++){
		double sum=0;
		for(int i=0; i<row; i++){
		   sum+=arr[i][j];
		}
		tempmean[j] = sum/(row-1);
	}
}


void calSigma1(double **arr1, int row, double **tempsigma){

	for(int i =0; i<n; i++){
		for(int k =0; k<n; k++){
			double sum=0;
			for(int j =0; j<row; j++){
			sum+=(arr1[j][i]-mu1[i])*(arr1[j][k]-mu1[i]);
			}
			tempsigma[i][k] = sum/(row-1);
			//cout<<sum<<endl;
		}
	}

}
void calSigma2(double **arr1, int row, double **tempsigma){

	for(int i =0; i<n; i++){
		for(int k =0; k<n; k++){
			double sum=0;
			for(int j =0; j<row; j++){
			sum+=(arr1[j][i]-mu2[i])*(arr1[j][k]-mu2[i]);
			}
			tempsigma[i][k] = sum/(row-1);
			//cout<<sum<<endl;
		}
	}

}
void para1(){

	mu1 = new double[n];
	sigma1 = new double*[n];
	for(int i=0; i<n; i++){
		sigma1[i] = new double[n];
	}
	calMean(train1, mtrain, mu1);
	calSigma1(train1, mtrain, sigma1);
}
void para2(){
	mu2 = new double[n];
	sigma2 = new double*[n];
	for(int i=0; i<n; i++){
		sigma2[i] = new double[n];
	}
	calMean(train2, mtrain, mu2);
	calSigma2(train2, mtrain, sigma2);
}
void trainModel(){
	para1();
	//cout<<"j";
	para2();
}

void findToy(){
	int i, j, extra;
	sigma1inv = new double*[n];
	for(int i=0; i<n; i++){
		sigma1inv[i] = new double[n];
	}
	sigma2inv = new double*[n];
	for(int i=0; i<n; i++){
		sigma2inv[i] = new double[n];
	}




  	gsl_matrix * malt1 = gsl_matrix_alloc (n,n);
  	gsl_permutation *perm = gsl_permutation_alloc(n);

	for (i = 0; i < n; i++)
	   for (j = 0; j < n; j++){
	   	gsl_matrix_set (malt1, i, j, sigma1[i][j]);
	   }

	gsl_linalg_LU_decomp(malt1, perm, &extra);
	//cout<<gsl_matrix_get(malt1, 1, 1);
	gsl_matrix *inv = gsl_matrix_alloc(n, n);

	det1 = gsl_linalg_LU_det (malt1, extra);

    gsl_linalg_LU_invert(malt1, perm, inv);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sigma1inv[i][j] = gsl_matrix_get(inv, i, j);
            //printf("%f ", element);
        }
        //printf("\n");
    }
    // printArr(sigma1inv, n,n);
    //malt1 will be changed
    for (i = 0; i < n; i++)
	   for (j = 0; j < n; j++){
	   	gsl_matrix_set (malt1, i, j, sigma2[i][j]);
	}

   // cout<<endl;
	//for second
	gsl_linalg_LU_decomp(malt1, perm, &extra);

	det2 = gsl_linalg_LU_det (malt1, extra);
    gsl_linalg_LU_invert(malt1, perm, inv);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            sigma2inv[i][j] = gsl_matrix_get(inv, i, j);
            //printf("%f ", element);
        }
        //printf("\n");
    }

    //printArr(sigma2inv, n,n);

    gsl_matrix_free(inv);
    gsl_permutation_free(perm);


}
double likelihood1(double **arr1, int index){
	double sum=0, xx,yy, left, right;
	for(int i =0; i<n; i++){
		xx = arr1[index][i];
		left = xx-mu1[i];

		for(int k =0; k<n; k++){
			yy = arr1[index][k];
			right = yy-mu1[k];
			sum+=left*sigma1inv[i][k]*right;

		}

	}
	sum = exp(-sum);
	return sum/det1;
}

double likelihood2(double **arr1, int index){
	double sum=0, xx,yy, left, right;
	for(int i =0; i<n; i++){
		xx = arr1[index][i];
		left = xx-mu2[i];

		for(int k =0; k<n; k++){
			yy = arr1[index][k];
			right = yy-mu2[k];
			sum+=left*sigma2inv[i][k]*right;
         //cout<<sum<<endl;
		}

		//cout<<sum<<endl;
	}
	sum = exp(-sum);
	return sum/det2;
}

void takeInput1Large(){
	ifstream path("class1_train.txt");
    string line;
    double *raw_arr = new double[70000*50];
	int ii=0, zz=0;
    while(getline(path, line))
    {
        stringstream ss(line);
        double number;
        while(ss >> number)
            raw_arr[ii++] = number;
        m++;
        //cout<<m;
    }

    if(m==70000){
    	n = 50;
	}
	else{
		n=10;
	}
	///cout<<m<<endl<<m;
	mtrain = m;

	train_long = new double*[2*mtrain];
	for(int i=0; i<2*mtrain; i++){
		train_long[i] = new double[n];

	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train_long[i][j] = raw_arr[zz++];
		}
	}

	delete[] raw_arr;
}

void takeInput2Large(){
	ifstream path("class2_train.txt");
    string line;
    double *raw_arr = new double[70000*50];
	int ii=0, zz=0;
    while(getline(path, line))
    {
        stringstream ss(line);
        double number;
        while(ss >> number)
            raw_arr[ii++] = number;
        m++;
        //cout<<m;
    }



	for(int i=mtrain; i<2*mtrain; i++){
		for(int j = 0; j<n; j++){
			train_long[i][j] = raw_arr[zz++];
		}
	}

	delete[] raw_arr;
}



void calSigma(double **arr1, int row, double **tempsigma){

	for(int i =0; i<n; i++){
		for(int k =0; k<n; k++){
			double sum=0;
			for(int j =0; j<row; j++){
			sum+=(arr1[j][i]-mu_long[i])*(arr1[j][k]-mu_long[i]);
			}
			tempsigma[i][k] = sum/(row-1);
			//cout<<sum<<endl;
		}
	}

}
void takeTest(string fname){
    ifstream path(fname);
    string line;
    double *raw_arr = new double[200000*50];
	int ii=0, zz=0;
    while(getline(path, line))
    {
        stringstream ss(line);
        double number;
        while(ss >> number)
            raw_arr[ii++] = number;
        mtest++;
    }

	test_long = new double*[mtest];
	for(int i=0; i<mtest; i++){
		test_long[i] = new double[n];
	}

	for(int i=0; i<mtest; i++){
		for(int j = 0; j<n; j++){
			test_long[i][j] = raw_arr[zz++];
		}
	}

	delete[] raw_arr;
}
void reduceData(){
    mu_long = new double[n];
	sigma_long = new double*[n];
	for(int i=0; i<n; i++){
		sigma_long[i] = new double[n];
	}
	calMean(train_long, 2*mtrain, mu_long);
	calSigma(train_long, 2*mtrain, sigma_long);

	gsl_matrix * gsl_sigma = gsl_matrix_alloc (n,n);
	for (int i = 0; i < n; i++)
	   for (int j = 0; j < n; j++){
	   	gsl_matrix_set (gsl_sigma, i, j, sigma_long[i][j]);
	   }


	gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv (gsl_sigma, eval, evec, w);
    gsl_eigen_symmv_free (w);

    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    double top_val[2][n];
	for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < n; ++j) {
             top_val[i][j] = gsl_matrix_get(evec, j, i);
            //printf("%f ", element);
        }
    }

    train1 = new double*[mtrain];
    train2 = new double*[mtrain];
    test = new double*[mtest];
	for(int i=0; i<mtrain; i++){
		train1[i] = new double[2];
		train2[i] = new double[2];
		test[i] = new double[2];

	}


	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<2; j++){
                for(int k=0; k<n; k++){
                    train1[i][j] += (train_long[i][k]-mu_long[k])*top_val[j][k];

                }

		}

	}
	for(int i=0; i<2; i++){
		for(int j = 0; j<mtrain; j++){
                for(int k=0; k<n; k++)
                    train2[j][i] += top_val[i][k]*(train_long[j+mtrain][k]-mu_long[k]);

		}

	}
	for(int i=0; i<mtest; i++){
		for(int j = 0; j<2; j++){
                for(int k=0; k<n; k++){
                    test[i][j] += (test_long[i][k]-mu_long[k])*top_val[j][k];

                }

		}

	}

    n=2;

}


void predictOutput(){
	double likely1, likely2;

	for(int i=0; i<mtest; i++){
		likely1 = likelihood1(test, i);
		likely2 = likelihood2(test, i);

		if(likely1-likely2>0){
                op<<1<<' ';
		}
		else{
            op<<2<<' ';
		}
	}

}

int main(int args, char **argv)
{
    takeInput1Large();
	takeInput2Large();
    string testfile = argv[1];
    takeTest(testfile);
	reduceData();
    trainModel();
    findToy();

    predictOutput();
    //predictClassTrain();
    return 0;
}

