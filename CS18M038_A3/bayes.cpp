#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

using namespace std;
int m=0,n=0, mtrain, mval, rowind=0, colind=0, mtest=0;
double **train1, **val1, **train2, **val2, **sigma1, **sigma2,**sigma1inv, **sigma2inv, *mu1, *mu2, accuracy1, accuracy2;
double det1, det2;
double **test;
//ofstream gammano1 ("gamma_no1.txt");
//ofstream gammano2 ("gamma_no2.txt");
ofstream op ("output.txt");
void printArr(double **M, int mm, int nn){
	for(int i =0; i<mm; i++){
		for(int j=0; j<nn; j++){
			cout<<M[i][j]<<" ";
		}
		cout<<'\n';
	}
}
void takeInput1(){
	ifstream path("dist_3/dist3_txt/class1_train.txt");
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
	mtrain = 0.8*m;
	mval = 0.2*m;
	train1 = new double*[mtrain];
	for(int i=0; i<mtrain; i++){
		train1[i] = new double[n];
	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train1[i][j] = raw_arr[zz++];
		}
	}

	val1 = new double*[mval];
	for(int i=0; i<mval; i++){
		val1[i] = new double[n];
	}

	for(int i=0; i<mval; i++){
		for(int j = 0; j<n; j++){
			val1[i][j] = raw_arr[zz++];
		}
	}


	delete[] raw_arr;
}
void takeInput2(){
	ifstream path("dist_3/dist3_txt/class2_train.txt");
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

	train2 = new double*[mtrain];
	for(int i=0; i<mtrain; i++){
		train2[i] = new double[n];
	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train2[i][j] = raw_arr[zz++];
		}
	}

	val2 = new double*[mval];
	for(int i=0; i<mval; i++){
		val2[i] = new double[n];
	}

	for(int i=0; i<mval; i++){
		for(int j = 0; j<n; j++){
			val2[i][j] = raw_arr[zz++];
		}
	}


	delete[] raw_arr;
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

    cout<<endl;
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
/*
void predictClassTrain(){
	double likely1, likely2;
	//cout<<"predict";
	int predict11=0, predict12=0, predict21=0, predict22=0;
	for(int i=0; i<mtrain; i++){
		likely1 = likelihood1(train1, i);
		likely2 = likelihood2(train1, i);

		if(likely1-likely2>0){

			predict11++;
			gammano1<<1<<' ';

		}
		else{
            predict12++;
            gammano1<<2<<' ';

		}
	}
    //cout<<predict11<<" "<<predict12<<" "<<predict21<<" "<<predict22<<endl;
	for(int i=0; i<mtrain  ; i++){
		likely1 = likelihood1(train2, i);
		likely2 = likelihood2(train2, i);
        //cout<<likely1<<" "<<likely2<<endl;
		if(likely1-likely2>0){

			predict21++;
			gammano1<<1<<' ';

		}
		else{

			predict22++;
			gammano1<<2<<' ';
		}
	}
	cout<<predict11<<" "<<predict12<<" "<<predict21<<" "<<predict22<<endl;
	gammano1<<"\n";
	accuracy1 = (double)(predict11+predict22)/(2*mtrain);
    gammano1<<accuracy1*100;
    gammano1<<"\n";
    //gammano1.close();
}

void predictClassVal(){
	double likely1, likely2;
	//cout<<"predict";
	int predict11=0, predict12=0, predict21=0, predict22=0;
	for(int i=0; i<mval; i++){
		likely1 = likelihood1(val1, i);
		likely2 = likelihood2(val1, i);
        //cout<<likely1<<" "<<likely2<<endl;
		if(likely1-likely2>0){
                gammano1<<1<<' ';
			predict11++;

		}
		else{
            gammano1<<2<<' ';
			predict12++;
		}
	}

	for(int i=0; i<mval; i++){
		likely1 = likelihood1(val2, i);
		likely2 = likelihood2(val2, i);
        //cout<<likely1<<" "<<likely2<<endl;
		if(likely1-likely2>0){
            gammano1<<1<<' ';
			predict21++;

		}
		else{
            gammano1<<2<<' ';
			predict22++;
		}
	}
	gammano1<<"\n";
	accuracy2 = (double)(predict11+predict22)/(2*mval);
    gammano1<<accuracy2*100;
    gammano1<<"\n";

}*/







void takeInput1Large(){
	ifstream path("dist3_txt/class1_train.txt");
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
	mtrain = 0.8*m;
	mval = 0.2*m;
	train1 = new double*[mtrain];
	for(int i=0; i<mtrain; i++){
		train1[i] = new double[n];
	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train1[i][j] = raw_arr[zz++];
		}
	}

	val1 = new double*[mval];
	for(int i=0; i<mval; i++){
		val1[i] = new double[n];
	}

	for(int i=0; i<mval; i++){
		for(int j = 0; j<n; j++){
			val1[i][j] = raw_arr[zz++];
		}
	}


	delete[] raw_arr;
}
void takeInput2Large(){
	ifstream path("dist3_txt/class2_train.txt");
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

	train2 = new double*[mtrain];
	for(int i=0; i<mtrain; i++){
		train2[i] = new double[n];
	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train2[i][j] = raw_arr[zz++];
		}
	}

	val2 = new double*[mval];
	for(int i=0; i<mval; i++){
		val2[i] = new double[n];
	}

	for(int i=0; i<mval; i++){
		for(int j = 0; j<n; j++){
			val2[i][j] = raw_arr[zz++];
		}
	}


	delete[] raw_arr;
}


int clearvar(){
m=0,n=0;
return 0;
}

void takeInputClass1(){
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
	train1 = new double*[mtrain];
	for(int i=0; i<mtrain; i++){
		train1[i] = new double[n];
	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train1[i][j] = raw_arr[zz++];
		}
	}

	delete[] raw_arr;
}


void takeInputClass2(){
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

    }

	train2 = new double*[mtrain];
	for(int i=0; i<mtrain; i++){
		train2[i] = new double[n];
	}

	for(int i=0; i<mtrain; i++){
		for(int j = 0; j<n; j++){
			train2[i][j] = raw_arr[zz++];
		}
	}

	delete[] raw_arr;
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

	test = new double*[mtest];
	for(int i=0; i<mtest; i++){
		test[i] = new double[n];
	}

	for(int i=0; i<mtest; i++){
		for(int j = 0; j<n; j++){
			test[i][j] = raw_arr[zz++];
		}
	}

	delete[] raw_arr;
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
    takeInputClass1();
	takeInputClass2();
    trainModel();
    findToy();
    string testfile = argv[1];
    takeTest(testfile);
    predictOutput();
    return 0;
}
