#include<time.h>
#include <fstream>
#include <cctype>
#include <math.h>
using namespace std;
clock_t start1, end1, start2, end2, start3, end3;
double extime1, extime2, extime3, fNorm, **A, **B, **C, m, n, p;
int size;

void addition(double **a, double **b, int size,double **c){
    int i,j;       
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            c[i][j] = a[i][j] + b[i][j];   
        }
    }
}

void subtraction(double **a,double **b,int size,double **c){
    int i,j;
    for(i=0;i<size;i++){
                for(j=0;j<size;j++){
                        c[i][j]= a[i][j] - b[i][j];
                }
        }
}

void  multiply(double **c,double **d,int size,int size2,double **C){
    if(size == 1){   
        C[0][0] = c[0][0] *d[0][0];   
    }
    else {
        int i,j;
        int tempsize =size/2;
        double **c11 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            c11[i]= new double[tempsize];
        }
        double **c12 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            c12[i]= new double[tempsize];
        }
        double **c21 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            c21[i]= new double[tempsize];
        }
        double **c22 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            c22[i]= new double[tempsize];
        }
        double **d11 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            d11[i]= new double[tempsize];
        }
        double **d12 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            d12[i]= new double[tempsize];
        }
        double **d21 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            d21[i]= new double[tempsize];
        }
        double **d22 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            d22[i]= new double[tempsize];
        }
        double **m1 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m1[i]= new double[tempsize];
        }
        double **m2 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m2[i]= new double[tempsize];
        }
        double **m3 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m3[i]= new double[tempsize];
        }
        double **m4 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m4[i]= new double[tempsize];
        }
        double **m5 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m5[i]= new double[tempsize];
        }
        double **m6 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m6[i]= new double[tempsize];
        }
        double **m7 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            m7[i]= new double[tempsize];
        }
        for(i=0;i<tempsize;i++){
            for(j=0;j<tempsize;j++){
                c11[i][j]=c[i][j];
                c12[i][j]=c[i][j+tempsize];
                c21[i][j]=c[i+tempsize][j];
                c22[i][j]=c[i+tempsize][j+tempsize];                   
                d11[i][j]=d[i][j];
                d12[i][j]=d[i][j+tempsize];
                d21[i][j]=d[i+tempsize][j];
                d22[i][j]=d[i+tempsize][j+tempsize];
            }
        }
        double **temp1 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp1[i]= new double[tempsize];
        }
        double **temp2 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp2[i]= new double[tempsize];
        }
       
        addition(c11,c22,tempsize,temp1);
        addition(d11,d22,tempsize,temp2);
        multiply(temp1,temp2,tempsize,size,m1);
        delete temp1;
        delete temp2;
       
        double **temp3 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp3[i]= new double[tempsize];
        }       
        addition(c21,c22,tempsize,temp3);
        multiply(temp3,d11,tempsize,size,m2);
        delete temp3;


        double **temp4 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp4[i]= new double[tempsize];
        }
        subtraction(d12,d22,tempsize,temp4);
        multiply(c11,temp4,tempsize,size,m3);
        delete temp4;


        double **temp5 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp5[i]= new double[tempsize];
        }
        subtraction(d21,d11,tempsize,temp5);
        multiply(c22,temp5,tempsize,size,m4);
        delete temp5;


        double **temp6 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp6[i]= new double[tempsize];
        }
        addition(c11,c12,tempsize,temp6);
        multiply(temp6,d22,tempsize,size,m5);
        delete temp6;

        double **temp7 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp7[i]= new double[tempsize];
        }   
        double **temp8 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp8[i]= new double[tempsize];
        }
        subtraction(c21,c11,tempsize,temp7);
        addition(d11,d12,tempsize,temp8);
        multiply(temp7,temp8,tempsize,size,m6);
        delete temp7;   
        delete temp8;

        double **temp9 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp9[i]= new double[tempsize];
        }   
        double **temp10 = new double*[tempsize];
        for(i=0;i<tempsize;i++){
            temp10[i]= new double[tempsize];
        }       
        subtraction(c12,c22,tempsize,temp9);
        addition(d21,d22,tempsize,temp10);
        multiply(temp9,temp10,tempsize,size,m7);
        delete temp9;
        delete temp10;   
       

        double **temporary1 = new double*[tempsize];
        double **temporary2 = new double*[tempsize];
        double **temporary3 = new double*[tempsize];
        double **temporary4 = new double*[tempsize];
        double **temporary5 = new double*[tempsize];
        double **temporary6 = new double*[tempsize];
        double **temporary7 = new double*[tempsize];
        double **temporary8 = new double*[tempsize];
        
        for(i=0;i<tempsize;i++){
            temporary1[i]= new double[tempsize];
			temporary2[i]= new double[tempsize];
			temporary3[i]= new double[tempsize];
			temporary4[i]= new double[tempsize];
			temporary5[i]= new double[tempsize];
			temporary6[i]= new double[tempsize];
			temporary7[i]= new double[tempsize];
			temporary8[i]= new double[tempsize];
			
        }
        
        
        addition(m1,m7,tempsize,temporary1);
        subtraction(m4,m5,tempsize,temporary2);
        addition(temporary1,temporary2,tempsize,temporary3);    //c11
           
        addition(m3,m5,tempsize,temporary4);//c12   
        addition(m2,m4,tempsize,temporary5);//c21
       
        addition(m3,m6,tempsize,temporary6);
        subtraction(m1,m2,tempsize,temporary7);
       
        addition(temporary6,temporary7,tempsize,temporary8);//c22
       
        int a=0;
        int b=0;
        int c=0;   
        int d=0;
        int e=0;
        int tempsize2= 2*tempsize;
        
        for(i=0;i<tempsize2;i++){
            for(j=0;j<tempsize2;j++){
                if(j>=0 && j<tempsize && i>=0 && i<tempsize){
                    C[i][j] = temporary3[i][j];
                }
                if(j>=tempsize && j<tempsize2 && i>=0 && i<tempsize){
                    a=j-tempsize;
                    C[i][j] = temporary4[i][a];
                }
                if(j>=0 && j<tempsize && i>= tempsize && i < tempsize2){
                    c=i-tempsize;
                    C[i][j] = temporary5[c][j];
                }
                if(j>=tempsize && j< tempsize2 && i>= tempsize && i< tempsize2 ){
                    d=i-tempsize;
                    e=j-tempsize;
                    C[i][j] =temporary8[d][e];
                }
            }   
        }
      
 	for(int i=0; i<tempsize; i++){
		delete[] (m1[i], m2[i], m2[i], m3[i], m4[i], m5[i], m6[i], m7[i], temporary1[i], temporary2[i], temporary3[i], temporary4[i], temporary5[i], temporary6[i],  
     temporary7[i], temporary8[i], c11[i], c12[i], c21[i], c22[i], d11[i], d12[i], d21[i], d22[i]);
		
	}
	
    delete [] (m1, m2, m2, m3, m4, m5, m6, m7, temporary1, temporary2, temporary3, temporary4, temporary5, temporary6,  
     temporary7, temporary8, c11, c12, c21, c22, d11, d12, d21, d22);   
    }   
}
void resetM(){
	for(int i=0; i<size; i++){
		delete[] C[i];
		
	}
	delete [] C;
	C = new double*[size];
	for(int i=0; i<size; i++){
		C[i] = new double[size]();
	}
	
}

void fileInputA(ifstream &cin){

	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			cin>>A[i][j];
		}
	}
	
}
void fileInputB(ifstream &cin){

	for(int i=0; i<n; i++){
		for(int j=0; j<p; j++){
			cin>>B[i][j];
		}
	}
}
void fileOutput(ofstream &cout){
	for(int i=0; i<m; i++){
		for(int j=0; j<p; j++){
			cout<<C[i][j]<<" ";
		}	
		cout<<endl;
	}
	cout<<endl;
}
void fileOutputVar(ofstream &cout){
	cout<<extime1;
	cout<<endl;
	cout<<extime2;
	cout<<endl;
	cout<<extime3;
	cout<<endl;
	cout<<fNorm;
	cout<<endl;
}
void strassenMM(){

    multiply(A, B,size,size,C);
	
}
void iterativeMM(){
	for(int i=0; i<m; i++){
		for(int j=0; j<p; j++){
			for(int k=0; k<n; k++){
			C[i][j] += A[i][k]*B[k][j];
			}
		}	
	}
}
void multiplyPartition(int Ai, int Aj, int Bi, int Bj, int Ci, int Cj){
	for(int i=Ai; i<Aj; i++){
		for(int j=Ci; j<Cj; j++){
			for(int k=Bi; k<Bj; k++){
			C[i][j] += A[i][k]*B[k][j];
			}
		}	
	}
}
void partitionMM(){
	int m1, m2, n1, n2, p1, p2;
	m1 = m/2;
	m2 = m - m/2;
	n1 = n/2;
	n2 = n - n/2;
	p1 = p/2;
	p2 = p - p/2;
	multiplyPartition(0, m1, 0, n1, 0, p1);
	multiplyPartition(0, m1, n1, n1+n2, 0, p1);
	multiplyPartition(0, m1, 0, n1, p1, p1+p2);
	multiplyPartition(0, m1, n1, n1+n2, p1, p1+p2);
	
	multiplyPartition(m1, m1+m2, 0, n1, 0, p1);
	multiplyPartition(m1, m1+m2, n1, n1+n2, 0, p1);
	multiplyPartition(m1, m1+m2, 0, n1, p1, p1+p2);
	multiplyPartition(m1, m1+m2, n1, n1+n2, p1, p1+p2);
	
}

double findFrobenius()
{
	ifstream cinA("matC_N.txt"), cinB("matC_RP.txt");
	
	double sqsum = 0.0;
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<p;j++)
		{
			double x , y ;
			 cinA>>x;
			 cinB>>y;
			 
			 x=x-y ;
			sqsum +=(x*x); 
		}
	}
	sqsum = sqrt(sqsum);

	cinA.close();
	cinB.close();
	
	return sqsum;			
}

int main(){
	int bitm = 1;
	
    
    ifstream cinA("matA.txt"), cinB("matB.txt");
	ofstream cout1("matC_N.txt"),cout2("matC_P.txt"), cout3("matC_RP.txt"), cout4("output.txt");
	cinA>>m>>n;
	cinB>>n>>p;
	size = (m>(n>p?n:p)?m:(n>p?n:p));
    while(bitm<size){
    	bitm<<=1;
	}
    size = bitm;
	A = new double*[size];
	B = new double*[size];
	C = new double*[size];		
	for(int i=0; i<size; i++){
		A[i] = new double[size]();
	B[i] = new double[size]();
	C[i] = new double[size]();
	}
	
  	fileInputA(cinA);
	fileInputB(cinB);
	
	start1 = clock();
  	iterativeMM();
  	end1 = clock();
  	
	fileOutput(cout1);
	
	resetM();
	
	start2 = clock();
	partitionMM();
	end2 = clock();
	
	fileOutput(cout2);
	
	start3 = clock();
	strassenMM();
	end3 = clock();
	
	fileOutput(cout3);
	
	extime1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
	extime2 = ((double) (end2 - start2)) / CLOCKS_PER_SEC;
	extime3 = ((double) (end3 - start3)) / CLOCKS_PER_SEC;
	fNorm = findFrobenius();
	fileOutputVar(cout4);
	cinA.close();
	cinB.close();
   	cout1.close();
   	cout2.close();
   	cout3.close();
   	cout4.close();
	return 0;
}



















