#include<iostream> 
#include<sstream>
#include<time.h>
#include <fstream>
using namespace std; 

clock_t start1, end1, start2, end2, start3, end3;
double extime1, extime2, extime3, fNorm, **A, **B, **C, m, n, p;
string finex = "";
long long **bracket;

// Function for printing the optimal 
// parenthesization of a matrix chain product 
string printAscii(long long n){
	ostringstream ss;
	ss<<n;
	return ss.str();
}
void optimalSubroutine(long long i, long long j, long long n, 
                      long long &name) 
{ 
    // If only one matrix left in current segment 
    
    if (i == j) 
    { 
        finex+=printAscii(++name); 
        return; 
    } 
  
    finex+= "("; 
  
    optimalSubroutine(i, bracket[i][j], n, 
                     name); 
    finex+= ","; 
    optimalSubroutine(bracket[i][j] + 1, j, 
                     n, name);
   
    finex+= ")";  
} 
  
long long optimalMatrixPairing(long long p[], long long n) 
{ 
   
    long long m[n][n]; 
  
  
    for (long long i=1; i<n; i++) 
        m[i][i] = 0; 
  
     
    for (long long L=2; L<n; L++) 
    { 
        for (long long i=1; i<n-L+1; i++) 
        { 
            long long j = i+L-1; 
            m[i][j] = LLONG_MAX; 
            for (long long k=i; k<=j-1; k++) 
            { 
                
                long long q = m[i][k] + m[k+1][j] + p[i-1]*p[k]*p[j]; 
                if (q < m[i][j]) 
                { 
                    m[i][j] = q; 
                    bracket[i][j] = k; 
                } 
            } 
        } 
    } 
  
    long long name = 0; 
  
    
    optimalSubroutine(1, n-1, n, name);
	return m[1][n-1]; 
    //coutfile << "nOptimal Cost is : " << m[1][n-1]<<endl<<endl; 
} 

long long nScalarMul(long long arr[], long long n){
	long long ans=0;
	for(long long i=1; i<n-1; i++){
		ans+=arr[0]*arr[i]*arr[i+1];
	}
	return ans;
}

long long lrpair(long long arr[], long long n){
	
	long long ans=0;
	if(n==1)return 0;
	long long ind[n];
	for(long long i=0; i<n; i++){
		ind[i] = arr[i];
	}
	long long num = n,single = 0, odd;
	
	for(long long j=0;j<n/2+1; j+=1){
		long long i,z=0;
		for(i=0; i<num-2;i+=2){
			ans+=ind[i]*ind[i+1]*ind[i+2];
			ind[z++]=ind[i];
			}
		ind[z++]=ind[i++];
		if(num%2==0)
		ind[z++]=ind[i];
		num=z;
	}
	return ans;
}



long long time_nScalarMul(long long arr[], long long n){
	long long ans=0;
	long long** dyn1,**dyn2, **dyn3;
	
	for(long long i=1; i<n-1; i++){
		dyn1 = new long long*[arr[0]];
		dyn2 = new long long*[arr[i]];
		dyn3 = new long long*[arr[0]];
		for(long long j=0; j<arr[0]; j++){
			dyn1[j] = new long long[arr[i]]();
			dyn3[j] = new long long[arr[i+1]]();
		}
		for(long long j=0; j<arr[i]; j++){
		
			dyn2[j] = new long long[arr[i+1]]();
		}
		
		for(long long ii=0; ii<arr[0]; ii++){
		for(long long j=0; j<arr[i+1]; j++){
			for(long long k=0; k<arr[i]; k++){
				dyn3[ii][j] += dyn1[ii][k]*dyn2[k][j];
				}
			}	
		}
		ans+=arr[0]*arr[i]*arr[i+1];
		
		for(long long j=0; j<arr[i]; j++){
		
			delete[] dyn2[j];
			
		}
		for(long long j=0; j<arr[0]; j++){
			delete[] dyn3[j];
			delete[] dyn1[j];
		}
		delete[] dyn1;
		delete[] dyn2;
		delete[] dyn3;
		
	}
	return ans;
}

long long time_lrpair(long long arr[], long long n){	
	long long ans=0;
	if(n==1)return 0;
	long long ind[n];
	for(long long i=0; i<n; i++){
		ind[i] = arr[i];
	}
	long long num = n,single = 0, odd;
	
	long long** dyn1,**dyn2, **dyn3;
	for(long long j=0;j<n/2+1; j+=1){
		long long i,z=0;
		for(i=0; i<num-2;i+=2){
			
		dyn1 = new long long*[ind[i]];
		dyn2 = new long long*[ind[i+1]];
		dyn3 = new long long*[ind[i]];
		for(long long j=0; j<ind[i]; j++){
			dyn1[j] = new long long[ind[i+1]]();
			dyn3[j] = new long long[ind[i+2]]();
		}
		for(long long j=0; j<ind[i+1]; j++){
			dyn2[j] = new long long[ind[i+2]]();
		}
		
		for(long long ii=0; ii<ind[i]; ii++){
		for(long long j=0; j<ind[i+2]; j++){
			for(long long k=0; k<ind[i+1]; k++){
				dyn3[ii][j] += dyn1[ii][k]*dyn2[k][j];
				}
			}	
		}
		
			//ans+=ind[i]*ind[i+1]*ind[i+2];
			ind[z++]=ind[i];
			
		for(long long j=0; j<ind[i+1]; j++){
			delete[] dyn2[j];
		}
		for(long long j=0; j<ind[i]; j++){
			delete[] dyn3[j];
			delete[] dyn1[j];
		}
		delete[] dyn1;
		delete[] dyn2;
		delete[] dyn3;
		}
		ind[z++]=ind[i++];
		if(num%2==0)
		ind[z++]=ind[i];
		num=z;
	}
	return ans;
}
 
void timeOptHelper(long long i, long long j, long long n, long long arr[]) 
{ 
    // If only one matrix left in current segment 
   
    if (i == j) 
        return; 
    
    timeOptHelper(i, bracket[i][j], n, arr);
	
    timeOptHelper(bracket[i][j] + 1, j, n, arr);
    
     
	long long** dyn1,**dyn2, **dyn3;
	
    dyn1 = new long long*[arr[i-1]];
		dyn2 = new long long*[arr[i]];
		dyn3 = new long long*[arr[i-1]];
		for(long long jj=0; jj<arr[i-1]; jj++){
			dyn1[jj] = new long long[arr[i]]();
			dyn3[jj] = new long long[arr[j]]();
		}
		for(long long jj=0; jj<arr[i]; jj++){
			dyn2[jj] = new long long[arr[j]]();
		}
	
		for(long long ii=0; ii<arr[i-1]; ii++){
		for(long long jj=0; jj<arr[j]; jj++){
			for(long long k=0; k<arr[i]; k++){
				dyn3[ii][jj] += dyn1[ii][k]*dyn2[k][jj];
				}
			}	
		}
		
		for(long long jj=0; jj<arr[i]; jj++){
			delete[] dyn2[jj];
		}
		for(long long jj=0; jj<arr[i-1]; jj++){
			delete[] dyn3[jj];
			delete[] dyn1[jj];
		}
		delete[] dyn1;
		delete[] dyn2;
		delete[] dyn3;
	
    
    //cout<<i<<" "<<j<<endl; 
} 

void timeOpt(long long arr[], long long n){
	timeOptHelper(1, n-1, n, arr);
}
 
int main() 
{   ifstream cinfile("input.txt");
	ofstream coutfile("output.txt");
	long long n;
     cinfile>>n;
     n++;
     long long arr[n],countOpt;
    for(int i=0; i<n; i++){
     cinfile>>arr[i];	
	}
    bracket = new long long*[n];
    for(long long i=0; i<n; i++){
    	bracket[i] = new long long[n]();
	}
	
    countOpt = optimalMatrixPairing(arr, n); 
    coutfile<<finex<<endl;
    
//    //cout<<endl
    coutfile<<nScalarMul(arr, n)<<endl;
    coutfile<<lrpair(arr,n)<<endl;
    coutfile<<countOpt<<endl;
    
    start1 = clock();
    timeOpt(arr, n);
    end1 = clock();
    
    extime1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
   
    start2 = clock();
    time_nScalarMul(arr,n);
    end2 = clock();
    extime2 = ((double) (end2 - start2)) / CLOCKS_PER_SEC;
    
    start3 = clock();
    time_lrpair(arr, n);
    end3 = clock();
    extime3 = ((double) (end3 - start3)) / CLOCKS_PER_SEC;
    
	
    coutfile<<extime2<<endl;//time for left to right pairing
    coutfile<<extime3<<endl;//time for pairing
    coutfile<<extime1<<endl;//time for optimal matrix multiplication
    
	cinfile.close();
    coutfile.close();
    return 0; 
} 
