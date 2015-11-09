#include <iostream>
#include <math.h>
using std::cerr;
using namespace std;
#include <fstream>
using std::ofstream;
#include <cstdlib>       // for exit function
#include <ctime>


// parameter definitions
#define N 1024    // system size
#define p_range 100000   // no. of p values in range 0-1      
#define T 10000      // No. of time steps
#define rep 1



//====== function to initialize the matrix ===============//
// This function creates a matrix with entries 0's and 1's. Here x is the proportion of 1's

void create_random_matrix(int*U,float x) {  

	int i1, i2;
	
	for (i1 = 0; i1 < N; i1++) {
		for (i2 = 0; i2 < N; i2++) {
			float number = rand()/(float)RAND_MAX;
			if (number<x) U[i1*N+i2]=1;
			else U[i1*N+i2]=0;
		}
	}

}


//=======================average density=========================//
// Calculates the proportion of 1's in the given matrix
float density(int* A){
	int i,j;
	float sum=0; 
	for(i=0;i<N;++i){
		for(j=0;j<N;++j){
			sum+= A[i*N+j];
		}
	}
	float average= sum/(N*N);
return average;
}

//=========================average along rows======================//
// If you have a matrix with each row as different realization, this function will give the average over all the realizations.

void average_along_row(float* A,float*A1){
	int l,k;
	float sum=0;
	
	for(k=0;k<p_range;++k){
		sum=0;
		for(l=0;l<rep;++l){
			sum+=A[l*p_range+k];
		}
		A1[k]=sum/rep;
	}

}	

//====================select random site==============================//
void select_neighbor_of_site(int i ,int j , int*neighbor){
	int left,right,bottom,top ;
	int in = i ,jn = j;

	float test = rand()/(float)RAND_MAX ;
	
	if (i==0) top=N-1;
	else top = i-1;
	
	if (i==N-1) bottom = 0;
	else bottom = i+1 ;
	
	if (j==0) left = N-1;
	else left = j-1 ;
	
	if (j==N-1) right = 0;
	else right = j+1;
				
	if (test <= 0.25) in = top;
	else if ( test <= 0.5) in = bottom;
	else if ( test <=0.75) jn =left;
	else jn = right;
	
	neighbor[0] = in;
	neighbor[1] = jn;
	
}

//========select neighbor of pair=====================================//
int select_neighbor_of_pair(int in ,int jn, int i, int j){
	int left,right,top,bottom,leftn,rightn,topn,bottomn , neighbor_of_pair;
	
	if (i==0) top=N-1;                 //periodic boundary
	else top = i-1;
	
	if (i==N-1) bottom = 0;
	else bottom = i+1 ;
	
	if (j==0) left = N-1;
	else left = j-1 ;
	
	if (j==N-1) right = 0;
	else right = j+1;
	
	if (in==0) topn=N-1;
	else topn = in-1;
	
	if (in==N-1) bottomn = 0;
	else bottomn = in+1 ;
	
	if (jn==0) leftn = N-1;
	else leftn = jn-1 ;
	
	if (jn==N-1) rightn = 0;
	else rightn = jn+1;
	
	int nn[6] ,c=0;
	
	if ((top*N +j) != (in*N+jn)) {
		nn[c]=top*N + j;
		c+=1;
	}
	if ((bottom*N + j) != (in*N+jn)) {
		nn[c]=bottom*N  + j;
		c+=1;
	}
	if ((i*N +right) != (in*N+jn)) {
		nn[c]= i*N + right;
		c+=1;
	}
	if ((i*N +left) != (in*N+jn)) {
		nn[c] = i*N + left;
		c+=1;
	}
	if ((topn*N +jn) != (i*N+j)) {
		nn[c]=topn*N + jn;
		c+=1;
	}
	if ((bottomn*N +jn) != (i*N+j)) {
		nn[c]=bottomn*N  + jn;
		c+=1;
	}
	if ((in*N +rightn) != (i*N+j)) {
		nn[c]= in*N + rightn;
		c+=1;
	}
	if ((in*N +leftn) != (i*N+j)) {
		nn[c] = in*N + leftn;
		c+=1;
	}
	
	float test =rand()/(float)RAND_MAX ;
	
	if (test <=(0.1666)) neighbor_of_pair= nn[0];
	else if ( test <= (2*0.1666)) neighbor_of_pair= nn[1];
	else if ( test <= (3*0.1666)) neighbor_of_pair= nn[2];
	else if ( test <= (4*0.1666)) neighbor_of_pair= nn[3];
	else if ( test <= (5*0.1666)) neighbor_of_pair= nn[4];
	else neighbor_of_pair = nn[5];
	

return neighbor_of_pair;

} 

//-------------------coarse grained variance-------------------//
float cg_var (int* A, int n, float mean) {				//n is the dimension of coarse graining i.e. dimension of the submatrix
	int i,j,k,l,count = 0;
	float reduced_mean[N*N/(n*n)] ;
	for (i=0; i<N-n+1; i+=n) {
		for (j=0; j<N-n+1; j+=n) {	
			float sum = 0;
			for (k=0; k<n; ++k) {
				for (l=0; l<n ; ++l) {
					sum+= (A[(i+k)*N+(j+l)]);
				}
			}
			
			reduced_mean[count] = sum/(float)(n*n);         //mean of the sub-matrix
			count+=1;
		}
	}
	float sum_sq=0;
	for (i=0; i<N*N/(n*n); i++) {
		sum_sq+= pow((reduced_mean[i]-mean),2.0);
	}
	float variance = sum_sq/(float)(N*N/(n*n)); 
	//float variance = mean_sq - (mean*mean) ;
return variance;
}

//-------------------coarse grained skewness-------------------//
float cg_skew (int* A, int n, float mean, float variance) {				//n is the dimension of coarse graining i.e. dimension of the submatrix
	int i,j,k,l,count = 0;
	float reduced_mean[N*N/(n*n)] ;
	float skewness;
	for (i=0; i<N-n+1; i+=n) {
		for (j=0; j<N-n+1; j+=n) {	
			float sum = 0;
			for (k=0; k<n; ++k) {
				for (l=0; l<n ; ++l) {
					
					sum+= A[(i+k)*N+(j+l)] ;
				}
			}
			
			reduced_mean[count] = sum/(float)(n*n);
			count+=1;
		}
	}
	
	float sum_cube=0;
	for (i=0; i<N*N/n/n; i++) {
		sum_cube+= reduced_mean[i]*reduced_mean[i]*reduced_mean[i];
	}
	float mean_cube = sum_cube/(float)(N*N/(n*n));
	float pwvar = pow(variance,1.5f);
	if (pwvar == 0) skewness = 0;
	else skewness = ( mean_cube - 3.0f*mean*variance - mean*mean*mean )/ pwvar;
	
return skewness;
}
			
		
////////////// main function //////////////////////////////////
int main(){
	
	srand(time(NULL));
	int x,l,t,i,j,z;
	float p[p_range], q=0.88;
	float skewness1, skewness2, skewness4, skewness8, variance1, variance2, variance4, variance8 , mean ,skewness_random1, skewness_random2, skewness_random4, skewness_random8, variance_random1, variance_random2, variance_random4, variance_random8;
	
	int* neighbor = new int[2];
	
	for(i=0;i<p_range;++i) {
		p[i]= i/(float)p_range;
	}
	 
	int*A = new int[N*N];
	int*B = new int[N*N];
	
	create_random_matrix(A,0.333);
	//float init = density(A);
	//cout<<init<<endl;
	for(x=32153;x>31499;x--){
			
			for(t=0;t<T;t++){
		
				for(z=0;z<N*N ; ++z){                // so that each site gets selected once on an average
				
					i = rand()%N;           // selecting one random site
					j = rand()%N;
					
					float test = rand()/(float)RAND_MAX;
					float test1 = rand()/(float)RAND_MAX;
					
					if (A[i*N+j]==1){     //if the site is occupied
						
						select_neighbor_of_site(i, j ,neighbor);    //look for a neighbor
						int in = neighbor[0] , jn = neighbor[1];
						
						if (A[in*N +jn]==0) {                     //if neighbor is empty
							if (test < p[x]) 
								A[in*N+jn]=1;                 //regular cp
							else A[i*N+j]=0;
						} 

						else {							  
							if (test < q){
								
								int neighbor_of_pair=select_neighbor_of_pair (in, jn, i, j);  //look for the neighbor of pair 
								A[neighbor_of_pair]=1;
							}
							else if (test1 < 1-p[x] )
								A[i*N+j]=0;	
						}
					}
					

				}	
				
			}
			mean = density(A);   
			
			
			cout<< p[x]<<"\t"<<mean<<endl;
			
			
			variance1  = cg_var (A,1,mean) ;
			variance2  = cg_var (A,2,mean) ;
			variance4  = cg_var (A,4,mean) ;  		
			variance8  = cg_var (A,8,mean) ;
			skewness1 = cg_skew (A,1,mean,variance1) ;
			skewness2 = cg_skew (A,2,mean,variance2) ;
			skewness4 = cg_skew (A,4,mean,variance4) ;
			skewness8 = cg_skew (A,8,mean,variance8) ;
			
			 
			
			create_random_matrix(B,mean);
			variance_random1  = cg_var (B,1,mean) ;
			variance_random2  = cg_var (B,2,mean) ;
			variance_random4  = cg_var (B,4,mean) ;  		
			variance_random8  = cg_var (B,8,mean) ;
			skewness_random1 = cg_skew (B,1,mean,variance_random1) ;
			skewness_random2 = cg_skew (B,2,mean,variance_random2) ;
			skewness_random4 = cg_skew (B,4,mean,variance_random4) ;
			skewness_random8 = cg_skew (B,8,mean,variance_random8) ;

	
	//====================saving data in a file======================//
			ofstream spatial_structure_fout;
			ofstream mean_cover_fout;
			
			ofstream variance1_fout;
			ofstream variance2_fout;
			ofstream variance4_fout;
			ofstream variance8_fout;
			ofstream variance_random1_fout;
			ofstream variance_random2_fout;
			ofstream variance_random4_fout;
			ofstream variance_random8_fout;
			
			ofstream skewness1_fout;
			ofstream skewness2_fout;
			ofstream skewness4_fout;
			ofstream skewness8_fout;
			ofstream skewness_random1_fout;
			ofstream skewness_random2_fout;
			ofstream skewness_random4_fout;
			ofstream skewness_random8_fout;
			
			spatial_structure_fout.open("tcp_snapshots_1024*1024_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			mean_cover_fout.open("tcp_mean_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			
			variance1_fout.open("tcp_variance_cg1_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			variance2_fout.open("tcp_variance_cg2_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			variance4_fout.open("tcp_variance_cg4_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			variance8_fout.open("tcp_variance_cg8_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			variance_random1_fout.open("tcp_variance_random_cg1_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			variance_random2_fout.open("tcp_variance_random_cg2_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			variance_random4_fout.open("tcp_variance_random_cg4_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			variance_random8_fout.open("tcp_variance_random_cg8_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			
			skewness1_fout.open("tcp_skewness_cg1_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			skewness2_fout.open("tcp_skewness_cg2_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			skewness4_fout.open("tcp_skewness_cg4_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			skewness8_fout.open("tcp_skewness_cg8_p=0.32513:-0.00001:0.31500_q0.88.dat",ios::app);
			skewness_random1_fout.open("tcp_skewness_random_cg1_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			skewness_random2_fout.open("tcp_skewness_random_cg2_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			skewness_random4_fout.open("tcp_skewness_random_cg4_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			skewness_random8_fout.open("tcp_skewness_random_cg8_p=0.32513:-0.00001:0.31500_q0.88_n1024.dat",ios::app);
			
			for(i=0;i<N;++i){
				for(j=0; j<N ;++j){
					spatial_structure_fout<< A[i*N+j]<< endl;
				}
			}
			
			mean_cover_fout << mean << endl;
			
			variance1_fout << variance1 << endl;
			variance2_fout << variance2 << endl;
			variance4_fout << variance4 << endl;
			variance8_fout << variance8 << endl;
			variance_random1_fout << variance_random1 << endl;
			variance_random2_fout << variance_random2 << endl;
			variance_random4_fout << variance_random4 << endl;
			variance_random8_fout << variance_random8 << endl;
			
			skewness1_fout << skewness1<< endl;
			skewness2_fout << skewness2<< endl;
			skewness4_fout << skewness4<< endl;
			skewness8_fout << skewness8<< endl;
			skewness_random1_fout << skewness_random1<< endl;
			skewness_random2_fout << skewness_random2<< endl;
			skewness_random4_fout << skewness_random4<< endl;
			skewness_random8_fout << skewness_random8<< endl;
			
			spatial_structure_fout.close();
			mean_cover_fout.close();

			variance1_fout.close();
			variance2_fout.close();
			variance4_fout.close();
			variance8_fout.close();
			variance_random1_fout.close();
			variance_random2_fout.close();
			variance_random4_fout.close();
			variance_random8_fout.close();
			
			skewness1_fout.close();
			skewness2_fout.close();
			skewness4_fout.close();
			skewness8_fout.close();
			skewness_random1_fout.close();
			skewness_random2_fout.close();
			skewness_random4_fout.close();
			skewness_random8_fout.close();
			
		create_random_matrix(A,mean);	
	} 
	
	
	delete[] A;
	delete[] B;
	delete[] neighbor;
	return 0;
}
