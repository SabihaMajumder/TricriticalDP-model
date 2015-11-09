#include <iostream>
using std::cerr;
using namespace std;
#include <fstream>
using std::ofstream;
#include <cstdlib>       // for exit function
#include <ctime>


// parameter definitions
#define N 128    // system size
#define dpdt 0.0005   // no. of p values in range 0-1      
#define T 400      // No. of time steps
#define rep 1


//====== function to initialize the matrix ===============//
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
	
	if (i==0) top=N-1;
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
float cg_var (int* A, int n) {				//n is the dimension of coarse graining i.e. dimension of the submatrix
	int i,j,k,l,count = 0;
	float var[N*N/n/n] ;
	for (i=0; i<N-n+1; i+=n) {
		for (j=0; j<N-n+1; j+=n) {	
			float sum = 0,sum_sq = 0;
			for (k=0; k<n; ++k) {
				for (l=0; l<n ; ++l) {
					sum+= A[(i+k)*N+(j+l)];
					sum_sq+= A[(i+k)*N+(j+l)]*A[(i+k)*N+(j+l)];
				}
			}
			float mean = sum/n/n;
			float mean_sq = sum_sq/n/n;
			var[count] = mean_sq - mean*mean;
			count+=1;
		}
	}
	float sum=0;
	for (i=0; i<N*N/n/n; i++) {
		sum+= var[i];
	}
	float variance = sum/(float)(N*N/n/n) ;
return variance;
}

//--------------correlation fn ----------------------------/
float correlation_lag1 (int*A){
	
	int i,j, right, bottom;
	float sum =0;
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
		
			if (i==N-1) bottom = 0;
			else bottom = i+1 ;
	
			if (j==N-1) right = 0;
			else right = j+1;
			 
			sum = sum + ( A[i*N+j] * A[bottom*N+j] + A[i*N+j]* A[i*N+right] );
		}
	}
	
	float average = sum/(float)(2*N*N);
	
return average;

}		
		
////////////// main function //////////////////////////////////
int main(){
	
	srand(time(NULL));
	int x,l,t,i,j,z;
	float p[T],q=0.0;
	int* neighbor = new int[2];
	float den,correlation,variance;

	
	p[0] = 0.71;

	int*A= new int[N*N]; 
			
		create_random_matrix(A,0.7);
		
		for(t=1;t<T;t++){
			
			p[t] = p[t-1]-dpdt;
				
			for(z=0;z<N*N ; ++z){                // so that each site gets selected once on an average
			
				i = rand()%N;           // selecting one random site
				j = rand()%N;
				
				float test = rand()/(float)RAND_MAX;
				float test1 = rand()/(float)RAND_MAX;
				
				if (A[i*N+j]==1){     //if the site is occupied
					
					select_neighbor_of_site(i, j ,neighbor);    //look for a neighbor
					int in = neighbor[0] , jn = neighbor[1];
					
					if (A[in*N +jn]==0) {                     //if neighbor is empty
						if (test < p[t]) 
							A[in*N+jn]=1;                 //regular cp
						else A[i*N+j]=0;
					} 
					else {							  
						if (test < q){
							
							int neighbor_of_pair=select_neighbor_of_pair (in, jn, i, j);  //look for the neighbor of pair 
							A[neighbor_of_pair]=1;
						}
						else if (test1 < 1-p[t] )
							A[i*N+j]=0;	
					}
				}
		
					
			}	
				
			den= density(A);
			correlation = correlation_lag1 (A);
			variance = cg_var (A,4);			
		
//============================saving data in a file===============================//
			ofstream spatial_structure_fout;
			ofstream mean_cover_fout;
			ofstream correlation_fout;
			ofstream variance_fout;
			
			spatial_structure_fout.open("tcp_spatial_data_cont.dat",ios::app);
			mean_cover_fout.open("tcp_mean_cover_cont.dat",ios::app);
			correlation_fout.open("tcp_correlation_cont.dat",ios::app);
			variance_fout.open("tcp_variance_cont.dat",ios::app);
			
			for(i=0;i<N;++i){
				for(j=0; j<N ;++j){
					spatial_structure_fout<< A[i*N+j]<< endl;
				}
			}
			mean_cover_fout << den << endl;
			correlation_fout << correlation -den*den << endl;
			variance_fout << variance << endl;
			
			spatial_structure_fout.close();
			mean_cover_fout.close();
			correlation_fout.close();
			variance_fout.close();
	
	}
			
	
	delete [] A;
	//delete [] A1;
	//delete [] den;
	//delete [] correlation;
	//delete [] variance;
	delete [] neighbor;
	//delete [] neighbor_pair;
	return 0;
}
