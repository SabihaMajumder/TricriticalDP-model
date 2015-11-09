#include <iostream>
using std::cerr;
using namespace std;
#include <fstream>
using std::ofstream;
#include <cstdlib>       // for exit function
#include <ctime>


// parameter definitions
#define N 1024    // system size
#define range 400   // no. of p values in range 0-1      



//====== function to initialize the matrix ===============//
void create_random_matrix(int*U,float x) {      // x is the initial density of the system

	int i1, i2;
	
	for (i1 = 0; i1 < N; i1++) {
		for (i2 = 0; i2 < N; i2++) {
			float number = rand()/(float)RAND_MAX;
			if (number<x) U[i1*N+i2]=1;
			else U[i1*N+i2]=0;
		}
	}

}




		
////////////// main function //////////////////////////////////
int main(){
	

	int*A = new int[N*N];
	float mean;

	ifstream fin;
	fin.open ("tcp_mean_p=0.999:-0.001:0.6_q0_n1024.dat");
	

	for (int p=0 ; p<range; p++){
	  	fin >> mean;
	
		create_random_matrix(A,mean);    //creates random matrix given mean density
		

		ofstream outdata; 
		outdata.open("random_corresponding_to_tcp_p=0.999:-0.001:0.6_q0_n1024.dat",ios::app); 
		
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				outdata<<A[i*N+j]<<"\n";
			}
		}
		outdata.close(); 
		
	}

	
	fin.close();
	
	
	
	delete [] A;
	return 0;
}
