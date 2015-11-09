#include <iostream>
using namespace std;
#include <cstdlib> 
#include <ctime>




int* create_random_matrix(int N,int x,int y) {  
	// X,Y is the dimension , x is the lower value of elements and y is the number of values. eg. in Kefi model there are -1,0,1 .
 	//Here x=-1 and y=3
	int L=N+2 ;
	int* U = new int[L*L];
	int i1, i2;
 
	srand (time(NULL));
	for (i1 = 0; i1 < L; i1++) {
		for (i2 = 0; i2 < L; i2++) {
			U[i1*L+i2]=rand()%y + x;
		
			U[0*L +i2] = U[(L-2)*L + i2];      // periodic boundary
			U[(L-1)*L + i2] = U[1*L + i2];
			U[i1*L + 0] = U[i1*L + (L-2)];   //periodic boundary
			U[i1*L + (L-1)] = U[i1*L+ 1];
		}
	}

return (int*) U;
}
   
 int main(){
 int X=4;;
 int* A= create_random_matrix( X,-1,3);
 int i1, i2;
 

 for (i1 = 0; i1 < X+2; i1++) {
  cout<<endl;
  for (i2 = 0; i2 < X+2; i2++) {
  cout<< A[i1*(X+2)+i2]<<"\t";
  }
  }
  delete A;
 return 0;} 
