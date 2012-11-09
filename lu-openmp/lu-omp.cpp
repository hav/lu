/*	filename: LU-omp.cpp 
	Parallel LU Factorization implementation in OpenMP
*/
// Parallel assignment.cpp : Defines the entry point for the console application.


//#include "stdafx.h"	/* used in VS */
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <omp.h>	/* OpenMP header */
//#include <time.h> 	/* used in VS timer */
#include <sys/time.h> 	/* used in *nix timer */

using namespace std;

int main(int argc, char* argv[]) {

  if(argc < 3) {
    fprintf(stderr, "Usage: ./lu-omp <size-of-linear-system> <num-of-threads>\n\ne.g: ./lu-omp 900 8\n\n");
    exit(1);
  }

  for (int repeat=0; repeat<5; repeat++)
  {	
  /* Timestamp variables */
  time_t start,end;
  double start_time, end_time;
  struct timeval tz;
  struct timezone tx;
 
  /* Set number of thread omp. now: Khawarizmi = 8 logical processors */
  int threadx = atoi(argv[2]);
  omp_set_nested(0);
  omp_set_dynamic(1);
  omp_set_num_threads(threadx);
  //-------------------------------------------------------------- Initiate Matrix A and Matrix B
  int n = atoi(argv[1]); /* enable n in terminal, Otherwise use below */	 
  // int n = 900; //----------------- Dimension n

  // printf("n:%d\n\n", n);
  float *x = new float [n]; //--- initiation matrix X
  float *y = new float [n]; //--- initiation matrix y
  float sum1 =0, sum2=0, sum3=0, sum4=0, sum=0;

  float **a; //------------------- initiate Matrix A
  a = new float *[n];
  
  for (int i = 0; i<n; i++) {	
    a[i] = new float [n];
  }

  float *b= new float [n]; //------ initiate Matrix B

  for (int i=0; i<(n); i++ ) {
    for (int j=0; j<(n); j++ ) {
      if (i==j)
      a[i][j] = n;
      else
      a[i][j] = 1;
      //cout<<a[i][j]<<" "; 
    }
    
  b [i]= (2*n)-1;
  // cout<<" = "<<b [i]<<endl ;
  }

  //-------------------------------------------------------------- Upper A matrix initiation

  /* Windows timer start */
  int start_s=clock();
  
  /* *nix timer start */
  gettimeofday(&tz, &tx);
  start_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;

  float **upperA ;
  upperA= new float *[n];

  for (int i = 0; i<n; i++) {
    upperA[i] = new float [n];
  } 

  //cout<<"\n\n\n";
  for (int i=0; i<(n); i++ ) { 
    for (int j=0; j<(n); j++ ) {
      if (i>j)
      upperA[i][j] = 0;      
    } 
  }

  //------------------------------------------------ initiate lower A

  float ** lowerA;
  lowerA = new float *[n];
  //cout<<"\n\n\n";
  for (int i = 0; i<n; i++) {
    lowerA[i] = new float [n];
  }

  for (int i=0; i<(n); i++ ) { 
    for (int j=0; j<(n); j++ ) {
      if (i==j)
      lowerA[i][j] = 1;
      else if (i<j)
      lowerA[i][j] = 0;
    }
  
  }

 //------------------------------calculate upper A
 for (int i=0; i<(n); i++ ) {
   #pragma omp parallel for private ( sum1)
   for(int j=i; j<(n); j++ ) {
   sum1=0;
    
   #pragma omp parallel for  reduction (+: sum1)
   for (int e=0; e<(n); e++ ){if (e!=i) {
     sum1=sum1 + (lowerA[i][e]* upperA[e][i]);}
     //cout <<upperA[i][j]<<" "; 	
   }

   sum1=a[i][j]-sum1;
   upperA[i][j]=sum1;
   //cout<<upperA[i][j]<<" "; //<<" (sum "<<sum1<<" )"
   sum1=0;
    
   }
  
   if(i+1<n) { 
   #pragma omp parallel for private(sum1)
   for(int j=i+1; j<n; j++) {
     sum1=0;
     #pragma omp parallel for  reduction (+: sum1)
     for (int e=0; e<(n); e++ ){if (e!=i) {
       sum1=sum1 + (lowerA[i][e]* upperA[e][i]);}
       //cout <<upperA[i][j]<<" "; 
     } 
      
     lowerA[j][i]=(a[j][i]-sum1)/upperA[i][i];	
     sum1=0;
      
    }
    
   }  
  }// end for	

  //---------------------------calculate lower B
  //-------------------------------------------------------------- Calculate for y[n]
  //cout<<"\n\n";

  double reduct=0;
  for (int i = 0; i<n; i ++) {
    #pragma omp parallel for reduction(+:sum)
    for (int j=0; j<n; j++) {
      if (i!=j)
      sum = sum +(lowerA[i][j]*y[j]);
    }
 
    y[i] = (b[i]- sum)/ lowerA[i][i];
    sum = 0;
  }

  //-------------------------------------------------------------- Calculate for x[n], the final answer

  for (int i = (n-1); i>-1; i --) {
    #pragma omp parallel for  reduction (+: sum)
    for (int j=(n-1); j>-1; j--) {
      if (i!=j)
	sum = sum +(upperA[i][j]*x[j]);
    }

    x[i] = (y[i]- sum)/ upperA[i][i];
    sum = 0;
  }

  /* time stop - Windows */
  int stop_s=clock();
  cout<<"\n test number :"<<repeat+1<<endl;
  cout << "WINDOWS Execute time: " <<(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 <<" miliseconds"<< endl;
  
  /* time stop - *nix */
  gettimeofday(&tz, &tx);
  end_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;  
  /* cal */
  printf("Execution time: %f seconds.\n", end_time - start_time);
  printf("Execution time: %f miliseconds.\n\n", (end_time - start_time) * 1000.0);
  //printf("Execution time: %f microseconds.\n", (end_time - start_time) * 1000000.0);

  //DISABLE ALL PRINT ALL COUT
  /*
  for (int i = 0; i<n; i ++) {
    cout<<"X"<<(i+1)<<" = "<<x[i]<<endl;
  }

  cout<<"\n\n";
  for (int i=0; i<n; i++) {
    cout<<y[i]<<endl;
  }

  cout<<"\n\n";
  for (int i=0; i<n; i++) {
    cout<<x[i]<<endl;
  }

  cout<<"lowerA :\n";
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) cout <<lowerA[i][j]<<" ";
    cout<<endl;
  }

  cout<<"upperA :\n";
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      cout <<upperA[i][j]<<" ";
    }
    
  cout<<endl;
  }

  */
  //cout<<"answer for x is in progress.. please wait, the coding is very confusing, hahaha \n\n";
  }

  //system("PAUSE");	/* used in VS */
  //return EXIT_SUCCESS;/* used in VS */
}


