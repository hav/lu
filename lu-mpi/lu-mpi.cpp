// Parallel assignment.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
//#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
//#include <stdio.h>
//#include <math.h>
#include <iomanip>
//#include <omp.h>

/* timer header */
#include <time.h>
#include <sys/time.h> /* *nix time */

/* MPI Header */
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{

if(argc < 2) {
	fprintf(stderr, "Usage: mpirun -np <num-of-process> lu-mpi <size-of-linear-system>\n\ne.g: mpirun -np 16 lu-mpi 900\n\n");
	exit(1);
}

int mysize,myrank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &mysize);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//for (int repeat=0; repeat<5; repeat++)
//{
MPI_Status   status;	

time_t start,end;
double start_time, end_time;
struct timeval tz;
struct timezone tx;
//-------------------------------------------------------------- Initiate Matrix A and Matrix B	 
	int n = atoi(argv[1]); // parse size of linear sys
	// int n = 900; //----------------- Dimension n
	// printf("n:%d\n\n", n);

	 float *x = new float [n]; //--- initiation matrix X
	 float *y = new float [n]; //--- initiation matrix y
	 float sum1 =0, sum2=0, sum3=0, sum4=0, sum=0;

	 float **a; //------------------- initiate Matrix A
a = new float *[n];
for (int i = 0; i<n; i++)
{	
	a[i] = new float [n];
}

	float *b= new float [n]; //------ initiate Matrix B

for (int i=0; i<(n); i++ )

  {
	for (int j=0; j<(n); j++ )
	{
		if (i==j)
		a[i][j] = n;
		else
		a[i][j] = 1;

	//	cout<<a[i][j]<<" "; 
	}	
  
	b [i]= (2*n)-1;
   // cout<<" = "<<b [i]<<endl ;
  }

//-------------------------------------------------------------- Upper A matrix initiation

/* timer start */
int start_s=clock();
gettimeofday(&tz, &tx);
start_time = (double)tz.tv_sec + (double)tz.tv_usec / 1000000.0;

float **upperA ;
upperA= new float *[n];

for (int i = 0; i<n; i++)
{	
	upperA[i] = new float [n];
} 

//cout<<"\n\n\n";
for (int i=0; i<(n); i++ )

  { 
	for (int j=0; j<(n); j++ )
	{
		
		if (i>j)
		upperA[i][j] = 0;
		
	}	
  
  }
//------------------------------------------------ initiate lower A

float ** lowerA;
lowerA = new float *[n];


//cout<<"\n\n\n";
for (int i = 0; i<n; i++)
{	
	lowerA[i] = new float [n];
}

for (int i=0; i<(n); i++ )

  { 
	for (int j=0; j<(n); j++ )
	{
		if (i==j)
		lowerA[i][j] = 1;
		else if (i<j)
		lowerA[i][j] = 0;
		}	
  
  }	


int BlockSize, MyBlockSize, MyStartPoint;

if(n%mysize == 0)
	{	BlockSize = n/mysize;
		MyBlockSize = n/mysize;
		MyStartPoint = MyBlockSize * myrank;
//cout<<"myrank : "<<myrank<<"  My starting point : "<<MyStartPoint<<"  MyBlocksize :"<<MyBlockSize<<endl;
	}
else
	{
		if(myrank+1==mysize)
		{
		BlockSize= n%mysize + n/mysize;
		MyBlockSize = n%mysize + n/mysize;
		MyStartPoint = (n/mysize) * myrank;
//cout<<"myrank : "<<myrank<<"  My starting point : "<<MyStartPoint<<"  MyBlocksize :"<<MyBlockSize<<endl;	
		}
		
		else
		{
		BlockSize = n/mysize;
		MyBlockSize=n/mysize;
		MyStartPoint = MyBlockSize * myrank;
//		cout<<"myrank : "<<myrank<<"  My starting point : "<<MyStartPoint<<"  MyBlocksize :"<<MyBlockSize<<endl;
		}
	}

//for (int repeat=0; repe
sum1=0;
MPI_Barrier(MPI_COMM_WORLD);
//------------------------------calculate upper 
//

for (int i=0; i<(n); i++ )

  { 

	  
//	  #pragma omp parallel for private ( sum1)
//	  for(int j=i; j<(n); j=j+mysize )
//	  {  
		float Buff[MyBlockSize];
		float readBuff[n];
		int countBuff=0;  
		  if(myrank==0 )
		  {
	
		for(int j =0; j<MyBlockSize; j++)
		{

		  if(j>=i)
			{
			for (int e=0; e<(n); e++ ){if (e!=i)
				{sum1=sum1 + (lowerA[i][e]* upperA[e][i]);}
				}

			sum1=a[i][j]-sum1;
			upperA[i][j]=sum1;
			readBuff[j]=upperA[i][j];
	//		printf("sum :%f",sum1);
	//		cout<<upperA[i][j]<<endl;
			sum1=0;
			//cout<<<<upperA[i][j]<<endl;
			}
		  else if (j<i)
		  {upperA[i][j]=0;
		  readBuff[j]=0;}
	//	  if (i==1)
	//  {cout<<"RANK 0 "<<upperA[i][j]<<" for j :"<<j<<" with blocksize : "<<MyBlockSize<<endl;}
		}
			int Point=MyBlockSize;
			for(int e = 1; e < mysize; e++)
			{
				if (e==mysize-1)
				{//	cout<<"receiving from"<<e<<endl;
					int Temp = n%mysize + n/mysize;
				MPI_Recv(&readBuff[Point],Temp,MPI_FLOAT,e,e,MPI_COMM_WORLD,&status);}
				else
				{
				//	cout<<"receiving from"<<e<<endl;
				MPI_Recv(&readBuff[Point],BlockSize,MPI_FLOAT,e,e,MPI_COMM_WORLD,&status);
				Point= Point+BlockSize;
				}
			}
			
	//		for(int e = 0;e<n;e++)
	//		{
	//	        cout<<readBuff[e]<<" ";
	//		}
	//		cout<<endl;
			MPI_Bcast(&readBuff[0],n,MPI_FLOAT,0,MPI_COMM_WORLD);
			for(int e =0;e<n;e++)
			{upperA[i][e]=readBuff[e];
			readBuff[e]=0;//---------reinitialize 
			}

		  }




		
		  else 
{
		//	  cout<<"hehehe its me"<<myrank<<endl;
		for(int j = MyStartPoint; j<MyStartPoint+MyBlockSize; j++)
	{  
		
		
		  if(j>=i)
		{
		  	for (int e=0; e<(n); e++ )
			{if (e!=i)
			  {
				sum1=sum1 + (lowerA[i][e]* upperA[e][i]);
				}
			}
				sum1=a[i][j]-sum1;
				Buff[countBuff]=sum1;
	//			cout<<Buff[countBuff]<<endl;
				countBuff++;
						         } 
		       	else
		{    
		       	countBuff++; }
		sum1=0;
		       }
	
MPI_Send(Buff,MyBlockSize,MPI_FLOAT,0,myrank,MPI_COMM_WORLD);	
MPI_Bcast(readBuff,n,MPI_FLOAT,0,MPI_COMM_WORLD);

for(int e=0;e<n;e++)
{upperA[i][e]=readBuff[e];
//cout<<readBuff[e];
readBuff[e]=0;}
countBuff=0;
//cout<<endl;
}





		  //broadcast here, with interval the size of processor, (int j=i; j<n; j=j+ProcSize)i
		  //use for loop and MPI SEND to send everyone
		  sum1=0;
		//  #pragma omp parallel for  reduction (+: sum1)
			
	//cout<<upperA[i][j]<<" "; //<<" (sum "<<sum1<<" )"
	//receive all upper A here.
	  sum1=0;

	
//}
//MPI_Barrier(MPI_COMM_WORLD);
	  if(i+1<n)
		{ 
//	#pragma omp parallel for private(sum1)	
	//	for(int j=i+1; j<n; j++)
	//	  {

			if(myrank==0) //----------------master rank 0
			  { for(int j=0; j<MyBlockSize; j++)
				  { 
					  if(j>i)
					{
				  	for (int e=0; e<(n); e++ )
					{
					if (e!=i)
						{sum1=sum1 + (lowerA[i][e]* upperA[e][i]);}
						//cout <<upperA[i][j]<<" "; 

					}
			 

					lowerA[j][i]=(a[j][i]-sum1)/upperA[i][i];
				        readBuff[j]=lowerA[j][i];	
					sum1=0;
					}
					  else  if(i==j)
					  {readBuff[j]=1;}
					  else if(j<i)
					  {readBuff[j]=0;}
				}
                        //------------------------------------------ lower gathring
		        int Point=MyBlockSize;
			for(int e = 1; e < mysize; e++)
			{
				if (e==mysize-1)
				{	int Temp = n%mysize + n/mysize;
				MPI_Recv(&readBuff[Point],Temp,MPI_FLOAT,e,e,MPI_COMM_WORLD,&status);
			//	cout<<"receiving from"<<e<<"lower part"<<endl
				}

				else
				{
				//	cout<<"receiving from"<<e<<"lower part"<<endl;
				MPI_Recv(&readBuff[Point],BlockSize,MPI_FLOAT,e,e,MPI_COMM_WORLD,&status);
				Point= Point+BlockSize;
			//	cout<<readBuff[e];
				}
			}
	//		for(int e = 0;e<n;e++)
	//		{
	//		cout<<readBuff[e]<<" ";
	//		}
	//			cout<<endl<<endl;
			MPI_Bcast(&readBuff[0],n,MPI_FLOAT,0,MPI_COMM_WORLD);
			for(int e =0;e<n;e++)
			{lowerA[e][i]=readBuff[e];}


		       	}
			 
		       	else //-----------------------slave
			  {					 
				  for (int j = MyStartPoint; j<MyStartPoint+MyBlockSize;j++)
				  { if(j>i)
					  {
				 	 for (int e=0; e<(n); e++ )
					{
					if (e!=i)
					{sum1=sum1 + (lowerA[i][e]* upperA[e][i]);}
					//cout <<upperA[i][j]<<" "; 

					}
			 

					lowerA[j][i]=(a[j][i]-sum1)/upperA[i][i];
					Buff[countBuff]= lowerA[j][i];
			//		cout<<Buff[countBuff]<<endl;
					countBuff++;
					sum1=0;
					}
				  else  if(i==j)
					  {Buff[countBuff]=1;
					  countBuff++;}
					  else if(j<i)
					  {Buff[countBuff]=0;
					  countBuff++;}
				  }
	
MPI_Send(Buff,MyBlockSize,MPI_FLOAT,0,myrank,MPI_COMM_WORLD);	
			  MPI_Bcast(readBuff,n,MPI_FLOAT,0,MPI_COMM_WORLD);
			  for(int e=0;e<n;e++)
			  {lowerA[e][i]=readBuff[e];
			// cout<<readBuff[e];
			  
			  }
			//  cout<<endl;
			  countBuff=0;
			}

			  sum1=0;

			
	//	}
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
  }	

//---------------------------calculate lower B


//-------------------------------------------------------------- Calculate for y[n]
//cout<<"\n\n";
MPI_Barrier(MPI_COMM_WORLD);
if(myrank ==0)
{
double reduct=0;
for (int i = 0; i<n; i ++) 
{	

// #pragma omp parallel for reduction(+:sum)
	for (int j=0; j<n; j++)
	{ 	

		if (i!=j)
	
		sum = sum +(lowerA[i][j]*y[j]);
	}
 
	
y[i] = (b[i]- sum)/ lowerA[i][i];

 sum = 0;


}

//-------------------------------------------------------------- Calculate for x[n], the final answer




for (int i = (n-1); i>-1; i --)
{	
//#pragma omp parallel for  reduction (+: sum)
	for (int j=(n-1); j>-1; j--)	
	{ 

		if (i!=j)
		sum = sum +(upperA[i][j]*x[j]);
	}

x[i] = (y[i]- sum)/ upperA[i][i];

sum = 0;
}

/* time stop */

        int stop_s=clock();

        cout << "Execute time: " <<(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 <<"miliseconds"<< endl;
        
gettimeofday(&tz, &tx);	
end_time = (double)tz.tv_sec + (double)tz.tv_usec / 1000000.0;
printf("*nix execution time: %f seconds.\n", end_time - start_time);
printf("*nix execution time: %f miliseconds.\n", (end_time - start_time) * 1000.0);

        

 //DISABLE ALL PRINT ALL COUT

//for (int i = 0; i<n; i ++)
//{
//	cout<<"X"<<(i+1)<<" = "<<x[i]<<endl;


//}

}
/*


cout<<"\n\n";
for (int i=0; i<n; i++)
{
	cout<<y[i]<<endl;
	
}

cout<<"\n\n";
for (int i=0; i<10; i++)
{
	cout<<x[i]<<endl;
	
}

cout<<"lowerA :\n";
for (int i=0; i<n; i++)
{
	for (int j=0; j<n; j++)
	{
		cout <<lowerA[i][j]<<" ";

	}

	cout<<endl;
}

cout<<"upperA :\n";
for (int i=0; i<n; i++)
{
	for (int j=0; j<n; j++)
	{
		cout <<upperA[i][j]<<" ";

	}

	cout<<endl;
}
*/
//cout<<"answer for x is in progress.. please wait, the coding is very confusing, hahaha \n\n";


//	}
MPI_Finalize();
	//system("PAUSE");
    //return EXIT_SUCCESS;


}


