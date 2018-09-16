#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int N = 100;

double *x = new double[N];
double *y = new double[N];

double x_0 = 0.8;
double y_0 = 0.5;

int main()
{	
	
	ofstream file1;
	file1.open("datos1.dat");
	
	x[0]=x_0;
	y[0]=y_0;
	
	double a = 1.4;
	double b = 0.3;	

	file1 << 0 << " " << x[0] << " " << y[0] << std::endl;	
	for(int i=0; i<N;i++)
	{
		x[i+1] = y[i]+1-a*pow(x[i],2);
		y[i+1] = b*x[i];
		file1 << i+1 << " " << x[i+1] << " " << y[i+1] << std::endl;
				
	}	

	ofstream file2;
	file2.open("datos2.dat");
	
	double epsilon = 0.001;
	
	x[0]=x_0+epsilon;
	y[0]=y_0+epsilon;


	file2 << 0 << " " << x[0] << " " << y[0] << std::endl;	
	for(int i=0; i<N;i++)
	{
		x[i+1] = y[i]+1-a*pow(x[i],2);
		y[i+1] = b*x[i];
		if(i == 10 || i == 50 || i == 100 || i == 1000 || i == 10000)
		{
			file2 << i << " " << x[i] << " " << y[i] << std::endl;
		}		
	}
		
	return 0;
}


