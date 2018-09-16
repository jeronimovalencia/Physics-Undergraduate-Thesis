#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

double t_min = 0.0;
double t_max= 1000.0;

int N = 100000;

double h = (t_max-t_min)/N;
	
double *x = new double[N];
double *y = new double[N];

double x_0 = 0.8;
double y_0 = 0.0;

double mu=1.0;

double derivadaX(double t, double x, double y);
double derivadaY(double t, double x, double y);
void temporalEvolution(double t_min, double t_max, int N, double inicial[]);

int main()
{	
	
	ofstream file1;
	file1.open("datosRK1.dat");
	
	x[0]=x_0;
	y[0]=y_0;

	double inicial[2] = {x_0,y_0};

	double t = t_min;

	temporalEvolution(t_min, t_max, N, inicial);	
		
	for(int i=0; i<N;i++)
	{
		file1 << t << " " << x[i] << " " << y[i]<< std::endl;	
		t += h;		
	}

	ofstream file2;
	file2.open("datosRK2.dat");
	
	x_0=2.0;
	y_0=0.0;

	double inicial2[2] = {x_0,y_0};

	t = t_min;

	temporalEvolution(t_min, t_max, N, inicial2);	
		
	for(int i=0; i<N;i++)
	{
		file2 << t << " " << x[i] << " " << y[i]<< std::endl;	
		t += h;		
	}

	ofstream file3;
	file3.open("datosRK3.dat");
	
	x_0=-2.0;
	y_0=0.0;

	double inicial3[2] = {x_0,y_0};

	t = t_min;

	temporalEvolution(t_min, t_max, N, inicial3);	
		
	for(int i=0; i<N;i++)
	{
		file3 << t << " " << x[i] << " " << y[i]<< std::endl;	
		t += h;		
	}	

	ofstream file4;
	file4.open("datosRK4.dat");
	
	x_0=-1.0;
	y_0=0.6;

	double inicial4[2] = {x_0,y_0};

	t = t_min;

	temporalEvolution(t_min, t_max, N, inicial4);	
		
	for(int i=0; i<N;i++)
	{
		file4 << t << " " << x[i] << " " << y[i]<< std::endl;	
		t += h;		
	}

	ofstream file5;
	file5.open("datosRK5.dat");
	
	x_0=-2.0;
	y_0=-2.0;

	double inicial5[2] = {x_0,y_0};

	t = t_min;

	temporalEvolution(t_min, t_max, N, inicial5);	
		
	for(int i=0; i<N;i++)
	{
		file5 << t << " " << x[i] << " " << y[i]<< std::endl;	
		t += h;		
	}

	ofstream file6;
	file6.open("datosRK6.dat");
	
	x_0=0.05;
	y_0=0.0;

	double inicial6[2] = {x_0,y_0};

	t = t_min;

	temporalEvolution(t_min, t_max, N, inicial6);	
		
	for(int i=0; i<N;i++)
	{
		file6 << t << " " << x[i] << " " << y[i]<< std::endl;	
		t += h;		
	}
		
	return 0;
}

void temporalEvolution(double t_min, double t_max, int N, double inicial[])
{	
	x[0]=inicial[0];
	y[0]=inicial[1];
	double t = t_min;		
	for(int i=1;i<N;i++)
	{	
		double k1X,k2X,k3X,k4X,k1Y,k2Y,k3Y,k4Y = 0;
		k1X = h * derivadaX(t,x[i-1],y[i-1]);
		k1Y = h * derivadaY(t,x[i-1],y[i-1]);

		k2X = h*derivadaX(t+0.5*h,x[i-1]+0.5*k1X,y[i-1]+0.5*k1Y);
		k2Y = h*derivadaY(t+0.5*h,x[i-1]+0.5*k1X,y[i-1]+0.5*k1Y);

		k3X = h*derivadaX(t+0.5*h,x[i-1]+0.5*k2X,y[i-1]+0.5*k2Y);
		k3Y = h*derivadaY(t+0.5*h,x[i-1]+0.5*k2X,y[i-1]+0.5*k2Y);

		k4X = h*derivadaX(t+h, x[i-1]+k3X, y[i-1]+k3Y);
		k4Y = h*derivadaY(t+h, x[i-1]+k3X, y[i-1]+k3Y);

		x[i] = x[i-1] + (1.0/6.0)*(k1X+2.0*k2X+2.0*k3X+k4X);
		y[i] = y[i-1] + (1.0/6.0)*(k1Y+2.0*k2Y+2.0*k3Y+k4Y);
	
		t = t+h;
	}	
}

double derivadaX(double t, double x, double y)
{
	return y;
}


double derivadaY(double t, double x, double y)
{
	return x-pow(x,3)-mu*y*(pow(y,2)-pow(x,2)+0.5*pow(x,4));
}






