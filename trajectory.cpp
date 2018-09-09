#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

double t_min = 0.0;
double t_max= 100.0;

int N = 10000;

double h = (t_max-t_min)/N;
	
double *x = new double[N];
double *y = new double[N];
double *z = new double[N];

double x_0 = 10.0;
double y_0 = 5.0;
double z_0 = -2.0;

double sigma = 10.0;
double b = 8.0/3.0;

double derivadaX(double t, double x, double y, double z, double gam);
double derivadaY(double t, double x, double y, double z, double gam);
double derivadaZ(double t, double x, double y, double z, double gam);
void temporalEvolution(double t_min, double t_max, int N, double inicial[], double gam);

int main()
{	
	
	ofstream file1;
	file1.open("datosRK.dat");
	
	x[0]=x_0;
	y[0]=y_0;
	z[0]=z_0;

	double inicial[3] = {x_0,y_0,z_0};

	double t = t_min;

	double gam;
	
	cout << "Enter a value for gamma: "; 	
	cin >> gam;

	temporalEvolution(t_min, t_max, N, inicial, gam);	
		
	for(int i=0; i<N;i++)
	{
		file1 << t << " " << x[i] << " " << y[i] << " " << z[i] << std::endl;	
		t += h;		
	}	
		
	return 0;
}

void temporalEvolution(double t_min, double t_max, int N, double inicial[], double b)
{	
	x[0]=inicial[0];
	y[0]=inicial[1];
	z[0]=inicial[2];
	double t = t_min;		
	for(int i=1;i<N;i++)
	{	
		double k1X,k2X,k3X,k4X,k1Y,k2Y,k3Y,k4Y,k1Z,k2Z,k3Z,k4Z = 0;
		k1X = h * derivadaX(t,x[i-1],y[i-1],z[i-1],b);
		k1Y = h * derivadaY(t,x[i-1],y[i-1],z[i-1],b);
		k1Z = h * derivadaZ(t,x[i-1],y[i-1],z[i-1],b);

		k2X = h*derivadaX(t+0.5*h,x[i-1]+0.5*k1X,y[i-1]+0.5*k1Y,z[i-1]+0.5*k1Z,b);
		k2Y = h*derivadaY(t+0.5*h,x[i-1]+0.5*k1X,y[i-1]+0.5*k1Y,z[i-1]+0.5*k1Z,b);
		k2Z = h*derivadaZ(t+0.5*h,x[i-1]+0.5*k1X,y[i-1]+0.5*k1Y,z[i-1]+0.5*k1Z,b);

		k3X = h*derivadaX(t+0.5*h,x[i-1]+0.5*k2X,y[i-1]+0.5*k2Y,z[i-1]+0.5*k2Z,b);
		k3Y = h*derivadaY(t+0.5*h,x[i-1]+0.5*k2X,y[i-1]+0.5*k2Y,z[i-1]+0.5*k2Z,b);
		k3Z = h*derivadaZ(t+0.5*h,x[i-1]+0.5*k2X,y[i-1]+0.5*k2Y,z[i-1]+0.5*k2Z,b);

		k4X = h*derivadaX(t+h, x[i-1]+k3X, y[i-1]+k3Y, z[i-1]+k3Z,b);
		k4Y = h*derivadaY(t+h, x[i-1]+k3X, y[i-1]+k3Y, z[i-1]+k3Z,b);
		k4Z = h*derivadaZ(t+h, x[i-1]+k3X, y[i-1]+k3Y, z[i-1]+k3Z,b);

		x[i] = x[i-1] + (1.0/6.0)*(k1X+2.0*k2X+2.0*k3X+k4X);
		y[i] = y[i-1] + (1.0/6.0)*(k1Y+2.0*k2Y+2.0*k3Y+k4Y);
		z[i] = z[i-1] + (1.0/6.0)*(k1Z+2.0*k2Z+2.0*k3Z+k4Z);
	
		t = t+h;
	}	
}

double derivadaX(double t, double x, double y, double z, double gam)
{
	return -sigma*(x-y);
}


double derivadaY(double t, double x, double y, double z, double gam)
{
	return (gam-z)*x-y;
}


double derivadaZ(double t, double x, double y, double z, double gam)
{
	return x*y - b*z;
}




