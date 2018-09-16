#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

double x_f = 0;
double y_f = 0;
double z_f = 0;

double x_base = 0;
double y_base = 0;
double z_base = 0;

double *col1 = new double[3];
double *col2 = new double[3];;
double *col3 = new double[3];;

double sigma = 10.0;
double b = 8.0/3.0;


// Evolucion temporal
int K = 1000000;
double t_min = 0.0;
double t_max= 0.005*K;
int N = K;

double h = (t_max-t_min)/N;
	
double *x = new double[K];
double *y = new double[K];
double *z = new double[K];


//Funciones

double derivadaX(double t, double x, double y, double z, double gam);
double derivadaY(double t, double x, double y, double z, double gam);
double derivadaZ(double t, double x, double y, double z, double gam);

double derivadaXFV(double t, double x, double y, double z, double x_1, double y_1, double z_1, double gam);
double derivadaYFV(double t, double x, double y, double z, double x_1, double y_1, double z_1, double gam);
double derivadaZFV(double t, double x, double y, double z, double x_1, double y_1, double z_1, double gam);

void temporalEvolution(double t_min, double t_max, int N, double inicial[], double gam);
double * temporalEvolutionFirstVariation(double t_min, double t_max, int N, double inicialV[], double inicial[], double gam);
double * LyapunovExponent(double t_min, double t_max, int N, int n, double inicial[], double base1[], double base2[], double base3[], double gam);

int main()
{	
	// Punto inicial

	double L = 1;

	double base1[3] = {1.0/sqrt(3*L),0.0/sqrt(3*L),0.0/sqrt(3*L)};
	double base2[3] = {0.0/sqrt(3*L),1.0/sqrt(3*L),0.0/sqrt(3*L)};
	double base3[3] = {0.0/sqrt(3*L),0.0/sqrt(3*L),1.0/sqrt(3*L)};
	
	ofstream file1;
	file1.open("datos3.dat");		
	

	double x_0 = 10.0;
	double y_0 = 10.0;
	double z_0 = 30.0;
	double initial[3] = {x_0,y_0,z_0}; 
	double gam = 0;
	
	cout << "Value for r: ";
	cin >> gam;
	
	double* lambda_spectrum = new double[3];
	lambda_spectrum = LyapunovExponent(t_min, t_max, N, K, initial, base1, base2, base3, gam);
	
	cout << "Lyapunov exponents of canonical Lorenz system" << endl;
	cout << lambda_spectrum[0] << " " << lambda_spectrum[1] << " " << lambda_spectrum[2] << endl;
	cout << "Sum of spectrum is " << lambda_spectrum[0] + lambda_spectrum[1] + lambda_spectrum[2] << " and it should be " << -(1+b+sigma) << endl;

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

double * temporalEvolutionFirstVariation(double t_min, double t_max, int N, double inicialV[], double inicial[], double b)
{
	double dx=0;
	double dy=0;
	double dz=0;
	
	dx=inicialV[0];
	dy=inicialV[1];
	dz=inicialV[2];

	double t=0;

	double k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1z,k2z,k3z,k4z   = 0;
	k1x = h * derivadaXFV(t,dx,dy,dz,inicial[0],inicial[1],inicial[2],b);
	k1y = h * derivadaYFV(t,dx,dy,dz,inicial[0],inicial[1],inicial[2],b);
	k1z = h * derivadaZFV(t,dx,dy,dz,inicial[0],inicial[1],inicial[2],b);
			
	k2x = h*derivadaXFV(t+0.5*h,dx+0.5*k1x,dy+0.5*k1y,dz+0.5*k1z,inicial[0],inicial[1],inicial[2],b);
	k2y = h*derivadaYFV(t+0.5*h,dx+0.5*k1x,dy+0.5*k1y,dz+0.5*k1z,inicial[0],inicial[1],inicial[2],b);
	k2z = h*derivadaZFV(t+0.5*h,dx+0.5*k1x,dy+0.5*k1y,dz+0.5*k1z,inicial[0],inicial[1],inicial[2],b);

	k3x = h*derivadaXFV(t+0.5*h,dx+0.5*k2x,dy+0.5*k2y,dz+0.5*k2z,inicial[0],inicial[1],inicial[2],b);
	k3y = h*derivadaYFV(t+0.5*h,dx+0.5*k2x,dy+0.5*k2y,dz+0.5*k2z,inicial[0],inicial[1],inicial[2],b);
	k3z = h*derivadaZFV(t+0.5*h,dx+0.5*k2x,dy+0.5*k2y,dz+0.5*k2z,inicial[0],inicial[1],inicial[2],b);

	k4x = h*derivadaXFV(t+h, dx+k3x, dy+k3y, dz+k3z,inicial[0],inicial[1],inicial[2],b);
	k4y = h*derivadaYFV(t+h, dx+k3x, dy+k3y, dz+k3z,inicial[0],inicial[1],inicial[2],b);
	k4z = h*derivadaZFV(t+h, dx+k3x, dy+k3y, dz+k3z,inicial[0],inicial[1],inicial[2],b);	

	double* ret = new double[3];
	ret[0] = dx + (1.0/6.0)*(k1x+2.0*k2x+2.0*k3x+k4x);
	ret[1] = dy + (1.0/6.0)*(k1y+2.0*k2y+2.0*k3y+k4y);
	ret[2] = dz + (1.0/6.0)*(k1z+2.0*k2z+2.0*k3z+k4z);
	return ret;
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

double derivadaXFV(double t, double x, double y, double z, double x_1, double y_1, double z_1, double gam)
{
	return -sigma*x+sigma*y;
}

double derivadaYFV(double t, double x, double y, double z, double x_1, double y_1, double z_1, double gam)
{
	return (gam-z_1)*x - y - x_1*z;
}

double derivadaZFV(double t, double x, double y, double z, double x_1, double y_1, double z_1, double gam)
{
	return y_1*x+x_1*y-b*z;
}


double * LyapunovExponent(double t_min, double t_max, int N, int n, double inicial[], double base1[], double base2[], double base3[], double b)
{
	double x_base1 = base1[0];
	double y_base1 = base1[1];
	double z_base1 = base1[2];

	double x_base2 = base2[0];
	double y_base2 = base2[1];
	double z_base2 = base2[2];	

	double x_base3 = base3[0];
	double y_base3 = base3[1];
	double z_base3 = base3[2];
	
	double sum1=0;
	double sum2=0;
	double sum3=0;

	double x_actual=inicial[0];
	double y_actual=inicial[1];
	double z_actual=inicial[2];
	
	int T = K;

	double points[3] = {x_actual, y_actual, z_actual};
	temporalEvolution(t_min, t_max, T, points, b);
	
	double points2[3] = {x[T-1], y[T-1], z[T-1]};
	temporalEvolution(t_min, t_max, N, points2, b);
		
	for(int j=0; j<n; j++)
	{
		x_actual = x[j];
		y_actual = y[j];
		z_actual = z[j];
			
		double point[3] = {x_actual,y_actual,z_actual};
		double tang1[3] = {x_base1,y_base1,z_base1};
		double tang2[3] = {x_base2,y_base2,z_base2};
		double tang3[3] = {x_base3,y_base3,z_base3};
		
		double *evolvedBasis1 = new double[3];
		double *evolvedBasis2 = new double[3];
		double *evolvedBasis3 = new double[3];
	
		evolvedBasis1 = temporalEvolutionFirstVariation(t_min, t_max, N, tang1, point,b);

		evolvedBasis2 = temporalEvolutionFirstVariation(t_min, t_max, N, tang2, point,b);

		evolvedBasis3 = temporalEvolutionFirstVariation(t_min, t_max, N, tang3, point,b);
		
		x_base1 = evolvedBasis1[0];
		y_base1 = evolvedBasis1[1];
		z_base1 = evolvedBasis1[2];

		x_base2 = evolvedBasis2[0];
		y_base2 = evolvedBasis2[1];
		z_base2 = evolvedBasis2[2];

		x_base3 = evolvedBasis3[0];
		y_base3 = evolvedBasis3[1];
		z_base3 = evolvedBasis3[2];
		
		double norm = sqrt(pow(x_base1,2)+pow(y_base1,2)+pow(z_base1,2));	
		
		sum1 += log(norm);	
	
		x_base1 = (x_base1/norm);
		y_base1 = (y_base1/norm);
		z_base1 = (z_base1/norm);

		//Gram-Schmidt

		double OneDotTwo = x_base1*x_base2 + y_base1*y_base2 + z_base1*z_base2;
	
		x_base2	-= OneDotTwo * x_base1;
		y_base2	-= OneDotTwo * y_base1;
		z_base2	-= OneDotTwo * z_base1;	

		norm = sqrt(pow(x_base2,2)+pow(y_base2,2)+pow(z_base2,2));	
	
		sum2 += log(norm);	
	
		x_base2 = (x_base2/norm);
		y_base2 = (y_base2/norm);
		z_base2 = (z_base2/norm);

		//Gram-Schmidt

		double OneDotThree = x_base1*x_base3 + y_base1*y_base3 + z_base1*z_base3;
	
		double TwoDotThree = x_base3*x_base2 + y_base3*y_base2 + z_base3*z_base2;
	
		x_base3	-= OneDotThree * x_base1 + TwoDotThree * x_base2;
		y_base3	-= OneDotThree * y_base1 + TwoDotThree * y_base2;
		z_base3	-= OneDotThree * z_base1 + TwoDotThree * z_base2;	

		norm = sqrt(pow(x_base3,2)+pow(y_base3,2)+pow(z_base3,2));	
		
		sum3 += log(norm);	
	
		x_base3 = (x_base3/norm);
		y_base3 = (y_base3/norm);
		z_base3 = (z_base3/norm);


	}
	double cant = n;
	
	double *spectrum = new double[3];
	spectrum[0] = sum1/(cant*h);
	spectrum[1] = sum2/(cant*h);
	spectrum[2] = sum3/(cant*h);
	
	return spectrum;
}









