#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

double epsilon = 0.0;
double t_min = 0.0;
double t_max= 100.0;

//Número de Partículas
const int P = 500;

//Pasos de tiempo
const int N = 5000; 
double h = (t_max-t_min)/N;

//Campo externo
double Hx0 = -3.0;
double Hy0 = -3.0;
double Hx = 0.0;
double Hy = 0.0;

vector<double> theta_anterior(P);
vector<double> p_anterior(P);
	
vector<double> theta_actual(P);
vector<double> p_actual(P);

double derivadaT(int N, int i, double epsilon, vector<double> theta, vector<double> p);
double derivadaP(int N, int i, double epsilon, vector<double> theta, vector<double> p);

double getElementT(int i, int j);
double getElementP(int i, int j);
void temporalEvolution(double t_min, double t_max, int N, double epsilon);
double Mod2Pi(double num);

int main()
{	
	//Condiciones iniciales
	double cont=4.0;
	for(int i=0; i<P; i++)
	{
		double angle = 2*M_PI*i/(P*cont);
		
		theta_anterior[i] = 2*M_PI*rand()/RAND_MAX;
		p_anterior[i]= (rand() % 5 - 5/2.0)/5;
	}
	
	ofstream file1;
	file1.open("datosRK.dat");	
	
	cout << "Enter a value for epsilon: ";
	cin >> epsilon;

	double t = t_min;

	for(int i=0; i<P; i++)
	{
		file1 << getElementT(i,0) << " ";
	}	
	file1 << endl;

	temporalEvolution(t_min, t_max, N, epsilon);

	for(int i=0; i<P; i++)
	{
		file1 << getElementT(i,N-1) << " ";
	}	
	file1 << endl;
				
	return 0;
}

void temporalEvolution(double t_min, double t_max, int N, double epsilon)
{
	for(int i=1;i<N;i++)
	{	
		p_anterior = p_actual;		
	
		vector<double> k1T(P);
		vector<double> k2T(P);
		vector<double> k3T(P);
		vector<double> k4T(P);
		vector<double> k1P(P);
		vector<double> k2P(P);
		vector<double> k3P(P);
		vector<double> k4P(P);
		vector<double> rowT(P);
		vector<double> rowP(P);

		if(i==0)
		{
			Hx = Hx0;			
			Hy = Hy0;	
		}
		
		if(i==int(N/10))
		{
			Hx = 0.0;			
			Hy = 0.0;	
		}

		for(int k=0; k<P;k++)
		{
			rowT[k] = getElementT(k,i-1);
			rowP[k] = getElementP(k,i-1) ;
		}	
		for(int j=0; j<P;j++)			
		{
			k1T[j] = h*derivadaT(N,j,epsilon,rowT,rowP);
			k1P[j] = h*derivadaP(N,j,epsilon,rowT,rowP);
		}
		for(int k=0; k<P;k++)
		{
			rowT[k] = getElementT(k,i-1)+0.5*k1T[k];
			rowP[k] = getElementP(k,i-1)+0.5*k1P[k] ;
		}
		for(int j=0; j<P;j++)			
		{
			k2T[j] = h*derivadaT(N,j,epsilon,rowT,rowP);
			k2P[j] = h*derivadaP(N,j,epsilon,rowT,rowP);
		}
		for(int k=0; k<P;k++)
		{
			rowT[k] = getElementT(k,i-1)+0.5*k2T[k];
			rowP[k] = getElementP(k,i-1)+0.5*k2P[k] ;
		}
		for(int j=0; j<P;j++)			
		{
			k3T[j] = h*derivadaT(N,j,epsilon,rowT,rowP);
			k3P[j] = h*derivadaP(N,j,epsilon,rowT,rowP);
		}
		for(int k=0; k<P;k++)
		{
			rowT[k] = getElementT(k,i-1)+k3T[k];
			rowP[k] = getElementP(k,i-1)+k3P[k] ;
		}
		for(int j=0; j<P;j++)			
		{
			k4T[j] = h*derivadaT(N,j,epsilon,rowT,rowP);
			k4P[j] = h*derivadaP(N,j,epsilon,rowT,rowP);
		}

		for(int k=0; k<P; k++)
		{
			theta_actual[k] = Mod2Pi(getElementT(k,i-1) + (1.0/6.0)*(k1T[k]+2.0*k2T[k]+2.0*k3T[k]+k4T[k]));
			p_actual[k] = getElementP(k,i-1) + (1.0/6.0)*(k1P[k]+2.0*k2P[k]+2.0*k3P[k]+k4P[k]);
		}	

		theta_anterior = theta_actual;
	
	}	
}

double Mod2Pi(double num)
{
	int cont=0;
	while(num - 2*M_PI > 0)
	{
		num -= 2*M_PI;
		cont++;
	}
	return num - 2*cont*M_PI;
}

double derivadaP(int N, int i, double epsilon, vector<double> theta, vector<double> p)
{	
	double sum=0;
	for(int j=0; j<P; j++)
	{
		sum -= sin(theta[i]-theta[j]) - Hx*sin(theta[i]) + Hy*cos(theta[i]);
	}
	return (epsilon/(2*P))*(sum);
}
double derivadaT(int N, int i, double epsilon, vector<double> theta, vector<double> p)
{
	return p[i];
}

double getElementT(int i, int j)
{
	return theta_anterior[i];
}

double getElementP(int i, int j)
{
	return p_anterior[i];
}

