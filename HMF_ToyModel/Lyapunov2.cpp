#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

double epsilon = 1.0;
double t_min = 0.0;
double t_max= 20.0;

//Número de Partículas
const int P = 100;

//Pasos de tiempo
const int N = 50000; 

double initialEnergy = 0.0;

double h = 0.008;

double LyapunovExponent = 0.0;

//External field
double Hx0 = 0.0;
double Hy0 = 0.0;

double Hx = Hx0;
double Hy = Hy0;

vector<double> theta_anterior(P);
vector<double> p_anterior(P);
	
vector<double> theta_actual(P);
vector<double> p_actual(P);

vector<double> basis(2*P);

double derivadaT(int N, int i, double epsilon, vector<double> theta, vector<double> p);
double derivadaP(int N, int i, double epsilon, vector<double> theta, vector<double> p);

double derivadaTFV(int N, int i, double epsilon, vector<double> theta, vector<double> p, vector<double> tangent);
double derivadaPFV(int N, int i, double epsilon, vector<double> theta, vector<double> p, vector<double> tangent);

double Hamiltonian();
double HamiltonianAnterior();
double HamiltonianPotential();

double getElementT(int i, int j);
double getElementP(int i, int j);
void temporalEvolution(double t_min, double t_max, int N, double epsilon);
vector<double> temporalEvolutionFV(double t_min, double t_max, int N, double epsilon, vector<double> basis);
double Mod2Pi(double num);

double Mx();
double My();


int main()
{	
	//Condiciones iniciales

	double p_0 = 1.0;
	ofstream file1;
	file1.open("datosLyap.dat");

	ofstream file2;
	file2.open("datosLyap2.dat");		
	

	double cont=3;
	int cant = 40;
	double zeroEnergy=0.0;
	for(int j=1; j<cant; j++)
	{
		//Condición inicial para los ángulos
		for(int i=0; i<P; i++)
		{
			double angle = 2*M_PI*i/(P*cont);
			//theta_anterior[i] = angle;
			theta_anterior[i] = 2*M_PI*rand()/(RAND_MAX*cont);
			p_anterior[i]= (rand() % j - j/2.0)/5;	
			
			basis[i]=1.0/sqrt(2*P);
			basis[i+P] = 0.0/sqrt(2*P);
		}

		if(j==25){cont=1;}
		
		LyapunovExponent = 0.0;
		double energy = HamiltonianAnterior();
	
		temporalEvolution(t_min, t_max, N, epsilon);

		file1 << energy/P << " " << sqrt(pow(Mx(),2)+pow(My(),2)) << endl;
		cout << "Magnetization: " << energy/P << " " << sqrt(pow(Mx(),2)+pow(My(),2)) << endl;

		file2 << energy/P << " " << LyapunovExponent << endl;
		cout << "LCE: " << energy/P << " " << LyapunovExponent << endl;
	}				
	return 0;
}

vector<double> temporalEvolutionFV(double t_min, double t_max, int N, double epsilon, vector<double> basis)
{	
	vector<double> basis_anterior(2*P);
	vector<double> basis_actual(2*P);
		
	basis_anterior = basis;		
	
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


	for(int k=0; k<P;k++)
	{
		rowT[k] = theta_actual[k];
		rowP[k] = p_actual[k];
	}	
	for(int j=0; j<P;j++)			
	{
		k1T[j] = h*derivadaTFV(N,j,epsilon,rowT,rowP,basis_anterior);
		k1P[j] = h*derivadaPFV(N,j,epsilon,rowT,rowP,basis_anterior);
	}
	for(int k=0; k<P;k++)
	{
		basis_anterior[k] = basis_anterior[k]+0.5*k1T[k];
		basis_anterior[k+P] = basis_anterior[k+P]+0.5*k1P[k] ;
	}
	for(int j=0; j<P;j++)			
	{
		k2T[j] = h*derivadaTFV(N,j,epsilon,rowT,rowP,basis_anterior);
		k2P[j] = h*derivadaPFV(N,j,epsilon,rowT,rowP,basis_anterior);
	}
	for(int k=0; k<P;k++)
	{
		basis_anterior[k] = basis_anterior[k]+0.5*k2T[k];
		basis_anterior[k+P] = basis_anterior[k+P]+0.5*k2P[k] ;
	}
	for(int j=0; j<P;j++)			
	{
		k3T[j] = h*derivadaTFV(N,j,epsilon,rowT,rowP,basis_anterior);
		k3P[j] = h*derivadaPFV(N,j,epsilon,rowT,rowP,basis_anterior);
	}
	for(int k=0; k<P;k++)
	{
		basis_anterior[k] = basis_anterior[k]+k3T[k];
		basis_anterior[k+P] = basis_anterior[k+P]+k3P[k] ;
	}
	for(int j=0; j<P;j++)			
	{
		k4T[j] = h*derivadaTFV(N,j,epsilon,rowT,rowP,basis_anterior);
		k4P[j] = h*derivadaPFV(N,j,epsilon,rowT,rowP,basis_anterior);
	}
	for(int k=0; k<P; k++)
	{
		basis_actual[k] = basis[k] + (1.0/6.0)*(k1T[k]+2.0*k2T[k]+2.0*k3T[k]+k4T[k]);
		basis_actual[k+P] = basis[k+P] + (1.0/6.0)*(k1P[k]+2.0*k2P[k]+2.0*k3P[k]+k4P[k]);
	}		
	
	return basis_actual;
}

void temporalEvolution(double t_min, double t_max, int N, double epsilon)
{	
	for(int i=1;i<N;i++)
	{			
	
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

		if(i==int(N/4))
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
		p_anterior = p_actual;

		//Lyapunov exponents
		
		basis = temporalEvolutionFV(t_min, t_max, N, epsilon, basis);
		
		double norm = 0;
		for(int l=0; l<2*P; l++)
		{
			norm += pow(basis[l],2); 		
		}
		norm = sqrt(norm);
		
		LyapunovExponent += log(norm);
		
		for(int l=0; l<2*P; l++)
		{
			basis[l] = basis[l]/norm; 		
		}
	}
	
	LyapunovExponent /= N*h;
	
	Hx = Hx0;
	Hy = Hy0;
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

double Mx()
{
	double sum=0.0;
	
	for(int i=0; i<P; i++)
	{
		sum += cos(theta_anterior[i]);		
	}
	return sum/P; 
}


double My()
{
	double sum=0.0;
	
	for(int i=0; i<P; i++)
	{
		sum += sin(theta_anterior[i]);		
	}
	return sum/P; 
}


double HamiltonianAnterior()
{
	double ret = 0.0;
	for(int i=0; i<P; i++)
	{
		ret += pow(p_anterior[i],2)/2;
		double sum=0.0;
		for(int j=0; j<P; j++)
		{
			sum += (epsilon/(2.0*P))*(1-cos(theta_anterior[i]-theta_anterior[j]));
		}
		ret += sum;
	}
	return ret;
}

double HamiltonianPotential()
{
	double ret = 0.0;
	for(int i=0; i<P; i++)
	{
		double sum=0.0;
		for(int j=0; j<P; j++)
		{
			sum += (epsilon/(2.0*P))*(1-cos(theta_anterior[i]-theta_anterior[j]));
		}
		ret += sum;
	}
	return ret;
}

double Hamiltonian()
{
	double ret = 0.0;
	for(int i=0; i<P; i++)
	{
		ret += pow(p_actual[i],2)/2;
		double sum=0.0;
		for(int j=0; j<P; j++)
		{
			sum += (epsilon/(2.0*P))*(1-cos(theta_actual[i]-theta_actual[j]));
		}
		ret += sum;
	}
	return ret;
}

double derivadaP(int N, int i, double epsilon, vector<double> theta, vector<double> p)
{	
	return -sin(theta[i])*Mx() + cos(theta[i])*My() + Hx*sin(theta[i]) - Hy*cos(theta[i]);
}
double derivadaT(int N, int i, double epsilon, vector<double> theta, vector<double> p)
{
	return p[i];
}

double derivadaTFV(int N, int i, double epsilon, vector<double> theta, vector<double> p, vector<double> tangent)
{
	return tangent[i+P];
}
double derivadaPFV(int N, int i, double epsilon, vector<double> theta, vector<double> p, vector<double> tangent)
{
	double sum=0;
	for(int j=0; j<P && i!=j; j++)
	{
		sum += ((1.0/P)*cos(theta[i]-theta[j]))*(tangent[j]);
	}
	sum += (-cos(theta[i])*Mx()-sin(theta[i])*My() + 1.0/P )*(tangent[i]);
	return (epsilon)*(sum);
}

double getElementT(int i, int j)
{
	return theta_anterior[i];
}

double getElementP(int i, int j)
{
	return p_anterior[i];
}

