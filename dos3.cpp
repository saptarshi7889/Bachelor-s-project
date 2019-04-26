#include <iostream>
#include <complex>
#include <cmath>
#include<math.h>
#include<cmath>
#include<Eigen/Dense>
#include<fstream>
#include<random>

#define pi 3.14159

using namespace std;
using namespace Eigen;

double deltaf(double x,double sigma)
	{
		double a,b;
		a=exp(-(0.5*x*x)/(sigma*sigma));
		b=a/(sqrt(2.0*pi*sigma*sigma));
		return b;
	}
	
int main()
	{
	//defining the hamiltonian matrix
	int n=1000,h=1,i,j,k,count,countl=0;;
	double y, ymax=20, sigma=0.01,sum,sum1[40001];
	for(i=0;i<40001;i++)
		sum1[i]=0.0;
	ofstream myfile;
	myfile.open("dos3.dat",ios::app);
	MatrixXcd KE_ham(n,n), Pot_ham(n,n), ham(n,n);
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<double> unif(0,3);
	for(k=0;k<75;k++)
	{
		for(i=0;i<n;i++)
			{
			for(j=0;j<n;j++)
				{
				if(i==j)
					{
					KE_ham(i,j)=0.0;
					Pot_ham(i,j)=unif(gen);
					}
				else if(fabs(i-j)==1.0||(int)fabs(i-j)==(n-1))
					{
					KE_ham(i,j)=1.0;
					Pot_ham(i,j)=0.0;
					}
				else
					{
					KE_ham(i,j)=Pot_ham(i,j)=0.0;
					}
				}
			}
		ham=KE_ham+Pot_ham;
		
		//diagonalising the hamiltonian
		SelfAdjointEigenSolver<MatrixXcd> es(ham);
		complex<double> lambda[n];
		MatrixXcd evec(n,n);
		for(i=0; i<n;i++)
			{
			lambda[i]=es.eigenvalues()[i];
			//cout<<"the numerical eigenvalue is \t"<<lambda[i]<<endl; //eigenvalues obtained
			}
		
		//proceeding to calculate the DoS
		count=0;
		for(y=-ymax;y<=ymax;y=y+0.001)
			{
				sum=0.0;
				for(i=0;i<n;i++)
					{
						if(fabs(y-real(lambda[i]))<(60.0*sigma))
							{
								sum=sum+deltaf((y-real(lambda[i])),sigma);
							}
					}
				sum1[count]=sum1[count]+sum/n;
				count++;
			}
			countl++;
			cout<<countl<<endl;
	}
		count=0;
		for(y=-ymax;y<=ymax;y=y+0.001)
		{
		myfile<<y<<"\t"<<sum1[count]/75.0<<endl;
		count++;
		}
		myfile.close();
						
	}