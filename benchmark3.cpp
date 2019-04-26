#include <iostream>
#include <complex>
#include <cmath>
#include<math.h>
#include<Eigen/Dense>
#include<fstream>

#define pi 3.14159

using namespace std;
using namespace Eigen;

int main()
	{
	//defining the hamiltonian matrix
	int n=6,i,j,k;
	double dtheta,R;
	R=n/(2.0*pi);
	dtheta=(2.0*pi)/n;
	MatrixXcd KE_ham(n,n);
	MatrixXcd Pot_ham(n,n);
	MatrixXcd ham(n,n);
	for(i=0;i<n;i++)
		{
		for(j=0;j<n;j++)
			{
			if(i==j)
				{
				KE_ham(i,j)=-2.0;
				Pot_ham(i,j)=0.0;
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

	ham=(-KE_ham)/(2.0*dtheta*dtheta*R*R)+Pot_ham;
	//cout<<KE_ham<<endl<<ham<<endl;
	//diagonalising the hamiltonian
	SelfAdjointEigenSolver<MatrixXcd> es(ham);
	complex<double> lambda[n];
	MatrixXcd evec(n,n);
	cout<<"the numerical eigenvalues are \n";
	for(i=0; i<n;i++)
		{
		lambda[i]=es.eigenvalues()[i];
		cout<<lambda[i]<<endl;
		}
	double energy[n];
	cout<<"the theoretical eigenvalues are \n";
	for(i=0;i<n;i++)
		{
			energy[i]= -(cos(i*dtheta)-1)/(dtheta*dtheta*R*R);
			cout<<energy[i]<<endl;
		}
	evec=es.eigenvectors();
	cout<<"the eigenvectors are \n"<<evec<<endl;
	complex<double> theta[n];
	for(i=0;i<n;i++)
		{
			theta[i]=dtheta*i;
			//cout<<theta[i]<<endl;
		}
		MatrixXcd eigenf(n,n);
		complex<double> temp;
		for(i=0;i<n;i++)
			{
				for(j=0;j<n;j++)
					{
						temp=1i*j;
						eigenf(i,j)=(exp(theta[j]*temp))/sqrt(n);
					}
			}
			//cout<<"the primary eigenfunctions are \n"<<eigenf<<endl;
			for(i=0;i<n;i++)
				{
					eigenf(i,1)=cos(theta[i]*1.0);
					eigenf(i,5)=sin(theta[i]*1.0);
					eigenf(i,2)=cos(theta[i]*2.0);
					eigenf(i,4)=sin(theta[i]*2.0);
				}
				//cout<<"the theoretical eigenvectors are \n"<<eigenf<<endl;
			VectorXcd temp1(n);
			VectorXcd temp2(n);
		MatrixXcd identity1(n,n);		
		for(i=0;i<n;i++)
			{
				for(j=0;j<n;j++)
					{
						temp1=eigenf.col(i);
						temp2=eigenf.col(j);
						identity1(i,j)=temp1.dot(temp2);
					}
			}		
		//cout<<"dot products of the numerically computed eigenvectors \n"<<identity1<<endl;
		MatrixXcd eigenf_ortho(n,n);
		eigenf_ortho=eigenf.householderQr().householderQ();
		cout<<"orthonormalized eigenfunctions \n"<<eigenf_ortho<<endl;
		for(i=0;i<n;i++)
			{
				for(j=0;j<n;j++)
					{
						temp1=eigenf_ortho.col(i);
						temp2=eigenf_ortho.col(j);
						identity1(i,j)=temp1.dot(temp2);
					}
			}	
		//cout<<"dot products of the orthonormalized theoretical eigenvectors \n"<<identity1<<endl;						
		
	}		