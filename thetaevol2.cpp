#include <iostream>
#include <complex>
#include <cmath>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<Eigen/Dense>
#include<fstream>
#include<random>

#define pi 3.14159

using namespace std;
using namespace Eigen;

int main()
	{
		//defining the Hamiltonian
		int n=1000,i,j,k,count,count1=0,m;
		double p;
		double dtheta,R;
		double fin_time=500;
		MatrixXd occ_prob(n,1);
		cout<<"input the site no. of the initial wavefunction \n";
		cin>>m;
		R=n/(2.0*pi);
		dtheta=(2*pi)/n;
		MatrixXcd KE_ham(n,n), Pot_ham(n,n), ham(n,n);
		double theta[n];
		for(i=0;i<n;i++)
			{
				if(i<=n/2)
				theta[i]=dtheta*i;
				else if(i>n/2)
				theta[i]=-theta[n-i];
				//cout<<"theta for "<<i+1<<" th site is "<<theta[i]<<endl;				
			}
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> unif(-5,5);
		for(i=0;i<n;i++)
		occ_prob(i)=0.0;
		ofstream myfile1;
		myfile1.open("thetaevolt500.dat",ios::app);
		for(k=0;k<30;k++)
					{
						for(i=0;i<n;i++)
							{
								for(j=0;j<n;j++)
									{
										if(i==j)
											{
											KE_ham(i,j)=-2.0;
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
							//cout<<"the potential hamiltonian at "<<k<<" th turn is \n"<<Pot_ham<<endl;
						ham=(-KE_ham)/(2.0*dtheta*dtheta*R*R)+Pot_ham;//hamiltonian defined
						//diagonalising the hamiltonian
						SelfAdjointEigenSolver<MatrixXcd> es(ham);
						complex<double> lambda[n];
						MatrixXcd evec(n,n);
						for(i=0; i<n;i++)
							{
							lambda[i]=es.eigenvalues()[i];
							//cout<<"the eigenvalue is \t"<<lambda[i]<<endl;
							}
						evec=es.eigenvectors();
						//checking for the time evolution
							VectorXcd coeff(n);
							coeff=evec.row(m-1);
							//finding out the time evolution
							MatrixXcd time_evol(n,1);
							MatrixXcd temp(n,1);
							complex<double> t;
							count=0;
								t=-1i*fin_time;
								for(i=0;i<n;i++)
								time_evol(i)=0.0;	
								for(i=0;i<n;i++)
									{
										temp=evec.col(i);
										time_evol=time_evol+coeff(i)*exp(lambda[i]*t)*temp;
									}
								for(i=0;i<n;i++)
									{
										occ_prob(i)=occ_prob(i)+abs(time_evol(i))*abs(time_evol(i));
									}
						   cout<<count1<<endl;
						   count1++;
					}
					//taking the outputs in file
					for(i=(n/2)+1;i<n;i++)
						{
							myfile1<<R*theta[i]<<"\t"<<occ_prob(i)/30.0<<"\n";
						}
					for(i=0;i<=n/2;i++)
						{
							myfile1<<R*theta[i]<<"\t"<<occ_prob(i)/30.0<<"\n";
						}
						
					myfile1.close();	
	}			
	