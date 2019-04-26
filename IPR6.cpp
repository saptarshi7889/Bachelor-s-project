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
//defining the hamiltonian matrix
	int n,i,j,k;
	ofstream myfile;
	myfile.open("IPR.dat",ios::app);
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> unif(-5,5);
	//myfile<<"set xlabel 'L'"<<endl<<"set ylabel 'IPR'"<<endl<<"unset autoscale y"<<endl<<"plot 'IPR.dat' u 1:2 w l \n";
	for(n=100;n<=4000;n=n+300)
		{
			double IPR_avg=0.0;
			for(k=0;k<5;k++)
				{
					double dtheta,R;
					R=n/(2.0*pi);
					dtheta=(2.0*pi)/n;
					MatrixXcd KE_ham(n,n);
					MatrixXcd Pot_ham(n,n);
					MatrixXcd ham(n,n);
					srand(time(NULL));
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
					//cout<<"the randomized potential hamiltonian is \n"<<Pot_ham<<endl;
					ham=(-KE_ham)/(2.0*dtheta*dtheta*R*R)+Pot_ham;
				
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
				//calculating the IPR of a random band
				int k=n/100;
				int evec_level_max=55*k-1;
				int evec_level_min=54*k-1;
				double sum,IPR=0.0;
				for(i=evec_level_min;i<=evec_level_max;i++)
				{
					sum=0.0;
					for(j=0;j<n;j++)
					{
					sum=sum+pow(abs(evec(j,i)),4);
					}
					IPR=IPR+(1.0/sum);
				}
				IPR_avg=IPR_avg+(IPR/k);
			}
				cout<<"The system size is "<<n<<endl;
				double IPR_avg_avg;
				IPR_avg_avg=IPR_avg/5.0;
				cout<<"average IPR of the 55th band is "<< IPR_avg_avg<<endl;
				//taking the outputs in a file
				myfile<<n<<"\t"<<IPR_avg_avg<<endl;
		}
		myfile.close();
	}