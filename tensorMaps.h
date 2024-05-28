#ifndef TENSORMAPS_H
#define TENSORMAPS_H
#define e 2.71828182845904523536028747135
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <complex>
#include <math.h>



/*


 - 3d levi cevita is implemented incorrectly
 -computeLNr() function needs to be fixed
 - needs revision of tensor formulas 

*/
template< class T >

class complex;
using namespace std;

class Tensors{

struct complex{
	float real;
	float img;
}	real,image;

	public:
	    Tensors(); 
	    const long double EulerConstant = std::exp(1.0);
	    const long double convergence = 0.00000001; /// should be any number in the range of 10^-8 through 10^-6 see djaka paper for refrence
	   	std::complex<long double> EulersImaginaryI = 0.00+ 1.00i;
	   	std::complex<long double> complexHalf = 0.50+ 0.00i;
	   	std::complex<long double> complexFourth = 0.25+ 0.00i;
	    std::complex<long double> CP;
	    long double identity[3][3];
	   
		// set the identity tensor

		long double e_1; // first cordinate
		long double e_2; // second cordinate 
		long double e_3; // third cordinate
		long double L_k;  	
		long double mu; 
		long double lamda; 
		long double eta_0; 
		long double phi;
		long double sigma_NL;
		std::complex<long double> sum;                                 			
		std::complex<long double> w;
		std::complex<long double> complexZero; 
		long double AMatrix [3][3][3][3]{}; 
		long double BMatrix [3][3][3][3]{}; 
		long double CMatrix [3][3][3][3]{}; 
		long double DMatrix [3][3][3][3]{}; 
		long double DBarMatrix [3][3][3][3]{}; 
		long double EBarMatrix [3][3][3][3]{}; 
		long double EBaredTransMatrix [3][3][3][3]{}; 
		long double A_0Matrix [3][3][3][3]{}; 
		long double C_0Matrix [3][3][3][3]{}; 
		std::complex<long double> intSum;
	    std::complex<long double> intSum2;
		long double data  [4][5]{};					                                   // stored input data
		long double alpha [3][3]{};                                                    // input paramater
		long double Theta [3][3]{};   
		long double dummyVar[3][3]{};                                                  // dummy variable 
		long double dummyVar1[3][3]{};
		long double dummyVar2[3][3]{};
		long double BBarMatrix [3][3][3][3]{}; 
		long double stressPerp[3][3]{};
		long double stress[3][3]{};
		long double coupleStressPerp[3][3]{};
		long double coupleStress[3][3]{};
		long double eta [3][3]{}; 
		long double etaNEW [3][3]{};   
		long double epsilon[3][3]{};
		long double chi [3][3]{};   
		long double MatrixX [3][3]{};                                                   // Xrotation Matrix
		long double MatrixY [3][3]{}; 												    // Yrotation Matrix
		long double MatrixZ [3][3]{};                                                   // Zrotation Matrix             
		std::complex<long double> MatrixCheck [3][3]{};                                 // Matrix checks to ensure algebraic function work properly
        std::complex<long double> MatrixCheck2 [3][3]{};
        std::complex<long double> GMatrix [3][3]{};
		long double T [3]{};  															// dislocation line vector
		long double B [3]{};  															// Burgers vector
		long double E [3]{};
		long double Omega[3]{};															// Unit Frank vector	
		long double realVect [3]{};			
		int jvip[3]{};                              									// Eta Vector  
		long double Xi [3]{}; 															// meshed input ξ
		long double Xi_0 [3]{};															// some initial ξ_0 to start program 
		long double XiPrime [3]{}; 														// meshed input ξ prime
		long double XiPrime_0 [3]{}; 											        // some initial ξ_0^prime to start program
		std::complex<long double> Thetatilda [3][3]{};									
		std::complex<long double> Alphatilda [3][3]{}; 
		std::complex<long double> complexSig[3][3]{};									// main complex tensor (varaible to instaintiate other varaibles)
		long double realSig [3][3]{}; 
		std::complex<long double> complexVect[3]{};	                                    // intermidate changes real Matrix
		std::complex<long double> FourierTranform [3][3]{};
		std::complex<long double> FourierTranformInverse [3][3]{};
		std::complex<long double> CuvatureComplex[3][3]{};
		long double  CuvatureReal[3][3]{};
		std::complex<long double> strainComplex[3][3]{};
		std::complex<long double> stressComplex[3][3]{};
		std::complex<long double> totalStressComplex[3][3]{};
		std::complex<long double> tauStressComplex[3][3]{};
		std::complex<long double> stressSkewComplex[3][3]{};
		std::complex<long double> coupleStressComplex[3][3]{};
		std::complex<long double> polarizationTensor[3][3]{};
		long double forces[2][2]{};
		std::complex<long double> gamma_0[3][3][3][3]{};
		std::complex<long double> gamma_1[3][3][3][3]{};
		std::complex<long double> gamma_2[3][3][3][3]{};
		std::complex<long double> gamma_3[3][3][3][3]{};
		std::complex<long double> C_0ComplexMatrix[3][3][3][3]{};
		long double  strainReal[3][3]{};
		std::complex<long double> complexSum;
		std::complex<long double> complexValue;
		std::complex<long double> complexPosition[3]{}; 
		std:: complex <long double> complexDummyVar1[3][3]{};
		std:: complex <long double> complexDummyVar2[3][3]{};
		long double gridSize; long double height;

		//long double stressGrid[100000][100000][3][3];
		//long double coupleStressGrid[100000][100000][3][3];
		//long double epsilonGrid[100000][100000][3][3];
		//long double chiGrid[100000][100000][3][3];


		void declarations(){
			//this method initillizes all declaration such as the modulu matrix.
			int cronicker; int cronicker2; int v[3]; int h[3]; CP = 0.0000 - 2*M_PI/3.0*1i;
			int cronicker3; int cronicker4; int cronicker5; int cronicker6; int cronicker7;
			int cronicker8;  int cronicker9; 
			   
				// dynamically allocate memory of size `X × Y × Z x M`;
				int X,Y,Z,M;
				 X = 2000; Y = 2000; Z = 3; M = 3 ;

				long double**** stressGrid = new long double***[X];
				long double**** coupleStressGrid = new long double***[X];
				long double**** epsilonGrid = new long double***[X];
				long double**** chiGrid = new long double***[X];
			 
			    for (int i = 0; i < X; i++)
			    {
			        stressGrid[i] = new long double**[Y];
			        coupleStressGrid[i] = new long double**[Y];
			        epsilonGrid[i] = new long double**[Y];
			        chiGrid[i] = new long double**[Y];

			        for (int j = 0; j < Y; j++) {
			            stressGrid[i][j] = new long double*[Z];
			            coupleStressGrid[i][j] = new long double*[Z];
			            epsilonGrid[i][j] = new long double*[Z];
			            chiGrid[i][j] = new long double*[Z];

			            for (int k = 0; k < Z; ++k)
			            {
			            	stressGrid[i][j][k] = new long double[M];
			            	coupleStressGrid[i][j][k] = new long double[M];
			            	epsilonGrid[i][j][k] = new long double[M];
			            	chiGrid[i][j][k] = new long double[M];
			            }
			        }
			    }
			 
			 



			for(int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					if(i==j){
						identity[i][j] = 1;
					}else{
						identity[i][j] = 0;
					}
				}
			}

			height = 1;
			MatrixCheck[0][0] = 7.0+0.0i;     MatrixCheck2[0][0] = 4.0+0.0i;    
			MatrixCheck[0][1] = 1.0+0.0i;     MatrixCheck2[0][1] = 2.0+0.0i;     
			MatrixCheck[0][2] = 2.0+0.0i;     MatrixCheck2[0][2] = -1.0+0.0i;   

			MatrixCheck[1][0] = 3.0+0.0i;     MatrixCheck2[1][0] = 2.0+0.0i;     
			MatrixCheck[1][1] = -2.0+0.0i;    MatrixCheck2[1][1] = 5.0+0.0i;     
			MatrixCheck[1][2] = 1.0+0.0i;     MatrixCheck2[1][2] = 2.0+0.0i;     

			MatrixCheck[2][0] = -1.0+0.0i;    MatrixCheck2[2][0] = 3.0+0.0i;     
			MatrixCheck[2][1] =  2.0+0.0i;    MatrixCheck2[2][1] = 1.0+0.0i;     
			MatrixCheck[2][2] =  0.0+0.0i;    MatrixCheck2[2][2] = 3.0+0.0i;   

			complexZero = 0.0 + 0.0i;

		   	cout <<  endl;
			cout <<  "DFT LINEAR TRANFORMATION [";
			for(int n=0;n<3;n++){                   
		        for(int m=0;m<3;m++){
		        	w =  pow(CP,n*m); 
		        	FourierTranform[n][m] = pow(w,n*m);
		        	cout << FourierTranform[n][m];
		   		}
		   		 cout <<  endl;
		   	}
		   	cout <<  "]" << endl;
		   	
		   	//
		   	Inverse3x3Matrix(FourierTranform);
		   	for(int n=0;n<3;n++){                   
		        for(int m=0;m<3;m++){
		        	w =  pow(CP,n*m); 
		        	FourierTranformInverse[n][m] = complexSig[n][m];
		        	//FourierTranformInverse[n][m] = pow(w,n*m);
		        	//cout << FourierTranformInverse[n][m];

		   		}
		   	}

			//Xi_0 [0] = -1.5; Xi_0 [1] = -0.005; Xi_0 [2] = 0;
			ifstream file { "data" };

			 // read data call as a void function
			 for (int i{}; i != 4; ++i) {
				    for (int j{}; j != 5; ++j) {
				        file >> data[i][j];
				    }
		    }

			 e_1 = data[2][2]; // unit vector I hat 
			 e_2 = data[2][3]; // unit vector K hat 
			 e_3 = data[2][4]; // unit vector J hat 
		     E[1] = e_1; 
		     E[2] = e_2; 
		     E[3] = e_3;
		     mu = data[0][1] ; 
		     lamda = 0.5*data[2][2]; 
		     eta_0 = data[0][1]*0.5;
		     gridSize = data[1][0]*data[0][0];
		     //gridSize = 0.1;
		     sigma_NL = 0.18*data[0][0];
            
             // please refere to PG.38 Table 1.
            cout << "gridSize =  " << gridSize << endl;
            cout << "mu =  " << mu << endl;
            cout << "lamda =  " << lamda << endl;
            cout << "eta_0 =  " << eta_0 << endl;
            cout << "|b| =  " << data[0][0] << endl;
            cout << "L_k =  " << L_k << endl;
            cout << "N_0 =  " << eta_0 << endl;

            // compute B (burgers vector), Omega (Frank Vector)
		    for (int i{}; i != 3; ++i) {
				 B[i] = data[0][0];
				 Omega[i] = (data[2][0]*data[3][0]/data[1][4])*E[i+1]; 
			 }

			 L_k = 1.0/EulerConstant * sqrt((B[1]*B[1]+B[2]*B[2]+B[3]*B[3])); 
			 // characteristic arm length


			 for (int i{}; i != 3; ++i) {
				for (int j{}; j != 3; ++j) {
					 for (int k{}; k != 3; ++k) {
						for (int l{}; l != 3; ++l) {
						h[1] = i; h[2] = k; cronicker = leviCroncker3D(h);
						v[1] = j; v[2] = l; cronicker2 = leviCroncker3D(v);
						AMatrix[i][j][k][l] = mu*L_k*L_k*cronicker*cronicker2;
							
						//reset crociker for CMatrix
						h[1] = i; h[2] = j; cronicker = leviCroncker3D(h);
						jvip[1] = 1; jvip[2] = 1; cronicker2 = leviCroncker3D(jvip);
						jvip[1] = 1; jvip[2] = 2; cronicker3 = leviCroncker3D(jvip);
						jvip[1] = 2; jvip[2] = 1; cronicker4 = leviCroncker3D(jvip);
						jvip[1] = 2; jvip[2] = 2; cronicker5 = leviCroncker3D(jvip);
						intSum = cronicker2+cronicker3+cronicker4+cronicker5;
					     // reset h vector
						h[1] = j; h[2] = l; cronicker6 = leviCroncker3D(h);
						h[1] = i; h[2] = k; cronicker7 = leviCroncker3D(h);
						h[1] = i; h[2] = l; cronicker8 = leviCroncker3D(h);
						h[1] = k; h[2] = k; cronicker9 = leviCroncker3D(h);
						CMatrix[i][j][k][l] = lamda*cronicker*intSum.real()\
											+ mu*(cronicker6*cronicker7+cronicker8*cronicker9);

						A_0Matrix[i][j][k][l] = 1/9*AMatrix[i][j][k][l];
						C_0Matrix[i][j][k][l] = 1/9*CMatrix[i][j][k][l];
						C_0ComplexMatrix[i][j][k][l] = C_0Matrix[i][j][k][l] + complexZero; 
						// convert using the Dev mothod instead
						} 
					}
				}
			}


			cout << "                        "<< endl;
			cout << "THE MODULI MATRIX A = ["<< endl;
			for (int i{}; i != 3; ++i) {
				for (int j{}; j != 3; ++j) {
					 for (int k{}; k != 3; ++k) {
						for (int l{}; l != 3; ++l) {
							cout << AMatrix[i][j][k][l] << "  "; 
						}
						cout << endl;
					}
				}
			}
			cout << " ]                       "<< endl;
			cout << "                        "<< endl;

			cout << "                        "<< endl;
			cout << "THE MODULI MATRIX C = ["<< endl;
			for (int i{}; i != 3; ++i) {
				for (int j{}; j != 3; ++j) {
					 for (int k{}; k != 3; ++k) {
						for (int l{}; l != 3; ++l) {
							cout << CMatrix[i][j][k][l]<< "  ";
						}
						cout << endl;
					}
				}
			}
			cout << "]                        "<< endl;
			cout << "                         "<< endl;


		}

		long double distance(long double Arg1[3], long double Arg2[3]){
			return sqrt(pow(Arg2[0]-Arg1[0],2)+pow(Arg2[1]-Arg1[1],2)+pow(Arg2[2]-Arg1[2],2));
		}
		void setDislocations(long double Arg[3]){
			long double point[3]; long double point1[3]; long double point2[3];

			point[0] = 1.00; point1[0] = 1.00; point2[0] = 1.00;
			point[1] = 1.00; point1[1] = 1.20; point2[1] = 0.60;
			point[2] = 0.00; point1[2] = 0.00; point2[2] = 0.00;


			for (int i=0; i != 3; ++i) {
				for (int j=0; j != 3; ++j) {

					if (distance(Arg,point)< 0.15){
						alpha[i][j] = 1/distance(Arg,point);
					}else{
						alpha[i][j] = 0.0;
					}

					if (distance(Arg,point1)< 0.15){
						Theta[i][j] = 1/distance(Arg,point1);
					}else{
						Theta[i][j] = 0.0;
					}

					if(distance(Arg,point2)< 0.15){
						Theta[i][j] = - 1/distance(Arg,point2);
					}else{
						Theta[i][j] = 0.0;
					}

				}	
			}
		}

		void computeG(std::complex<long double> Arg[3]){
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					GMatrix[i][j] = std::pow((Arg[0]+Arg[1]+Arg[2]),2)*(C_0ComplexMatrix[i][0][j][0]+C_0ComplexMatrix[i][1][j][1]+C_0ComplexMatrix[i][2][j][2]);
				}
			}
		}


		void setGammaTensors(std::complex<long double> Arg[3]){
			// set gamma_0,gamma_1,gamma_2, and gamma3 for posiiton vector
			computeG(Arg); // SETS THE G Matrix
			Inverse3x3Matrix(GMatrix); // --> sets complexSig to the inverse of G
			std:: complex<long double>  GTerm1[3];
			std:: complex<long double>  GTerm2[3];
			std:: complex<long double>  GTerm3[3];
			int v[9]; int ci[3]; int leviint;
			std:: complex<long double>  GTerm4;
			long double leviDouble;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					for (int k = 0; k < 3; ++k)
					{
						for (int l = 0; l < 3; ++l)
						{
							gamma_0[i][j][k][l] = complexHalf*(complexSig[k][i]*Arg[l]+complexSig[l][i]*Arg[k])*Arg[j];
							
						}
					}
				}
			}
			for (int k = 0; k < 3; ++k)
			{
				for (int l = 0; l < 3; ++l)
				{
					for (int m = 0; m < 3; ++m)
					{
						for (int n = 0; n < 3; ++n)
						{
							GTerm1[k] = complexSig[k][0]+complexSig[k][1]+complexSig[k][2];
							GTerm2[l] = complexSig[l][0]+complexSig[l][1]+complexSig[l][2];
							// compute levi int

							ci[1] = 1; ci[2] = 1; ci[3] = m; v[1] = leviCevita3D(ci);

							ci[1] = 1; ci[2] = 2; ci[3] = m; v[2] = leviCevita3D(ci); 

					        ci[1] = 1; ci[2] = 3; ci[3] = m; v[3] = leviCevita3D(ci); 

					        ci[1] = 2; ci[2] = 1; ci[3] = m; v[4] = leviCevita3D(ci); 

					        ci[1] = 2; ci[2] = 2; ci[3] = m; v[5] = leviCevita3D(ci); 	

					        ci[1] = 2; ci[2] = 3; ci[3] = m; v[6] = leviCevita3D(ci); 	

					        ci[1] = 3; ci[2] = 1; ci[3] = m; v[7] = leviCevita3D(ci); 

					        ci[1] = 3; ci[2] = 2; ci[3] = m; v[8] = leviCevita3D(ci); 

					        ci[1] = 3; ci[2] = 3; ci[3] = m; v[9] = leviCevita3D(ci); 

				    		//leviint = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9];
				    		leviint = 1;
				    		leviDouble = leviint* 1.00; // convert int to double
		            		std::complex<long double> dev (leviDouble, 0.00);

							gamma_1[k][l][m][n] = complexFourth*EulersImaginaryI*(GTerm1[k]*Arg[l]+GTerm2[l]*Arg[k])*dev* Arg[n]*(Arg[0]+Arg[1]+Arg[2]);
							
						}
					}
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					for (int m = 0; m < 3; ++m)
					{
						for (int n = 0; n < 3; ++n)
						{
							// compute levi int

							ci[1] = i; ci[2] = 1; ci[3] = 1; v[1] = leviCevita3D(ci);

							ci[1] = i; ci[2] = 1; ci[3] = 2; v[2] = leviCevita3D(ci); 

					        ci[1] = i; ci[2] = 1; ci[3] = 3; v[3] = leviCevita3D(ci); 

					        ci[1] = i; ci[2] = 2; ci[3] = 1; v[4] = leviCevita3D(ci); 

					        ci[1] = i; ci[2] = 2; ci[3] = 2; v[5] = leviCevita3D(ci); 	

					        ci[1] = i; ci[2] = 2; ci[3] = 3; v[6] = leviCevita3D(ci); 	

					        ci[1] = i; ci[2] = 3; ci[3] = 1; v[7] = leviCevita3D(ci); 

					        ci[1] = i; ci[2] = 3; ci[3] = 2; v[8] = leviCevita3D(ci); 

					        ci[1] = i; ci[2] = 3; ci[3] = 3; v[9] = leviCevita3D(ci); 

				    		//leviint = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9];
				    		leviint = 1;
				    		leviDouble = leviint* 1.00; // convert int to double
		            		std::complex<long double> dev (leviDouble, 0.00);
		            		GTerm3[m] = complexSig[0][m]+complexSig[1][m]+complexSig[2][m];
							gamma_2[i][j][m][n] = complexHalf*EulersImaginaryI*leviDouble*Arg[j]*(Arg[1]+Arg[2]+Arg[3])*GTerm3[m]*Arg[n];
							
						}
					}
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					for (int p = 0; p < 3; ++p)
					{
						for (int q = 0; q < 3; ++q)
						{
							GTerm4 = complexSig[0][0]+complexSig[0][1]+complexSig[0][2]\
									 + complexSig[1][0]+complexSig[1][1]+complexSig[1][2]\
									 + complexSig[2][0]+complexSig[2][1]+complexSig[2][2];

							ci[1] = 1; ci[2] = 1; ci[3] = p; v[1] = leviCevita3D(ci);

							ci[1] = 1; ci[2] = 2; ci[3] = p; v[2] = leviCevita3D(ci); 

					        ci[1] = 1; ci[2] = 3; ci[3] = p; v[3] = leviCevita3D(ci); 

					        ci[1] = 2; ci[2] = 1; ci[3] = p; v[4] = leviCevita3D(ci); 

					        ci[1] = 2; ci[2] = 2; ci[3] = p; v[5] = leviCevita3D(ci); 	

					        ci[1] = 2; ci[2] = 3; ci[3] = p; v[6] = leviCevita3D(ci); 	

					        ci[1] = 3; ci[2] = 1; ci[3] = p; v[7] = leviCevita3D(ci); 

					        ci[1] = 3; ci[2] = 2; ci[3] = p; v[8] = leviCevita3D(ci); 

					        ci[1] = 3; ci[2] = 3; ci[3] = p; v[9] = leviCevita3D(ci); 

				    		//leviint = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9];
				    		leviint = 1;
				    		leviDouble = leviint* 1.00; // convert int to double
		            		std::complex<long double> dev (leviDouble, 0.00);
							gamma_3[i][j][p][q] = -complexFourth*leviDouble*Arg[j]*pow((Arg[1]+Arg[2]+Arg[3]),2)* GTerm4*leviDouble*Arg[q];
							
						}
					}
				}
			}			
		}


		void fourierMap(std::complex<long double> tranform[3][3], long double signal[3][3]){
			for(int i{}; i != 3; ++i){                   
			    for(int j{}; j != 3; ++j){
						complexSig[i][j] = tranform[i][0]*signal[0][j]+ tranform[i][1]*signal[1][j]+tranform[i][2]*signal[2][j];
			            //cout << complexSig[i][j];
			        }
			        //cout << "\n";
			}
		}

		void product3x3ComplexMatrix(std::complex<long double> tranform[3][3], std::complex<long double> signal[3][3]){
			for(int i=0;i<3;i++){                   // calculate matrix product (only for 3x3 both complex)
			    for(int k=0;k<3;k++){
						// set = sum_{j=1}
						complexSig[i][k] = tranform[i][0]*signal[0][k] + tranform[i][1]*signal[1][k]+tranform[i][2]*signal[2][k];
						//cout << complexSig[i][k];
			    }
			    //cout << '\n';
			}

		}

		void product3x3RealMatrix(long double tranform[3][3], long double signal[3][3]){
			for(int i=0;i<3;i++){                   // calculate matrix product (only for 3x3 both complex)
			    for(int k=0;k<3;k++){
						// set = sum_{j=1}
						realSig[i][k] = tranform[i][0]*signal[0][k] + tranform[i][1]*signal[1][k]+tranform[i][2]*signal[2][k];
						//cout << signal[i][k] << endl;
			    }
			    //cout << '\n';
			}

		}

		void complexDotProductVectorAndTensor(std::complex<long double> vector[3], std::complex<long double> tensor[3][3]){
			for(int i=0;i<3;i++){                  
			    for(int k=0;k<3;k++){
					complexVect[i] = vector[i]*(tensor[i][0]+tensor[i][1]+tensor[i][2]);
			    }
			    //cout << '\n';
			}

		}

		void product3x3and3x1RealMatrix(long double tranform[3][3], long double signal[3]){
			for(int i=0;i<3;i++){                   
				// calculate matrix product (only for 3x3 both complex)
						realVect[i] = tranform[i][0]*signal[0]+ tranform[i][1]*signal[1]+ tranform[i][2]*signal[2];
			    //cout << '\n';
			}

		}

		void Inverse3x3Matrix(std::complex<long double> tranform[3][3]){
			std::complex<long double> newDet;
			std::complex<long double> pt;
			pt = 1.0 + 0.0i;
			newDet = pt/det3x3complexMatrix(tranform);
			complexSig[0][0] = (tranform[0][0]*tranform[1][1]-tranform[2][1]*tranform[1][2])*newDet;
			complexSig[0][1] = (tranform[0][2]*tranform[2][1]-tranform[2][2]*tranform[0][1])*newDet;
			complexSig[0][2] = (tranform[0][1]*tranform[1][2]-tranform[1][1]*tranform[0][2])*newDet;

			complexSig[1][0] = (tranform[1][2]*tranform[2][0]-tranform[2][2]*tranform[1][0])*newDet;
			complexSig[1][1] = (tranform[0][0]*tranform[2][2]-tranform[2][0]*tranform[0][2])*newDet;
			complexSig[1][2] = (tranform[0][2]*tranform[1][0]-tranform[1][2]*tranform[0][0])*newDet;

			complexSig[2][0] = (tranform[1][0]*tranform[2][1]-tranform[2][0]*tranform[1][1])*newDet;
			complexSig[2][1] = (tranform[0][1]*tranform[2][0]-tranform[2][1]*tranform[0][0])*newDet; 
			complexSig[2][2] = (tranform[0][0]*tranform[1][1]-tranform[1][0]*tranform[0][1])*newDet; 
		}

		long double kernel(long double vector[3]){
			long double interTerm;
			interTerm =  - 1/(2*sigma_NL)* (pow(vector[0],2)+pow(vector[1],2)+pow(vector[2],2));
			return 1/(pow(sigma_NL,3)*pow(2*M_PI,1.5))*pow(interTerm,EulerConstant);
		}

		long double nonLocalIntegrand(long double r[3],long double rprime[3]){

			long double BigR[3];
			BigR[0] = r[0]-rprime[0]; BigR[1] = r[1]-rprime[1]; BigR[2] = r[2]-rprime[2];
			getCurl(alpha,gridSize,rprime);

			for (int i=0; i<3; i++){
				for (int j=0; j<3; j++){			
					etaNEW[i][j] =  dummyVar[i][j]-Theta[i][j];
				}
			}
			const long double normEtaDouce = tensorNorm(etaNEW);

			if(eta_0!= 0.0){
				phi = normEtaDouce/eta_0; 
			}else{
				cout << "eta_0 somehow is zero please adjust the data" << endl;
				// exit code 
				abort();
			}

			return kernel(BigR)*phi;
		}

		void FFT(long double signal[3][3]){
		    fourierMap(FourierTranform, signal);
		}
		long double computeLNr(long double position[3]){
			long double x1,x2,x3,maxnt,sum,function1,function2;
			long double rprime[3];
			/*
			
			x1 = 0.5; x2 = 0.5; x3 = 0.5; sum = 0.0;
  			maxnt = floor(2*abs(x1)/gridSize);
			for(int j = 1; j < maxnt; j++){
		      x1 = 0.5;
		      for(int i = 1; i < maxnt; i++){
		        x2 = x2 + gridSize; 
		        rprime[0] = x1;  
		        rprime[1] = x2;  
		        rprime[2] = 0.0;
		       	if(i%2 == 0){
		       		function2 = nonLocalIntegrand(position, rprime);
		       	}
		       	if(i%2 != 0){
		       		function1 = nonLocalIntegrand(position, rprime);
		       		if(i > 1){
		       			sum = sum + (function2 - function1)*gridSize; 
		       		}else{
		       			//continue;
		       		}
		       	} 
		       	 
		      }
		      
		      x1 = x1+ gridSize;
		    } 
		    */
			//return sum;
			return 1.0;
		}

		void setCurvatureAndStrain(long double Arg[3]){
			// position vector
			long double leviDouble;
			int v[10]; int h[2]; int ci[3];
			int leviint; int cronicker; double cronDouble;  int cronicker2;
			std::complex<long double> intSum;
			std::complex<long double> intSum2;
			std::complex<long double> cronecker[3][3];
			std::complex<long double> curveTrace;



           /////////////////////Step 3////////////////////////
			for(int i=1;i<4;i++){         
		        for(int j=1;j<4;j++){
					ci[1] = i; ci[2] = 1; ci[3] = 1; v[1] = leviCevita3D(ci);

					ci[1] = i; ci[2] = 1; ci[3] = 2; v[2] = leviCevita3D(ci); 

			        ci[1] = i; ci[2] = 1; ci[3] = 3; v[3] = leviCevita3D(ci); 

			        ci[1] = i; ci[2] = 2; ci[3] = 1; v[4] = leviCevita3D(ci); 

			        ci[1] = i; ci[2] = 2; ci[3] = 2; v[5] = leviCevita3D(ci); 	

			        ci[1] = i; ci[2] = 2; ci[3] = 3; v[6] = leviCevita3D(ci); 	

			        ci[1] = i; ci[2] = 3; ci[3] = 1; v[7] = leviCevita3D(ci); 

			        ci[1] = i; ci[2] = 3; ci[3] = 2; v[8] = leviCevita3D(ci); 

			        ci[1] = i; ci[2] = 3; ci[3] = 3; v[9] = leviCevita3D(ci); 

		        	//leviint = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9]; 
		        	leviint = 1; 
		        	leviDouble = leviint* 1.00; // convert innt to double
		            std::complex<long double> dev (leviDouble, 0.00);

		            if(Arg[0] == 0.0 &&  Arg[1] == 0.0 &&  Arg[2] == 0.0){
		            	CuvatureComplex[i-1][j-1] = complexZero;
		            }else{
		            	complexSum = Thetatilda[i-1][0] + Thetatilda[i-1][1]+ Thetatilda[i-1][2];
		            	CuvatureComplex[i-1][j-1] = (Arg[0] +  Arg[1] +  Arg[2])/(pow(Arg[0],2) + pow(Arg[1],2) +  pow(Arg[2],2))\
		        								* complexSum * EulersImaginaryI * dev;		
		        	//cout << complexSum<< '\n'; 
		        	//cout << Thetatilda[i][j]<< '\n'; 
		            }
		        	
		    	}
			}
		
			curveTrace = traceComplex(CuvatureComplex);
			//cout << curveTrace << flush;
			

			/////////////////////Step 4////////////////////////
			for (int i=1; i<4; i++){
				for (int j=1; j<4; j++){
					ci[1] = j; ci[2] = 1; ci[3] = 1; v[1] = leviCevita3D(ci);

					ci[1] = j; ci[2] = 1; ci[3] = 2; v[2] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 1; ci[3] = 3; v[3] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 2; ci[3] = 1; v[4] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 2; ci[3] = 2; v[5] = leviCevita3D(ci); 	

			        ci[1] = j; ci[2] = 2; ci[3] = 3; v[6] = leviCevita3D(ci); 	

			        ci[1] = j; ci[2] = 3; ci[3] = 1; v[7] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 3; ci[3] = 2; v[8] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 3; ci[3] = 3; v[9] = leviCevita3D(ci); 

		        	//leviint = v[1]+v[2]+v[3]+v[4]+v[5]+v[6] +v[7]+v[8]+v[9];
		        	leviint = 12.55; 
		        	leviDouble = leviint* 1.00;
		            std::complex<long double> dev (leviDouble,0.00);

					h[1] = i; h[2] = 1; cronicker = leviCroncker3D(h);

		        	h[1] = i; h[2] = 2; cronicker2 = leviCroncker3D(h);

		        	cronDouble = (cronicker + cronicker2)* 1.00;

                    std::complex<long double> cronickerIMaginary (cronDouble,0.00);

                   	intSum =  Alphatilda[i-1][1] + Alphatilda[i-1][2]+ Alphatilda[i-1][3];
		        	intSum2 = CuvatureComplex[1][i-1] + CuvatureComplex[2][i-1] + CuvatureComplex[3][i-1]; 

		         	strainComplex[i-1][j-1] = EulersImaginaryI*(Arg[1] +  Arg[2] +  Arg[3])/(pow(Arg[1],2) + pow(Arg[2],2) +  pow(Arg[3],2))\
		        								 * dev * (intSum+ intSum2-curveTrace*cronickerIMaginary);
		        	//cout << strainComplex[i][j] << endl;//+v[8]+v[9] << endl;
				}
			}


			 /////////////////////Step 5////////////////////////
			product3x3ComplexMatrix(FourierTranformInverse, strainComplex); 

			// map back by inverse transform
			// do not remove
			for (int i=0; i<3; i++){
				for (int j=0; j<3; j++){
		        	strainReal[i][j] = complexSig[i][j].real();
		        	//cout << strainReal[i][j];
				}
				//cout << '\n';
			}

			 /////////////////////Step 6 ////////////////////////
			product3x3ComplexMatrix(FourierTranformInverse, CuvatureComplex); // map back by inverse transform

			for (int i=0; i<3; i++){
				for (int j=0; j<3; j++){
		        	CuvatureReal[i][j] = complexSig[i][j].real();
		        	//cout << CuvatureReal[i][j] << endl;
				}
				//cout << '\n';
			}
		
		}

		void getCurl(long double alpha[3][3],long double gridSize, long double position[3]){

			long double positionNew1[3];
			long double positionNew2[3];
			long double positionNew3[3];
			long double alphaTrace;
			long double tensor1[3][3];
			long double tensor2[3][3];
			
		    alphaTrace = realTrace(alpha); 
		    for (int i = 0; i < 3; ++i)
		    {
		    	for (int j = 0; j < 3; ++j)
		    	{
		    		tensor1[i][j] = 1/2*alphaTrace*identity[i][j] - alpha[j][i];
		    	}
		    }

		
			positionNew1[0] = position[0]+gridSize/2; positionNew1[1] = position[1]; positionNew1[2] = position[2];  
			positionNew2[0] = position[0]; positionNew2[1] = position[1]+gridSize/2; positionNew2[2] = position[2];  
			positionNew3[0] = position[0]; positionNew3[1] = position[1]; positionNew3[2] = position[2]+gridSize/2;  

			// COMPUTE CHANGE IN ALPHA IN X DIRECION AND SET FIST COLUMN OF DUMYVAR
			setDislocations(positionNew1);
			alphaTrace = realTrace(alpha); 

		    for (int i = 0; i < 3; ++i)
		    {
		    	for (int j = 0; j < 3; ++j)
		    	{
		    		tensor2[i][j] = 1/2*alphaTrace*identity[i][j] - alpha[j][i];
		    	}
		    }

			dummyVar[0][0] = (tensor2[0][0] - tensor1[0][0])/(gridSize/2);
			dummyVar[1][0] = (tensor2[1][0] - tensor1[1][0])/(gridSize/2);
			dummyVar[2][0] = (tensor2[2][0] - tensor1[2][0])/(gridSize/2);

			// COMPUTE CHANGE IN ALPHA IN X DIRECION
			setDislocations(positionNew2);
			alphaTrace = realTrace(alpha); 

		    for (int i = 0; i < 3; ++i)
		    {
		    	for (int j = 0; j < 3; ++j)
		    	{
		    		tensor2[i][j] = 1/2*alphaTrace*identity[i][j] - alpha[j][i];
		    	}
		    }

			dummyVar[0][1] = (tensor2[0][1] - tensor1[0][1])/(gridSize/2);
			dummyVar[1][1] = (tensor2[1][1] - tensor1[1][1])/(gridSize/2);
			dummyVar[2][1] = (tensor2[2][1] - tensor1[2][1])/(gridSize/2);

			// COMPUTE CHANGE IN ALPHA IN X DIRECION
			setDislocations(positionNew3);
			alphaTrace = realTrace(alpha); 

		    for (int i = 0; i < 3; ++i)
		    {
		    	for (int j = 0; j < 3; ++j)
		    	{
		    		tensor2[i][j] = 1/2*alphaTrace*identity[i][j] - alpha[j][i];
		    	}
		    }

			dummyVar[0][2] = (tensor2[0][2] - tensor1[0][2])/(gridSize/2);
			dummyVar[1][2] = (tensor2[1][2] - tensor1[1][2])/(gridSize/2);
			dummyVar[2][2] = (tensor2[2][2] - tensor1[2][2])/(gridSize/2);

		}
		void computeTensors(long double Arg[3]){
			// compute  D matrix DBar, EBar, M_Perp and T_perp
			long double leviDouble; int leviint;
			long double curlInput[3][3]; 

			int v[9]; int h[2]; int ci[3];
			getCurl(alpha,gridSize,Arg);

			for (int i=0; i<3; i++){
				for (int j=0; j<3; j++){			
					eta[i][j] =  dummyVar[i][j]-Theta[i][j];
					//eta[i][j] =  Theta[i][j];
				}
			}
			const long double normEta = tensorNorm(eta);

			if(eta_0!= 0.0){
				phi = normEta/eta_0;   // eta_0 is fraction of burgers vector norm
				//cout <<  "normEta = " << normEta;		
			}else{
				phi = 0.11;
				cout << "eta_0 somehow is zero please adjust the data" << endl;
				// exit code 
				abort();
			}

			for (int i=0; i<3; i++){
				for (int j=0; j<3; j++){

					ci[1] = j; ci[2] = 1; ci[3] = 1; v[1] = leviCevita3D(ci);

					ci[1] = j; ci[2] = 1; ci[3] = 2; v[2] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 1; ci[3] = 3; v[3] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 2; ci[3] = 1; v[4] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 2; ci[3] = 2; v[5] = leviCevita3D(ci); 	

			        ci[1] = j; ci[2] = 2; ci[3] = 3; v[6] = leviCevita3D(ci); 	

			        ci[1] = j; ci[2] = 3; ci[3] = 1; v[7] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 3; ci[3] = 2; v[8] = leviCevita3D(ci); 

			        ci[1] = j; ci[2] = 3; ci[3] = 3; v[9] = leviCevita3D(ci); 

				    //leviint = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9];
				    leviint = 1;
				    leviDouble = leviint* 1.00;


				    // fix D bar to be non-local
					for (int k=0; k<3; k++){
						for (int l=0; l<3; l++){
				        	DMatrix[i][j][k][l] = phi*CMatrix[i][j][k][l];
				        	DBarMatrix[i][j][k][l] = -leviDouble*computeLNr(Arg)*(CMatrix[i][j][1][1]+CMatrix[i][j][1][2]+\
				        				CMatrix[i][j][2][1] + CMatrix[i][j][2][2]);
				        	//cout << DBarMatrix[i][j][k][l];
						}
					//cout << '\n';
					}
				}
			}

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					for (int k = 0; k < 3; ++k){
						for (int l = 0; l < 3; ++l){
							BBarMatrix[i][j][k][l] = DBarMatrix[k][l][i][j];
							//cout << BBarMatrix[i][j][k][l];
						}
					}
				}
				//cout << "\n";
			}	

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					for (int k = 0; k < 3; ++k){
						for (int l = 0; l < 3; ++l){
							EBarMatrix[i][j][k][l] = 0.5*(DBarMatrix[i][j][k][l]+BBarMatrix[k][l][i][j]);
						}
					}
				}
				//cout << "\n";
			}


			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					for (int k = 0; k < 3; ++k){
						for (int l = 0; l < 3; ++l){
							EBaredTransMatrix [i][j][k][l] = EBarMatrix[k][l][i][j];
						}
					}
				}
				//cout << "\n";
			}

			//cout << strainReal[1][1] << flush;
			realInnerProduct4x2(AMatrix,CuvatureReal); // sets realsig param to inner product

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					dummyVar1[i][j] = realSig[i][j]; // simply make function to return vector to reduce this complexity in future
					//cout << dummyVar1[i][j] ;
				}
			}

			realInnerProduct4x2(EBarMatrix,strainReal); // sets realsig param to inner product
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					dummyVar2[i][j] = realSig[i][j];
					//cout << dummyVar2[i][j] ;
				}
			}

			// set stress perp in step 7
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					stressPerp[i][j] = dummyVar1[i][j] + dummyVar2[i][j];
				}
			}

			realInnerProduct4x2(CMatrix,strainReal); // sets realsig param to inner product

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					dummyVar1[i][j] = realSig[i][j]; // simply make function to return array to reduce this complexity in future
					//cout << dummyVar1[i][j] ;
				}
			}

			realInnerProduct4x2(EBaredTransMatrix,CuvatureReal); // sets realsig param to inner product
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					dummyVar2[i][j] = realSig[i][j];
					//cout << dummyVar2[i][j] ;
				}
			}

			// set stress perp in step 8
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					coupleStressPerp[i][j] = dummyVar1[i][j] + dummyVar2[i][j];
				}
			}



			// sets epsilon in step 9
			realInnerProduct4x2(C_0Matrix,stressPerp);

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					epsilon[i][j] = realSig[i][j];
				}
			}


					// sets chi in step 10
			realInnerProduct4x2(A_0Matrix, coupleStressPerp);

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					chi[i][j] = realSig[i][j];
				}
			}

			// sets the initial stress in step 11
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					stress[i][j] = dummyVar1[i][j]+dummyVar2[i][j]+ stressPerp[i][j];
				}
			}


			realInnerProduct4x2(CMatrix, chi);

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					dummyVar1[i][j] = realSig[i][j];
				}
			}

			realInnerProduct4x2(EBaredTransMatrix,epsilon);

			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					dummyVar2[i][j] = realSig[i][j];
				}
			}
					
			// sets the initial couple stress in step 12
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					coupleStress[i][j] = dummyVar1[i][j]+dummyVar2[i][j]+coupleStressPerp[i][j];
				}
			}


			forces[0][0] = coupleStress[2][1]*Theta[2][2];    // F_1^theta = M_32 theta _33
			forces[0][1] = - coupleStress[2][0]*Theta[2][2];  // F_2^theta = - M_31 theta _33 
			forces[1][0] = stress[0][1]*alpha[0][2];          // F_1^alpha = T_12 alpha _13
			forces[1][1] = - - stress[0][1]*alpha[0][2];      // F_2^alpha = - T_12 alpha _13
			// , 


					// convert position vector to complex vector and call iterate tensors




					// need to convert by using std::complex<long double> dev (leviDouble, 0.00);
					//for (int i = 0; i < 3; ++i)
					//{
					//	for (int j = 0; j < 3; ++j)
					//	{

		            		//std::complex<long double> dev (position[i][j], 0.00);
							//complexPosition[i][j] = 
					//	}
					//}
					//complexPosition[0].real() = Arg[0] 
					//iterateOverTensors(Arg, epsilon, chi, stress, coupleStress);
					     //-- needs to do step 13-25 and complete the algo

		}


		std:: complex <long double> djackaTest(std:: complex <long double> position[3], std::complex<long double> tensor[3][3]){
			// need to return the norm of (dot product of position and tensor)
			// set complex vector 
			std::complex <long double> complexNorm;
			std::complex <long double> stressAtZero = 0.00+0.1i;
			complexDotProductVectorAndTensor(position,tensor); // resets complex
			complexNorm = complexVectorNorm(complexVect)/stressAtZero; 

			// stressAtZero needs to be calucalted by reseting alpha and setting stressAtZero
			// only need the stress one time because it is constant non need keep resetting the stress at every point
			// comput T_n(0) instead of using constant
			// need to return something different
			return complexNorm;
		}

		// need norm of complex 2nd order tensor
		void iterateOverTensors(std:: complex <long double> position[3], long double epsilon[3][3], long double chi[3][3], long double stress[3][3], long double coupleStress[3][3]){
			// this function will be do the iteration process between step 13 and step 25
			std:: complex <long double> e_n;
			std:: complex < long double> epsilonComplex[3][3];
			std:: complex < long double> chiComplex[3][3];

			e_n = djackaTest(position,stressComplex); 
			//while(e_n.real()> convergence){


				// reset e_n after doing the rest
			  //e_n = djackaTest(position,stressComplex);

			//}    
			// make an increment right here


			while(e_n.real() > convergence){ 
				FFT(stress);
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						stressComplex[i][j] = complexSig[i][j];
					}
				}

				FFT(coupleStress); 

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						coupleStressComplex[i][j] = complexSig[i][j];
					}
				}


				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						stressSkewComplex[i][j] =  complexHalf* EulersImaginaryI*(position[0]+position[1]+position[2])*coupleStressComplex[i][j];
						totalStressComplex[i][j] = stressComplex[i][j]+stressSkewComplex[i][j];
					}
				}

				realInnerProduct4x2(C_0Matrix,epsilon);
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						complexSig[i][j] = realSig[i][j] + complexZero ;//convert to compplex sig
						tauStressComplex[i][j] = stressComplex[i][j] - complexSig[i][j];
					}
				}

				// compute polarization tensor step 19
				realInnerProduct4x2(A_0Matrix,chi);
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						complexSig[i][j] = realSig[i][j] + complexZero ;//convert to compplex
						polarizationTensor[i][j] = coupleStressComplex[i][j] - complexSig[i][j]; 
					}
				}

				// step 20-25 
				complexInnerProduct4x2(gamma_0,tauStressComplex);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						complexDummyVar1[i][j] = - complexSig[i][j];
					}
				}

				complexInnerProduct4x2(gamma_1,polarizationTensor);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						complexDummyVar2[i][j] = - complexSig[i][j];
						epsilonComplex[i][j] = complexDummyVar1[i][j]+complexDummyVar2[i][j];
					}
				}

				complexInnerProduct4x2(gamma_2,tauStressComplex);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						complexDummyVar1[i][j] = - complexSig[i][j];
					}
				}

				complexInnerProduct4x2(gamma_3,polarizationTensor);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						complexDummyVar2[i][j] = - complexSig[i][j];
						chiComplex[i][j] = complexDummyVar2[i][j]+complexDummyVar2[i][j];
					}
				}
				product3x3ComplexMatrix(FourierTranformInverse, epsilonComplex); 

				for (int i = 0; i < 2; ++i)
				{
					for (int j = 0; j < 2; ++j)
					{
						epsilon[i][j] = complexSig[i][j].real();
					}
				}

				product3x3ComplexMatrix(FourierTranformInverse, chiComplex); 

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						chi[i][j] = complexSig[i][j].real();
					}
				}


				// compute polarization tensor step 19
				realInnerProduct4x2(CMatrix,epsilon);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						dummyVar1[i][j] = realSig[i][j];
					}
				}
				realInnerProduct4x2(EBarMatrix,chi);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						dummyVar2[i][j] = realSig[i][j];
						stress[i][j] = dummyVar1[i][j]+dummyVar2[i][j]+stressPerp[i][j];
					}
				}

				realInnerProduct4x2(AMatrix,chi);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						dummyVar1[i][j] = realSig[i][j];
					}
				}

				realInnerProduct4x2(EBaredTransMatrix,epsilon);

				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						dummyVar2[i][j] = realSig[i][j];
						coupleStress[i][j] = dummyVar2[i][j]+dummyVar2[i][j]+coupleStressPerp[i][j];
					}
				}
				e_n = djackaTest(position,stressComplex); 
			}

		}

			// levi cevita symbol
		int leviCevita3D(int vect[3]){
			if (vect[1] == 1 && vect[2] == 2 && vect[3] == 3){
				return  1;
			}
			if (vect[1] == 3 && vect[2] == 1 && vect[3] == 2){
				return  1;
			}
			if (vect[1] == 2 && vect[2] == 3 && vect[3] == 1){
				return  1;
			}
			if (vect[1] == 2 && vect[2] == 1 && vect[3] == 3){
				return -1;
			}
			if (vect[1] == 3 && vect[2] == 2 && vect[3] == 1){
				return -1;
			}
			if (vect[1] == 1 && vect[2] == 3 && vect[3] == 2){
				return -1;
			}else{
				return  0;
			}
		}
		
		// cronickerDelta function
		int leviCroncker3D(int vect[2]){
			if (vect[1] == vect[2]){
				return 1;
			}else{
				return 0;
			}
		}
		
		// determinant for complex 3x3 matrices
		std::complex<long double> det3x3complexMatrix(std::complex<long double> matrix[3][3]){
			// det complex sigconst 
			std::complex<long double> detiminant;
			detiminant = matrix[0][0]*(matrix[1][1]* matrix[2][2] -  matrix[1][2]* matrix[2][1]) -\
						 matrix[0][1]*(matrix[1][0]* matrix[2][2] -  matrix[1][2]* matrix[2][0]) +\
						 matrix[0][2]*(matrix[1][0]* matrix[2][1] -  matrix[1][1]* matrix[2][0]);
			return detiminant;
		}

		// 3d cronicker delta function
		void Croncker3D(int vect[3][3]){
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					if(i==j){
						vect[i][j] = 0;
					}else{

					}	
				}
			}
		}

		long double tensorNorm(long double vect[3][3]){
			return sqrt(pow(vect[1][1],2) + pow(vect[1][2],2) + pow(vect[1][3],2)\
				+  pow(vect[2][1],2)+pow(vect[2][2],2)+pow(vect[2][3],2)\
				+ pow(vect[3][1],2)+pow(vect[3][2],2)+pow(vect[3][3],2));
		}

		//std::sqrt(std::complex<double>(real, imaginary));
		std::complex<long double> complexVectorNorm(std::complex<long double> vect[3]){
			return std::sqrt(std::pow(vect[0],2) + std::pow(vect[1],2) + std::pow(vect[2],2));
		}

		std::complex<long double> traceComplex(std::complex<long double> vect[3][3]){
			return vect[1][1] + vect[2][2] + vect[3][3];
		}
		long double realTrace(long double vect[3][3]){
			return vect[1][1] + vect[2][2] + vect[3][3];
		}


	   void realInnerProduct4x2(long double Tensor1[3][3][3][3], long double Tensor2[3][3]){
	   	    // compute the inner product of 4th order tensor and 2nd order tensor and set to realSig;
	   		static long double v[9];
	   		for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					v[1] = Tensor1[i][j][1][1]*Tensor2[1][1]; //cout << v[1];
					v[2] = Tensor1[i][j][1][2]*Tensor2[1][2]; //cout << v[2];
					v[3] = Tensor1[i][j][1][3]*Tensor2[1][3]; //cout << v[3];
					v[4] = Tensor1[i][j][2][1]*Tensor2[2][1]; //cout << v[4];
					v[5] = Tensor1[i][j][2][2]*Tensor2[2][2]; //cout << v[5];
					v[6] = Tensor1[i][j][2][3]*Tensor2[2][3]; //cout << v[6];
					v[7] = Tensor1[i][j][3][1]*Tensor2[3][1]; //cout << v[7];
					v[8] = Tensor1[i][j][3][2]*Tensor2[3][2]; //cout << v[8];
					v[9] = Tensor1[i][j][3][3]*Tensor2[3][3]; //cout << v[9];
				    realSig[i][j] = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9];
				}	
			}

		}


	   void complexInnerProduct4x2(std::complex<long double> Tensor1[3][3][3][3], std::complex<long double> Tensor2[3][3]){
	   	    // compute the inner product of 4th order tensor and 2nd order tensor and set to realSig;
	   		std::complex<long double> v[9];
	   		for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					v[1] = Tensor1[i][j][1][1]*Tensor2[1][1]; 
					v[2] = Tensor1[i][j][1][2]*Tensor2[1][2]; 
					v[3] = Tensor1[i][j][1][3]*Tensor2[1][3];
					v[4] = Tensor1[i][j][2][1]*Tensor2[2][1]; 
					v[5] = Tensor1[i][j][2][2]*Tensor2[2][2]; 
					v[6] = Tensor1[i][j][2][3]*Tensor2[2][3];
					v[7] = Tensor1[i][j][3][1]*Tensor2[3][1]; 
					v[8] = Tensor1[i][j][3][2]*Tensor2[3][2]; 
					v[9] = Tensor1[i][j][3][3]*Tensor2[3][3];
				    complexSig[i][j] = v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]+v[9];
				}	
			}

		}

		/*
		void PrintMyArray(long double Arg[3][3]){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; i++){
					cout << Arg[i][j];
				}
				cout << '\n' << flush;
			}
		}*/

};

#endif
