/*
    -- Main reference implementation "Nonlocal elasticity tensors in dislocation and disclination cores"
	  
    Authors to main reference: 
            Capolungo, Laurent
		        Gbemou, kodjovi
		        Taupin, Vincent
		        Fressengeas, Claude
    Info:        
		      Provided by the author(s) and the Los Alamos National Laboratory (2017-04-14).
		      To be published in: Journal of the Mechanics and Physics of Solids
		      DOI to publisher's version: 10.1016/j.jmps.2017.01.003
		      Permalink to record: http://permalink.lanl.gov/object/view?what=info:lanl-repo/lareport/LA-UR-17-20035

	  --> terminal compilation with g++ (--> = enter key):
		--> cd <current directory> --> g++ dislocationMain.cpp tensorMaps.cpp --> ./a.exe
    -- interfaces have extension (.h) and Implementations hav extension (.cpp).
    -- see data for data[i][j] for parameter decloration and instalization
    -- alpha must not be zero as it denotes the curl of the plastic destortion tensor U (See Eq.(15))
    -- T_i must be a unit vector as well as e_1,e_2, and e_3 --> add check method before runnning algo
*/

/*

///
*/
#include <fstream>
#include <string>
#include <vector>
#include "tensorMaps.h"
#include <complex>
#include <math.h>
#include <iostream>
template< class T >
class complex;
using namespace std;
const double small_inc = 0.0061; // controlls mesh size
int maxnt;
long double height;
long double xx; 
long double yy; 
long double zz; 
long double angle; 
long double angle_0;
long double angleINC;
long double cubeLength;


int main(){
  // always use double qutoes with c++ in string otherwise it is a character
  cout << "Welcome to Roys Dislocation Discritization Algorithm" << endl;
  cout << "                                                    " << endl;
  cout << "======== ====       ====   =======   ==============  ============      ==========   ==============              " << endl;
  cout << "  ||    ||/\\     //||     //   /\\       ||        ||              //               ||     ||                " << endl;
  cout << "  ||    || /\\   // ||     ||    ||       ||        ||             ||                ||     ||                "  << endl;
  cout << "  ||    ||  /\\ //  ||      /\\           ||        ||==========   ||                ||=====||                " << endl;
  cout << "  ||    ||   =====  ||        /\\         ||        ||             /\\               ||     ||                " << endl;
  cout << "  ||    ||          ||    ||   //         ||        ||              /\\              ||     ||                " << endl;
  cout << "======  ====      =====   ======          ||        =============    ===========   =============              " << endl;
  cout << "                                                                                                               " << endl;

     


  Tensors crystals; // calling constructor to gain access to tensor library
  crystals.declarations(); // declares all parameter in code.
  ofstream fileStrain,fileCoupleStress,fileStress, filePositionVector,fileThetaTensor,\
             fileAlphaTensor,fileDTensor,fileBTensor,fileDensityTensor, fileForces;
  fileStrain.open ("strain_data");
  fileStress.open ("stress_data");
  fileCoupleStress.open ("coupleStress_data");
  filePositionVector.open ("positionVector_data");
  fileThetaTensor.open ("theta_data");
  fileAlphaTensor.open ("alpha_data");
  fileDTensor.open("D_data");
  fileBTensor.open("B_data");
  fileDensityTensor.open("varphi_data");
  fileForces.open("forces_data");

  // + small_inc*tt; 
  //cubeLength = 1;
  //xx = -cubeLength; yy = -cubeLength; zz = 0.0;
  //maxnt = floor(2*cubeLength/crystals.gridSize);

//cubeLength = 1;
xx = 0.5; yy = 0.5; zz = 0.0;
maxnt = floor(2*abs(xx)/crystals.gridSize);


    // MAIN lOOP compute tensors at each point in the grid

    for(int j=1; j < maxnt; j++){
      yy = 0.5;
     for(int i=1; i < maxnt;i++){
        yy = yy + crystals.gridSize;  
        // start curvature meshed at ∀ξ 
        // ξ and ξ' should never be equal (for non local contribution)

        crystals.Xi[0] =  xx; 
        crystals.Xi[1] =  yy; 
        crystals.Xi[2] =  0.0;

      // sets the Nye tensor α and rotation tensor θ at each point ξ
      crystals.setDislocations(crystals.Xi); 

      fileAlphaTensor << crystals.alpha[0][0] << " "  << crystals.alpha[0][1]     << " "  << crystals.alpha[0][2] << " "  <<\
                crystals.alpha[1][0] << " "    << crystals.alpha[1][1]   << " "    << crystals.alpha[1][2] << " "  <<\
                crystals.alpha[2][0] << " "    << crystals.alpha[2][1]   << " "    << crystals.alpha[2][2] << " "  << endl;

      

      fileThetaTensor << crystals.Theta[0][0] << " "  << crystals.Theta[0][1]     << " "  << crystals.Theta[0][2] << " "  <<\
                crystals.Theta[1][0] << " "    << crystals.Theta[1][1]   << " "    << crystals.Theta[1][2] << " "  <<\
                crystals.Theta[2][0] << " "    << crystals.Theta[2][1]   << " "    << crystals.Theta[2][2] << " "  << endl;



      crystals.FFT(crystals.alpha); // sets complex sig to fourier mapping

      // step 1
      for(int i=0;i<3;i++){                   
        for(int j=0;j<3;j++){
          crystals.Alphatilda[i][j] = crystals.complexSig[i][j];
                 //cout << crystals.Alphatilda[i][j]; 
        }
               //cout << '\n'; 
      }
      crystals.FFT(crystals.Theta);

      // step 2
      for(int i=0;i<3;i++){                   
        for(int j=0;j<3;j++){
          crystals.Thetatilda[i][j] = crystals.complexSig[i][j];
                //cout << crystals.Thetatilda[i][j]; 
        }
              //cout << '\n';
      }
    
      
      crystals.setCurvatureAndStrain(crystals.Xi);  //set curvature at each point in V (step 3-6)
      crystals.computeTensors(crystals.Xi); // completes steps 7-25

      filePositionVector <<  crystals.Xi[0] << " " <<  crystals.Xi[1] << " " << crystals.Xi[2] << " " << endl;

      fileStrain << crystals.strainReal[0][0] << " "  << crystals.strainReal[0][1]   << " "  << crystals.strainReal[0][2] << " "  <<\
                crystals.strainReal[1][0] << " "  << crystals.strainReal[1][1]   << " "    << crystals.strainReal[1][2] << " "  <<\
                crystals.strainReal[2][0] << " "  << crystals.strainReal[2][1]   << " "  << crystals.strainReal[2][2] << " "  << endl;


      fileDTensor << crystals.DMatrix[0][1][2][0] << " "  << endl;

      fileBTensor << crystals.BMatrix[2][0][0][0] << " "   << endl;

      fileForces << crystals.forces[0][0] << " " << crystals.forces[0][1] << " " << crystals.forces[1][0] << " " << crystals.forces[1][1] << " " << endl;

      fileCoupleStress << crystals.coupleStress[0][0] << " "  << crystals.coupleStress[0][1]     << " "  << crystals.coupleStress[0][2] << " "  <<\
                crystals.coupleStress[1][0] << " "    << crystals.coupleStress[1][1]   << " "    << crystals.coupleStress[1][2] << " "  <<\
                crystals.coupleStress[2][0] << " "    << crystals.coupleStress[2][1]   << " "    << crystals.coupleStress[2][2] << " "  << endl;

      fileStress << crystals.stress[0][0] << " "  << crystals.stress[0][1]     << " "  << crystals.stress[0][2] << " "  <<\
                crystals.stress[1][0] << " "    << crystals.stress[1][1]   << " "    << crystals.stress[1][2] << " "  <<\
                crystals.stress[2][0] << " "    << crystals.stress[2][1]   << " "    << crystals.stress[2][2] << " "  << endl;
    }
    xx = xx+ crystals.gridSize;
  }  
    
    // close files
    fileStrain.close();     
    fileCoupleStress.close();
    fileThetaTensor.close();
    filePositionVector.close();
    fileAlphaTensor.close();
    fileStress.close();
    fileDTensor.close();
    fileBTensor.close();
    fileDensityTensor.close();
    fileForces.close();
   return 0;
}

      
      //angle_0 = 0.0; xx = 1.0*i; angleINC = (2*M_PI/maxnt)*tt; angle = angle_0 + angleINC;
    /*
      // define rotation vectors in e_1, e_2, e_3 
      crystals.MatrixX[0][0] = crystals.e_1; crystals.MatrixZ[0][0] = cos(angle);     crystals.MatrixZ[0][0] = cos(angle);    
      crystals.MatrixX[0][1] = 0.0;          crystals.MatrixZ[0][1] = 0.0;            crystals.MatrixZ[0][1] = 0.0;  
      crystals.MatrixX[0][2] = 0.0;          crystals.MatrixZ[0][2] = sin(angle);     crystals.MatrixZ[0][2] = -sin(angle);      

      crystals.MatrixX[1][0] = 0.0;          crystals.MatrixZ[1][0] = 0.0;            crystals.MatrixZ[1][0] = cos(angle);  
      crystals.MatrixX[1][1] = cos(angle);   crystals.MatrixZ[1][1] = crystals.e_2;   crystals.MatrixZ[1][1] = -sin(angle);  
      crystals.MatrixX[1][2] = -sin(angle);   crystals.MatrixZ[1][2] = 0.0;            crystals.MatrixZ[1][2] = 0.0;  

      crystals.MatrixX[2][0] = 0.0;          crystals.MatrixZ[2][0] = -sin(angle);    crystals.MatrixZ[2][0] = 0.0;   
      crystals.MatrixX[2][1] = sin(angle);   crystals.MatrixZ[2][1] = 0.0;            crystals.MatrixZ[2][1] = 0.0;   
      crystals.MatrixX[2][2] = cos(angle);   crystals.MatrixZ[2][2] = cos(angle);     crystals.MatrixZ[2][2] = crystals.e_2;  
      */

      //crystals.product3x3RealMatrix(crystals.MatrixZ, crystals.MatrixZ);
      //crystals.product3x3RealMatrix(crystals.MatrixX, crystals.realSig);
      //crystals.product3x3RealMatrix(crystals.MatrixX, crystals.realSig);
