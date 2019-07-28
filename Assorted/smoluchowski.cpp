#include <iostream> // allows program print to screen
using std::cout;
using std::cin;
using std::endl;
using std::showpoint;
using std::fixed;
using std::left;
using std::ios;

#include <string>
using std :: string;
//using std :: to_string; 

#include <sstream>
using std::stringstream;

#include <iomanip>
using std::setw;
using std::setprecision;

#include<fstream> //file stream
using std::ofstream; //output file stream
using std::ifstream; // input file stream

#include<sstream>

#include<cmath>
using std::pow;

#include<cstdlib>
using std::rand;
using std::srand;

#include<ctime>
using std::time;

//Input parameters
static int tsteps   = 2;
static int max_agg  = 3; 
static int npart   = 1;
static double dt    = 0.001;
static double k     = 0.1; 

//Method declarations
void solveSmoluchowski(); 
void writeSolution(int); 


//Dependent variables
static int* num_agg; 
static int* num_agg_new; 
static double* form; 
static double* decay; 


/*****************************************************************************************************************************************
 *                  MAIN METHOD
 ****************************************************************************************************************************************/

int main()
{
    int t;
    int j; 

    cout <<"Hello Im in the main program"<< endl; 
 
    num_agg        = new int[max_agg]; 
    num_agg_new    = new int[max_agg]; 
    form            = new double[max_agg]; 
    decay           = new double[max_agg]; 

    num_agg[1] = npart; 
    num_agg_new[1] = 0; 
    form[1] = 0.0; 
    decay[1] = 0.0; 
        
    for (j=2; j<max_agg; j++){
        num_agg[j]      = 0; 
        num_agg_new[j]  = 0; 
        form[j]         = 0.0; 
        decay[j]        = 0.0; 
               
   }

    for (t=1; t<=tsteps; t++) {
    
    solveSmoluchowski(); 
    writeSolution(t);
     
    }
   delete[] num_agg; 
   delete[] num_agg_new; 
   delete[] form;
   delete[] decay; 

    return 0; 
}
/********************************************************************************************************************************************
 *              END MAIN METHOD 
 * ******************************************************************************************************************************************/


/********************************************************************************************************************************************
 *              FUNCTIONS
 * *****************************************************************************************************************************************/

void solveSmoluchowski(){
    int i;
    int j; 

    cout <<"Hello Im in the Smol" << endl;

   for (j=1; j<max_agg; j++){

        //formation
        form[j] = 0; 
 
        for (i=1; i<j; i++) {
            form[j] = form[j] + num_agg[i]*num_agg[j-i]; 
            cout << j << "  "<< i << "   "<< num_agg[j-i] << endl; 
            }
        //

        //decay
        decay[j] = 0;
 
        for (i=1; i<max_agg; i++){
           decay[j] = decay[j] -num_agg[j]*num_agg[i]; 
//
        }
  //      //

        num_agg_new[j] = num_agg[j] + 0.5*k*form[j]*dt + k*decay[j]*dt;
    
    }
    
    for (j=1; j<max_agg; j++){
        num_agg[j] = num_agg_new[j]; 
    }

}

void writeSolution(int tstep){
    stringstream ss;
    ss << "traj" << tstep <<".dat";
    
    ofstream outputSolution(ss.str().c_str(), ios::out);
    int j;
    for(j=1;j<max_agg;j++){
        outputSolution << j << "     "<< num_agg[j]<< endl;
    } 
}

















