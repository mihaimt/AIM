// Navier-Stokes solver for simulationof a one-dimensional arc within a gas insulated
// busbar
// Author: Deep Ray
// Date  : 9 November , 2015

#include <iostream>
#include "fv.h"


using namespace std;

map<string,double> constants;

void process_command_line (int argc, char* argv[], int& ifile);


int main (int argc, char* argv[])
{
   cout << "=================================================================================\n";
   cout << "   Starting 1D NSE solver for 1D arc within a gas insulated busbar\n";   
   cout << "   ---  Author: Deep Ray \n";
   cout << "=================================================================================\n";
   
   int ifile;
   process_command_line (argc, argv, ifile);
   
   FV fv(argv[ifile]);
   fv.run();
   
   return 0;
}
