#ifndef __IC_H__
#define __IC_H__

#include <iostream>
#include <string>
#include <vector>
#include "fparser.h"
#include "ext_constants.h"

//------------------------------------------------------------------------------
// Class to store initial condition functions
//------------------------------------------------------------------------------
class InitialCondition
{
   public:
      void    add (std::string, std::string);
      void value (const double& p, std::vector<double> & prim);

   private:
      FParser density;
      FParser velocity;
      FParser pressure;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string variable, std::string fun)
{
   if(variable == "density")
      density.FParse (fun);
   else if(variable == "velocity")
      velocity.FParse (fun);
   else if(variable == "pressure")
      pressure.FParse (fun);
   else
   {
      std::cout << "InitialCondition::add: Unknown variable " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
void InitialCondition::value (const double& p, std::vector<double>& prim)
{
   prim.resize(NVAR);

   double t = 0.0;
   double vals[] = {p, t};
   prim[0]     = density.Eval (vals);
   prim[1]     = velocity.Eval (vals);
   prim[2]     = pressure.Eval (vals);

}

#endif
