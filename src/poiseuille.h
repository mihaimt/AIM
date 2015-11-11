#ifndef __POISEUILLE_H__
#define __POISEUILLE_H__

#include <iostream>
#include <string>
#include <vector>
#include "fparser.h"
#include "ext_constants.h"

//------------------------------------------------------------------------------
// Class to store Poiseuille correction function
//------------------------------------------------------------------------------
class Poiseuille
{
   public:
      void    add (std::string, std::string);
      void value (const double& p, const double& t, const double& vel, const double& mu ,std::vector<double>& cor_val);

   private:
      FParser momentum_cor_coef;
      FParser energy_cor_coef;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void Poiseuille::add (std::string variable, std::string fun)
{
   if(variable == "momentum_cor_coef")
      momentum_cor_coef.FParse (fun);
   else if(variable == "energy_cor_coef")
      energy_cor_coef.FParse (fun);
   else
   {
      std::cout << "Sources::add: Unknown variable " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate source at a given point
//------------------------------------------------------------------------------
inline
void Poiseuille::value (const double& p, const double& t, const double& vel, const double& mu, std::vector<double>& cor_val)
{
   cor_val.resize(2);
   double vals[] = {p, t};
   cor_val[0]     = momentum_cor_coef.Eval (vals)*vel*mu;
   cor_val[1]     = energy_cor_coef.Eval (vals)*vel*vel*mu;

}

#endif
