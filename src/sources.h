#ifndef __SOURCES_H__
#define __SOURCES_H__

#include <iostream>
#include <string>
#include <vector>
#include "fparser.h"
#include "ext_constants.h"

//------------------------------------------------------------------------------
// Class to store source functions
//------------------------------------------------------------------------------
class Sources
{
   public:
      void    add (std::string, std::string);
      void value (const double& p, const double& t, std::vector<double>& source_val);

   private:
      FParser mass_source;
      FParser energy_source;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void Sources::add (std::string variable, std::string fun)
{
   if(variable == "mass_source")
      mass_source.FParse (fun);
   else if(variable == "energy_source")
      energy_source.FParse (fun);
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
void Sources::value (const double& p, const double& t, std::vector<double>& source_val)
{
   double vals[] = {p, t};
   source_val[0]     = mass_source.Eval (vals);
   source_val[1]     = energy_source.Eval (vals);

}

#endif
