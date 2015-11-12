#ifndef __MISC_FUNC_H__
#define __MISC_FUNC_H__

#include <cmath>

#define SIGN(a) (((a)<0) ? -1:1)

using namespace std;

//------------------------------------------------------------------------------
// Minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a,
               const double& b,
               const double& c)
{
   double result;
   
   if (a*b > 0.0 && b*c > 0.0)
   {
      result  = min( min(fabs(a), fabs(b)), fabs(c) );
      result *= SIGN(a);
   }
   else
      result = 0.0;
   
   return result;
   
}

//------------------------------------------------------------------------------
// vanleer limiter
//------------------------------------------------------------------------------
double vanleer (const double& a,
                const double& b)
{
   double du;

   if(fabs(a*b) > 0.0)
      du = (SIGN(a)+SIGN(b))*fabs(a)*fabs(b)/(fabs(a) + fabs(b));
   else
      du = 0.0;

   return du;
}
//------------------------------------------------------------------------------
// Reconstruct left state of right face
//------------------------------------------------------------------------------
vector<double> muscl_limited (const vector<double>& ul,
                      const vector<double>& uc,
                      const vector<double>& ur)
{
   unsigned int n = ul.size();
   vector<double> result (n);
   double dul, duc, dur;
   const double beta = 2.0;
   
   for(unsigned int i=0; i<n; ++i)
   {
      dul = uc[i] - ul[i];
      dur = ur[i] - uc[i];
      duc = (ur[i] - ul[i])/2.0;
      //result[i] = uc[i] + 0.5 * vanleer (dul, dur);
      result[i] = uc[i] + 0.5 * minmod (beta*dul, duc, beta*dur);
   }
   
   return result;
}

#endif