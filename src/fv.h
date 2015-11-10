#ifndef __FV_H__
#define __FV_H__

//------------------------------------------------------------------------------
// Main problem class
//------------------------------------------------------------------------------

#include <vector>
#include <fstream>
#include "ext_constants.h"
#include "parameter.h"

class FV
{
   public:
      FV(char* file)
      {
         param.input_file = file;
         param.read();
         ng_cell = 2;
      };
      void run ();
   
   private:
      Parameter param;
      
      void initialize ();
      void compute_dt ();
      void con_to_prim ();
      void prim_to_con ();
      void compute_face_derivatives ();
      void reconstruct (const unsigned int face,
                        std::vector<double>& left,
                        std::vector<double>& right) const;
      
      void roe_flux(const std::vector<double>& left,
                         const std::vector<double>& right,
                         std::vector<double>& flux) const;
      void roe_fixed_flux(const std::vector<double>& left,
                         const std::vector<double>& right,
                         std::vector<double>& flux) const;                   
      void rusanov_flux ( const std::vector<double>& left,
                         const std::vector<double>& right,
                         std::vector<double>& flux) const;                             
      void num_flux (const int&  f,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const double&         tau,
                     const double&         q,
                     std::vector<double>&) const;              
      void compute_face_values ();
      void compute_residual ();
      void compute_residual_norm ();
      void update_solution (const unsigned int rk);                                        
      void output ();
      void apply_bc();
      double temperature(const std::vector<double>& prim) const;
      double enthalpy(const std::vector<double>& prim) const;
      //double viscosity (const double T) const;
    
      double dt;
      std::vector<double> xc;
      std::vector<double> xf;
      std::vector<double> dx;
   
      std::vector< std::vector<double> > primitive;
      std::vector< std::vector<double> > primitive_old;
      std::vector< std::vector<double> > residual;
      std::vector< std::vector<double> > conserved;
      std::vector< std::vector<double> > conserved_old;
      std::vector< std::vector<double> > conl, conr;
      std::vector<double> res_norm;
      std::vector<double> tau_f,  q_f;
      double tau,q;
      
      unsigned int ng_cell, n_cell, n_face;
      unsigned int counter;
      
      double Cp;
      
};

#endif