#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "reader.h"
#include "ext_constants.h"
#include "ic.h"

class Parameter
{
   public:
      char* input_file;
      std::string mesh_file;
      std::string ic_file;
      
      bool mesh_file_read, ic_file_read;

      unsigned int n_rks;
      unsigned int max_iter;
      double cfl;
      double final_time;
      double min_residue;
      double xmin;
      double xmax;
      unsigned int n_cell;
      
      Boundary bc;
      
      double gamma;
      double gas_const;
      double Pr;
      double Cp;
      double T_0, T_ref, mu_ref; // constants for sutherland law
      double omega; // exponent in power law viscosity
      double max_pressure;
      
      FlowModel model;
      FluxScheme flux_scheme;
      MuModel mu_model;
      ReconstructionScheme reconstruct_scheme;
      TimeIntegrationScheme time_scheme;
      InitialCondition initial_condition;

      unsigned int write_frequency;
      void read ();
      double viscosity (const double T) const;

   private:
      void read_grid (Reader&);
      void read_constants (Reader&);
      void read_numeric (Reader&);
      void read_material (Reader&);
      void read_initial_condition (Reader&);
      void read_output (Reader&);
};

#endif
