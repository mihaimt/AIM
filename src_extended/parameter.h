#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "reader.h"
#include "ext_constants.h"
#include "ic.h"
#include "sources.h"
#include "poiseuille.h"

class Parameter
{
   public:
      char* input_file;
      std::string mesh_file;
      std::string ic_file;
      std::string source_file;
      
      bool mesh_file_read, ic_file_read, source_file_read;

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
      Sources sources;
      Poiseuille poiseuille_cor;
      std::string OUTPUT_PATH;
      double probe_pos;
      double cut_off_pre;
      double atm_pressure;
      Valve valve;

      double output_dt;
      void read ();
      double viscosity (const double T) const;

   private:
      void read_grid (Reader&);
      void read_source (Reader&);
      void read_poiseuille (Reader&);
      void read_constants (Reader&);
      void read_numeric (Reader&);
      void read_material (Reader&);
      void read_initial_condition (Reader&);
      void read_output (Reader&);
};

#endif
