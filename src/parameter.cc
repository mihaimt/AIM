#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include "parameter.h"

using namespace std;

extern map<string,double> constants;

//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void Parameter::read ()
{
   cout << "Reading input file " << input_file << endl;
   Reader fin(input_file);

   read_grid (fin);
   read_numeric (fin);
   read_material (fin);
   read_constants (fin);
   read_initial_condition (fin);
   read_source (fin);
   read_poiseuille (fin);
   read_output (fin);
}

//------------------------------------------------------------------------------
// Read grid section
//------------------------------------------------------------------------------
void Parameter::read_grid (Reader &fin)
{
   cout << "  Reading grid section\n";

   string input;

   fin.begin_section ("grid");

   fin.entry ("read_from_file");
   fin >> input;
   
   if(input == "yes")
   {
      fin.entry ("file_name");
      fin >> mesh_file;
      mesh_file_read = true;
   }
   else if(input == "no")
   {
      mesh_file_read = false;
      fin.entry ("x_min");
      fin >> xmin;
   
      fin.entry ("x_max");
      fin >> xmax;
   
      assert(xmin < xmax);
   
      fin.entry ("n_cell");
      fin >> n_cell;
   
      assert(n_cell > 0);
   }
   else
   {
      cout << "   Error: unknown option " << input << " for grid: read_from_file"<<endl;
      cout << "   Allowed values: yes, no"<<endl;
      exit (0);
   }
   fin.entry("probe");
   fin >> probe_pos;
   
   assert(probe_pos >= xmin);
   assert(probe_pos <= xmax);
   
   fin.end_section ();
   
}

//------------------------------------------------------------------------------
// Read numeric section
//------------------------------------------------------------------------------
void Parameter::read_numeric (Reader &fin)
{
   cout << "  Reading numeric section\n";

   string input;

   fin.begin_section ("numeric");

   fin.entry ("time_scheme");
   fin >> input;
   if(input=="rk1")
   {
      n_rks = 1;
      time_scheme = rk1;
   }   
   else if(input=="ssprk3")
   {
      n_rks = 3;
      time_scheme = ssprk3;
   } 
   else if(input=="jameson_rk4")
   {
      n_rks = 4;
      time_scheme = jameson_rk4;
   } 
   else
   {
      cout << "   Error: unknown time_scheme " << input << endl;
      exit (0);  
   }

   fin.entry ("cfl");
   fin >> cfl;
   assert (cfl > 0.0);

   fin.entry ("max_iter");
   fin >> max_iter;
   assert (max_iter > 0);

   fin.entry ("final_time");
   fin >> final_time;
   assert (final_time >= 0.0);

   fin.entry ("min_residue");
   fin >> min_residue;
   assert (min_residue > 0.0);

   fin.entry("reconstruct");
   fin >> input;
   if(input == "first")
      reconstruct_scheme = first;
   else if(input == "second")   
      reconstruct_scheme = muscl;
   else
   {
      cout << "   Error: unknown reconstruction method " << input << endl;
      exit (0);   
   }   
          
   fin.entry("bc_type");
   fin >> input;
   if(input=="wall")
   {
      bc = wall;
   }   
   else if(input=="periodic")
   {
      bc = periodic;
   } 
   else if(input=="neumann")
   {
      bc = neumann;
   } 
   else
   {
      cout << "   Error: unknown boundary type " << input << endl;
      exit (0);  
   }      
   
   fin.entry("cut_off_pre");
   fin >> cut_off_pre;

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read material section
//------------------------------------------------------------------------------
void Parameter::read_material (Reader &fin)
{
   cout << "  Reading material section\n";

   string input;

   fin.begin_section ("material");

   fin.entry ("gamma");
   fin >> gamma;
   assert (gamma > 1.0);

   fin.entry ("gas_const");
   fin >> gas_const;
   assert (gas_const > 0.0);

   fin.begin_section ("viscosity");
   fin.entry ("model");
   fin >> input;
   if(input == "constant")
   {
      mu_model = mu_constant;

      fin.entry ("mu_ref");
      fin >> mu_ref;
      assert (mu_ref >= 0.0);
   }
   else if(input == "sutherland")
   {
      mu_model = mu_sutherland;
      fin.entry ("mu_ref");
      fin >> mu_ref;
      fin.entry ("T_ref");
      fin >> T_ref;
      fin.entry ("T_0");
      fin >> T_0;
      assert (mu_ref > 0.0);
      assert (T_ref > 0.0);
      assert (T_0 > 0.0);
   }
   else if(input == "power")
   {
      mu_model = mu_power;
      fin.entry ("mu_ref");
      fin >> mu_ref;
      fin.entry ("T_ref");
      fin >> T_ref;
      fin.entry ("omega");
      fin >> omega;
      assert (mu_ref > 0.0);
      assert (T_ref > 0.0);
      assert (omega > 0.0);
   }
   else if(input == "mu_arc")
   {
      mu_model = mu_arc;
      // fin.entry ("mu_ref");
//       fin >> mu_ref;
//       fin.entry ("T_ref");
//       fin >> T_ref;
//       fin.entry ("omega");
//       fin >> omega;
//       assert (mu_ref > 0.0);
//       assert (T_ref > 0.0);
//       assert (omega > 0.0);
   }
   else
   {
      cout << "   Error: unknown viscosity type " << input << endl;
      exit (0);
   }
   fin.end_section ();

   fin.entry ("prandtl");
   fin >> Pr;
   assert (Pr > 0.0);

   fin.entry ("model");
   fin >> input;
   if(input == "euler")
      model = euler;
   else if(input == "ns")
      model = ns;
   else
   {
      cout << "   Error: unknown flow model " << input << endl;
      exit (0);
   }

   fin.entry ("flux");
   fin >> input;
   if(input == "roe")
      flux_scheme = roe;  
   else if(input == "roe_fixed")
      flux_scheme = roe_fixed;   
   else if(input == "rusanov")
      flux_scheme = rusanov;
   else   
   {
      cout << "   Error:: unknown flux scheme: " << input << endl;
      exit (0);
   }
   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read constants section
//------------------------------------------------------------------------------
void Parameter::read_constants (Reader &fin)
{
   cout << "  Reading constants section\n";

   string input;
   double value;

   fin.begin_section ("constants");

   while (!fin.eos())
   {
      fin >> input;
      fin >> value;
      cout << setw(16) << input << setw(16) << value << endl;
      constants.insert ( pair<string,double>(input, value) );
   }

}

//------------------------------------------------------------------------------
// Read initial condition
//------------------------------------------------------------------------------
void Parameter::read_initial_condition (Reader &fin)
{
   cout << "  Reading initial condition section\n";

   string input;

   fin.begin_section ("initial_condition");

   fin.entry ("read_from_file");
   fin >> input;
   
   if(input == "yes")
   {
      fin.entry("file_name");
      ic_file_read = true;
      fin >> ic_file;
   }
   else if(input == "no")
   {
      ic_file_read = false;
      
      fin.entry("density");
      fin.getline (input);
      initial_condition.add ("density", input);

      fin.entry ("velocity");
      fin.getline (input);
      initial_condition.add ("velocity", input);

      fin.entry ("pressure");
      fin.getline (input);
      initial_condition.add ("pressure", input);
   }
   else
   {
      cout << "   Error: unknown option " << input << " for initial_condition: read_from_file"<<endl;
      cout << "   Allowed values: yes, no"<<endl;
      exit (0);
   }
   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read sources
//------------------------------------------------------------------------------
void Parameter::read_source (Reader &fin)
{
   cout << "  Reading source section\n";

   string input;

   fin.begin_section ("source");
      
   fin.entry("mass_source");
   fin.getline (input);
   sources.add ("mass_source", input);

   fin.entry("energy_source");
   fin.getline (input);
   sources.add ("energy_source", input);
   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read sources
//------------------------------------------------------------------------------
void Parameter::read_poiseuille (Reader &fin)
{
   cout << "  Reading Poiseuille correction section\n";

   string input;

   fin.begin_section ("poiseuille");
      
   fin.entry("momentum_cor_coef");
   fin.getline (input);
   poiseuille_cor.add ("momentum_cor_coef", input);

   fin.entry("energy_cor_coef");
   fin.getline (input);
   poiseuille_cor.add ("energy_cor_coef", input);
   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read output section
//------------------------------------------------------------------------------
void Parameter::read_output (Reader &fin)
{
   cout << "  Reading output section\n";

   string input;

   fin.begin_section ("output");
   fin.entry ("output_dt");
   fin >> output_dt;
   assert (output_dt > 0);
   
   fin.entry ("output_path");
   fin >> OUTPUT_PATH;

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Viscosity coefficient
//------------------------------------------------------------------------------
double Parameter::viscosity (const double T) const
{
   switch (mu_model)
   {
      case mu_constant:
         return mu_ref;

      case mu_sutherland:
         return mu_ref * pow(T/T_ref, 1.5) * (T_ref + T_0) / (T + T_0);

      case mu_power:
         return mu_ref * pow(T/T_ref, omega);
         
      case mu_arc:
         return sqrt(gamma*T)*1.5e-6;   

      default:
         cout << "    viscosity: unknown model " << mu_model << endl;
         abort ();
   }
}
