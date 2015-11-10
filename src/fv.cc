#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include "fv.h"
#include "misc_func.h"


using namespace std;

//------------------------------------------------------------------------------
// Initializing problem
//------------------------------------------------------------------------------
void FV::initialize ()
{
   // Make grid
   cout<<"\n  Initializing mesh ...\n";
   if(param.mesh_file_read)
   {
      ifstream fi;
      fi.open (param.mesh_file.c_str());
      assert (fi.is_open());
      fi >> n_face;
      n_cell = param.n_cell;
      n_cell = n_face - 1;
      xf.resize(n_face);
      xc.resize(n_cell);
      dx.resize(n_cell);
      for(unsigned int i=0; i<n_face; ++i)
      {
         fi >> xf[i];
         if(i > 0)
         {
            dx[i-1] = xf[i] - xf[i-1];
            xc[i-1] = 0.5*(xf[i] + xf[i-1]);
         }   
      }   
      fi.close ();
   }
   else
   {
      n_cell = param.n_cell;
      n_face = n_cell + 1;
      double dx_uni = (param.xmax - param.xmin) / n_cell;
      xc.resize (n_cell);
      xf.resize (n_face);
      dx.resize(n_cell);
      for(unsigned int i=0; i<n_face; ++i)
      {
         xf[i] = param.xmin + i * dx_uni;
         if(i > 0)
         {
            dx[i-1] = dx_uni;
            xc[i-1] = 0.5*(xf[i] + xf[i-1]);
         }   
      }   
   }
   
   // Initial condition

   cout << "  Setting initial condition ...\n";
   primitive.resize (n_cell + 2*ng_cell, vector<double>(NVAR));
   residual.resize (n_cell, vector<double>(NVAR));
   conserved.resize (n_cell, vector<double>(NVAR));
   conserved_old.resize (n_cell, vector<double>(NVAR));
   tau_f.resize  (n_face);
   q_f.resize  (n_face);
   res_norm.resize (NVAR);
   if(param.ic_file_read)
   {
      ifstream fi;
      fi.open (param.ic_file.c_str());
      assert (fi.is_open());
      for(unsigned int i=0; i<n_cell; ++i)
         fi >> primitive[i][0]
            >> primitive[i][1]
            >> primitive[i][2];
      fi.close ();
   }
   else
   {
      for(unsigned int i=0; i<n_cell; ++i)
      {  
         param.initial_condition.value (xc[i],primitive[i+ng_cell]);
         assert (primitive[i+ng_cell][0] > 0.0);
         assert (primitive[i+ng_cell][2] > 0.0);
      }   
   }
   Cp  = (param.gamma * param.gas_const / (param.gamma - 1.0));
   apply_bc();
}

//------------------------------------------------------------------------------
// Apply Boundary condition
//------------------------------------------------------------------------------
void FV::apply_bc()
{
   if(param.bc == wall)
   {
      for(unsigned int i = 0; i< ng_cell; ++i)
      {
          primitive[i][0] = primitive[2*ng_cell - i - 1][0];
          primitive[i][1] = -primitive[2*ng_cell - i - 1][1];
          primitive[i][2] = primitive[2*ng_cell - i - 1][2];
          
          primitive[n_cell + ng_cell + i][0] = primitive[n_cell + ng_cell - i - 1][0];
          primitive[n_cell + ng_cell + i][1] = -primitive[n_cell + ng_cell - i - 1][1];
          primitive[n_cell + ng_cell + i][2] = primitive[n_cell + ng_cell - i - 1][2];
      }    
   }
   else if(param.bc == periodic)
   {
      for(unsigned int i = 0; i< ng_cell; ++i)
      {
          primitive[i][0] = primitive[n_cell - ng_cell + i - 1][0];
          primitive[i][1] = primitive[n_cell - ng_cell + i - 1][1];
          primitive[i][2] = primitive[n_cell - ng_cell + i - 1][2];
          
          primitive[n_cell + ng_cell + i][0] = primitive[ng_cell + i][0];
          primitive[n_cell + ng_cell + i][1] = primitive[ng_cell + i][1];
          primitive[n_cell + ng_cell + i][2] = primitive[ng_cell + i][2];
      }    
   }
   else if(param.bc == neumann)
   {
      for(unsigned int i = 0; i< ng_cell; ++i)
      {
          primitive[i][0] = primitive[ng_cell][0];
          primitive[i][1] = primitive[ng_cell][1];
          primitive[i][2] = primitive[ng_cell][2];
          
          primitive[n_cell + ng_cell + i][0] = primitive[n_cell-1][0];
          primitive[n_cell + ng_cell + i][1] = primitive[n_cell-1][1];
          primitive[n_cell + ng_cell + i][2] = primitive[n_cell-1][2];
      }    
   }
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FV::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      double speed = fabs(primitive[i+ng_cell][1]) + 
                     sqrt(param.gamma * primitive[i+ng_cell][2] / primitive[i+ng_cell][0]);
      dt = min (dt, dx[i]/speed);
      dt = min (dt,dx[i]*dx[i]/param.mu_ref);
   }

   dt *= param.cfl;
   
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
void FV::con_to_prim ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
      primitive[i+ng_cell][0] = conserved[i][0];
      primitive[i+ng_cell][1] = conserved[i][1]/conserved[i][0];
      primitive[i+ng_cell][2] = (param.gamma-1.0) * (conserved[i][2] - 
                           0.5 * pow(conserved[i][1], 2.0) / conserved[i][0]);
   }
   apply_bc();
}

//------------------------------------------------------------------------------
// Convert primitive to conserved
//------------------------------------------------------------------------------
void FV::prim_to_con ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
      conserved[i][0] = primitive[i+ng_cell][0];
      conserved[i][1] = primitive[i+ng_cell][0]*primitive[i+ng_cell][1];
      conserved[i][2] = 0.5*primitive[i+ng_cell][0]*pow(primitive[i+ng_cell][1],2.0) 
                        + primitive[i+ng_cell][2]/(param.gamma-1.0);              
   }
}

//------------------------------------------------------------------------------
// Compute derivatives at faces for tau and q
// Face centered derivatives
//------------------------------------------------------------------------------
void FV::compute_face_derivatives ()
{   
   double T_left, T_right, T, mu, k, dx_f;
   if(param.bc == periodic)
   {
      T_left  = temperature(primitive[ng_cell - 1]);
      T_right = temperature(primitive[ng_cell]);
      T   = 0.5 * (T_left + T_right);
      mu  = param.viscosity (T);
      k   = mu * Cp / param.Pr;
      dx_f = xc[0] - xf[0] + xf[n_face-1] - xc[n_cell-1];
      tau_f[0] = (4.0 * mu / 3.0 ) * 
                 (primitive[ng_cell][1] - primitive[ng_cell - 1][1]) / dx_f;
      q_f[0]   = -k * (T_right - T_left) / dx_f;  
      tau_f[n_face-1] = tau_f[0];
      q_f[n_face-1] = q_f[0];                      
   }
   else
   {
      tau_f [0] = 0.0;
      q_f   [0] = 0.0;
      tau_f[n_face-1] = tau_f[0];
      q_f[n_face-1] = q_f[0]; 
   }   
   
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      double T_left  = temperature(primitive[ng_cell + i - 1]);
      double T_right = temperature(primitive[ng_cell + i]);
      double T   = 0.5 * (T_left + T_right);
      double mu  = param.viscosity (T);
      double k   = mu * Cp / param.Pr;
      double dx_f = xc[i] - xc[i-1];

      // face centered values
      tau_f[i] = (4.0 * mu / 3.0 ) * 
                 (primitive[ng_cell + i][1] - primitive[ng_cell + i - 1][1]) / dx_f;
      q_f[i]   = -k * (T_right - T_left) / dx_f;
   }

}

// ------------------------------------------------------------------------------
// Compute tau and q ( at face centers generally)
// ------------------------------------------------------------------------------
// void FV::compute_tau_q (vector<double>& left,
//                                 vector<double>& right)
// {
//   double T_left  = temperature(left);
//   double T_right = temperature(right);
//   double T   = 0.5 * (T_left + T_right);
//   double mu  = viscosity (T);
//   double k   = mu * Cp / param.Pr;
// 
//   tau = (4.0 * mu / 3.0 ) * (right[1] - left[1]) / dx;
//   q= -k * (T_right - T_left) / dx;
//   
// }

//------------------------------------------------------------------------------
// Reconstruct left/right state at a face
//------------------------------------------------------------------------------
void FV::reconstruct (const unsigned int face,
                      vector<double>& left,
                      vector<double>& right) const
{
   if(param.reconstruct_scheme == first)
   {
      left = primitive[ng_cell+face-1];
      right= primitive[ng_cell+face];
   }
   else if(param.reconstruct_scheme == muscl)
   {
      left  = muscl_limited (primitive[ng_cell+face-2], primitive[ng_cell+face-1], primitive[ng_cell+face]);
      right = muscl_limited (primitive[ng_cell+face+1], primitive[ng_cell+face], primitive[ng_cell+face-1]);
   }
}



//------------------------------------------------------------------------------
// Numerical flux function
// Original Roe Scheme
//------------------------------------------------------------------------------
void FV::roe_flux(const vector<double>& left,
                  const vector<double>& right,
                  vector<double>& flux) const
{   
   double fl  = sqrt(left[0]);
   double fr  = sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   double Tl  = temperature (left);
   double Tr  = temperature (right);

   double Hl  =param.gamma*param.gas_const*Tl/(param.gamma-1.0) + 0.5*pow(left[1],2);
   double Hr  =param.gamma*param.gas_const*Tr/(param.gamma-1.0) + 0.5*pow(right[1],2);
   
   double El  = left[2]/(param.gamma-1.0) + 0.5*left[0]*pow(left[1],2);
   double Er  = right[2]/(param.gamma-1.0) + 0.5*right[0]*pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * pow(left[1],2) + 
                  right[2] + right[0] * pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add conservative dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = sqrt((param.gamma-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   
   double Lambda[] = { fabs(u-a), fabs(u), fabs(u+a)};
   
   double dU[] = {
      right[0] - left[0], right[0]*ur - left[0]*ul, Er - El
   };
   
   double aa[3];
   aa[1] = (param.gamma-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
   
}

//------------------------------------------------------------------------------
// Roe Scheme with entropy fix
//------------------------------------------------------------------------------
void FV::roe_fixed_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const
{   
   double fl  = sqrt(left[0]);
   double fr  = sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   double Tl  = temperature (left);
   double Tr  = temperature (right);

   double Hl  =param.gamma*param.gas_const*Tl/(param.gamma-1.0) + 0.5*pow(left[1],2);
   double Hr  =param.gamma*param.gas_const*Tr/(param.gamma-1.0) + 0.5*pow(right[1],2);
   
   double El  = left[2]/(param.gamma-1.0) + 0.5*left[0]*pow(left[1],2);
   double Er  = right[2]/(param.gamma-1.0) + 0.5*right[0]*pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * pow(left[1],2) + 
                  right[2] + right[0] * pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add conservative dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = sqrt((param.gamma-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   
   // Entropy fix
   double delta = 0.1*a;
   double Lambda1, Lambda3;
   if (fabs(u-a) > delta)
		Lambda1 = fabs(u-a);
   else
        Lambda1 = (pow(u-a,2) + pow(delta,2))/(2*delta);
   if (fabs(u+a) > delta)
		Lambda3 = fabs(u+a);
   else
        Lambda3 = (pow(u+a,2) + pow(delta,2))/(2*delta);
   	
   double Lambda[] = { Lambda1, fabs(u), Lambda3};
   
   double dU[] = {
      right[0] - left[0], right[0]*ur - left[0]*ul, Er - El
   };
   
   double aa[3];
   aa[1] = (param.gamma-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
   
}

//------------------------------------------------------------------------------
// Numerical flux function
// Rusanov Scheme
//------------------------------------------------------------------------------
void FV::rusanov_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const
{   
   double fl  = sqrt(left[0]);
   double fr  = sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   double Tl  = temperature (left);
   double Tr  = temperature (right);

   double Hl  = param.gamma*param.gas_const*Tl/(param.gamma-1.0) + 0.5*pow(left[1],2);
   double Hr  = param.gamma*param.gas_const*Tr/(param.gamma-1.0) + 0.5*pow(right[1],2);
   
   double El  = left[2]/(param.gamma-1.0) + 0.5*left[0]*pow(left[1],2);
   double Er  = right[2]/(param.gamma-1.0) + 0.5*right[0]*pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * pow(left[1],2) + 
                  right[2] + right[0] * pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add Rusanov dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = sqrt((param.gamma-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   
   double MaxLambda = fabs(u) + a;
   double Lambda[] = { MaxLambda,MaxLambda,MaxLambda};
   
   double dU[] = {
      right[0] - left[0], right[0]*ur - left[0]*ul, Er - El
   };
   
   double aa[3];
   aa[1] = (param.gamma-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
   
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FV::num_flux(const int&            f,
                  const vector<double>& left,
                  const vector<double>& right,
                  const double&         tau,
                  const double&         q,
                  vector<double>&       flux) const
{
   vector<double> flux_pos(NVAR);
   vector<double> flux_neg(NVAR);
   
   switch(param.flux_scheme)
   {
         case roe:
            roe_flux(left, right, flux);
            break;
         case roe_fixed:
            roe_fixed_flux(left, right, flux);
            break; 
         case rusanov:
            rusanov_flux(left, right, flux);
            break;       
             
   }
   
   double u = 0.5 * (left[1] + right[1]);
   // viscous flux
   flux[1] -= tau;
   flux[2] -= (u * tau - q);
     
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FV::compute_residual ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {   for(unsigned int j=0; j<NVAR; ++j)
         residual[i][j] = 0.0;
   }
   
   vector<double> flux (NVAR);
   vector<double> left (NVAR);
   vector<double> right(NVAR);
   
   // Flux through left boundary
   reconstruct (0, left, right);
   num_flux (0, left, right,tau_f[0], q_f[0], flux);
   for(unsigned int j=0; j<NVAR; ++j)
      residual[0][j] -= flux[j];
      
   // Flux through interior faces
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      reconstruct (i, left, right);
      num_flux (i, left, right, tau_f[i], q_f[i], flux);
      for(unsigned int j=0; j<NVAR; ++j)
      {
         residual[i-1][j] += flux[j];
         residual[i][j]   -= flux[j];           
      }
   }

   // Flux through right boundary
   reconstruct (n_face-1, left, right);
   num_flux (n_face-1, left, right,tau_f[n_face-1], q_f[n_face-1], flux);
   for(unsigned int j=0; j<NVAR; ++j)
       residual[n_cell-1][j] += flux[j];   
    
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FV::update_solution (const unsigned int rk)
{
   if (param.time_scheme == rk1 || param.time_scheme == ssprk3)
      for(unsigned int i=0; i<n_cell; ++i)
         for(unsigned int j=0; j<NVAR; ++j)
            conserved[i][j] = arks[rk] * conserved_old[i][j] +
                              brks[rk] * (conserved[i][j] - (dt/dx[i]) * residual[i][j]);
          
   else if (param.time_scheme == jameson_rk4)
      for(unsigned int i=0; i<n_cell; ++i)
         for(unsigned int j=0; j<NVAR; ++j)
            conserved[i][j] = conserved_old[i][j] +
                              jameson_rks[rk] * (- (dt/dx[i]) * residual[i][j]);   
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
void FV::output ()
{
   
   string filename = "sol";
   string extension = ".dat";

   string precount;
   if     (counter <= 9)    precount = "000";
   else if(counter <= 99)   precount = "00";
   else if(counter <= 999)  precount = "0";
   else if(counter <= 9999) precount = "";
   else
   {
      cout << "Writer::output: counter is too large !!!\n";
   }
   
   stringstream ss;
   ss <<precount<<counter;
   filename += ss.str();
   filename +=extension;
   filename = param.OUTPUT_PATH + "/" + filename;
    
   cout<<"Saving solutions in "<< filename<<endl;   
 
   ofstream fo(filename.c_str());
   for(unsigned int i=0; i<n_cell; ++i)
   {
      fo << xc[i] << " " 
         << primitive[i+ng_cell][0] << " " 
         << primitive[i+ng_cell][1] << " "
         << primitive[i+ng_cell][2] << endl;
      
   }
   fo.close ();
   counter++;
}


//------------------------------------------------------------------------------
// compute norm of residual
//------------------------------------------------------------------------------
void FV::compute_residual_norm ()
{
   for(unsigned int j=0; j<NVAR; ++j)
      res_norm[j] = 0;
   
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<NVAR; ++j)
         res_norm[j] += pow(residual[i][j],2);
   
   for(unsigned int j=0; j<NVAR; ++j)
   {
      res_norm[j] /= n_cell;
      res_norm[j] = sqrt(res_norm[j]);
   }
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FV::run ()
{
   initialize();
   prim_to_con ();
   counter = 0;
   output();
   double time = 0.0;
   unsigned int iter = 0;
   while (time < param.final_time && iter < param.max_iter)
   {
      conserved_old = conserved;
      compute_dt ();
      if(time+dt > param.final_time) dt = param.final_time - time;
      for(unsigned int rk=0; rk<param.n_rks; ++rk)
      {
         compute_face_derivatives ();
         compute_residual ();
         poiseulle_cor(); 
         source_terms(time);
         update_solution (rk);
         con_to_prim ();
      }
      time += dt;
      ++iter;
      compute_residual_norm();

      if(iter % param.write_frequency == 0 || time == param.final_time)
      {
         cout << "Iter = " << iter << " Time = " << time << endl;
         cout << "Residual norm = ";
         cout << res_norm[0] << " " << res_norm[1] << " " << res_norm[2] << endl;
         
         output ();
      }
   }
   

}

//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
double FV::temperature(const vector<double>& prim) const
{
   return prim[2] / (param.gas_const * prim[0]);
}

//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
double FV::enthalpy(const vector<double>& prim) const
{
   return param.gamma * prim[2] / (prim[0] * (param.gamma-1.0)) + 
          0.5 * pow(prim[1], 2);
}

// ------------------------------------------------------------------------------
// Hagen-Poiseulle correction
// ------------------------------------------------------------------------------
void FV::poiseulle_cor()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
     double T = temperature(primitive[ng_cell +i]);
     double mu  = param.viscosity (T);
     residual[i][1] += dx[i] * 8*mu*M_PI*primitive[ng_cell +i][1];
   }  
}

// ------------------------------------------------------------------------------
// Hagen-Poiseulle correction
// ------------------------------------------------------------------------------
void FV::source_terms(const double time)
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
     vector<double> source_val;
     source_val.resize(2);
     param.sources.value (xc[i],time,source_val);
     residual[i][0] += dx[i] * source_val[0];
     residual[i][2] += dx[i] * source_val[1];
   }  
}