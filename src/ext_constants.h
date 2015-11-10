#ifndef __EXT_CONSTANTS_H__
#define __EXT_CONSTANTS_H__

#define   NVAR   3

// Coefficients for RK-schemes
static const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};
static const double jameson_rks[] = {1.0/4.0,1.0/3.0,1.0/2.0, 1.0};

enum ReconstructionScheme {first, muscl};
enum FluxScheme {roe, roe_fixed, rusanov};
enum TimeIntegrationScheme {rk1, ssprk3, jameson_rk4};
enum FlowModel {euler, ns};
enum MuModel {mu_constant, mu_sutherland, mu_power, mu_arc};
enum Boundary {wall, periodic, neumann};

#endif
