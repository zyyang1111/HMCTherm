#ifndef PDNCONFIG_H
#define PDNCONFIG_H

static double R_PCB = 0.031e-3; // [ohm]
static double R_PKG = 0.333e-3; // [ohm]
static double R_GRID = 285e-3; // [ohm]
static double R_BUMP = 40e-3; // [ohm]
static double R_PDNTSV = 114.7e-3; // [ohm]

static double C_PCB = 240e-6; // [F]
static double L_PCB = 21e-12; // [H]
static double C_PKG = 26e-6; // [F]
static double L_PKG = 120e-12; // [H]
static double L_GRID = 22e-12; // [H]
static double L_BUMP = 72e-12; // [H]
static double C_BUMP = 0.65e-12; // [F]
static double C_BLK = 7.8e-12; // [F]
static double C_DEN = 33e-5; // [F/m^2];
static double C_PDNTSV = 0.65e-12; // [F]
static double L_PDNTSV = 32e-12; // [H]

#endif