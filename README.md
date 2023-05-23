# Mutz, J., J. S. Thaler, T. A. Ugine, B. D. Inouye, and N. Underwood.
# Predator densities alter the influence of non-consumptive effects on the population dynamics of an agricultural pest.
# Guide to supplemental code and simulation output.
---

## R Scripts.
**01_model-functions.R**   
The code in this script file defines:   
(A) parameters as described in Tables 1 and S1;   
(B) hatching simulations to calculate number of cannibal and non-cannibal larvae 
    emerging from each clutch;   
\(C\) functions for CPB densities at each life stage;   
(D) a function to iterate the model across one CPB generation (N_g=1,t=0 -> N_g=2,t=0).   

**02_simulations.R**   
The code in this script file:   
(A) defines a function for simulating the model across one generation for given values
    of initial CPB density, P (predator density), and mc (mortality of cannibal larvae
    upon predator encounter);   
(B) simulates the model across two full generations (N_g=0,t=0 to N_g=2,t=0) for a range
    of values of P, mc, and prop. of PE adults in initial CPB population, with and without
    NCEs (output provided in '05_simulation-output.RData');   
\(C\) uses simulation output to calculate: (i) absolute difference in cumulative CPB density with
    and without NCEs, (ii) relative difference in mean population growth rate with and without
    NCEs, and (iii) difference in effects of predator density (during first and second CPB 
    generations) on cumulative density and mean population growth rate with and without NCEs.    

**03_sensitivity.R**   
The code in this script file:   
(A) defines a function for simulating the model across one generation for given values
    of P, tPE, cPE, and lambdaPE;   
(B) simulates the model across two full generations (N_g=0,t=0 to N_g=2,t=0) for a range
    of values of P, tPE, cPE, lambdaPE (output provided in '06_sensitivity-output.RData');   
\(C\) uses simulation output to estimate partial rank correlation coefficients (PRCCs) for
    tPE, cPE, and lambdaPE for cumulative CPB density and mean population growth rate
    across values of P (predator density).   

## R Workspaces.
**04_clutch-sizes-and-hatching.RData**   
Contains the following objects:   

| Object | Description |  
| ------ | ----------- |  
| cs.freq.NE | Vector of length n.max with frequencies of clutch sizes 1...n.max for NE mothers |  
| cs.freq.PE | Vector of length n.max with frequencies of clutch sizes 1...n.max for PE mothers |  
| array.cs.L_0.05_0.18 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.05, tNE = 0.02} x *c* = {cPE = 0.18, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.05_0.18 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.05, tNE = 0.02} x *c* = {cPE = 0.18, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.05_0.26 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.05, tNE = 0.02} x *c* = {cPE = 0.26, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.05_0.26 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.05, tNE = 0.02} x *c* = {cPE = 0.26, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.05_0.4 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.05, tNE = 0.02} x *c* = {cPE = 0.4, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.05_0.4 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.05, tNE = 0.02} x *c* = {cPE = 0.4, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.15_0.18 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.15, tNE = 0.02} x *c* = {cPE = 0.18, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.15_0.18 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.15, tNE = 0.02} x *c* = {cPE = 0.18, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.15_0.26 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.15, tNE = 0.02} x *c* = {cPE = 0.26, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |   
| array.cs.C_0.15_0.26 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.15, tNE = 0.02} x *c* = {cPE = 0.26, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.15_0.4 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.15, tNE = 0.02} x *c* = {cPE = 0.4, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.15_0.4 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.15, tNE = 0.02} x *c* = {cPE = 0.4, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.4_0.18 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.4, tNE = 0.02} x *c* = {cPE = 0.18, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |   
| array.cs.C_0.4_0.18 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.4, tNE = 0.02} x *c* = {cPE = 0.18, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.4_0.26 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.4, tNE = 0.02} x *c* = {cPE = 0.26, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.4_0.26 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.4, tNE = 0.02} x *c* = {cPE = 0.26, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.L_0.4_0.4 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.4, tNE = 0.02} x *c* = {cPE = 0.4, cNE = 0.18}, containing prop. of eggs expected to transition into non-cannibal larvae (*h_L*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  
| array.cs.C_0.4_0.4 | Array of dimensions *v* = 1...n.max x *h_x* = 1...n.max x *t* = {tPE = 0.4, tNE = 0.02} x *c* = {cPE = 0.4, cNE = 0.18}, containing prop. of eggs expected to transition into cannibal larvae (*h_C*) given clutch size (*v*), prop. trophic eggs (*t*), and prob. of cannibalism (*c*) |  


**05_simulation-output.RData**   
Contains dataframes 'sim_set_g0.2' and 'sim_set_g0.2_noNCEs', with simulation output across values of P and mc, with and without NCEs, respectively.   

| Column name | Description |   
| ----------- | ----------- |   
| PE.prop | Prop. of adults in g = 0 that have previously been exposed to predators. |   
| ml.val | Mortality rate of non-cannibal larvae. |   
| mc.val | Mortality rate of cannibal larvae. |   
| P.g0 | Predator density during g = 0. |   
| N0.g0 | Initial adult CPB density at t = 0, g = 0. |   
| q1_PE.g0 | Number of clutches laid by PE adults at t = 1, g = 0. |   
| q1_NE.g0 | Number of clutches laid by NE adults at t = 1, g = 0. |   
| n2_L.g0 | Number of non-cannibal larvae at t = 2, g = 0. |   
| n2_C.g0 | Number of cannibal larvae at t = 2, g = 0. |   
| n3_PE.g0 | Number of PE pupae at t = 3, g = 0. |   
| n3_NE.g0 | Number of NE pupae at t = 3, g = 0. |   
| n0_PE.g1 | Number of PE adults at t = 0, g = 1. |   
| n0_NE.g1 | Number of NE adults at t = 0, g = 1. |   
| P.g1 | Predator density during g = 1. |   
| N0.g1 | Initial adult CPB density at t = 0, g = 1 |   
| q1_PE.g1 | Number of clutches laid by PE adults at t = 1, g = 1. |   
| q1_NE.g1 | Number of clutches laid by NE adults at t = 1, g = 1. |   
| n2_L.g1 | Number of non-cannibal larvae at t = 2, g = 1. |   
| n2_C.g1 | Number of cannibal larvae at t = 2, g = 1. |   
| n3_PE.g1 | Number of PE pupae at t = 3, g = 1. |   
| n3_NE.g1 | Number of NE pupae at t = 3, g = 1. |   
| n0_PE.g2 | Number of PE adults at t = 0, g = 2. |   
| n0_NE.g2 | Number of NE adults at t = 0, g = 2. |   
| N0.g2 | Total number of adults at t = 0, g = 2. |   
| n0.g1.sum | Total number of adults at t = 0, g = 1 (n0_PE.g1 + n0_NE.g1). |   
| cumul | Cumulative CPB density (n0.g1.sum + N0.g2) |   
| mean.growth | Geometric mean growth rate ( sqrt((N0.g2/N0.g1)*(n0.g1.sum/N0.g0)) ) |   


**06_sensitivity-output.RData**   
Contains dataframe 'sens_set_g0.2', with simulation output across values of tPE, cPE, and lamPE.   

| Object | Description |   
| ------ | ----------- |   
| PE.prop | Prop. of adults in g = 0 that have previously been exposed to predators. |   
| tPE.val | Prop. of each clutch that are non-viable trophic eggs for clutches laid by PE mothers. |   
| cPE.val | Prob. that newly hatched larva eats an egg in clutches laid by PE mothers. |   
| mc.val | Mortality rate of cannibal larvae. |   
| lamPE.val | Density-independent per capita fecundity (no. clutches) for PE mothers. |   
| P.g0 | Predator density during g = 0. |   
| N0.g0 | Initial adult CPB density at t = 0, g = 0. |   
| q1_PE.g0 | Number of clutches laid by PE adults at t = 1, g = 0. |   
| q1_NE.g0 | Number of clutches laid by NE adults at t = 1, g = 0. |     
| n2_L.g0 | Number of non-cannibal larvae at t = 2, g = 0. |   
| n2_C.g0 | Number of cannibal larvae at t = 2, g = 0. |   
| n3_PE.g0 | Number of PE pupae at t = 3, g = 0. |   
| n3_NE.g0 | Number of NE pupae at t = 3, g = 0. |   
| n0_PE.g1 | Number of PE adults at t = 0, g = 1. |   
| n0_NE.g1 | Number of NE adults at t = 0, g = 1. |   
| P.g1 | Predator density during g = 1. |   
| N0.g1 | Initial adult CPB density at t = 0, g = 1. |   
| q1_PE.g1 | Number of clutches laid by PE adults at t = 1, g = 1. |   
| q1_NE.g1 | Number of clutches laid by NE adults at t = 1, g = 1. |   
| n2_L.g1 | Number of non-cannibal larvae at t = 2, g = 1. |   
| n2_C.g1 | Number of cannibal larvae at t = 2, g = 1. |   
| n3_PE.g1 | Number of PE pupae at t = 3, g = 1. |   
| n3_NE.g1 | Number of NE pupae at t = 3, g = 1. |   
| n0_PE.g2 | Number of PE adults at t = 0, g = 2. |   
| n0_NE.g2 | Number of NE adults at t = 0, g = 2. |   
| N0.g2 | Total number of adults at t = 0, g = 2. |   
| n0.g1.sum | Total number of adults at t = 0, g = 1 (n0_PE.g1 + n0_NE.g1). |   
| cumul | Cumulative CPB density (n0.g1.sum + N0.g2) |   
| mean.growth | Geometric mean growth rate ( sqrt((N0.g2/N0.g1)*(n0.g1.sum/N0.g0)) ) |   
