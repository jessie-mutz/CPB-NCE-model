# Mutz, J., J.S. Thaler, T.A. Ugine, B.D. Inouye, and N. Underwood. Predator densities 
#  alter the influence of non-consumptive effects on the population dynamics of 
#  an agricultural pest.
# Supplemental code: 03_sensitivity

# Code written and annotated by JM.
# The following code:
# (A) defines a function for simulating the model across one generation for given values
#     of P, tPE, cPE, and lambdaPE;
# (B) simulates the model across two full generations (N_g=0,t=0 to N_g=2,t=0) for a range
#     of values of P, tPE, cPE, lambdaPE (output provided in '06_sensitivity-output.RData');
# (C) uses simulation output to estimate partial rank correlation coefficients (PRCCs) for
#     tPE, cPE, and lambdaPE for cumulative CPB density and mean population growth rate
#     at P_g=0 = P_g=1 = 0 and P_g=0 = P_g=1 = 8.
# *Abbr: NCEs = non-consumptive effects, PE = predator-exposed, NE = not predator-exposed*

# Load required packages.
require(sensitivity)

# Load required functions and data.
source("01_model-functions.R")
load("04_clutch-sizes-and-hatching.RData")


# (A) Function for simulating the model across one generation for given values of P, tPE,
#     cPE, and lambdaPE.================================================================

# vals.df = dataframe with parameter values and initial N values to use
# ref.df = dataframe with parameter values to use
# PE.prop = prop. of adults that are PE
# NE.prop = prop. of adults that are NE
# n0_PE = number of PE adults
# n0_NE = number of NE adults
# Uses arrays for prop. of eggs expected to transition into cannibal larvae (h_C; 'array.cs.C') and
#  non-cannibal larvae (h_L; 'array.cs.L') given clutch size (v), prop. trophic eggs (t), and 
#  prob. of cannibalism (c). Arrays were built using code in '01_model-functions.R' for values of 
#  tPE = {0.05,0.15,0.4} and cPE = {0.18,0.26.0.4}.

CPB_sens <- function(vals.df, ref.df, PE.prop, NE.prop, n0_PE, n0_NE, start.gen = 0) {
  out.df <- data.frame(matrix(nrow = 0, ncol = dim(vals.df)[2]+11))
  params.temp <- p.list

  for ( row in 1:length(ref.df[,1]) ) {    #for each parameter combination (i.e., row in ref.df)
    params.temp["mc"] <- ref.df$mc.val[row]   #assign mc value from ref.df
    params.temp["lamPE"] <- ref.df$lamPE.val[row]   #assign lamPE value from ref.df

    array.cs.L <- get(paste0("array.cs.L_", ref.df$tPE.val[row], "_", ref.df$cPE.val[row]))
    array.cs.C <- get(paste0("array.cs.C_", ref.df$tPE.val[row], "_", ref.df$cPE.val[row]))

    temp.out.bind <- data.frame(matrix(nrow = 0, ncol = 11))

    for ( pred in 1:length(preds) ) {  #run 'CPB_iter' for each value in 'preds' vector
      vals.df.temp <- vals.df[ which( vals.df$tPE.val == ref.df$tPE.val[row] &
                                      vals.df$cPE.val == ref.df$cPE.val[row] &
                                      vals.df$mc.val == ref.df$mc.val[row] &
                                      vals.df$lamPE.val == ref.df$lamPE.val[row] ),]
      temp.out <- cbind(vals.df.temp, CPB_iter(n0_PE = PE.prop*vals.df.temp[,n0_PE], 
                                               n0_NE = NE.prop*vals.df.temp[,n0_NE], 
                                               params = params.temp, array.cs_L = array.cs.L, 
                                               array.cs_C = array.cs.C, P = preds[pred], 
                                               start.gen = start.gen))
      temp.out.bind <- rbind(temp.out.bind, temp.out)
    } #close predator density loop

  out.df <- rbind(out.df, temp.out.bind)
  cat(row, "\n")
  } #close parameter value loop
  return(out.df)
} #close function


# (B) Simulating the model across two full generations, for a range of P, tPE, cPE,
#     and lamPE.=======================================================================

# Set values for initial adult CPB density (N_g=0,t=0).
Ns <- c(0.5)
# Set values for predator density.
preds <- seq(0, 8, by = 1)

# Create dataframe with parameter values and N for which to simulate CPB model ('vals.df')
vals <- data.frame(expand.grid(
 tPE.val = c(0.05, 0.15, 0.4),
 cPE.val = c(0.18, 0.26, 0.4),
 mc.val = c(0.62),
 lamPE.val = c(6.5, 7),
 N0 =  Ns ))

# Create dataframe to use for reference with parameter values ('ref.df')
ref <- data.frame(expand.grid(
 tPE.val = c(0.05, 0.15, 0.4),
 cPE.val = c(0.18, 0.26, 0.4),
 mc.val = c(0.62),
 lamPE.val = c(6.5, 7) ))

# With NCEs: Simulate the model across one generation using initial conditions and 
#  parameter values in CPB.sens.vals and CPB.sens.ref (above), for prop. adults in 
#  initial CPB population (PE.prop) = 0, 0.5, and 1. 
for ( PE.prop.temp in 0:2 ) {
  sens_out <- cbind(PE.prop = PE.prop.temp/2, 
                   CPB_sens(vals.df = vals, ref.df = ref, PE.prop = PE.prop.temp/2, 
                            NE.prop = 1 - PE.prop.temp/2, n0_PE = "N0", n0_NE = "N0"))
  assign(paste0("sens_out_", PE.prop.temp), sens_out)
}
sens_set_g0.1 <- rbind(sens_out_0, sens_out_1, sens_out_2)

# Add second generation -- PE and NE proportions are reduced from 1 to 0.75 to represent 
#  the ~25% of beetles that go directly into diapause after the first full generation.
sens_set_g0.2 <- CPB_sens(vals.df = sens_set_g0.1, ref.df = ref,
 PE.prop = 0.75, NE.prop = 0.75, n0_PE = "n0_PE.g1", n0_NE = "n0_NE.g1", start.gen = 1)

sens_set_g0.2$n0.g1.sum <- with(sens_set_g0.2, n0_NE.g1 + n0_PE.g1)

# Calculate cumulative beetle density (adults) and mean population growth rate across the 
#  season (two full CPB generations) for the output with and without NCEs.
sens_set_g0.2$cumul <- with( sens_set_g0.2, N0.g1 + N0.g2 )
sens_set_g0.2$mean.growth <- with(  sens_set_g0.2, sqrt((N0.g2/N0.g1)*(N0.g1/N0.g0)) )


# (C) Use simulation output to estimate partial rank correlation coefficients (PRCCs) for
#     tPE, cPE, and lambdaPE for cumulative CPB density and mean population growth rate
#     across values of P (predator density).================================================

# Function for estimating PRCCs across values of P (predator density).
# df = dataframe with simulation output
# response = column of simulation output with response variable of interest
# params.to.est = parameters for which to estimate PRCC

PRCC.func <- function(df, response, params.to.est) {
  out.df <- data.frame(matrix(nrow = 0, ncol = 7))
  for ( P.dens in 1:length(P.dens.vec) ) {
    temp <- subset(df, P.g0 == P.dens.vec[P.dens] & P.g1 == P.dens.vec[P.dens]) 
    sens1 <- pcc(temp[,params.to.est], temp[,response], nboot = 1000, rank = TRUE)
    out.temp1 <- cbind(param = params.to.est, P = P.dens.vec[P.dens], sens1$PRCC)
    out.df <- rbind(out.df, out.temp1)
  }
  return(out.df)
}

P.dens.vec <- c(0,8)
params.to.est <- c("tPE.val", "cPE.val", "lamPE.val")

PRCC_cumul <- PRCC.func( df = subset(sens_set_g0.2, PE.prop == 0.5), response = 'cumul',
                         params.to.est = params.to.est)
PRCC_growth <- PRCC.func( df = subset(sens_set_g0.2, PE.prop == 0.5), response = 'mean.growth',
                          params.to.est = params.to.est)

colnames(PRCC_cumul) <- colnames(PRCC_growth) <- c("param", "P", "original", "bias", "SE", "min.ci", "max.ci")


