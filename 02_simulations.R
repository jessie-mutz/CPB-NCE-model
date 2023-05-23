# Mutz, J., J.S. Thaler, T.A. Ugine, B.D. Inouye, and N. Underwood. Predator densities 
#  alter the influence of non-consumptive effects on the population dynamics of 
#  an agricultural pest.
# Supplemental code: 02_simulations

# Code written and annotated by JM.
# The following code:
# (A) defines a function for simulating the model across one generation for given values
#     of initial CPB density, P (predator density), and mc (mortality of cannibal larvae
#     upon predator encounter);
# (B) simulates the model across two full generations (N_g=0,t=0 to N_g=2,t=0) for a range
#     of values of P, mc, and prop. of PE adults in initial CPB population, with and without
#     NCEs (output provided in '05_simulation-output.RData');
# (C) uses simulation output to calculate: (i) absolute difference in cumulative CPB density with
#     and without NCEs, (ii) relative difference in mean population growth rate with and without
#     NCEs, and (iii) difference in effects of predator density (during first and second CPB 
#     generations) on cumulative density and mean population growth rate with and without NCEs. 
# *Abbr: NCEs = non-consumptive effects, PE = predator-exposed, NE = not predator-exposed*

source("01_model-functions.R")
load("04_clutch-sizes-and-hatching.RData")

# (A) Function for simulating the model for given values of initial CPB density, 
#     P (predator density), and mc (mortality of cannibal larvae upon predator 
#     encounter).======================================================================

# vals.df = dataframe with parameter values and initial N values to use
# ref.df = dataframe with parameter values to use
# PE.prop = prop. of adults that are PE
# NE.prop = prop. of adults that are NE
# n0_PE = number of PE adults
# n0_NE = number of NE adults

CPB_sim <- function(vals.df, ref.df, PE.prop, NE.prop, n0_PE, n0_NE, start.gen = 0) {
  out.df <- data.frame(matrix(nrow = 0, ncol = dim(vals.df)[2]+11))
  params.temp <- p.list

  for ( row in 1:length(ref.df[,1]) ) {    #for each parameter combination (i.e., row in ref.df)
    params.temp["mc"] <- ref.df$mc.val[row]   #assign mc value from ref.df
    params.temp["ml"] <- ref.df$ml.val[row]   #assign mc value from ref.df
    temp.out.bind <- data.frame(matrix(nrow = 0, ncol = 11))

    for ( pred in 1:length(preds) ) {  #run 'CPB_iter' for each value in 'preds' vector
      vals.df.temp <- vals.df[ which(vals.df$mc.val == ref.df$mc.val[row] & vals.df$ml.val == ref.df$ml.val[row]),]
      temp.out <- cbind( vals.df.temp[,-dim(vals.df.temp)[2]],
                         CPB_iter(n0_PE = PE.prop*vals.df.temp[,n0_PE], n0_NE = NE.prop*vals.df.temp[,n0_NE], 
                                  params = params.temp, P = preds[pred], start.gen = start.gen) )
      temp.out.bind <- rbind(temp.out.bind, temp.out)
    } #close predator density loop

  out.df <- rbind(out.df, temp.out.bind)
  cat(row, "\n")
  } #close parameter value loop
  return(out.df)
} #close function


# (B) Simulating the model across two full generations, for a range of P, mc,
#     and PE.prop .============================================================

# Set values for initial adult CPB density (N_g=0,t=0).
Ns <- c(0.5)
# Set values for predator density.
preds <- seq(0, 8, by = 0.2)

# Create dataframe with values of mc and N for which to simulate CPB model ('vals.df')
vals <- data.frame(expand.grid(
  ml.val = c(0.71),
  mc.val = c(0.35, 0.53, 0.62),
  N0 =  Ns))
# Create dataframe to use for reference with values of mc ('ref.df')
ref <- data.frame(expand.grid(
  ml.val = c(0.71),
  mc.val = c(0.35, 0.53, 0.62) ))

# With NCEs: Simulate the model across one generation using initial conditions and 
#  parameter values in CPB.sens.vals and CPB.sens.ref (above), for prop. adults in 
#  initial CPB population (PE.prop) = 0, 0.5, and 1. 
for ( PE.prop.temp in 0:2 ) {
  sim_out <- cbind(PE.prop = PE.prop.temp/2, 
                   CPB_sim(vals.df = vals, ref.df = ref, PE.prop = PE.prop.temp/2, 
                           NE.prop = 1 - PE.prop.temp/2, n0_PE = "N0", n0_NE = "N0"))
  assign(paste0("sim_out_", PE.prop.temp), sim_out)
}
sim_set_g0.1 <- rbind(sim_out_0, sim_out_1, sim_out_2)

# Add second generation -- PE and NE proportions are reduced from 1 to 0.75 to represent 
#  the ~25% of beetles that go directly into diapause after the first full generation.
sim_set_g0.2 <- CPB_sim(vals.df = sim_set_g0.1, ref.df = ref,
 PE.prop = 0.75, NE.prop = 0.75, n0_PE = "n0_PE.g1", n0_NE = "n0_NE.g1", start.gen = 1)

# Without NCEs: Assuming no effects of previous predator exposure -- for the first generation,
#  this is just PE.prop = 0. Add second generation, assigning all adults to the "NE" category.
sim_set_g0.2_noNCEs <- CPB_sim(vals.df = subset(sim_set_g0.1, PE.prop == 0), ref.df = ref,
 PE.prop = 0, NE.prop = 0.75, n0_PE = "N0.g1", n0_NE = "N0.g1", start.gen = 1)

sim_set_g0.2$n0.g1.sum <- with(sim_set_g0.2, n0_NE.g1 + n0_PE.g1)
sim_set_g0.2_noNCEs$n0.g1.sum <- with(sim_set_g0.2_noNCEs, n0_NE.g1 + n0_PE.g1)

# Calculate cumulative beetle density (adults) and mean population growth rate across the 
#  season (two full CPB generations) for the output with and without NCEs.
sim_set_g0.2$cumul <- with( sim_set_g0.2, n0.g1.sum + N0.g2 )
sim_set_g0.2$mean.growth <- with(  sim_set_g0.2, sqrt((N0.g2/N0.g1)*(n0.g1.sum/N0.g0)) )

sim_set_g0.2_noNCEs$cumul <- with( sim_set_g0.2_noNCEs, n0.g1.sum + N0.g2 )
sim_set_g0.2_noNCEs$mean.growth <- with(  sim_set_g0.2_noNCEs, sqrt((N0.g2/N0.g1)*(n0.g1.sum/N0.g0)) )


# (C) Using simulation output to calculate differences in cumulative CPB density
#     and mean population growth rate, with and without NCEs.===========================

# Dataframe for differences in cumulative CPB density and mean population growth rate
#  between simulations with and without NCEs.
diff.df <- subset(sim_set_g0.2, PE.prop == 1)
diff.df$cumul.noNCEs <- with( sim_set_g0.2_noNCEs, cumul )
diff.df$diff.cumul <- with( diff.df, cumul - cumul.noNCEs )
diff.df$growth.noNCEs <- with( sim_set_g0.2_noNCEs, mean.growth )
diff.df$rel.growth <- with( diff.df, round ( mean.growth / growth.noNCEs, 3 ) )  

# Define function to calculate the effect of P.g0 (predator density during generation g = 0) 
#  and P.g1 (predator density during generation g = 1) on cumulative CPB density and mean
#  population growth rate.
P_effect <- function(sim_set) {

  df <- get('sim_set')

  for ( i in 1:dim(df)[1] ) {
    temp <- df[i,]
    df[i,c('cumul_P.g0','growth_P.g0')] <- 
      df[which( with(df, mc.val == temp[,'mc.val'] & PE.prop == temp[,'PE.prop'] & 
                P.g0 == 0 & P.g1 == temp[,'P.g1'])), c('cumul', 'mean.growth')]
    df[i,c('cumul_P.g1','growth_P.g1')] <- 
      df[which( with(df, mc.val == temp[,'mc.val'] & PE.prop == temp[,'PE.prop'] & 
                P.g0 == temp[,'P.g0'] & P.g1 == 0)), c('cumul', 'mean.growth')]
  }

  df$P.g0_effe_cumul <- with( df, (cumul - cumul_P.g0) / cumul_P.g0 )
  df$P.g1_effe_cumul <- with( df, (cumul - cumul_P.g1) / cumul_P.g1 )

  df$P.g0_effe_growth <- with( df, (mean.growth - growth_P.g0) / growth_P.g0 )
  df$P.g1_effe_growth <- with( df, (mean.growth - growth_P.g1) / growth_P.g1 )

  return(df)

}

# Apply the function 'P_effect' to calculate effects of P.g0 and P.g1 in simulations with
#  and without NCEs.
P.effect.df <- P_effect(sim_set_g0.2)
P.effect.df_noNCEs <- P_effect(sim_set_g0.2_noNCEs)

# Dataframe for differences in the effect of P.g0 and P.g1 on cumulative CPB density and mean 
#  population growth rate between simulations with and without NCEs.
P.effect.df <- subset(P.effect.df, PE.prop == 1)

P.effect.df$P.g0_effe_cumul.noNCEs <- with( P.effect.df_noNCEs, P.g0_effe_cumul )
P.effect.df$P.g1_effe_cumul.noNCEs <- with( P.effect.df_noNCEs, P.g1_effe_cumul )
P.effect.df$P.g0_effe_growth.noNCEs <- with( P.effect.df_noNCEs, P.g0_effe_growth )
P.effect.df$P.g1_effe_growth.noNCEs <- with( P.effect.df_noNCEs, P.g1_effe_growth )

P.effect.df$P.g0_effe_cumul_rel <- with( P.effect.df, (P.g0_effe_cumul - P.g0_effe_cumul.noNCEs) / P.g0_effe_cumul )
P.effect.df$P.g1_effe_cumul_rel <- with( P.effect.df, (P.g1_effe_cumul - P.g1_effe_cumul.noNCEs) / P.g1_effe_cumul )
P.effect.df$P.g0_effe_growth_rel <- with( P.effect.df, (P.g0_effe_growth - P.g0_effe_growth.noNCEs) / P.g0_effe_growth )
P.effect.df$P.g1_effe_growth_rel <- with( P.effect.df, (P.g1_effe_growth - P.g1_effe_growth.noNCEs) / P.g1_effe_growth )



