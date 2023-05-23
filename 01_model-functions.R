# Mutz, J., J.S. Thaler, T.A. Ugine, B.D. Inouye, and N. Underwood. Predator densities 
#  alter the influence of non-consumptive effects on the population dynamics of 
#  an agricultural pest.
# Supplemental code: 01_model-functions

# Code written and annotated by JM and BDI.
# The following code defines:
# (A) parameters as described in Tables 1 and S1;
# (B) hatching simulations to calculate number of cannibal and non-cannibal larvae 
#     emerging from each clutch;
# (C) functions for CPB densities at each life stage;
# (D) a function to iterate the model across one CPB generation (N_g=1,t=0 -> N_g=2,t=0).
# *Abbr: NCEs = non-consumptive effects, PE = predator-exposed, NE = not predator-exposed*

# clear environment
rm(list = ls(all = TRUE))

# Set working directory
#setwd("/Users/")

# Load required data.
load("04_clutch-sizes-and-hatching.RData")


# (A) Setting parameters.=============================================================

#Define parameter vector.
p.list <- c(
 lamPE = 6.5,      # DI reproductive rate (no. clutches) for PE adults
 lamNE = 7,        # DI reproductive rate (no. clutches) for NE adults
 b = 6,            # DD reproduction term
 n.max = 239,      # maximum clutch size
 reps = 100,       # number of reps for determining clutch size -> larvae transitions
 tPE = 0.15,       # prop. trophic eggs in clutches laid by PE mothers
 tNE = 0.02,       # prop. trophic eggs in clutches laid by NE mothers
 cPE = 0.26,       # prob. of hatchling eating a sib (in PE-laid clutches)
 cNE = 0.18,       # prob. of hatchling eating a sib (in NE-laid clutches)
 s = 0.8,          # baseline probability of larval survival during each instar
 a = 0.007,        # predator area searched per unit time (m^2)
 T = 28,           # search time (hours)
 h = 0.62,         # predator handling time (hours)
 mc = 0.62,        # mortality of cannibal larvae upon predator encounter
 ml = 0.71,        # mortality of regular larvae upon predator encounter
 ma = 0.2,         # mortality of adults upon predator encounter
 pup.surv = 0.7    # survival from pupa to adult
)


## (B) Clutch-level hatching simulation.===============================================

# First, define the function that will assign fates (regular larva, cannibal, trophic egg,
#  or dead) of eggs in a clutch.
# E = total eggs in a clutch
# V = viable eggs, T = Trophic eggs, E = V+T
# t = proportion of clutch that is inviable trophic eggs
# c = probability that a newly hatched larvae eats an egg
# L = number of Larvae walking away at the of the process.
# C = number of Cannibal larvae walking away at the end of the process. (L+C) <= E

clutchfate = function(E, c, t)  {
  V = ceiling(E*(1 - t))  ; Ttemp = E - V; Vtemp = V  ; L = 0  ;  C = 0 
  
  while(Vtemp > 0){
    Vtemp <- Vtemp - 1   #a viable egg is assigned during each run	
    random.num.1 <- runif(1)
    if(random.num.1 >= c) L <- L + 1    #if an egg not eaten, make an L
      if((Ttemp > 0 & random.num.1 < c) | (Vtemp > 1 & random.num.1 < c)) {   
        C <- C + 1                                 #otherwise, an egg is eaten and the eater becomes a C 
        Vprop <- (Vtemp / (Vtemp + Ttemp))         #proportion of available eggs that are V
        random.num.2 <- runif(1)
        Vtemp <- ifelse(random.num.2 <= Vprop, Vtemp - 1, Vtemp)      #there is one less V
        Ttemp <- ifelse(random.num.2 > Vprop, Ttemp - 1, Ttemp)       #or one less T,   
        Ttemp <- ifelse(Ttemp < 0, 0, Ttemp)          
      } #close if loop
  } #close while loop

  if (Ttemp > 0) {                                   #if there are still trophic eggs left,
    LT.temp <- ceiling( Ttemp*( L / (L + C) ))       #they get eaten by L and C in proportion to their abundances
    C <- ifelse( LT.temp <= L, C + LT.temp, C + L)   #L who eat remaining trophic eggs become Cs
    L <- ifelse( LT.temp <= L, L - LT.temp, 0)       #...meaning fewer Ls at the end of the day
  }

  return(c(L, C, Ttemp))
} #close function

## to visualize how the proportion of trophic eggs affects clutch fates
tvec = seq(0, .75, by = 0.05)  # range of "t" values to plot
L.to.plot = rep(0, length(tvec))
C.to.plot = rep(0, length(tvec))

for (i in 1: length(tvec)) {	#clutchfate for each t value
output = clutchfate(40, 0.3, tvec[i])
L.to.plot[i] = output[1]
C.to.plot[i] = output[2]
}

## plot it
plot(C.to.plot ~ tvec, ylim = c(0,30), xlab = "proportion trophic eggs", ylab = "C(black), L (red), total (blue)")
points(L.to.plot ~ tvec, col = "red")
points((L.to.plot + C.to.plot) ~ tvec, col = "dark blue", type = "b")

## NOT RUN -- output provided in 'array.cs.L_XXX_XXX' dataframe in '04_clutch-sizes-and-hatchlings.RData'.

# Next, build an array that: 
# (1) distributes eggs laid by NE and PE females into clutches based on empirically measured 
#     clutch size distributions, and
# (2) generates probability distributions for the number of regular larvae and cannibal larvae
#     emerging from clutches of a given size for defined values of t (prop. trophic eggs) and
#     c (prob. that a hatchling eats an egg) -- this comes from "clutchfate" simulations

#t.vals <- c(p.list["tPE"], p.list["tNE"])	#set values of t for which to generate prob. dist.
#c.vals <- c(p.list["cPE"], p.list["cNE"])	#set values of c for which to generate prob. dist.

# Create empty arrays for regular larvae (L) and cannibal larvae (C) with dimensions of range 
#  of clutch sizes x no. larvae emerging x values of t x values of c.
#array.cs.L <- array.cs.C <- array(NA, dim = c(p.list["n.max"], p.list["n.max"], length(t.vals), length(c.vals)),
#                                  dimnames = list(v = seq(1, p.list["n.max"]), h_x = seq(1, p.list["n.max"]),
#                                  t = c("tPE", "tNE"), c = c("cPE", "cNE")))


#for (m in 1:p.list["n.max"]) {        #for each clutch size (1...n.max)
#  for (p in 1:length(t.vals)) {       #and given values of t
#    for (q in 1:length(c.vals)) {     #and c
#      L.out <- rep(NA, p.list["reps"]); C.out <- rep(NA, p.list["reps"]); #empty vectors for each rep 
#      E = m; t = t.vals[p]; c = c.vals[q]; #assign values of E, t, and c for the 'clutchfate' function

# Run the 'clutchfate' function a bunch of times for the assigned E, c, and t values 
#  and save the output values of L and C (no. of regular and cannibal larvae).
#      for (i in 1:p.list["reps"]) {
#        clutch.temp <- clutchfate(E, c, t)
#        L.out[i] <- clutch.temp[1]; C.out[i] <- clutch.temp[2]
#	}

# Populate the array by counting up the number of times that L or C was a given value 
#  (starting with L or C = 1) across reps and dividing by the number of reps to get prob.
#      array.cs.L[m, , p, q] <- tabulate(L.out, nbins = p.list["n.max"])/p.list["reps"]
#      array.cs.C[m, , p, q] <- tabulate(C.out, nbins = p.list["n.max"])/p.list["reps"]
#}}}


## (C) Functions for CPB densities at each life stage.====================================

# In all functions below,
# P = predator density
# params = parameter vector

# Number of clutches laid by PE mothers (q1_PE) and NE mothers (q1_NE) (Eq. 1-2 in main text).
# n0_PE = number of PE adults (t = 0)
# n0_NE = number of NE adults (t = 0)

clutches <- function( n0_PE, n0_NE, params ) {
  q1_NE = (n0_NE*params["lamNE"]*params["b"])/(params["b"] + n0_PE + n0_NE) 
  q1_PE = (n0_PE*params["lamPE"]*params["b"])/(params["b"] + n0_PE + n0_NE)
  return(cbind(q1_PE, q1_NE))
}

# Number of regular or cannibal larvae emerging from PE or NE clutches.
# array = array that specifies the prop. of eggs expected to transition into larvae
#         (regular or cannibal) given clutch size, prop. trophic eggs (tPE and tNE), 
#         and prob. of cannibalism (cPE and cNE)
# vec = vector with number of clutches of each clutch size, 1-n.max
# XE = use dimensions of the array associated with PE or NE parameter values?

larv_trans <- function( array, vec, params, XE ) {
  u <-  t( array[,,paste0('t', XE),paste0('c', XE)] ) %*% vec*seq(1, params["n.max"] ) 
  return( sum(u) )
}

# Total number of regular larvae (n2_L) and cannibal larvae (n2_C) emerging after the 
#  hatching process (Eq. 3 in main text).
# q1_NE = number of clutches laid by NE mothers (t = 1)
# q1_PE = number of clutches laid by PE mothers (t = 1)
# cs.freq_NE = vector of clutch size frequencies for NE mothers
# cs.freq_PE = vector of clutch size frequencies for PE mothers
# array.cs_L = array that specifies the prop. of eggs expected to transition into regular 
#              larvae given clutch size, prop. trophic eggs (tPE and tNE), and prob. of 
#              cannibalism (cPE and cNE)
# array.cs_C = array that specifies the prop. of eggs expected to transition into cannibal
#              larvae given clutch size, prop. trophic eggs (tPE and tNE), and prob. of 
#              cannibalism (cPE and cNE)
# using 'larv_trans' function

larvae <- function( q1_PE, q1_NE, cs.freq_PE, cs.freq_NE, params, array.cs_L, array.cs_C ) {
# Vectors of the number of clutches of each size, [1,n.max], for PE and NE mothers (one column
#  for each value of q1_PE or q1_NE -- allows for vector values of q1_PE and q1_NE)
  cs.PE = matrix(q1_PE, nrow = length(cs.freq_PE), ncol = length(q1_PE), byrow = TRUE)*cs.freq_PE
  cs.NE = matrix(q1_NE, nrow = length(cs.freq_NE), ncol = length(q1_NE), byrow = TRUE)*cs.freq_NE
# Apply 'larv_func' to matrix columns to calculate the total number of regular and cannibal
#  larvae across clutch sizes and mother types (NE and PE)
  n2_L <- apply(cs.NE, 2, function(x) larv_trans(array.cs_L, x, params, "NE")) + 
          apply(cs.PE, 2, function(x) larv_trans(array.cs_L, x, params, "PE"))
  n2_C <- apply(cs.NE, 2, function(x) larv_trans(array.cs_C, x, params, "NE")) + 
          apply(cs.PE, 2, function(x) larv_trans(array.cs_C, x, params, "PE"))
  return(cbind(n2_L, n2_C))
}


# Functions for density-dependent predation of larvae:
# Prop. of larvae detected by predators during one instar (prop_detect; Eq. S5b in Appendix A) and 
#  larval density after one instar of predation of PE (larv_surv_PE) and NE (larv_surv_NE) larvae 
#  (see Eq. 4 in main text).
# n2ins_PE = number of PE larvae during previous instar
# n2ins_NE = number of NE larvae during previous instar
# n2ins_L = number of regular (i.e., non-cannibal) larvae during previous instar
# n2ins_C = number of cannibal larvae during previous instar
# m = mortality rate

prop_detect <- function( n2ins_L, n2ins_C, P, params ) {
  temp <- ( params["a"]*params["T"]*P ) / 
          ( 1 + params["a"]*params["h"]*(params["mc"]*n2ins_C + params["ml"]*n2ins_L) )
  sapply( temp, function(x) { min(x, 1) } )
}

larv_surv_PE <- function( n2ins_PE, n2ins_NE, n2ins_L, n2ins_C, P, m, params ) {
  params["s"] * (n2ins_NE + n2ins_PE) * prop_detect(n2ins_L = n2ins_L, n2ins_C = n2ins_C, P = P, params = params) * (1 - m) +
  params["s"] * n2ins_PE * (1 - prop_detect(n2ins_L = n2ins_L, n2ins_C = n2ins_C, P = P, params = params))
}

larv_surv_NE <- function( n2ins_PE, n2ins_NE, n2ins_L, n2ins_C, P, m, params ) {
  params["s"] * n2ins_NE * (1 - prop_detect(n2ins_L = n2ins_L, n2ins_C = n2ins_C, P = P, params = params))
}

# Number of PE pupae (n3_PE) and NE pupae (n3_NE) after predation across four larval instars
#  (see Eq. 4 in main text).
# n2_L = total no. of regular (i.e., non-cannibal) larvae emerging from clutches (t = 2)
# n2_C = total no. of cannibal larvae emerging from clutches (t = 2)

pupae_DD <- function( n2_L, n2_C, P, params ) {

# 1st to 2nd instar
  ins2_PE.C <- larv_surv_PE( n2ins_NE = n2_C, n2ins_PE = 0, n2ins_C = n2_C, n2ins_L = n2_L, 
                             P = P, m = params["mc"], params )
  ins2_NE.C <- larv_surv_NE( n2ins_NE = n2_C, n2ins_PE = 0, n2ins_C = n2_C, n2ins_L = n2_L, 
                             P = P, m = params["mc"], params )
  ins2_PE.L <- larv_surv_PE( n2ins_NE = n2_L, n2ins_PE = 0, n2ins_C = n2_C, n2ins_L = n2_L, 
                             P = P, m = params["ml"], params )
  ins2_NE.L <- larv_surv_NE( n2ins_NE = n2_L, n2ins_PE = 0, n2ins_C = n2_C, n2ins_L = n2_L, 
                             P = P, m = params["ml"], params )

# 2nd to 3rd instar
  ins3_PE.C <- larv_surv_PE( n2ins_NE = ins2_NE.C, n2ins_PE = ins2_PE.C, n2ins_C = (ins2_PE.C + ins2_NE.C), 
                             n2ins_L = (ins2_PE.L + ins2_NE.L), P = P, m = params["mc"], params)
  ins3_NE.C <- larv_surv_NE( n2ins_NE = ins2_NE.C, n2ins_PE = ins2_PE.C, n2ins_C = (ins2_PE.C + ins2_NE.C),
                             n2ins_L = (ins2_PE.L + ins2_NE.L), P = P, m = params["mc"], params)
  ins3_PE.L <- larv_surv_PE( n2ins_NE = ins2_NE.L, n2ins_PE = ins2_PE.L, n2ins_C = (ins2_PE.C + ins2_NE.C),
                             n2ins_L = (ins2_PE.L + ins2_NE.L), P = P, m = params["ml"], params)
  ins3_NE.L <- larv_surv_NE( n2ins_NE = ins2_NE.L, n2ins_PE = ins2_PE.L, n2ins_C = (ins2_PE.C + ins2_NE.C),
                             n2ins_L = (ins2_PE.L + ins2_NE.L), P = P, m = params["ml"], params)

# 3rd to 4th instar
  ins4_PE.C <- larv_surv_PE( n2ins_NE = ins3_NE.C, n2ins_PE = ins3_PE.C, n2ins_C = (ins3_PE.C + ins3_NE.C),
                             n2ins_L = (ins3_PE.L + ins3_NE.L), P = P, m = params["mc"], params)
  ins4_NE.C <- larv_surv_NE( n2ins_NE = ins3_NE.C, n2ins_PE = ins3_PE.C, n2ins_C = (ins3_PE.C + ins3_NE.C),
                             n2ins_L = (ins3_PE.L + ins3_NE.L), P = P, m = params["mc"], params)
  ins4_PE.L <- larv_surv_PE( n2ins_NE = ins3_NE.L, n2ins_PE = ins3_PE.L, n2ins_C = (ins3_PE.C + ins3_NE.C),
                             n2ins_L = (ins3_PE.L + ins3_NE.L), P = P, m = params["ml"], params)
  ins4_NE.L <- larv_surv_NE( n2ins_NE = ins3_NE.L, n2ins_PE = ins3_PE.L, n2ins_C = (ins3_PE.C + ins3_NE.C),
                             n2ins_L = (ins3_PE.L + ins3_NE.L), P = P, m = params["ml"], params)

# 4th instar to pupa
  pup_PE.C <- larv_surv_PE( n2ins_NE = ins4_NE.C, n2ins_PE = ins4_PE.C, n2ins_C = (ins4_PE.C + ins4_NE.C), 
                            n2ins_L = (ins4_PE.L + ins4_NE.L), P = P, m = params["mc"], params)
  pup_NE.C <- larv_surv_NE( n2ins_NE = ins4_NE.C, n2ins_PE = ins4_PE.C, n2ins_C = (ins4_PE.C + ins4_NE.C),
                            n2ins_L = (ins4_PE.L + ins4_NE.L), P = P, m = params["mc"], params)
  pup_PE.L <- larv_surv_PE( n2ins_NE = ins4_NE.L, n2ins_PE = ins4_PE.L, n2ins_C = (ins4_PE.C + ins4_NE.C),
                            n2ins_L = (ins4_PE.L + ins4_NE.L), P = P, m = params["ml"], params)
  pup_NE.L <- larv_surv_NE( n2ins_NE = ins4_NE.L, n2ins_PE = ins4_PE.L, n2ins_C = (ins4_PE.C + ins4_NE.C),
                            n2ins_L = (ins4_PE.L + ins4_NE.L), P = P, m = params["ml"], params)

  n3_PE <- pup_PE.C + pup_PE.L ; n3_NE <- pup_NE.C + pup_NE.L ;

  return( cbind( n3_PE, n3_NE ) )
}

# Prop. of adults detected by predators (Eq. S5b in Appendix A).
# n3 = density of pupae (t = 3)

prop_detect_ad <- function(n3, P, params) {
  temp <- (params["a"]*params["T"]*P) / 
          ( 1 + params["a"]*params["h"]*params["ma"]*n3 )
  min(temp, 1)
}

# Number of reproductive adults after predation (n0_PE_g1 and n0_NE_g1) (Eq. 5 in main text).
# n3_PE = density of PE pupae (t = 3)
# n3_NE = density of NE pupae (t = 3)

adults_DD <- function(n3_PE, n3_NE, P, params) {
  n0_PE_g1 <- params["pup.surv"] * (n3_NE + n3_PE) * prop_detect_ad(n3 = (n3_NE + n3_PE), P = P, params = params) * (1 - params["ma"]) +
              params["pup.surv"] * n3_PE * (1 - prop_detect_ad(n3 = (n3_NE + n3_PE), P = P, params = params))
  n0_NE_g1 <- params["pup.surv"] * n3_NE * ( 1 - prop_detect_ad(n3 = (n3_NE + n3_PE), P = P, params = params) )
  return( cbind(n0_PE_g1, n0_NE_g1) )
}


# (D) Function to iterate the model across one CPB generation.====================================

# Put it all together:
# Simulate across one generation, calculating values of
# - adult CPB density at time t
# - no. of clutches laid by PE and NE mothers (t + 1)
# - total no. regular and cannibal larvae (t + 2) 
# - no. of PE and NE pupae after four round of larval predation (t + 3)
# - no. of PE and NE adults surviving to reproduction (t + 4)
# using 'clutches', 'larvae', 'pupae_DD', and 'adults_DD' functions

CPB_iter <- function(n0_PE, n0_NE, params, array.cs_L = array.cs.L_0.15_0.26, array.cs_C = array.cs.C_0.15_0.26,
                     cs.freq_PE = cs.freq.PE, cs.freq_NE = cs.freq.NE, P, start.gen = 0) {

  clutches1 <- clutches(n0_PE = n0_PE, n0_NE = n0_NE, params = params)
  larvae1 <- larvae(q1_PE = clutches1[,'q1_PE'], q1_NE = clutches1[,'q1_NE'],
                    cs.freq_PE, cs.freq_NE, params, array.cs_L, array.cs_C)
  pupae1 <- pupae_DD(n2_L = larvae1[,'n2_L'], n2_C = larvae1[,'n2_C'], P = P, params = params)
  adults1 <- adults_DD(n3_PE = pupae1[,'n3_PE'], n3_NE = pupae1[,'n3_NE'], P = P, params = params)

  out <- data.frame(cbind(
           P = P, N0 = n0_PE + n0_NE, q1_PE = clutches1[,'q1_PE'], q1_NE = clutches1[,'q1_NE'],
           n2_L = larvae1[,'n2_L'], n2_C = larvae1[,'n2_C'], n3_PE = pupae1[,'n3_PE'], n3_NE = pupae1[,'n3_NE'],
           n0_PE = adults1[,'n0_PE_g1'], n0_NE = adults1[,'n0_NE_g1'], adults1[,'n0_PE_g1'] + adults1[,'n0_NE_g1'] ))
  colnames(out) <- c( sapply( colnames(out)[1:(length(colnames(out))-3)], 
                              function(x) paste0(x, ".g", start.gen) ),
                      sapply( colnames(out)[(length(colnames(out))-2):(length(colnames(out))-1)], 
                              function(x) paste0(x, ".g", start.gen + 1) ), paste0("N0.g", start.gen + 1) )
  return(out)
}



