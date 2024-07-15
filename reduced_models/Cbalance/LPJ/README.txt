
## This file contains info for understanding simulation code and how it generates the results in Bathiany et al., 2024, GCB.


############### naming conventions in simulations and the paper

##### 1. naming convention of the model versions (Carbon and population density) in the code:
# A0 = prescribed allocated fractions
# A1 = interactive allocation
# A2 = C pools per ind are prescribed (but not constant), no allocation
   # note: prescribing constant pools would lead to no variability in N

# N0 = prescribed N from LPJ, mort is loss term for pools
# N1 = interactive N
# N2 = est is used as loss term for C pool dynamics
# N3 = const N in time, uses function N0 with constant N
# N7 = prescribe N as an AR1 process, uses function N0


## naming in the paper (which presents a reduced set of simulations):
 # C stands for C dynamics (pools per ind)
 # N stands for population dynamics

# A1N1 = LPJ-CNi (interactive mortality rate); LPJ-CN (fixed mortality rate)
# A1N3 = LPJ-C; N is fixed (constant)
# A1N2 = CiNp (not used in revised paper; N is prescribed, C interactive)
# A2N1 = LPJ-N; C is prescribed


###### 2. naming of variants for the population dynamics part (establishment and adjustment)

# M0 = mortality is constant in space and time, mort_prescribed=0.015
# M1 = interactive mortality (like original LPJmL)
# M2 = mort is constant in time, but at each grid point set to the time mean
        # mort_prescribed=np.mean(mort_LPJ[:,cellind+1])
# M3 = prescribed from LPJ: 
        # mort_prescribed=mort_LPJ[year_M,cellind+1]                


# E0 = est is constant in space and time, est_prescribed=0.01
# E1 = interactive est (like original model)
# E2 =  fixed to local mean
        # est_prescribed=np.mean(est_LPJ[:,cellind+1])
# E3 = prescribed but shuffled in time (destroys correlations)
        # est_prescribed=est_LPJ[year_E,cellind+1]   

# A0 = adjustment off
# A1 = adjustment on

# note: for N dynamics, allocated fluxes are used by default, but not actually needed if mort is not interactive


## naming in the paper (which presents a reduced set of simulations):

# A and E have the same meaning as above
 # Er = randomised establishment
 # Ei = interactive establishment
 # A0 = no adjustment
 # A1 = with adjustment       
        

