#####################################################
##     Run Neutral Drift Stiff-Stalk Population    ##
#####################################################
##  Run neutral drift with sampling an initial     ##
##  population of 96 individuals. Those frequencies##
##  then are used to simulate drift in a population##
##  of 60 males and 120 females (see paper for     ##
##  details) over 15 generations. A final neutral  ##
##  drift simulation is run to sample 96           ##
##  individuals. This was submitted to a HPC and   ##
##  run 10 times with different random seeds due   ##
##  to memory constraints.                         ##
#####################################################

library(plyr)
#provide the path to this script; this is the work-horse in
# which neutral drift is simulated
source("DriftSimulator.R")

args <- commandArgs(trailingOnly = TRUE)

chunk <- args[1]
random_seed <- args[2]

num_gen <- 15

digits <- 2 ##number of digits to round to
step <- 0.01 ##interval between allele frequencies
sim <- 100000 ##number of simulations to run
total <- 180 ##total individuals/generation
males <- 60 ##total males/generation
females <- 120 ##total females/generation

set.seed(random_seed) ##insure we are getting a randomly-generated value

##create a vector of allele frequencies with step; exclude 0 and 1.
fs <- round(seq(step, 1 - step, by = step), digits)
rows <- length(fs)

##Simulate sampling 96 individuals with which to start
## the selection series (see paper for details).
##The sim frequencies generatedwill become the starting
## allele frequencies for the 15 generation drift simulation.
sample_96_gen0 <- sapply(
        fs, function(x) driftNeutral(
            x, total.pop = 96, no.males = 48,
            no.females = 48, sims = sim, cycles = 1
            )
    )
print(sample_96_gen0)
##Simulate neutral drift for 15 generations
gen15_cycles <- sapply(
    sample_96_gen0, function(x) driftNeutral(
        x, total.pop = total, no.males = males,
        no.females = females, sims = 1, cycles = num_gen
        )
    )

##Sample the final generation as 96 inviduals sampled after a final generation.
final_sample <- sapply(
    gen15_cycles, function(x) driftNeutral(
        x, total.pop = 96, no.males = 48,
        no.females = 48, sims = 1, cycles = 1
        )
    )

##change list to a matrix to write.
final_mat <- matrix(final_sample, nrow = rows, byrow = TRUE)

filename <- sprintf("Drift_Sim_SS_%s.txt", chunk)
write.table(final_mat, filename, row.name = FALSE, sep = "\t")

q()
n
