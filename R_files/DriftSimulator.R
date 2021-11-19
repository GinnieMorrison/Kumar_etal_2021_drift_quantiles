###################################################
###  Function for Calculating Drift             ###
###################################################

## Timothy M Beissinger
## Updated: 1-25-2012

###################################################
# This function will simulate drift for a neutral locus.
#
# Load the function into R with the command: source("/path/to/this/file")
#
#
# FUNCTION PARAMETERS:
# initial = Assumed allele frequency in the initial population
# total.pop = population size in each generation (number of progeny)
# no.males = Number of males contributing gametes
# no.females = Number of females contributing gametes
# cycles = Number of generations of random mating
# plot = T/F, whether or not a plot should be drawn
#        If T, the plot will sequentially add a line as each
#        simulation is performed.
# sims = number of simulations to perform
#
# The function will return a vector with one element for each simulation performed.
# Each element of the returned vector represents the final allele frequency for one
# simulation.

driftNeutral <- function(initial,total.pop,no.males=total.pop/2,no.females=total.pop/2,cycles=30,plot=F,sims=1,sleep=0){
  prog <- c(rep("A",initial*total.pop*2),rep("B",(1-initial)*total.pop*2))
  maleprog <- c()
  femaleprog <- c()
  frequency <- c()
  frequency[1] <- length(which(prog=="A"))/length(prog)
  for(i in 1:cycles){
    maleprog <- sample(prog,size=no.males,replace=T)  #sample male alleles from total population
    maleprog <- rep(maleprog,length.out=total.pop)
    femaleprog <- sample(prog,size=no.females,replace=F) #sample females from total population
    femaleprog <- rep(femaleprog,length.out=total.pop)
    prog <- c(maleprog,femaleprog)
    frequency[i+1] <- length(which(prog=="A"))/(2*total.pop)
  }
  dat <- matrix(nrow={cycles+1},ncol=2)
  dat[,1] <- c(0:cycles)
  dat[,2] <- frequency
  finalfrequencies <- c()
  finalfrequencies[1] <- frequency[cycles+1]
  if(plot==T) plot(dat,type="l",ylim=c(0,1),xlab="Generation", ylab="Frequency")
  if(sims > 1){
    for(j in 2:sims){
        Sys.sleep(sleep)
        prog <- c(rep("A",initial*total.pop*2),rep("B",(1-initial)*total.pop*2))
        maleprog <- c()
        femaleprog <- c()
        frequency <- c()
        frequency[1] <- length(which(prog=="A"))/length(prog)
        for(i in 1:cycles){
          maleprog <- sample(prog,size=no.males,replace=T)  #sample male alleles from total population
          maleprog <- rep(maleprog,length.out=total.pop)
          femaleprog <- sample(prog,size=no.females,replace=F) #sample females from total population
          femaleprog <- rep(femaleprog,length.out=total.pop)
          prog <- c(maleprog,femaleprog)
          frequency[i+1] <- length(which(prog=="A"))/(2*total.pop)
        }
        dat <- matrix(nrow={cycles+1},ncol=2)
        dat[,1] <- c(0:cycles)
        dat[,2] <- frequency
        if(plot==T)lines(dat, type="l")
        finalfrequencies[j] <- frequency[cycles+1]
    }
  }
  return(finalfrequencies)
}

###
###
###
### Note: This function is free to use and distribute.  Please provide appropriate credit when doing so.
