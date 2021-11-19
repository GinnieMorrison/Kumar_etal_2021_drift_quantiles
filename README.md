# Kumar_etal_2021_drift_quantiles
I worked on the scripts that simulated drift in a neutral version of the divergently-selected synthentic stiff-stalk population and the calculated the quantiles of each allele frequency change. The R scripts included here were included in Kumar et al. 2021 (https://academic.oup.com/pcp/article/62/7/1199/6279219?login=true). I was only using the `apply` family of R functions and had to chunk it all up to run on a HPC. In order to force myself to practice numpy, I've decided to try and translate those scripts into numpy and see if I can't run it all at once, on my local box! By translate, I hope to not literally do a line-by-line translation, but accomplish something more like how a good translater translates a book into a new language.
Note: DriftSimulator.R is a program created by Timothy M Beissinger; please credit him (as per the script) when using DriftSimulator.R. 
