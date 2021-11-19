#####################################################
##    Quantiles Drift Sim Stiff Stalk Selection    ## 
#####################################################
##  Combines the output files from all the neutral ##
##  drift simulation output, and  then applies     ##
##  quantiles at 1e-6 (0.999999) to 0.1 (0.9), by  ##
##  steps of powers of ten, of the simulated allele##
##  frequencies for each starting frequency.       ##
##  Returns a table of starting frequency and the  ##
##  frequency after all sampling events at the     ##
##  lower and upper quantiles.                     ##
#####################################################

step=0.01
digits=2

fs<-round(seq(step,1-step,by=step),digits)

##read in and combine all the output files from Run_Neutral_Drift_StiffStalk_Population.R
all_sims<-read.table('Drift_Sim_SS_1.txt', header=T,sep='\t')
for (i in c(2:10)){
  fname=sprintf('Drift_Sim_SS_%s.txt',i)
  holder=read.table(fname,header=T, sep='\t')
  cbind(all_sims, holder)
}

##Apply quantiles at the given probabilies to the simpulatate allele frequencies
q6<-apply(all_sims, 1, function(x) quantile(x, probs=c(0.000001,0.999999)))
q5<-apply(all_sims, 1, function(x) quantile(x, probs=c(0.00001,0.99999)))
q4<-apply(all_sims, 1, function(x) quantile(x, probs=c(0.0001,0.9999)))
q4<-apply(all_sims, 1, function(x) quantile(x, probs=c(0.001,0.999)))
q2<-apply(all_sims, 1, function(x) quantile(x, probs=c(0.01,0.99)))
q1<-apply(all_sims, 1, function(x) quantile(x, probs=c(0.1,0.9)))

compare_dataframe<-t(rbind(q1,q2,q3,q4,q5,q6)) ##create a data.frame of the upper and lower quantiles

compare_dataframe<-cbind(fs, compare_dataframe) ##add on the initial frequency data

##name columns for ease of understanding
colnames(compare_dataframe)<-c('Allele.Frequency.Gen0', 'Lower10', 'Upper90', 'Lower01', 'Upper99', 'Lower001', 'Upper999', 'Lower0001', 'Upper9999', 'Lower00001', 'Upper99999', 'Lower000001', 'Upper999999')

##write a tab-delimited file
write.table(compare_dataframe,"Quantiles_1milsim_15gen_01step_180tot_96sample_gen0_gen15_3.txt",row.name=F, sep= '\t')
q()
n
