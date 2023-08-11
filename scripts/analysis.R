library(mobr)
library(vegan)

dat <- read.csv('./data/Sweep_Species_r - Sheet1.csv')

head(dat)


# drop all samples where replicate is not known
dat <- subset(dat, !is.na(dat$rep))

dat$samp_id <- with(dat, paste(site, rep, sep = '-'))
              
comm <- with(dat, tapply(abu, list(samp_id, order_num),
                        sum, rm.na = TRUE))
comm

comm <- ifelse(is.na(comm), 0, comm)

# drop samples with no individuals
comm <- subset(comm, rowSums(comm) > 0)

hab <- ifelse(grepl('HH', rownames(comm)), 'wetland', 'upland')

sum(comm)

N <- rowSums(comm)
S <- rowSums(comm > 0)

par(mfrow=c(2,2))
hist(N)
hist(S)
boxplot(N ~ hab)
boxplot(S ~ hab)

# test if means different at the replicate scale
t.test(N ~ hab)
#Welch Two Sample t-test
#data:  N by hab
#t = -0.42447, df = 47.234, p-value = 0.6731
#alternative hypothesis: true difference in means between group upland and group wetland is not equal to 0
#95 percent confidence interval:
#  -19.52221  12.71862
#sample estimates:
#  mean in group upland mean in group wetland 
#62.3913               65.7931 
t.test(S ~ hab)
#Welch Two Sample t-test
#data:  S by hab
#t = 2.2668, df = 41.465, p-value = 0.02868
#alternative hypothesis: true difference in means between group upland and group wetland is not equal to 0
#95 percent confidence interval:
#  0.1235041 2.1343670
#sample estimates:
#  mean in group upland mean in group wetland 
#8.956522              7.827586 

bug_mob <- make_mob_in(comm, data.frame(hab))

pdf('./figs/up_vs_wet_ibr.pdf', width = 7*1.75, height = 5)
par(mfrow=c(1,3))
plot_rarefaction(bug_mob, "hab", 
                 ref_level = "upland", method = 'IBR',
                 scales = c('alpha', 'gamma', 'study'),
                 leg_loc = 'bottomright', log = 'xy')
dev.off()


SBR_wet <- rarefaction(comm[hab == 'wetland',], 'SBR')
SBR_up <- rarefaction(comm[hab == 'upland',], 'SBR')

pdf('./figs/up_vs_wet_sbr.pdf')
plot(SBR_wet, type='l', col="#78D3EC", lwd =2,
     xlab = '# of sweep net samples', 
     ylab = '# of taxonomic orders')
lines(SBR_up, col="#FFB3B5", lwd = 2)
legend('bottomright', c('upland', 'wetland'),
       col = c("#FFB3B5","#78D3EC"), lty = 1, lwd =2,
       bty = 'n')
dev.off()

anova(cca(comm ~ hab))

