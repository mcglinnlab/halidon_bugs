library(mobr)
library(vegan)

#' read in data
dat <- read.csv('./data/Sweep_Species_r - Sheet1.csv')

head(dat)

# orders
dat$order[1:13]
uni_orders <- c('Hymenoptera', 'Orthoptera', 
                'Coleoptera', 'Areinada',
                'Lepidoptera', 'Dipteran',
                'Odonata', 'Mantodea',
                'Phasmotodea', 'Hemiptera',
                'Araneae', 'Oppliones',
                'unknown')
com_orders <- c('wasps', 'grasshoppers',
                'beetles','mites',
                'moths/butterflys','flys',
                'dragonflys', 'matids',
                'stickbugs','bugs',
                'spiders','harvestmen')


#' drop all samples where replicate is not known
dat <- subset(dat, !is.na(dat$rep))

dat$samp_id <- with(dat, paste(site, rep, sep = '-'))
              
comm <- with(dat, tapply(abu, list(samp_id, order_num),
                        sum, rm.na = TRUE))
comm

comm <- ifelse(is.na(comm), 0, comm)

#' drop samples with no individuals
comm <- subset(comm, rowSums(comm) > 0)

hab <- ifelse(grepl('HH', rownames(comm)), 'wetland', 'upland')

sum(comm)

#' compute abundances and richness at sweep net scale
N <- rowSums(comm)
S <- rowSums(comm > 0)

par(mfrow=c(1,2))
hist(N)
hist(S)

#pdf('./figs/N_&_S_habitats_rep_scale.pdf', width = 7*1.5)
par(mfrow=c(1,2))
boxplot(N ~ hab, ylab = '# of individuals')
boxplot(S ~ hab, ylab = '# of taxonomic orders')
#dev.off()

n_avg <- tapply(N, hab, mean, na.rm = TRUE)
n_sd <- tapply(N, hab, sd, na.rm = TRUE)
n_n <- length(N)
n_se <- n_sd / sqrt(n_n)

s_avg <- tapply(S, hab, mean, na.rm = TRUE)
s_sd <- tapply(S, hab, sd, na.rm = TRUE)
s_n <- length(S)
s_se <- s_sd / sqrt(s_n)


s_ord_avg <- apply(comm, 2, function(x) tapply(x, hab, mean))
s_ord_sd <- apply(comm, 2, function(x) tapply(x, hab, sd))
s_ord_n <- nrow(comm)
s_ord_se <- s_ord_sd / sqrt(s_ord_n)


#pdf('./figs/n_s_orders.pdf', width = 7*2)
par(mfrow=c(1,3))
n_plt <- barplot(n_avg, width = 0.25, ylim = c(0, 80), 
                  ylab = 'Number of Individuals', 
                 col = c('lightgreen', 'lightblue'))
arrows(n_plt, n_avg - (n_se * 1.96), y1 = n_avg + (n_se * 1.96),
       angle = 90, code = 3, length = 0.1)
s_plt <- barplot(s_avg, width = 0.25, ylim = c(0, 10), 
                 ylab = 'Number of Taxonomic Orders', 
                 col = c('lightgreen', 'lightblue'))
arrows(s_plt, s_avg - (s_se * 1.96), y1 = s_avg + (s_se * 1.96),
       angle = 90, code = 3, length = 0.1)
s_ord_plt <- barplot(s_ord_avg, beside = TRUE, ylim = c(0, 25),
                     ylab = 'Number of Individuals',
                     col = c('lightgreen', 'lightblue'))
arrows(s_ord_plt, s_ord_avg - (s_ord_se * 1.96), y1 = s_ord_avg + (s_ord_se * 1.96),
       angle = 90, code = 3, length = 0.01)
#dev.off()


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

