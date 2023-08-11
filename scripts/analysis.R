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

t.test(N ~ hab)
t.test(S ~ hab)

bug_mob <- make_mob_in(comm, data.frame(hab))

par(mfrow=c(1,3))
plot_rarefaction(bug_mob, "hab", 
                 ref_level = "upland", method = 'IBR',
                 scales = c('alpha', 'gamma', 'study'),
                 leg_loc = 'bottomright', log = 'xy')



SBR_wet <- rarefaction(comm[hab == 'wetland',], 'SBR')
SBR_up <- rarefaction(comm[hab == 'upland',], 'SBR')

plot(SBR_wet, type='l')
lines(SBR_up, col='red')

anova(cca(comm ~ hab))
