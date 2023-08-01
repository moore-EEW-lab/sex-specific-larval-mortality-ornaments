setwd('/Users/michaelmoore/desktop/Working Directory')

library(car)
library(MASS)
library(lme4)
library(lmerTest)
library(ggplot2)
library(nlme)
library(phytools)
library(MCMCglmm)


# load data, etc
dat <- read.csv('drags.larval.sex.ratio.csv')
phylo <- read.tree('drags.larval.sex.ratio.phylo.tre')

# for MCMCglmm, need to create a new column called "animal" for matching the phylogeny up to the species-level identities
dat$animal <- dat$binom
dat$animal.matched <- match(dat$animal, phylo$tip.label)
dat[which(is.na(dat$animal.matched)),]

phylo.tips <- data.frame(phylo$tip.label)
phylo.tips$matched <- match(phylo.tips$phylo.tip.label, dat$binom)
phylo.tips[which(is.na(phylo.tips$matched)),]


### compare sex ratios at the end of the larval stage between species that have male wing coloration and species that do not
prior.1 = list(R = list(V = 1, nu = 0.0001), G = list(G1  = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu = 0.0001))) # extremely weak priors

ratio01 <- MCMCglmm(fixed = prop.male ~ male.wing.color, random = ~binom + animal, prior = prior.1, pedigree = phylo, node = 'tips', burnin = 200000, thin = 1000, nitt = 1000000, data = dat, pr = TRUE, verbose = FALSE) 
plot(ratio01) # looks okay
summary(ratio01) # credible intervals do not overlap 0

 # Iterations = 200001:999001
 # Thinning interval  = 1000
 # Sample size  = 800 

 # DIC: 903.6759 

 # G-structure:  ~binom

      # post.mean  l-95% CI u-95% CI eff.samp
# binom    0.7244 1.685e-05     3.78      800

               # ~animal

       # post.mean  l-95% CI u-95% CI eff.samp
# animal    0.3856 1.563e-05    2.189    684.2

 # R-structure:  ~units

      # post.mean l-95% CI u-95% CI eff.samp
# units      36.2    28.05    45.47    910.9

 # Location effects: prop.male ~ male.wing.color 

                 # post.mean l-95% CI u-95% CI eff.samp  pMCMC   
# (Intercept)         47.875   46.402   49.248    715.7 <0.001 **
# male.wing.colory    -2.949   -5.638   -0.674    865.6 0.0175 * 

### run analysis without the one unusually low value

# first, subset the data to exclude the weirdo
dat[which(dat$prop.male < 20), ]
dat$z.prop.male <- (dat$prop.male - mean(dat$prop.male))/sd(dat$prop.male)
dat$outlier <- ifelse(dat$z.prop.male <= -3, 'y', 'n')
dat.sub <- dat[which(dat$outlier == 'n'), ]

# then generate a phylogeny for the subsetted data
match_phylo <- match(phylo$tip.label, dat.sub$binom)
drop.species <- phylo$tip.label[is.na(match_phylo)==TRUE]
phylo.sub <- drop.tip(phylo, drop.species)
summary(phylo.sub)

# create new "animal" column to make sure points and phylogeny align
dat.sub$animal <- dat.sub$binom
dat.sub$rematched <- match(dat.sub$animal, phylo.sub$tip.label)
dat.sub[which(is.na(dat.sub$rematched)),]

# run analysis
prior.2 = list(R = list(V = 1, nu = 0.0001), G = list(G1  = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu = 0.0001))) # extremely weak prior

ratio01b <- MCMCglmm(fixed = prop.male ~ male.wing.color, random = ~binom + animal, prior = prior.2, pedigree = phylo.sub, node = 'tips', burnin = 200000, thin = 100, nitt = 1000000, data = dat.sub, pr = TRUE, verbose = FALSE) 
plot(ratio01b) # looks okay
summary(ratio01b) # credible intervals do not overlap 0


  # Iterations = 200001:999901
 # Thinning interval  = 100
 # Sample size  = 8000 

 # DIC: 874.3793 

 # G-structure:  ~binom

      # post.mean  l-95% CI u-95% CI eff.samp
# binom    0.4621 1.201e-05    2.672     3440

               # ~animal

       # post.mean  l-95% CI u-95% CI eff.samp
# animal    0.3783 7.363e-06    1.997     2383

 # R-structure:  ~units

      # post.mean l-95% CI u-95% CI eff.samp
# units     30.71    23.43     38.4     8000

 # Location effects: prop.male ~ male.wing.color 

                 # post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)        47.8700  46.5869  49.1231     8024 <1e-04 ***
# male.wing.colory   -2.1747  -4.5005  -0.1587     7016  0.055 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##### for species with known body sizes

# first, subset dataset and phylogeny
dat.m <- dat[which(!is.na(dat$male.size)), ]
match_phylo.m <- match(phylo$tip.label, dat.m$binom)
drop.species <- phylo$tip.label[is.na(match_phylo.m)==TRUE]
phylo.m <- drop.tip(phylo, drop.species)
summary(phylo.m)

# create new vector to match to phylogeny
dat.m$animal <- dat.m$binom
dat.m$rematched <- match(dat.m$animal, phylo.m$tip.label)
dat.m[which(is.na(dat.m$rematched)),]


#run analysis
prior.3 = list(R = list(V = 1, nu = 0.0001), G = list(G1  = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu = 0.0001))) # extremely week prior

ratio.m.01 <- MCMCglmm(fixed = prop.male ~ male.wing.color + male.size, random = ~binom + animal, prior = prior.3, pedigree = phylo.m, node = 'tips', burnin = 200000, thin = 1000, nitt = 1000000, data = dat.m, pr = TRUE, verbose = FALSE) 
plot(ratio.m.01) # looks okay
summary(ratio.m.01) # credible intervals do not overlap 0

 # Iterations = 200001:999001
 # Thinning interval  = 1000
 # Sample size  = 800 

 # DIC: 890.8829 

 # G-structure:  ~binom

      # post.mean  l-95% CI u-95% CI eff.samp
# binom    0.7055 1.023e-05    4.195      800

               # ~animal

       # post.mean  l-95% CI u-95% CI eff.samp
# animal     1.158 2.163e-05    6.313      675

 # R-structure:  ~units

      # post.mean l-95% CI u-95% CI eff.samp
# units        36    27.93    45.83      800

 # Location effects: prop.male ~ male.wing.color + male.size 

                 # post.mean l-95% CI u-95% CI eff.samp  pMCMC   
# (Intercept)       46.60466 41.09481 52.04011      800 <0.001 **
# male.wing.colory  -2.84340 -5.38536 -0.27035      800 0.0325 * 
# male.size          0.02456 -0.06448  0.11839      800 0.6325   


## ornamented species with known amounts of ornamentation

# first, subset dataset and phylogeny
dat.o <- dat[which(!is.na(dat$estim.col.amount)), ]
match_phylo.o <- match(phylo$tip.label, dat.o$binom)
drop.species <- phylo$tip.label[is.na(match_phylo.o)==TRUE]
phylo.o <- drop.tip(phylo, drop.species)
summary(phylo.o)

# 12 species, 28 estimates

# create new vector to match to phylogeny
dat.o$animal <- dat.o$binom
dat.o$rematched <- match(dat.o$animal, phylo.o$tip.label)
dat.o[which(is.na(dat.o$rematched)),]

# run analysis
prior.4 = list(R = list(V = 1, nu = 0.0001), G = list(G1  = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu = 0.0001))) # extremely week prior

ratio.o.01 <- MCMCglmm(fixed = prop.male ~ estim.col.amount, random = ~binom + animal, prior = prior.4, pedigree = phylo.o, node = 'tips', burnin = 200000, thin = 1000, nitt = 1000000, data = dat.o, pr = TRUE, verbose = FALSE) 
plot(ratio.o.01) # looks okay
summary(ratio.o.01) 

 # Iterations = 200001:999001
 # Thinning interval  = 1000
 # Sample size  = 800 

 # DIC: 171.735 

 # G-structure:  ~binom

      # post.mean  l-95% CI u-95% CI eff.samp
# binom     8.508 1.728e-05    62.47    647.5

               # ~animal

       # post.mean  l-95% CI u-95% CI eff.samp
# animal     103.1 6.588e-05      249    720.9

 # R-structure:  ~units

      # post.mean l-95% CI u-95% CI eff.samp
# units     19.95    7.972    36.87    917.8

 # Location effects: prop.male ~ estim.col.amount 

                 # post.mean l-95% CI u-95% CI eff.samp  pMCMC   
# (Intercept)         40.493   26.294   53.042      800 <0.001 **
# estim.col.amount    -8.453  -50.570   33.789      800   0.67   

