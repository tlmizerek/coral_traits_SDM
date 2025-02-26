rm(list=ls())
library(reshape2) #melt 
library(lme4) 
library(biomod2) 
library(ncdf4)

source('R/functions.R')

#traits
#corallite size / turbidity, sediment rejection increased with larger corallite size
#growth form / wave power, robust growth forms are more resilient to wave energy
#growth rate / wave power, if you grow fast, you can be more stable to energy
#calcification or growth rate / SST, higher growth rate does better in higher temp??


#symbiodinium diversity or presence of a particular clade
#symbiodinium density
#range?
#colony max size


CTDB <- read.csv("data/ctdb_20150525.csv", as.is=TRUE)
CTDB[is.na(CTDB) == T] = 0
CTDB[CTDB == "unknown"] <- NA

data_reshaped <- acast(CTDB, specie_name~trait_name, value.var="value", fun.aggregate=change_CTDB, fill=NA_character_)
#traits contains max values for traits with more than 1 value 
traits <- data.frame(data_reshaped, stringsAsFactors=FALSE)
traits$molefam <- CTDB$family_molecules[match(rownames(traits), CTDB$specie_name)]
traits$molefam[traits$molefam == "Incertae sedis"] <- NA

traits$morphfam <- CTDB$family_morphology[match(rownames(traits), CTDB$specie_name)]
traits$master_species <- rownames(traits) 


veron <- read.csv("data/Veron20150603.csv", header=TRUE) #presence/absence

test<-veron$Revised.CTDB.name
test2<-traits$master_species
test[test %in% test2 == F] 



traitsveron <- merge(x=veron, y=traits, by.x="Revised.CTDB.name", by.y="master_species")


# Grab a vector of species names, as characters
sppList <- as.character(unique(traitsveron$Revised.CTDB.name))
latlon <- paste0(as.character(traitsveron[,'Latitude...original']),":",as.character(traitsveron[, 'Longitude...original'])) # paste lat and long together into a character string
uniquelatlon = unique(latlon) # find all unique lat-lon locations

dat <- matrix(data = 0, nrow = length(uniquelatlon), ncol = 313) # make an empty matrix where number of rows is number of unique locations and number of columns is for species (306) + lat/long (2) + envt var (5). 
dat <- as.data.frame(dat) # convert to a data frame
names(dat)[1:2] = c("Latitude","Longitude")
names(dat)[3:7] = names(as.data.frame(traitsveron))[7:11] # grab names of envt vars from newpredictors matrix
names(dat)[8:ncol(dat)] = sppList # rename the rest of the columns using the names in 'cols' from above

for (i in 1:length(uniquelatlon)){
  dat$Latitude[i] = as.numeric(strsplit(uniquelatlon[i],":")[[1]][1])
  dat$Longitude[i] = as.numeric(strsplit(uniquelatlon[i],":")[[1]][2])
}

for (i in 1:nrow(dat)){
  lat = dat$Latitude[i]
  lon = dat$Longitude[i]
  rows1 = which(lat == traitsveron$Latitude...original)
  rows2 = which(lon == traitsveron$Longitude...original)
  matches = intersect(rows1,rows2) # find rows that are the same in both groups
  dat[i,3:7] = traitsveron[matches[1],7:11] # grab columns 3:6 from the 1st row in 'matches', since all rows in matches should have the same data in columns 3:6 (they're all the same site)
  
  # Now go through the contents of 'matches', find out what species is listed in each row of ALAtraits contained within 'matches', and put a 1 into the correct species column in 'dat' to denote species presence at the current site. 
  for (j in 1:length(matches)){
    # matches contains a vector of row numbers that match the current lat lon
    spp = as.character(traitsveron$Revised.CTDB.name[matches[j]]) # go to row of ALA that is stored in matches[j], grab the scientific name stored there
    # find matching column name in dat and put a 1 in it to show presence of this species
    sppMatch = which(spp == sppList)
    # sppMatch contains a number, corresponding to the entry in sppList that matches the species name we just extracted from sppList. Column numbers in dat don't correspond directly to index number in sppList, because there are 6 extract columns of data at the start of dat. So offset the sppMatch index value by 6 to find the correct column number in dat that corresponds to the current species. 
    sppMatch = sppMatch + 7
    dat[i,sppMatch] = 1 # insert a 1 to show this species is present at this site. 
  } # end of for loop j
  
} # end of for loop i

enviro <- dat #each lat/long listed in rows and environmental variable values and species presence/absence for each site

#########################
#DATA
enviro["site"] <- NA # creates the new column named "site" filled with "NA". Sites are unique lat/longs

#?
enviro$site <- c(1:225) 

PA <- melt(enviro[c(8:314)],id.vars="site") #extract species' presence absence data

colnames(PA)[2:3] <- c("species","present")

model.data <- merge(PA, enviro[-c(8:313)], by='site',all.x=T) #merge presence/absences with environmental data (site is the last column, hence 314 and not 315)

#rm(PA)

data <- merge(model.data,traits,by.y='master_species',by.x='species',all.x=T) #add species traits

colnames(data) <- c("species", "site", "present", "lat", "long", "turb", "highSST", "lowSST", "meanSST", "wave.power", "CalcRt", "Col", "Cor", "GF", "GR", "polyps", "skel", "skel.micro", "fam", "fammorph")

model.data <- subset(data, wave.power<= 2)


#These species have corallite size and calcification and diverse ranges
model.datasp <- model.data[ which((model.data$species=='Seriatopora hystrix') |
                                    (model.data$species=='Acropora valida') | 
                                    (model.data$species=='Pocillopora damicornis') | 
                                    (model.data$species=='Acropora hyacinthus') | 
                                    (model.data$species=='Acropora solitaryensis') |
                                    (model.data$species=='Porites annae') | 
                                    (model.data$species=='Porites cylindrica') |
                                    (model.data$species=='Porites lobata') | 
                                    (model.data$species=='Acropora muricata') | 
                                    (model.data$species=='Acropora pulchra') | 
                                    (model.data$species=='Acropora glauca') | 
                                    (model.data$species=='Acropora longicyathus') | 
                                    (model.data$species=='Porites lutea') | 
                                    (model.data$species=='Isopora palifera') | 
                                    (model.data$species=='Montipora verrucosa') | 
                                    (model.data$species=='Porites australiensis') |
                                    (model.data$species=='Acropora horrida') | 
                                    (model.data$species=='Acropora humilis') |
                                    (model.data$species=='Acropora aspera') |
                                    (model.data$species=='Pavona clavus') |
                                    (model.data$species=='Leptoria phrygia') |
                                    (model.data$species=='Gardineroseris planulata') |
                                    (model.data$species=='Goniastrea retiformis') |
                                    (model.data$species=='Pavona varians') |
                                    (model.data$species=='Dipsastraea pallida') |
                                    (model.data$species=='Oulophyllia crispa')),]
######################DIFF SPECIES######################
#model.datasp <- model.data[ which(
(model.data$species=='Acropora gemmifera') | #Solitaries
  (model.data$species=='Acropora hyacinthus') | #Japan
  (model.data$species=='Acropora microclados') | #Solitaries
  (model.data$species=='Acropora monticulosa') | #Solitaries
  (model.data$species=='Pocillopora damicornis') | #AB-brooder
  (model.data$species=='Isopora cuneata') | #AB-brooder
  (model.data$species=='Micromussa lordhowensis') | #AB - more abundandant in subtropics
  (model.data$species=='Paragoniastrea australensis') | #AB - more abundandant in subtropics
  (model.data$species=='Turbinaria radicalis') | #AB
  (model.data$species=='Turbinaria mesenterina') | #AB
  (model.data$species=='Plesiastrea versipora') |
  (model.data$species=='Acropora humilis') |
  (model.data$species=='Dipsastraea favus') |
  (model.data$species=='Fungia fungites') |
  (model.data$species=='Acropora glauca') |
  (model.data$species=='Turbinaria reniformis') |
  (model.data$species=='Platygyra sinensis') |
  (model.data$species=='Porites stephensoni') |
  (model.data$species=='Stylophora pistillata') |
  (model.data$species=='Mycedium elephantotus') |
  
  (model.data$species=='Acanthastrea echinata') |
  (model.data$species=='Isopora brueggemanni') |
  (model.data$species=='Acropora digitifera') |
  (model.data$species=='Acropora donei') |
  (model.data$species=='Acropora robusta') |
  (model.data$species=='Acropora yongei') |
  (model.data$species=='Coscinaraea exesa')|
  (model.data$species=='Coscinaraea columna') |
  (model.data$species=='Dipsastraea pallida') |
  (model.data$species=='Psammocora contigua') |
  (model.data$species=='Porites vaughani')),]




#select species
#species, p/a, all envt vars
model.covs <- data.frame(present=model.datasp$present,
                         species=model.datasp$species, 
                         turbraw=model.datasp$turb,
                         highSSTraw=model.datasp$highSST, #SST of the warmest month
                         lowSSTraw=model.datasp$lowSST, #SST of the coldest, ice free month
                         meanSSTraw=model.datasp$meanSST,
                         waveraw=model.datasp$wave.power) #combined contribution to the bed shear stress by waves, tides, wind and density-driven circulation


par(mfrow=c(2,3))
plot(density(na.omit(model.covs$turbraw)),main="Turbidity raw")
plot(density(na.omit(model.covs$highSSTraw)),main="High SST raw")
plot(density(na.omit(model.covs$lowSSTraw)),main="Low SST raw")
plot(density(na.omit(model.covs$meanSSTraw)),main="mean SST raw")
plot(density(na.omit(model.covs$waveraw)),main="Wave power raw")


pairs(~turbraw+highSSTraw+lowSSTraw+meanSSTraw+ waveraw, data=model.covs,
      main="Scatterplot Matrix of raw environmental variables")
log.turb <- log(model.covs$turbraw)
log.highSST <- log(model.covs$highSSTraw)
log.lowSST <- log(model.covs$lowSSTraw)
log.meanSST <- log(model.covs$meanSSTraw)
log.wavepower <- log(model.covs$waveraw+1)

model.covs$turb <- ((log.turb-mean(log.turb, na.rm=TRUE))/(2*sd(log.turb, na.rm=TRUE))) 
model.covs$highSST <- ((log.highSST-mean(log.highSST, na.rm=TRUE))/(2*sd(log.highSST, na.rm=TRUE))) 

model.covs$lowSST <- ((log.lowSST-mean(log.lowSST, na.rm=TRUE))/(2*sd(log.lowSST, na.rm=TRUE)))

model.covs$meanSST <- ((log.meanSST-mean(log.meanSST, na.rm=TRUE))/(2*sd(log.meanSST, na.rm=TRUE)))

model.covs$wavepower <- ((log.wavepower -mean(log.wavepower , na.rm=TRUE))/(2*sd(log.wavepower , na.rm=TRUE)))



par(mfrow=c(2,3))
plot(density(na.omit(model.covs$turb)),main="Turbidity")
plot(density(na.omit(model.covs$highSST)),main="High SST")
plot(density(na.omit(model.covs$lowSST)),main="Low SST")
plot(density(na.omit(model.covs$meanSST)),main="mean SST")
plot(density(na.omit(model.covs$wavepower)),main="Wave power")

pairs(~turb+highSST+lowSST+ meanSST+wavepower,data=model.covs,
      main="Scatterplot Matrix of transformed environmental variables")


#ADD trait data
model.covs$calc_rate <- as.numeric(as.vector(model.datasp$CalcRt))
model.covs$col_size <- as.numeric(as.vector(model.datasp$Col))
model.covs$cor_size <- as.numeric(as.vector(model.datasp$Cor))
model.covs$gf <- factor(model.datasp$GF)
model.covs$gr_rate <- as.numeric(as.vector(model.datasp$GR))
model.covs$polyps <- as.numeric(as.vector(model.datasp$polyps))
model.covs$skel_dens <- as.numeric(as.vector(model.datasp$skel))
model.covs$skel_micro <- as.numeric(as.vector(model.datasp$skel.micro))
model.covs$fam <- factor(model.datasp$fam)

par(mfrow=c(1,1))
plot(density(na.omit(model.covs$calc_rate)),main="Calc rate")
plot(density(na.omit(model.covs$col_size)),main="colony sz")
plot(density(na.omit(model.covs$cor_size)),main="corallite")
plot(density(na.omit(model.covs$gr_rate)),main="gr rate")
plot(density(na.omit(model.covs$polyps)),main="polyps")
plot(density(na.omit(model.covs$skel_dens)),main="skel_dens")

model.covs$gf <- classify(model.covs$gf,
                              c(branching_closed="branching", 
branching_open="branching", 
columnar="massive", 
corymbose="digitate", 
digitate="digitate", 
encrusting="encrusting", 
encrusting_long_uprights="encrusting", 
hispidose="branching",
laminar="tabular",
massive="massive" , 
tables_or_plates="tabular"))
model.covs$gf <- factor(model.covs$gf, levels = c("branching", "digitate","encrusting","tabular", "massive"))
plot(model.covs$gf)



###TURBIDITY (and corallite size)
modT <- glmer(present ~ turb + (1 +turb|species), data= model.covs,family=binomial(link=logit))
modT
coef(modT)
#species distributions of those with negative turb coef should generally be found where water is less turbid
plot(coef(modT))

#Species with larger corallite should be in more turbid water, they're better at expelling sediments and handling the lower light due to turbid water
modT1 <- glmer(present ~ turb + cor_size:turb+(1 +turb|species), data= model.covs,family=binomial(link=logit))
modT1
coef(modT1)


###WAVE POWER (and growth form)
modW <- glmer(present ~ wavepower + (1 +wavepower|species), data= model.covs,family=binomial(link=logit))
modW
coef(modW)
#species with positive coef should be more where wave power is greater (generally offshore/outer reef, higher lat, etc)
plot(coef(modW))

#Massive/encrusting more robust to wave energy compared to digitate, tabular can go either way
modW1 <- glmer(present ~ wavepower + gf:wavepower+(1 +wavepower|species), data= model.covs,family=binomial(link=logit))
modW1
coef(modW1)
ggCaterpillar(ranef(modW1, postVar=TRUE))
dotplot(ranef(modW1, postVar=TRUE))$species[1]
dotplot(ranef(modW1, postVar=TRUE))$species[2]
dotchart(fixef(modW1))
GFwave  <-  function() {
  dotchart(fixef(modW1))
}
to.pdf(GFwave(), 'output/GFwave .pdf', width=6, height=6)


###LOW TEMP (SST of the coldest ice free month)
modLT <- glmer(present ~ lowSST + (1 +lowSST|species), data= model.covs,family=binomial(link=logit))
modLT
coef(modLT)
#species with negative coef should be distributed in colder water (high lat)
plot(coef(modLT))

#The tolerance of test corals to low water temperature was closely related to their morphologies, with the branching corals being the most vulnerable to bleaching and death by separating the symbiotic polyps from their skeletons. Li, S., Yu, K. F., Shi, Q., Chen, T. R. and Zhao, M. X. (2009) Low water temperature tolerance and responding mode of scleractinian corals in Sanya Bay. Yingyong Shengtai Xuebao [Chinese Journal of Applied Ecology], 20 9: 2289-2295. 
modLT1 <- glmer(present ~ lowSST + gf:lowSST+(1 +lowSST|species), data= model.covs,family=binomial(link=logit))
modLT1
coef(modLT1)
ggCaterpillar(ranef(modLT1, postVar=TRUE))
dotplot(ranef(modLT1, postVar=TRUE))$species[1]
dotplot(ranef(modLT1, postVar=TRUE))$species[2]
dotchart(fixef(modLT1))

GFtemp  <-  function() {
  dotchart(fixef(modLT1))
}
to.pdf(GFtemp(), 'output/GFtemp .pdf', width=6, height=6)















mod <- glmer(present ~ turb+ MARsst+ var + wavepower+ 
                  gr_rate:turb + gr_rate:MARsst  +gr_rate:var +gr_rate:wavepower +
                  cor_size:turb +cor_size:MARsst  +cor_size:var  +cor_size:wavepower  +
                  gf:turb + gf:MARsst + gf:var + gf:wavepower+
                  (1 + turb+ MARsst+ var + wavepower|species),
                data= model.covs,family=binomial(link=logit))

ggCaterpillar(ranef(modW1, postVar=TRUE))
dotplot(ranef(modW1, postVar=TRUE))$species[1]
dotplot(ranef(modW1, postVar=TRUE))$species[2]
dotplot(ranef(modW1, postVar=TRUE))$species[3]
dotplot(ranef(mod, postVar=TRUE))$species[4]
dotplot(ranef(mod, postVar=TRUE))$species[5]
dotchart(fixef(modW1))
dotchart(coef(modW1))


mod.wp <- glmer(present ~ wavepower+
               gr_rate:wavepower +
               cor_size:wavepower  +
               col_size:wavepower +
               calc_rate:wavepower+
               gf:wavepower+
               (1 + wavepower|species),
             data= model.covs,family=binomial(link=logit))
dotchart(fixef(mod.wp)[2:8])

mod.turbSST <- glmer(present ~ turb+ MARsst+
               gr_rate:turb + gr_rate:MARsst+  
               cor_size:turb +cor_size:MARsst+  
               col_size:turb + col_size:MARsst+  
               calc_rate:turb + calc_rate:MARsst+
               gf:turb + gf:MARsst +
               (1 + turb+ MARsst|species),
             data= model.covs,family=binomial(link=logit))
ggCaterpillar(ranef(mod.turbSST, postVar=TRUE))

dotplot(ranef(mod.turbSST, postVar=TRUE))$species[1]
dotplot(ranef(mod.turbSST, postVar=TRUE))$species[2]
dotplot(ranef(mod.turbSST, postVar=TRUE))$species[3]
dotchart(fixef(mod.turbSST)[2:17])
dotchart(coef(mod), xlim=c(-4.5,4.5))
dotchart(coefplot2(mod)[4:9], xlim=c(-4.5,4.5))














FULLmodel <- lmer(present ~ SST + var +  Chl + TSM + PAR + skw +
                    abundance:SST + abundance:var +  abundance:Chl+abundance:TSM + abundance:PAR + abundance:skw +
                    Ldepth:SST + Ldepth:var +  LdepthChl+Ldepth:TSM +Ldepth:PAR + Ldepth:skw +
                    Rdepth:SST + Rdepth:var +  Rdepth:Chl+Rdepth:TSM +Rdepth:PAR + Rdepth:skw +
                    reef:SST + reef:var +  reef:Chl + reef:TSM +reef:PAR + reef:skw +
                    growth:SST + growth:var +  growth:Chl + growth:TSM +growth:PAR + growth:skw +
                    rangesize:SST + rangesize:var +  rangesize:Chl + rangesize:TSM +rangesize:PAR + rangesize:skw +
                    broodspawn:SST + broodspawn:var +  broodspawn:Chl + broodspawn:TSM +broodspawn:PAR + broodspawn:skw +
                    hermgono:SST + hermgono:var +  hermgono:Chl + hermgono:TSM +hermgono:PAR + hermgono:skw +
                    subattach:SST + subattach:var +  subattach:Chl + subattach:TSM +subattach:PAR + subattach:skw +
                    claritypref:SST + claritypref:var +  claritypref:Chl + claritypref:TSM +claritypref:PAR + claritypref:skw +
                    exppref:SST + exppref:var +  exppref:Chl + exppref:TSM +exppref:PAR + exppref:skw +
                    zoox:SST + zoox:var +  zoox:Chl + zoox:TSM +zoox:PAR + zoox:skw +
                    wave_exp:SST + wave_exp:var +  wave_exp:Chl + wave_exp:TSM +wave_exp:PAR + wave_exp:skw +
                    water_clarity:SST + water_clarity:var +  water_clarity:Chl + water_clarity:TSM +water_clarity:PAR + water_clarity:skw +
                    depthrange:SST + depthrange:var +  depthrange:Chl + depthrange:TSM +depthrange:PAR + depthrange:skw +
                    (1 + SST + var +  Chl + TSM + PAR + skw |species),
                  data= model.covsB,family=binomial(link=logit),
                  control=list(maxIter=10000)) 

model.temp <- lmer(present ~ SST + var + skwtr +
                     abundance:SST + abundance:var + abundance:skwtr + #correlation of response to trait
                     Ldepth:SST + Ldepth:var +  Ldepth:skwtr +
                     reef:SST + reef:var +  reef:skwtr +
                     broodspawn:SST + broodspawn:var +  broodspawn:skwtr +
                     hermgono:SST + hermgono:var +  hermgono:skwtr +
                     subattach:SST + subattach:var +  subattach:skwtr +
                     (1 + SST + var +  skwtr |species), #random slope of env within species with correlated intercept, #response to enviro variables allowed to vary by species
                   data= model.covsB,family=binomial(link=logit),
                   control=list(maxIter=10000)) 
#67% is explained by between species variation
#33% is explained by within species variation

(1 + SST + var +  skwtr |species)


modelVCP<-lmer(Feeding ~ (1|Ind_ID), data)
model.tempbasic <- lmer(present ~ (1 + SST + var +  skwtr |species),
                        data= model.covsB,family=binomial(link=logit),
                        control=list(maxIter=10000))

model.temp1 <- lmer(present ~ SST + var + skwtr +
                      abundance:SST + abundance:var + abundance:skwtr +
                      Ldepth:SST + Ldepth:var +  Ldepth:skwtr +
                      reef:SST + reef:var +  reef:skwtr +
                      hermgono:SST + hermgono:var +  hermgono:skwtr +
                      subattach:SST + subattach:var +  subattach:skwtr +
                      (1 + SST + var +  skwtr |species),
                    data= model.covsB,family=binomial(link=logit),
                    control=list(maxIter=10000)) 

varCorr1

model.red <- lmer(present ~ SST + var +  TSM + skw +
                    
                    abundance:SST + abundance:var +  abundance:TSM + abundance:skw +
                    Ldepth:SST + Ldepth:var + Ldepth:TSM +Ldepth:skw +
                    reef:SST + reef:var +  reef:TSM +reef:skw +
                    broodspawn:SST + broodspawn:var +  broodspawn:TSM +broodspawn:skw +
                    hermgono:SST + hermgono:var +  hermgono:TSM + hermgono:skw +
                    subattach:SST + subattach:var +  subattach:TSM +subattach:skw +
                    claritypref:SST + claritypref:var +  claritypref:TSM +claritypref:skw +
                    exppref:SST + exppref:var +  exppref:TSM +exppref:skw +
                    zoox:SST + zoox:var +  zoox:TSM +zoox:skw +
                    (1 + SST + var +  TSM + skw |species),
                  data= model.covsB,family=binomial(link=logit),
                  control=list(maxIter=10000)) 

model.red <- lmer(present ~ SST + var +  TSM + 
                    
                    
                    Ldepth:SST + Ldepth:var + Ldepth:TSM +
                    claritypref:SST + claritypref:var +  claritypref:TSM +
                    exppref:SST + exppref:var +  exppref:TSM +
                    zoox:SST + zoox:var +  zoox:TSM +
                    (1 + SST + var +  TSM|species),
                  data= model.covsB,family=binomial(link=logit),
                  control=list(maxIter=10000)) 




model.sstvar <- lmer(present ~ SST + var + 
                       Ldepth:SST + Ldepth:var +  
                       reef:SST + reef:var +  
                       broodspawn:SST + broodspawn:var + 
                       hermgono:SST + hermgono:var +  
                       subattach:SST + subattach:var +  
                       (1 + SST + var |species),
                     data= model.covsB,family=binomial(link=logit),
                     control=list(maxIter=10000)) 
#start2:53
ggCaterpillar(ranef(model3, postVar=TRUE))
#OR
qqmath(ranef(model2, postVar = TRUE), strip = FALSE)$species
randoms<-ranef(model2, postVar = TRUE)
qq <- attr(ranef(model2, postVar = TRUE)[[1]], "postVar")
rand.interc<-randoms$species
df<-data.frame(Intercepts=randoms$species[,1],
               sd.interc=2*sqrt(qq[,,1:length(qq)]),#####WHY ISNT THIS WORKING???
               lev.names=rownames(rand.interc))
df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])
library(ggplot2)
p <- ggplot(df,aes(lev.names,Intercepts,shape=lev.names))
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Levels") + ylab("")

#Final adjustments of plot
p <- p + theme(axis.text.x=element_text(size=rel(1.2)),
               axis.title.x=element_text(size=rel(1.3)),
               axis.text.y=element_text(size=rel(1.2)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())

#To put levels on y axis you just need to use coord_flip()
p <- p+ coord_flip()
print(p)



dotplot(ranef(mod.wp, postVar=TRUE))$species[1]
dotplot(ranef(mod.wp, postVar=TRUE))$species[2]
dotplot(ranef(mod.wp, postVar=TRUE))$species[3]
dotplot(ranef(mod, postVar=TRUE))$species[4]
dotchart(fixef(mod)[2:6], xlim=c(-1,1))
dotchart(coef(mod), xlim=c(-4.5,4.5))
dotchart(coefplot2(mod)[4:9], xlim=c(-4.5,4.5))
f<- fixef(model2)


residuals <- resid(model)
summary(residuals)
hist(residuals)


fixef$lower <- fixef$Est - fixef$SE
fixef$upper <- fixef$Est + fixef$SE
dotplot(Trait ~ Est| Env, data = fixef, 
        aspect = 1.5,
        xlab = "Trait Coefficient",
        xlim = c(-4, 4),
        scales = list(cex = .6),
        
        panel = function (x, y) {
          panel.xyplot(x, y, pch = 16, col = "black")
          panel.segments(fixef$lower, as.numeric(y),
                         fixef$upper, as.numeric(y), lty = 1, col = "black")})

dotplot(var.labels ~ sample.means, data = new.data,
        aspect = 1.5,
        xlim = c(3.3, 4.3),
        xlab = "Mean importance rating",
        panel = function (x, y) {
          panel.xyplot(x, y, pch = 16, col = "black")
          panel.segments(fixef$lower, as.numeric(y),
                         fixef$upper, as.numeric(y), lty = 1, col = "black")} )

library(ggplot2) 
pd <- position_dodge(width=0.4) 
g0 <- ggplot(model2,aes(x=age,y=pred,colour=Sex))+ 
  geom_point(position=pd) 
g0 + geom_linerange(aes(ymin=pred-2*SE,ymax=pred+2*SE), position=pd)

## prediction intervals 
g0 + geom_linerange(aes(ymin=pred-2*SE2,ymax=pred+2*SE2), position=pd)


dotplot(site ~ yield | variety, data=barley,
        >        groups=year, yield2=barley$yield2,
        >        col=c("gray", "black"),
        >        panel = panel.superpose,
        >        panel.groups = function(x, y, subscripts, yield2, col, ...) {
          >            panel.xyplot(x, y, col = col, ...)
          >            panel.segments(x, y, yield2[subscripts], y, col = col)
          >        })

mypanel = function(x,y,...){ 
  agg<-aggregate(y,list(x),function(x)c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))) 
  mns = agg[,2][,1] 
  sds = agg[,2][,2] 
  for(i in 1:nrow(agg))llines(c(i-.1,i+.1),rep(mns[i],2),lwd=3) 
  for(i in 1:nrow(agg)){llines(c(i-.1,i+.1),rep(mns[i] + 1.96 * sds[i],2)); 
                        llines(c(i-.1,i+.1),rep(mns[i] - 1.96 * sds[i],2))} 
  panel.dotplot(x,y,...) 
} 

dotplot(Env ~ Est | Trait, groups=Trait, data=fixef)




# JOSH

modelselspACRS2

fixef(model2)

sst <- seq(min(model.covs[,"SST"]), max(model.covs[,"SST"]), 0.05) #x variable
tsm <- seq(min(model.covs[,"TSM"]), max(model.covs[,"TSM"]), 0.05) #y variable

#1narrowdepth, clear
dd <- 1
cc <- 1
tt <- 0

#2narrowdepth, turbid
dd <- 1
cc <- 0
tt <- 1
#3narrowdepth, both
dd <- 1
cc <- 0
tt <- 0

#4broaddepth, clear
dd <- 0
cc <- 1
tt <- 0

#5broaddepth, turbid
dd <- 0
cc <- 0
tt <- 1

#6broaddepth, both
dd <- 0
cc <- 0
tt <- 0


mat <- matrix(NA, length(sst), length(tsm))
bf <- function(x) {exp(x) / (1 + exp(x))}

for (i in 1:length(sst)) {
  ssti <- sst[i]
  pred <- -5.4048 + 3.3705 * tsm + (-1.5140 * ssti) + (-3.7082 * tsm) * dd + (-0.3810 * ssti) * dd + (-2.2030 * tsm) * cc + (1.6040 * ssti) * cc + (3.6607 * tsm) * tt + (2.7954 * ssti) * tt
  mat[i, ] <- bf(pred)
}


colfunc <- colorRampPalette(c("yellow", "red"))
colfunc(20)
contour(mat)
image(mat, col = colfunc(20))

#predict(modelselspACRS2, data.frame(TSM=1, var=1, depth="s", water_clarity="clear"))

#Formula: present ~ TSM + var + depth:TSM + depth:var + water_clarity:TSM +      water_clarity:var + (1 + TSM + var | species) 



plotCI(x = means.depth, uiw = stderr.depth,lty = 2, xaxt ="n", ylim = c(-4.5,4.5), xlim = c(.5,2.5), gap = 0, pch=16, ylab="Mean ", xlab="Trait Coefficients", col="blue",  cex=1.2, main = "Traits Modulate responses to the environment")
Error: could not find function "plotCI"

plotCI(x = means.clear, uiw = stderr.clear, lty = 1, col="coral", xaxt ="n", ylim = c(-4.5,4.5),  gap = 0, cex=1.2, add = TRUE)

plotCI(x = means.Turbid, uiw = stderr.Turbid, lty = 1, col="black",xaxt ="n", ylim = c(-4.5,4.5),  gap = 0, add = TRUE)

#Draw the x-axis (omitted above)

axis(side = 1, at = 1:2, cex=1.2, labels = names(stderr.clear), cex = 0.7)

Error in as.graphicsAnnot(labels) : object 'stderr.clear' not found
# Add legend for male and female participants

legend(1.5,1, legend=c("clear","depth", "turbid"),lty=1:2)
Error in strwidth(legend, units = "user", cex = cex, font = text.font) : 
  
  means.depth <- c(-3.7082, -0.3810)
stderr.depth <- c(1.0115, 0.8690)

names(means.depth) <- c("TSM","SD")
names(stderr.depth) <- c("TSM","SD")

means.clear <- c(-2.2030, 1.6040)
stderr.clear <- c(0.9499, 1.0497)

names(means.clear) <- c("TSM","SD")
names(stderr.clear) <- c("TSM","SD")

means.Turbid<- c(3.6607, 2.7954)
stderr.Turbid <- c(1.4962, 2.2129)

names(means.Turbid) <- c("TSM","SD")
names(stderr.Turbid) <- c("TSM","SD")






