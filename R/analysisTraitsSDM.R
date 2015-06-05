rm(list=ls())
library(reshape2) #melt 
library(lme4) 
library(biomod2) 
library(ncdf4)

source('R/functions.R')

#traits
#corallite size / turbidity? sediment rejection (from Rachello dolmen)

#lower depth
#depth range
#growth rate

#symbiodinium diversity or presence of a particular clade
#symbiodinium density
#range?
#colony max size

 
CTDB <- read.csv("data/ctdb_20150525.csv", as.is=TRUE)
CTDB[is.na(CTDB) == T] = 0
CTDB[CTDB == "unknown"] <- NA

data_reshaped <- acast(CTDB, specie_name~trait_name, value.var="value", fun.aggregate=change_CTDB, fill=NA_character_)
#traits contains averages for traits with more than 1 value 
traits <- data.frame(data_reshaped, stringsAsFactors=FALSE)
traits$molefam <- CTDB$family_molecules[match(rownames(traits), CTDB$specie_name)]
traits$molefam[traits$molefam == "Incertae sedis"] <- NA

traits$morphfam <- CTDB$family_morphology[match(rownames(traits), CTDB$specie_name)]
traits$master_species <- rownames(traits) 


veron <- read.csv("data/Veron20150525.csv", header=TRUE) #presence/absence

test<-veron$Revised.CTDB.name
test2<-traits$master_species
test[test %in% test2 == F] 



traitsveron <- merge(x=veron, y=traits, by.x="Revised.CTDB.name", by.y="master_species")

PAtraits <- traitsveron[c("Revised.CTDB.name", "Latitude...original", "Longitude...original","Calcification.rate","Colony.maximum.diameter","Corallite.width.maximum","Growth.form.typical","Growth.rate" ,"Polyps.per.area","Skeletal.density","Skeletal.micro.density","molefam","morphfam")] #species, lat, long, traits

COOR <- as.matrix(PAtraits[,3:2])#long, lat
layers <- list.files(path="/Users/tmizerek/Documents/Biomod/BIOMOD2/Data_tif", pattern='tif', full.names=TRUE )
predictors <- stack(layers[c(54,48, 1, 40, 60, 28)]) 
#(stdev(temp variability), SST_p50, CHL, skewness, TSM, PAR) datum=WGS84
newpredictors <- extract(predictors, COOR)
newpredictors <- cbind(COOR,newpredictors)
newpredictors <- as.data.frame(newpredictors)

#Sea surface temperature derived variables were obtained from the second version of the coral reef temperature anomaly database (CoRTAD) [23]. This database contains global SST and related thermal stress metrics at an approximately 4-km resolution weekly from 1982 through 2008, derived from measurements from the Advanced Very High Resolution Radiometer onboard NOAA suite of polar orbiting satellites.


# Grab a vector of species names, as characters
sppList <- as.character(unique(PAtraits$Revised.CTDB.name))
latlon <- paste0(as.character(newpredictors[,'Latitude...original']),":",as.character(newpredictors[, 'Longitude...original'])) # paste lat and long together into a character string
uniquelatlon = unique(latlon) # find all unique lat-lon locations

dat <- matrix(data = 0, nrow = length(uniquelatlon), ncol = 315) # make an empty matrix where number of rows is number of unique locations and number of columns is for species (307) + lat/long (2) + envt var (6). 
dat <- as.data.frame(dat) # convert to a data frame
names(dat)[1:2] = c("Latitude","Longitude")
names(dat)[3:8] = names(as.data.frame(newpredictors))[3:8] # grab names of envt vars from newpredictors matrix
names(dat)[9:ncol(dat)] = sppList # rename the rest of the columns using the names in 'cols' from above

for (i in 1:length(uniquelatlon)){
  dat$Latitude[i] = as.numeric(strsplit(uniquelatlon[i],":")[[1]][1])
  dat$Longitude[i] = as.numeric(strsplit(uniquelatlon[i],":")[[1]][2])
}

for (i in 1:nrow(dat)){
  lat = dat$Latitude[i]
  lon = dat$Longitude[i]
  rows1 = which(lat == newpredictors$Latitude...original)
  rows2 = which(lon == newpredictors$Longitude...original)
  matches = intersect(rows1,rows2) # find rows that are the same in both groups
  dat[i,3:8] = newpredictors[matches[1],3:8] # grab columns 3:6 from the 1st row in 'matches', since all rows in matches should have the same data in columns 3:6 (they're all the same site)
  
  # Now go through the contents of 'matches', find out what species is listed in each row of ALAtraits contained within 'matches', and put a 1 into the correct species column in 'dat' to denote species presence at the current site. 
  for (j in 1:length(matches)){
    # matches contains a vector of row numbers that match the current lat lon
    spp = as.character(PAtraits$Revised.CTDB.name[matches[j]]) # go to row of ALA that is stored in matches[j], grab the scientific name stored there
    # find matching column name in dat and put a 1 in it to show presence of this species
    sppMatch = which(spp == sppList)
    # sppMatch contains a number, corresponding to the entry in sppList that matches the species name we just extracted from sppList. Column numbers in dat don't correspond directly to index number in sppList, because there are 6 extract columns of data at the start of dat. So offset the sppMatch index value by 6 to find the correct column number in dat that corresponds to the current species. 
    sppMatch = sppMatch + 8
    dat[i,sppMatch] = 1 # insert a 1 to show this species is present at this site. 
  } # end of for loop j
  
} # end of for loop i

enviro <- dat #each lat/long listed in rows and environmental variable values and species presence/absence for each plot across columns

######################### POLLOCK SCRIPT ###################################
#DATA
enviro["site"] <- NA # creates the new column named "site" filled with "NA". Sites are unique lat/longs

#?
enviro$site <- c(1:225) 

PA <- melt(enviro[c(9:316)],id.vars="site") #extract species' presence absence data

colnames(PA)[2:3] <- c("species","present")

model.data <- merge(PA, enviro[-c(9:315)], by='site',all.x=T) #merge presence/absences with environmental data (site is the last column, hence 315 and not 316)

#rm(PA)

model.data <- merge(model.data,traits,by.y='master_species',by.x='species',all.x=T) #add species traits

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
                          SSTraw=model.datasp$SST_p50, 
                          varraw=model.datasp$stdev,
                          skewnessraw=model.datasp$skewnes, 
                          Chlraw=model.datasp$Chl_50p, 
                          PARraw=model.datasp$PAR_50P, 
                          TSMraw=model.datasp$TSM_50P) 



par(mfrow=c(2,3))
plot(density(na.omit(model.covs$SST)),main="SST")
plot(density(na.omit(model.covs$TSMraw)),main="TSM raw")
plot(density(na.omit(model.covs$Chlraw)),main="Chl raw")
plot(density(na.omit(model.covs$PARraw)),main="PAR raw")
plot(density(na.omit(model.covs$skewnessraw)),main="skewness raw")
plot(density(na.omit(model.covs$varraw)),main="var raw")

pairs(~SST+turb+chl+ salmin + salmax, data=model.data,
      main="Scatterplot Matrix of raw environmental variables")
log.SST <- log(model.data$SST)
log.turb <- log(model.covs$TSMraw)
log.chl <- log(model.data$chl)
log.salmax <- log(model.data$salmax)
log.salmin <- log(model.data$salmin)


model.data$SST <- ((log.SST-mean(log.SST, na.rm=TRUE))/(sd(log.SST, na.rm=TRUE))) 

model.covs$turb <- ((log.turb-mean(log.turb, na.rm=TRUE))/(sd(log.turb, na.rm=TRUE)))

model.data$chl <- ((log.chl-mean(log.chl, na.rm=TRUE))/(sd(log.chl, na.rm=TRUE)))

model.data$salmax <- ((log.salmax -mean(log.salmax , na.rm=TRUE))/(sd(log.salmax , na.rm=TRUE)))

model.data$salmin <- ((log.salmin-mean(log.salmin, na.rm=TRUE))/(sd(log.salmin, na.rm=TRUE)))


pairs(~SST+turb+chl+ salmax + salmin,data=model.data,
      main="Scatterplot Matrix of transformed environmental variables")


pairs(~SSTraw+varraw+skewnessraw+Chlraw+ PARraw + TSMraw, data=model.covs,
      main="Scatterplot Matrix of raw environmental variables")
######################
log.SST <- log(model.covs$SSTraw)
log.var <- log(model.covs$varraw)
log.Chl <- log(model.covs$Chlraw)
log.PAR <- log(model.covs$PARraw)
log.TSM <- log(model.covs$TSMraw)
log.skw <- log(model.covs$skewnessraw)


# Rescale variables as necessary (e.g. SST logged, centered, scaled)


model.covs$SST <- ((log.SST-mean(log.SST, na.rm=TRUE))/(sd(log.SST, na.rm=TRUE))) 

model.covs$var <- ((log.var-mean(log.var, na.rm=TRUE))/(sd(log.var, na.rm=TRUE)))

model.covs$Chl <- ((log.Chl-mean(log.Chl, na.rm=TRUE))/(sd(log.Chl, na.rm=TRUE)))

model.covs$PAR <- ((log.PAR-mean(log.PAR, na.rm=TRUE))/(sd(log.PAR, na.rm=TRUE)))

model.covs$TSM <- ((log.TSM-mean(log.TSM, na.rm=TRUE))/(sd(log.TSM, na.rm=TRUE)))

model.covs$skw <- model.covs$skewnessraw
model.covs$skwtr <- ((log.skw-mean(log.skw, na.rm=TRUE))/(sd(log.skw, na.rm=TRUE)))

pairs(~SST+var+Chl+ PAR + TSM + skw,data=model.covs,
      main="Scatterplot Matrix of transformed environmental variables")


#


#ADD trait data
model.covs$calc_rate <- as.numeric(as.vector(model.datasp$Calcification.rate))
model.covs$col_size <- as.numeric(as.vector(model.datasp$Colony.maximum.diameter))
model.covs$cor_size <- as.numeric(as.vector(model.datasp$Corallite.width.maximum))
model.covs$gf <- factor(model.datasp$Growth.form.typical)
model.covs$gr_rate <- as.numeric(as.vector(model.datasp$Growth.rate))
model.covs$polyps <- as.numeric(as.vector(model.datasp$Polyps.per.area))
model.covs$skel_dens <- as.numeric(as.vector(model.datasp$Skeletal.density))
model.covs$skel_micro <- as.numeric(as.vector(model.datasp$Skeletal.micro.density))
model.covs$fam <- factor(model.datasp$molefam)

plot(model.covs$cor_size)



model <- glmer(present ~turb + (1 +  turb |species),
                  data= model.covs,family=binomial(link=logit)) 

model1 <- glmer(present ~SST + skel_dens:SST + 
                 (1 +  SST |species),
               data= model.covsB,family=binomial(link=logit)) 

model2 <- glmer(present ~ TSM+ SST+ 
                  skel_dens:TSM + skel_dens:SST + 
                  cor_size:TSM +cor_size:SST  +
                          (1 + TSM+SST|species),
                        data= model.covsB,family=binomial(link=logit),
                        control=list(maxIter=10000)) 

model3 <- glmer(present ~ TSM+ SST+ 
                  gr_rate:TSM + gr_rate:SST + 
                  cor_size:TSM +cor_size:SST  +
                  (1 + TSM+SST|species),
                data= model.covsB,family=binomial(link=logit)) 

model5 <- glmer(present ~ turb+ SST+ chl + salmin + salmax+
                  gr_rate:turb + gr_rate:SST + gr_rate:chl +gr_rate:salmin +gr_rate:salmax +
                  cor:turb +cor:SST  +cor:chl +cor:salmin  +cor:salmax  +
                  col:turb + col:SST + col:chl +col:salmin +col:salmax +
                  (1 + turb+ SST+ chl + salmin + salmax|species),
                data= model.data,family=binomial(link=logit))

model5 <- glmer(present ~ turb+ SST+ chl + salmin + salmax+
                  gr_rate:turb + gr_rate:SST + gr_rate:chl +gr_rate:salmin +gr_rate:salmax +
                  cor:turb +cor:SST  +cor:chl +cor:salmin  +cor:salmax  +
                  col:turb + col:SST + col:chl +col:salmin +col:salmax +
                  (1 + turb+ SST+ chl + salmin + salmax|species),
                data= model.data,family=binomial(link=logit))

col_size
ggCaterpillar(ranef(model4, postVar=TRUE))
















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
              


dotplot(ranef(model4, postVar=TRUE))$species[1]
dotplot(ranef(model4, postVar=TRUE))$species[2]
dotplot(ranef(model4, postVar=TRUE))$species[3]
dotplot(ranef(model4, postVar=TRUE))$species[4]
dotchart(fixef(model4)[2:6], xlim=c(-1,1))
dotchart(coef(model4), xlim=c(-4.5,4.5))
dotchart(coefplot2(model2)[4:9], xlim=c(-4.5,4.5))
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





