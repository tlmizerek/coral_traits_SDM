rm(list=ls())
library(reshape2) #melt 
library(lme4) 
library(biomod2) 
library(ncdf4)

source('R/functions.R')

#traits
#corallite size / turbidity, sediment rejection increased with larger corallite size
#growth form / wave power, robust growth forms are more resilient to wave energy
#growth rate / wave power, if you grow fast, you can be more stable in high wave energy
#calcification or growth rate / SST, higher growth rate does better in higher temp??
#sediment rejection less to do with branching corals because not as much sediment settles on branches

 
CTDB <- read.csv("data/ctdb_20150630.csv", as.is=TRUE)
CTDB[is.na(CTDB) == T] = 0
CTDB[CTDB == "unknown"] <- NA

data_reshaped <- acast(CTDB, specie_name~trait_name, value.var="value", fun.aggregate=change_CTDB, fill=NA_character_)
#traits contains averages for traits with more than 1 value 
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

PAtraits <- traitsveron[c("Revised.CTDB.name", "Latitude...original", "Longitude...original","Calcification.rate","Colony.maximum.diameter","Corallite.width.maximum","Growth.form.typical","Growth.rate" ,"Polyps.per.area","Depth.lower","Depth.upper" ,"Skeletal.density","Skeletal.micro.density","molefam","morphfam")] #species, lat, long, traits

#COOR <- as.matrix(PAtraits[,3:2])#long, lat
#layers <- list.files(path="/Users/tmizerek/Documents/Biomod/BIOMOD2/Data_tif", pattern='tif$', full.names=TRUE )
# predictors <- stack(layers[c(168, 37, 48,49, 60,61)]) #42 is low temperature??
# #(stdev(temp variability), SST_p50, CHL, skewness, TSM, PAR) datum=WGS84
# plot(predictors)
# newpredictors <- extract(predictors, COOR)
# newpredictors <- cbind(COOR,newpredictors)
# newpredictors <- as.data.frame(newpredictors)
# 
# library(raster)
# # get a vector of file names of the rasters with extent A. This code will return a vector of all .tif files in the given path. Change tif to whatever the extension is
# ff1 <- list.files('c:/path/to/directory/containing/the/rasters', patt='\\.tif$') # change the first argument to the path containing the first set of rasters
# # same for the rasters with extent B
# ff2 <- list.files('c:/path/to/directory/containing/the/rasters', patt='\\.tif$') # change the first argument to the path containing the second set of rasters
# # check the vectors to make sure they contain the file names you expected:
# ff1
# ff2
# 
# # stack them
# s1 <- stack(ff1)
# s2 <- stack(ff2)
# 
# layersnew <- layers[,c(28:167)]
# layersold <- layers[c(1:27, 168:203)]

### 
library(sp)
library(raster)
layers <- list.files(path="/Users/tmizerek/Documents/Biomod/BIOMOD2/Data_tif", pattern='tif$', full.names=TRUE)
layers1 <- layers[12:151]
layers2 <- layers[-(12:151)]
s1 <- stack(layers1)
s2 <- stack(layers2)
xy <- setNames(PAtraits[,3:2], c('lon', 'lat'))
d1 <- extract(s1, xy)
d2 <- extract(s2, xy)
ss <- cbind(xy, d1, d2)

# months <- format(as.Date(names(s1), format='par_%Y_%m_%d'), '%B %Y')
# winter_ind_split <- split(names(s1), list(format(as.Date(names(s1), format='par_%Y_%m_%d'), '%Y'),
#                       grepl('June|July', format(as.Date(names(s1), format='par_%Y_%m_%d'), '%B'))))
# winter_ind_split <- winter_ind_split[grep('TRUE', names(winter_ind_split))]
# winter_ind_split <- winter_ind_split[sapply(winter_ind_split, length) > 0]
# winter_dat <- setNames(lapply(seq_along(winter_ind_split), function(i) {
#   nm <- names(winter_ind_split)[i]
#   d <- ss[, winter_ind_split[[i]]]
#   rowMeans(d)
# }), paste0('par_junejuly_', sub('\\.TRUE', '', names(winter_ind_split))))
# ss <- cbind.data.frame(dat, winter_dat)
# ss <- cbind.data.frame(dat, par_all_years_mean_junejuly=rowMeans(as.data.frame(winter_dat)))

## Or just this one line...
#ss$winter_par <- rowMeans(ss[, grep('_06_|_07_', names(s1), val=T)])
## if only want a subset of years
ss$winter_par_early <- rowMeans(ss[, grep('9._06_|9._07_|2000_06_|2000_07_|2001_06_|2001_07_|2002_06_|2002_07_', names(s1), val=TRUE)]) #5 year average

ss<- ss[,c(1,2,143,160, 161, 166, 169,171)]
#Sea surface temperature derived variables were obtained from the second version of the coral reef temperature anomaly database (CoRTAD). This database contains global SST and related thermal stress metrics at an approximately 4-km resolution weekly from 1982 through 2008, derived from measurements from the Advanced Very High Resolution Radiometer onboard NOAA suite of polar orbiting satellites.


# Grab a vector of species names, as characters
sppList <- as.character(unique(PAtraits$Revised.CTDB.name))
latlon <- paste0(as.character(ss[,'lat']),":",as.character(ss[, 'lon'])) # paste lat and long together into a character string
uniquelatlon = unique(latlon) # find all unique lat-lon locations

dat <- matrix(data = 0, nrow = length(uniquelatlon), ncol = 314) # make an empty matrix where number of rows is number of unique locations and number of columns is for species (306) + lat/long (2) + envt var (6). 
dat <- as.data.frame(dat) # convert to a data frame
names(dat)[1:2] = c("Latitude","Longitude")
names(dat)[3:8] = names(as.data.frame(ss))[c(3:8)] # grab names of envt vars from newpredictors matrix
names(dat)[9:ncol(dat)] = sppList # rename the rest of the columns using the names in 'cols' from above

for (i in 1:length(uniquelatlon)){
  dat$Latitude[i] = as.numeric(strsplit(uniquelatlon[i],":")[[1]][1])
  dat$Longitude[i] = as.numeric(strsplit(uniquelatlon[i],":")[[1]][2])
}

for (i in 1:nrow(dat)){
  lat = dat$Latitude[i]
  lon = dat$Longitude[i]
  rows1 = which(lat == ss$lat)
  rows2 = which(lon == ss$lon)
  matches = intersect(rows1,rows2) # find rows that are the same in both groups
  dat[i,3:8] = ss[matches[1],3:8] # grab columns 3:6 from the 1st row in 'matches', since all rows in matches should have the same data in columns 3:6 (they're all the same site)
  
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

#####
#DATA
enviro["site"] <- NA # creates the new column named "site" filled with "NA". Sites are unique lat/longs

enviro$site <- c(1:225) 

PA <- melt(enviro[c(9:315)],id.vars="site") #extract species' presence absence data

colnames(PA)[2:3] <- c("species","present")

model.data <- merge(PA, enviro[-c(9:314)], by='site',all.x=T) #merge presence/absences with environmental data (site is the last column, hence 315 and not 316)

#rm(PA)

model.data <- merge(model.data,traits,by.y='master_species',by.x='species',all.x=T) #add species traits

#These have corallite size, lower depth, growth form
model.datasp <- model.data[ which((model.data$species=='Seriatopora hystrix') |  
                                    (model.data$species=='Seriatopora caliendrum') |                                 
                                    (model.data$species=='Acropora gemmifera') | #Solitaries   no colony size
                                    (model.data$species=='Acropora hyacinthus') | #Japan
                                    (model.data$species=='Acropora nasuta') |
                                    (model.data$species=='Acropora microclados') | #Solitaries
                                    (model.data$species=='Acropora monticulosa') | #Solitaries
                                    (model.data$species=='Pocillopora damicornis') | #AB-brooder
                                    (model.data$species=='Pocillopora eydouxi') | 
                                    (model.data$species=='Isopora cuneata') | #AB-brooder
                                    (model.data$species=='Isopora palifera') |
                                    (model.data$species=='Paragoniastrea australensis') | #AB - more abund in subtropics
                                    (model.data$species=='Turbinaria radicalis') | #AB
                                    (model.data$species=='Turbinaria mesenterina') | #AB
                                    (model.data$species=='Turbinaria frondens') | 
                                    (model.data$species=='Plesiastrea versipora') |
                                    (model.data$species=='Porites lichen') |
                                    (model.data$species=='Porites lobata') |
                                    (model.data$species=='Goniastrea lobata') |
                                    (model.data$species=='Goniastrea retiformis') |
                                    (model.data$species=='Paragoniastrea australiensis') |
                                    (model.data$species=='Galaxea fascicularis') |
                                    (model.data$species=='Porites vaughani')),]


aggregate(cbind(present) ~ species, data = model.datasp, sum)
aggregate(cbind(present) ~ site, data = model.datasp, sum)

#select species
#species, p/a, all envt vars

level.plot(model.datasp[,6], model.datasp[,5:4]) #Chl
level.plot(model.datasp[,7], model.datasp[,5:4]) #low temp
level.plot(model.datasp[,8], model.datasp[,5:4]) #avg temp
level.plot(model.datasp[,9], model.datasp[,5:4]) #temp var
level.plot(model.datasp[,10], model.datasp[,5:4]) #turb
level.plot(model.datasp[,11], model.datasp[,5:4]) #PAR

par(mfrow=c(2,3))
plot(density(na.omit(model.datasp$Chl_50p)),main="Chl")
plot(density(na.omit(model.datasp$SST_2p)),main="low temp")
plot(density(na.omit(model.datasp$SST_50p)),main="avg temp")
plot(density(na.omit(model.datasp$stdev)),main="temp var")
plot(density(na.omit(model.datasp$TSM_50P)),main="turb")
plot(density(na.omit(model.datasp$winter_par_early)),main="winter par")

pairs(~Chl_50p+SST_2p+SST_50p+ stdev + TSM_50P +winter_par_early, data=model.datasp,
      main="Scatterplot Matrix of raw environmental variables")
log.chl <- log(model.datasp$Chl_50p)
log.turb <- log(model.datasp$TSM_50p)
log.lowtemp <- log(model.datasp$SST_2p)
log.avgtemp <- log(model.datasp$SST_50p)
log.var <- log(model.datasp$stdev)
log.PAR <- log(model.datasp$winter_par_early)

model.datasp$PAR <- ((log.PAR-mean(log.PAR, na.rm=TRUE))/(sd(log.PAR, na.rm=TRUE)))

model.datasp$turb <- ((log.turb-mean(log.turb, na.rm=TRUE))/(sd(log.turb, na.rm=TRUE)))

model.datasp$chl <- ((log.chl-mean(log.chl, na.rm=TRUE))/(sd(log.chl, na.rm=TRUE)))

model.datasp$lowtemp <- ((log.lowtemp -mean(log.lowtemp, na.rm=TRUE))/(sd(log.lowtemp, na.rm=TRUE)))

model.datasp$var <- ((log.var-mean(log.var, na.rm=TRUE))/(sd(log.var, na.rm=TRUE)))

model.datasp$avgtemp <- ((log.avgtemp-mean(log.avgtemp, na.rm=TRUE))/(sd(log.avgtemp, na.rm=TRUE)))


######################



#ADD trait data
model.datasp$calc_rate <- as.numeric(as.vector(model.datasp$Calcification.rate))
model.datasp$col_size <- as.numeric(as.vector(model.datasp$Colony.maximum.diameter))
model.datasp$cor_size <- as.numeric(as.vector(model.datasp$Corallite.width.maximum))
model.datasp$gf <- factor(model.datasp$Growth.form.typical)
model.datasp$gr_rate <- as.numeric(as.vector(model.datasp$Growth.rate))
model.datasp$polyps <- as.numeric(as.vector(model.datasp$Polyps.per.area))
model.datasp$skel_dens <- as.numeric(as.vector(model.datasp$Skeletal.density))
model.datasp$skel_micro <- as.numeric(as.vector(model.datasp$Skeletal.micro.density))
model.datasp$fam <- factor(model.datasp$molefam)
model.datasp$lower <- as.numeric(as.vector(model.datasp$Depth.lower))
model.datasp$upper <- as.numeric(as.vector(model.datasp$Depth.upper))

plot(model.datasp$cor_size)
plot(model.datasp$lower)
plot(model.datasp$gr_rate)


#PAR
#turbidity
#temperature (low. variability?)
#############
#corallite size
#depth



###PAR
modP<- glmer(present ~PAR + (1 +  PAR |species), data= model.datasp,family=binomial(link=logit))
modP
coef(modP)
#species distributions of those with negative PAR coef should generally be found where less light
plot(coef(modP))

#Species that can be found deeper can be found where PAR is lower because they can handle low light 
modP1 <- glmer(present ~ PAR + lower:PAR+(1 +PAR|species), data= model.datasp,family=binomial(link=logit))
modP1
coef(modP1)
plot(coef(modP1))

#Species with larger corallites can be found where PAR is low because they have more exposure to access light
modP2 <- glmer(present ~ PAR + cor_size:PAR+(1 +PAR|species), data= model.datasp,family=binomial(link=logit))
modP2
coef(modP2)
plot(coef(modP2))



modP3 <- glmer(present ~ PAR + 
                 cor_size:PAR+
                 lower:PAR+(1 +PAR|species), data= model.datasp,family=binomial(link=logit))
modP3
coef(modP3)

modP4 <- glmer(present ~ PAR + cor_size:PAR+lower:PAR+col_size:PAR+(1 +PAR|species), data= model.datasp,family=binomial(link=logit))
modP4
coef(modP4)



###low temp
modSST<- glmer(present ~lowtemp + (1 +  lowtemp |species), data= model.datasp,family=binomial(link=logit))
modSST
coef(modSST)
#species distributions of those with negative SST coef should generally be found where temp lower
plot(coef(modSST))

#Species with faster growth rates can do better in colder water to outcompete macroalgae 
modSST1 <- glmer(present ~ lowtemp +
                 gr_rate:lowtemp+(1 +lowtemp|species), data= model.datasp,family=binomial(link=logit))
modSST1
coef(modSST1)
plot(coef(modSST1))

#Species with larger colony size may be better competitors in cold water
modSST2 <- glmer(present ~ PAR + col_size:PAR+(1 +PAR|species), data= model.datasp,family=binomial(link=logit))
modSST2
coef(modSST2)
plot(coef(modSST2))




###TURBIDITY (and corallite size)
modT <- glmer(present ~ turb + (1 +turb|species), data= model.datasp,family=binomial(link=logit))
modT
coef(modT)
#species distributions of those with negative turb coef should generally be found where water is less turbid
plot(coef(modT))

#Species with larger corallite should be in more turbid water, they're better at expelling sediments and handling the lower light due to turbid water
modT1 <- glmer(present ~ turb + cor_size:turb+(1 +turb|species), data= model.datasp,family=binomial(link=logit))
modT1
coef(modT1)

modT2 <- glmer(present ~ turb + cor_size:turb+col_size:turb+(1 +turb|species), data= model.datasp,family=binomial(link=logit))
modT2
coef(modT2)

#coral shape should vary with resistence to turbidity (sediment doesn't settle well on branching morphology but would settle on table/plate forms easily)
mod T3 <- glmer(present ~ turb + gf:turb+(1 +turb|species), data= model.datasp,family=binomial(link=logit))
modT3


modPt <- glmer(present ~ PAR + turb + cor_size:PAR+lower:PAR+cor_size:turb+lower:turb+(1 +PAR+turb|species), data= model.datasp,family=binomial(link=logit))
modPt
coef(modPt)
ggCaterpillar(ranef(modPt, postVar=TRUE))
dotplot(ranef(modPt, postVar=TRUE))$species[1]
dotplot(ranef(modPt, postVar=TRUE))$species[2]
dotchart(fixef(modPt))














model1 <- glmer(present ~SST + skel_dens:SST + 
                 (1 +  SST |species),
               data= model.dataspB,family=binomial(link=logit)) 

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





