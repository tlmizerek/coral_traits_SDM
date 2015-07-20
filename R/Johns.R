d <- readLines('http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.par.seawifs.php', encoding = 'UTF-8')
grep('\\.gz', d, val=T)
ff <- grep('\\.gz', grep('\\./data', unlist(strsplit(d, 'a href')), val=T), val=T)
ff <- paste0('http://orca.science.oregonstate.edu/data', gsub('.*\\./data|\">.*', '', ff))

download.file(ff[140], f <- tempfile(fileext = '.gz'), mode = 'wb')
untar(f, exdir='~/Desktop/coral_traits_SDM/data')


library(raster)
library(ncdf4)
library(rgdal)
raster('~/Downloads/par.1998091.hdf')
GDALinfo('~/Downloads/par.1998091.hdf')
getSds('~/Downloads/par.1998091.hdf')
gdalDrivers()

library(MODIS)
