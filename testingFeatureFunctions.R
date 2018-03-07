##### functions for feature building


library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(lubridate)
library(raster)
library(ncdf4)
library(rgdal)

ABoVE = readOGR("../../data/above_shp",
								"ABoVE_30mgrid_tiles_Final")

## spatial functions
getPxShp = function(px,py,ti){

	# get the tile
	Bh = as.numeric(substr(ti, 3,4))
	Bv = as.numeric(substr(ti, 6,7))

	theAh = Bh %/% 6
	theAv = Bv %/% 6

	theBh = Bh %% 6
	theBv = Bv %% 6

	tile = ABoVE[ABoVE@data$Ahh == theAh & ABoVE@data$Avv == theAv &
							 ABoVE@data$Bh == theBh  & ABoVE@data$Bv == theBv,]

						 ULx = extent(tile)@xmin
						 ULy = extent(tile)@ymax

						 # each pixel is 30m big. just the upper left is enough for this
						 pxUL = matrix(c(ULx + px*30, ULy - py*30), ncol=2)

						 pxUL_sp = SpatialPoints(coords = pxUL, proj4string = CRS(projection(tile)))

						 return(pxUL_sp)
}

# lat/lon picker
getLatLon = function(px,py, ti, latFlag) {

	pxUL_sp = getPxShp(px,py,ti)

	pxll = spTransform(pxUL_sp,CRS("+proj=longlat +datum=WGS84"))

	return(ifelse(latFlag,extent(pxll)@ymin,extent(pxll)@xmin))

}

llApplier = function(i,latFlag){
	out = try(getLatLon(coefsdt[i,px],coefsdt[i,py],coefsdt[i,tile],latFlag), silent=T)
	if(class(out)=="try-error"){
#		print(i) # ok there was an error before but now there isn't no clue why
		return(NULL)
	}
	return(out)
}

# ecozones!
akeco = readOGR(dsn = "../../data/ecoregions/alaska",
								layer = "akeco_lvl2")
caneco = readOGR(dsn = "../../data/ecoregions/canada",
								 layer = "can_ecozones")

# keep the ecozones nice and aligned
akeco = spTransform(akeco, CRS(projection(caneco)))

getEcozone = function(i){

	# get coordinates
	pxUL_sp = getPxShp(coefsdt[i,px], coefsdt[i,py], coefsdt[i,tile])

	# try intersecting with canada first
	ecozone = over(pxUL_sp, caneco)

	# try intersecting with alaska
	if(is.na(ecozone)){
		ecozone = as.character(over(pxUL_sp, akeco)$LEVEL_2)
	}else{
		ecozone = as.character(ecozone[1,"zone"])
	}
	# return the name of the intersected polygon
	return(ecozone)

}




getPredict = function(coefs_r){
	# get synthetic data
	nmcoefs = c(rep("robust_a0_",7),
							rep("robust_c1_",7),
							rep("robust_a1_",7),
							rep("robust_b1_",7),
							rep("robust_a2_",7),
							rep("robust_b2_",7),
							rep("robust_a3_",7),
							rep("robust_b3_",7))
	bands = c("blue", "green", "red", "nir", "swir1", "swir2", "bt")
	cb = data.table(nmcoefs=nmcoefs, bands=bands)
	cb[,coefs_nms:=paste0(nmcoefs,bands)]
	w = 2*pi/365.25
	days = coefs_r$start:coefs_r$end
	coefs = matrix(unlist(coefs_r[1,cb$coefs_nms,with=F]),
								 byrow=T,
								 ncol=7)

	#coefs=as.matrix(coefs)

	#rcoefs_names = paste0("robust_", coefs_names)
	#rcoefs = as.matrix(coefs_r[1,rcoefs_names,with=F])
	xmatr = cbind(rep(1,length(days)),
								days,
								cos(w*days),
								sin(w*days),
								cos(2*w*days),
								sin(2*w*days),
								cos(3*w*days),
								sin(3*w*days))


	synth = xmatr %*% coefs
	synth = synth/10000

	out = as.data.table(cbind.data.frame(days,synth))

	out[,date := ymd("0001-1-1") + days(days)]
	out[,yr:=year(date)]
	out[,doy:=yday(date)]

	setnames(out, c("days","blue","green","red","nir","swir1","swir2","bt", "date", "yr", "doy"))
	#out[,ndvi:=(nir-red)/(nir+red)]
	return(out)
	#   rsynth = xmatr %*% t(rcoefs)
}


getMonthlySynthetics = function(i){

	# get synthetic values
	test = getPredict(coefsdt[i,])

	if(i %% 1000 == 0){
		print(Sys.time())
		print(i)
	}

	# calc indices
	G = 2.5; L = 1; C1 = 6; C2 = 7.5
	tcb_b = 0.2043; tcb_g = 0.4158; tcb_r = 0.5524; tcb_n = 0.5741; tcb_s1 = 0.3124; tcb_s2 = 0.2303
	tcg_b =-0.1603; tcg_g =-0.2819; tcg_r =-0.4934; tcg_n = 0.7940; tcg_s1 = 0.0002; tcg_s2 = 0.1446
	tcw_b = 0.0315; tcw_g = 0.2021; tcw_r = 0.3102; tcw_n = 0.1594; tcw_s1 = 0.6806; tcw_s2 = 0.6109
	test[,ndvi:= (nir - red) / (nir + red)]
	test[,nbr := (nir - swir2) / (nir + swir2)]
	# test[,lswi := (nir - swir1) / (nir + swir1)]
	test[,evi := G * (nir - red) / (nir + C1*red - C2*blue + L)]
	test[,tcb := tcb_b*blue + tcb_g*green + tcb_r*red + tcb_n*nir + tcb_s1*swir1 + tcb_s2*swir2]
	test[,tcg := tcg_b*blue + tcg_g*green + tcg_r*red + tcg_n*nir + tcg_s1*swir1 + tcg_s2*swir2]
	test[,tcw := tcw_b*blue + tcw_g*green + tcw_r*red + tcw_n*nir + tcw_s1*swir1 + tcw_s2*swir2]

	# cleaning up bad data
	test[nbr< -1 | nbr>1, nbr:=NA]
	test[evi< -1 | evi>1, evi:=NA]
	test[ndvi< -1 | ndvi>1, ndvi:=NA]
	# test[lswi< -1 | lswi>1, lswi:=NA]
	test[tcb< -1 | tcb>1, tcb:=NA]
	test[tcw< -1 | tcw>1, tcw:=NA]
	test[tcg< -1 | tcg>1, tcg:=NA]

	# calc differences
	test[,nbrevi := nbr - evi]
	test[,tcwgd  := tcw - tcg]

	# some cleaning up
	# read netcdf files to get some Landsat data characteristics

	# get px, py
	# also need year start/end
	px = coefsdt[i,px]
	py = coefsdt[i,py]
	yr_start = coefsdt[i,yr_start]
	yr_end = coefsdt[i,yr_end]
	tile = coefsdt[i,tile]

	# read nc file
	#ncf_l <-paste0("../../../NC_FILES/", tile, "/output");  the slow version
	#ncf_ls <-list.files(ncf_l, full.names=T);
	#ncf <-ncf_ls[grep(paste0(py,"-|-",py,".nc"),ncf_ls)]
#testspeed = function(){
	itsTheFirstOne <-py%%2 == 0
	itsTheSecondOne <-py%%2 == 1
	filenm <-paste0(ifelse(itsTheFirstOne, paste0(py, '-', py+1), paste0(py-1,'-',py)), '.nc')
	ncf <- paste0('../../../NC_FILES/',tile,'/output/',filenm)

	if(itsTheFirstOne) gety = 1
	if(itsTheSecondOne) gety = 2

	nc = nc_open(ncf)

	# nc["dim"] gives a list 3 row matrix of sensor (5,7), some value, and yeardoy (yyyyddd)
	# there's five items in the list... bands, time, x, y, head
	# nc["dim"][["bands"]] gives some details that are mostly blank, useful thing seems to be $vals in 1:8
	# [["time"]] has no units :( but a length indicating the number of images used. vals for just a 1:length of
	# the dataset in the time dimension...
	# man there's like nothing here
	# nc["dim"]$dim$head contains the yeardoy matrix. use $vals. seriously none of the metadata is filled out!!
	yeardoy = nc["dim"]$dim$head$vals[3,]
	yeardoydt = data.table(i = 1:length(yeardoy),
												 yr = as.numeric(substr(yeardoy,1,4)),
												 doy = as.numeric(substr(yeardoy,5,7)))


	# the data
	# nc[["var"]]$data it's in
	# nc[["var"]]$data$dim[1] gives bands for vals
	# nc[["var"]]$data$dim[[2]] gives the object with the data I think, under vals, but it's just 1:length. poo

	# nc[["var"]][["data"]]

	# get actual data:
	# ncvar_get(nc, "data") is 4dim array where it's band x yeardoy x px x py
	# need to identify which py!!!!!!!!
	# # much much faster if we don't read in all the data first
	fmaskdat = try(ncvar_get(nc, "data", start = c(8,1,px,gety), count = c(-1, -1, 1, 1)),silent=T)
	if(class(fmaskdat) == "try-error"){
		print("nc error in:")
		print(i)
		snowmelt = 120
		snowfall = 300
		snowiness = 0.3
		cloudiness = 0.6
		nc_close(nc)
	}else{
		#test = ncvar_get(nc, "data")[8,,px,gety]

		yeardoydt[, fmask := fmaskdat]

		yeardoydt = yeardoydt[yr %in% yr_start:yr_end,]
		nc_close(nc)

		# get statistics on fmask frequency ... 3 is snow, 2 is cloud shadow, 4 is cloud, 255 is no data

		# throw out 255s - they don't overlap this pixel
		yeardoydt = yeardoydt[fmask != 255,]

		# how snowy is it when it's cloud free
		snowiness = nrow(yeardoydt[fmask==3,])/nrow(yeardoydt[fmask %in% c(0,1,3),])

		# how cloudy is it when it's snow free
		cloudiness = nrow(yeardoydt[fmask==4,])/nrow(yeardoydt[fmask %in% c(0,1,2,4),])

		# throw out clouds and cloud shadows
		yeardoydt = yeardoydt[!fmask %in% c(2,4),]

		# when does the snow melt, also considering over water
		# snowmelt = min(yeardoydt[fmask == 0, doy])
		snowmelt = as.integer(quantile(yeardoydt[fmask %in% c(1,0), doy], 0.01, na.rm=T))


		# when does the snow fall
		# snowfall = max(yeardoydt[fmask == 0, doy])
		snowfall = as.integer(quantile(yeardoydt[fmask %in% c(1,0), doy], 0.99, na.rm=T))


		# what if this snowmelt/snowfall value results in having just one value in a month? will cause problems
		if(snowmelt %in% yday(ymd(paste0('2011-',1:12,'-1')))) snowmelt = snowmelt - 2
		if(snowfall %in% yday(ymd(paste0('2011-',1:12,'-1')))-1) snowfall = snowfall + 2
	}
#}
	# filter on year and on doy
	if(!is.na(snowmelt) & !is.na(snowfall)){
		test = test[doy %in% snowmelt:snowfall,]
	}else{
		test = test[doy %in% 100:300,]
	}

	# let's keep colors for now
	#test[,c("blue", "red", "green", "nir", "swir1", "swir2"):=NULL]

	# band names
	bnames = names(test)
	bnames = bnames[!grepl("days|date|yr|doy",bnames)]

	# if nothing remains, throw it all out
	if(nrow(na.omit(test))==0)return(NULL)

	# overall trends
	test[, paste0("trnd_",bnames) := lapply(.SD, function(x)coefficients(lm(x ~ days))[2]), .SDcol=bnames]

	# first get mean/amplitude/min/max for each year, then average over all years, to effectively detrend
	# summary by month
	test[, paste0("m_",bnames) := lapply(.SD, mean, na.rm=T), .SDcol = bnames, by = yr]
	test[, paste0("m_",bnames) := lapply(.SD, mean, na.rm=T), .SDcol = paste0("m_",bnames)]

	# get rough differential?
	#getDiff = function(band, doy, molength){
	#subdt = data.table(band = band, doy = doy)
	#setkey(subdt, doy)
	#return((last(subdt$band) - first(subdt$band))/molength*30)
	#}

	# test[, paste0("diff_",bnames) := lapply(.SD, getDiff, doy = doy, molength=molength), .SDcol = bnames, by = c("mo", "yr")]
	# test[, paste0("diff_",bnames) := lapply(.SD, mean), .SDcol = paste0("diff_",bnames), by = mo]

	test[, paste0("q25_",bnames) := lapply(.SD, function(x)quantile(x, c(0.25),na.rm=T)), .SDcols=bnames, by = yr] # get yearly q25litudes
	test[, paste0("q25_",bnames) := lapply(.SD, mean), .SDcols=paste0("q25_",bnames)] # get yearly q25litudes

	test[, paste0("q75_",bnames) := lapply(.SD, function(x)quantile(x,c(0.75),na.rm=T)), .SDcols=bnames, by = yr] # get yearly q75litudes
	test[, paste0("q75_",bnames) := lapply(.SD, mean), .SDcols=paste0("q75_",bnames)] # get yearly q75litudes

	test[, paste0("amp_",bnames) := lapply(.SD, function(x)diff(range(x, na.rm=T))), .SDcols=bnames, by = yr] # get yearly amplitudes
	test[, paste0("amp_",bnames) := lapply(.SD, mean), .SDcols=paste0("amp_",bnames)] # get yearly amplitudes

	test[, paste0("min_",bnames) := lapply(.SD, min, na.rm=T), .SDcols=bnames, by = yr]
	test[, paste0("min_",bnames) := lapply(.SD, mean), .SDcols=paste0("min_",bnames)] # get yearly min

	test[, paste0("max_",bnames) := lapply(.SD, max, na.rm=T), .SDcols=bnames, by = yr]
	test[, paste0("max_",bnames) := lapply(.SD, mean), .SDcols=paste0("max_",bnames)] # get yearly max

	# clean up monthly values
	# mo = unique(test[,c("mo",paste0("m_",bnames), paste0("diff_",bnames)),with=F])
	ye = unique(test[,c(paste0("m_",bnames), paste0("trnd_",bnames),paste0("amp_",bnames), paste0("min_",bnames), paste0("q25_",bnames), paste0("q75_",bnames), paste0("max_",bnames)),with=F])
	# mo = dcast(mo, . ~ mo, value.var = c(paste0("m_",bnames), paste0("diff_",bnames)))
	# mo[,".":=NULL]

	# merging info
	ye[,yr_start := coefsdt[i,yr_start]]
	ye[,px := coefsdt[i,px]]
	ye[,py := coefsdt[i,py]]

	# and the fmask params
	ye[, c("snowiness", "cloudiness", "snowmelt", "snowfall") := .(snowiness, cloudiness, snowmelt, snowfall)]

	return(ye)
}


# assign epochs in order
assignEpoch =function(yr_start, yr_end){

	epoch = 0
	if(yr_start > 2010)epoch=2014
	if(yr_end < 1990)epoch=1986
	if(1990 %in% yr_start:yr_end)epoch=1990
	if(1994 %in% yr_start:yr_end)epoch=1994
	if(1998 %in% yr_start:yr_end)epoch=1998
	if(2002 %in% yr_start:yr_end)epoch=2002
	if(2006 %in% yr_start:yr_end)epoch=2006
	if(2010 %in% yr_start:yr_end)epoch=2010
	return(epoch)

}

## testing for NA filling
#library(feather)
#f = "../../data/ccdc_feathers/Bh04v06/yatsm_c9.npz.ft"
#coefsdt = as.data.table(read_feather(f))
#
## drop unnecessary vars
#dropnames = names(coefsdt)[grepl("^a0|^a1|^b1|^c1|^a2|^b2|^a3|^b3|magnitude", names(coefsdt))]
#coefsdt[,(dropnames):=NULL]
#
##print(Sys.time())
#print('read in coefficient data!')
#
## get spacetime info
#coefsdt[,yr_start:=year(days(start)+ymd("0001-1-1"))]
#coefsdt[,yr_end:=year(days(end)+ymd("0001-1-1"))]
#
## filter by epoch
#coefsdt[,epoch := assignEpoch(yr_start, yr_end), by=1:nrow(coefsdt)]
#coefsdt = coefsdt[epoch == thepoch,]
#if(nrow(coefsdt)==0){
#	  # ok so if there's nothing in this feather then skip... make sure to check for existence on featurelist
#	  next
#}
#coefsdt = coefsdt[yr_start <= 2010 & yr_end >= 2010,]
#print(paste0('there are ', nrow(coefsdt), ' rows!'))
#
## some spatial info
#coefsdt[,pixel:=paste(px,py,sep="-")]
#coefsdt[,tile:=ti]
#
## some change-based information
#coefsdt[,nbreaks:=(.N-1L),by=pixel]
#coefsdt[,seglength:=(end-start)]

getMissingPx = function(coefsdt){

	# determine which pixels are missing
	py = unique(coefsdt$py)
	px = 0:5999

	allPx = expand.grid(px, py)
	allPx = paste0(allPx$Var1, "-", allPx$Var2)

	missingPx = allPx[!allPx %in% coefsdt$pixel]
	return(missingPx)
}

#missingPx = getMissingPx(coefsdt)

fillGap = function(pixel, tile){

	px=as.numeric(strsplit(pixel,"-")[[1]][1])
	py=as.numeric(strsplit(pixel,"-")[[1]][2])

	# for the missing pixel, grab netCDF data

	itsTheFirstOne <-py%%2 == 0
	itsTheSecondOne <-py%%2 == 1
	filenm <-paste0(ifelse(itsTheFirstOne, paste0(py, '-', py+1), paste0(py-1,'-',py)), '.nc')
	ncf <- paste0('../../../NC_FILES/',tile,'/output/',filenm)

	if(itsTheFirstOne) gety = 1
	if(itsTheSecondOne) gety = 2

	nc = nc_open(ncf)

	# nc["dim"] gives a list 3 row matrix of sensor (5,7), some value, and yeardoy (yyyyddd)
	# there's five items in the list... bands, time, x, y, head
	# nc["dim"][["bands"]] gives some details that are mostly blank, useful thing seems to be $vals in 1:8
	# [["time"]] has no units :( but a length indicating the number of images used. vals for just a 1:length of
	# the dataset in the time dimension...
	# man there's like nothing here
	# nc["dim"]$dim$head contains the yeardoy matrix. use $vals. seriously none of the metadata is filled out!!
	yeardoy = nc["dim"]$dim$head$vals[3,]
	yeardoydt = data.table(i = 1:length(yeardoy),
												 yr = as.numeric(substr(yeardoy,1,4)),
												 doy = as.numeric(substr(yeardoy,5,7)))


	LSdat = try(ncvar_get(nc, "data", start = c(1,1,px,gety), count = c(-1, -1, 1,1)), silent=T)
	LSyddat = yeardoydt
	LSyddat[, blue := LSdat[1,]/10000]
	LSyddat[, green := LSdat[2,]/10000]
	LSyddat[, red := LSdat[3,]/10000]
	LSyddat[, nir := LSdat[4,]/10000]
	LSyddat[, swir1 := LSdat[5,]/10000]
	LSyddat[, swir2 := LSdat[6,]/10000]
	LSyddat[, bt := LSdat[7,]/10000]
	LSyddat[, fmask := LSdat[8,]]

	# get relevant LS data for pixel
	LSyddat = LSyddat[fmask != 255,]
	LSyddat = LSyddat[blue != 1.6 & green != 1.6 & red != 1.6 & nir != 1.6 & swir1 != 1.6 & swir2 != 1.6,] 
		
		# throw out 255s - they don't overlap this pixel
		LSyddat = LSyddat[fmask != 255,]

		# how snowy is it when it's cloud free
		snowiness = nrow(LSyddat[fmask==3,])/nrow(LSyddat[fmask %in% c(0,1,3),])

		# how cloudy is it when it's snow free
		cloudiness = nrow(LSyddat[fmask==4,])/nrow(LSyddat[fmask %in% c(0,1,2,4),])

		# throw out clouds and cloud shadows
		LSyddat = LSyddat[!fmask %in% c(2,4),]

		# when does the snow melt, also considering over water
		# snowmelt = min(LSyddat[fmask == 0, doy])
		snowmelt = as.integer(quantile(LSyddat[fmask %in% c(1,0), doy], 0.01, na.rm=T))

		# when does the snow fall
		# snowfall = max(LSyddat[fmask == 0, doy])
		snowfall = as.integer(quantile(LSyddat[fmask %in% c(1,0), doy], 0.99, na.rm=T))
	nc_close(nc)

	# calc indices
	G = 2.5; L = 1; C1 = 6; C2 = 7.5
	tcb_b = 0.2043; tcb_g = 0.4158; tcb_r = 0.5524; tcb_n = 0.5741; tcb_s1 = 0.3124; tcb_s2 = 0.2303
	tcg_b =-0.1603; tcg_g =-0.2819; tcg_r =-0.4934; tcg_n = 0.7940; tcg_s1 = 0.0002; tcg_s2 = 0.1446
	tcw_b = 0.0315; tcw_g = 0.2021; tcw_r = 0.3102; tcw_n = 0.1594; tcw_s1 = 0.6806; tcw_s2 = 0.6109
	LSyddat[,ndvi:= (nir - red) / (nir + red)]
	LSyddat[,nbr := (nir - swir2) / (nir + swir2)]
	# LSyddat[,lswi := (nir - swir1) / (nir + swir1)]
	LSyddat[,evi := G * (nir - red) / (nir + C1*red - C2*blue + L)]
	LSyddat[,tcb := tcb_b*blue + tcb_g*green + tcb_r*red + tcb_n*nir + tcb_s1*swir1 + tcb_s2*swir2]
	LSyddat[,tcg := tcg_b*blue + tcg_g*green + tcg_r*red + tcg_n*nir + tcg_s1*swir1 + tcg_s2*swir2]
	LSyddat[,tcw := tcw_b*blue + tcw_g*green + tcw_r*red + tcw_n*nir + tcw_s1*swir1 + tcw_s2*swir2]

	# cleaning up bad data
	LSyddat[nbr< -1, nbr:=-1]
	LSyddat[evi< -1, evi:=-1]
	LSyddat[ndvi< -1, ndvi:=-1]
  LSyddat[tcb< -1, tcb:=-1]
	LSyddat[tcw< -1, tcw:=-1]
	LSyddat[tcg< -1, tcg:=-1]

	LSyddat[nbr>1, nbr:=1]
	LSyddat[evi>1, evi:=1]
	LSyddat[ndvi>1, ndvi:=1]
  LSyddat[tcb>1, tcb:=1]
	LSyddat[tcw>1, tcw:=1]
	LSyddat[tcg>1, tcg:=1]

	# calc differences
	LSyddat[,nbrevi := nbr - evi]
	LSyddat[,tcwgd  := tcw - tcg]

	# calc metrics
	# if nothing remains, throw it all out
	if(nrow(na.omit(LSyddat))==0)return(NULL)
	bnames = names(LSyddat)[!names(LSyddat) %in% c("i", "yr", "doy", "fmask")]

	# overall trends
	LSyddat[, paste0("trnd_",bnames) := 0]

	# first get mean/amplitude/min/max for each year, then average over all years, to effectively detrend
	# summary by month
	LSyddat[, paste0("m_",bnames) := lapply(.SD, mean, na.rm=T), .SDcol = bnames, by = yr]
	LSyddat[, paste0("m_",bnames) := lapply(.SD, mean, na.rm=T), .SDcol = paste0("m_",bnames)]

	LSyddat[, paste0("q25_",bnames) := lapply(.SD, function(x)quantile(x, c(0.25),na.rm=T)), .SDcols=bnames, by = yr] # get yearly q25litudes
	LSyddat[, paste0("q25_",bnames) := lapply(.SD, mean), .SDcols=paste0("q25_",bnames)] # get yearly q25litudes

	LSyddat[, paste0("q75_",bnames) := lapply(.SD, function(x)quantile(x,c(0.75),na.rm=T)), .SDcols=bnames, by = yr] # get yearly q75litudes
	LSyddat[, paste0("q75_",bnames) := lapply(.SD, mean), .SDcols=paste0("q75_",bnames)] # get yearly q75litudes

	LSyddat[, paste0("amp_",bnames) := lapply(.SD, function(x)diff(range(x, na.rm=T))), .SDcols=bnames, by = yr] # get yearly amplitudes
	LSyddat[, paste0("amp_",bnames) := lapply(.SD, mean), .SDcols=paste0("amp_",bnames)] # get yearly amplitudes

	LSyddat[, paste0("min_",bnames) := lapply(.SD, min, na.rm=T), .SDcols=bnames, by = yr]
	LSyddat[, paste0("min_",bnames) := lapply(.SD, mean), .SDcols=paste0("min_",bnames)] # get yearly min

	LSyddat[, paste0("max_",bnames) := lapply(.SD, max, na.rm=T), .SDcols=bnames, by = yr]
	LSyddat[, paste0("max_",bnames) := lapply(.SD, mean), .SDcols=paste0("max_",bnames)] # get yearly max

	# clean up monthly values
	# mo = unique(LSyddat[,c("mo",paste0("m_",bnames), paste0("diff_",bnames)),with=F])
	ye = unique(LSyddat[,c(paste0("m_",bnames), paste0("trnd_",bnames),paste0("amp_",bnames), paste0("min_",bnames), paste0("q25_",bnames), paste0("q75_",bnames), paste0("max_",bnames)),with=F])

	# merging info                                                                                            
	ye[,yr_start := min(LSyddat$yr)]                                                                      
	ye[,px := as.numeric(strsplit(pixel,"-")[[1]][1])]                                                                                  
	ye[,py := as.numeric(strsplit(pixel,"-")[[1]][2])]                                                                                  
	                                                                                                          
	# and the fmask params                                                                                    
	ye[, c("snowiness", "cloudiness", "snowmelt", "snowfall") := .(snowiness, cloudiness, snowmelt, snowfall)]
	return(ye)
}
