####################
#
#		Some additional features besides the model coefficients
#
#			incl:
#				- num breaks per pixel
#				- length of segments
#				- latitude and longitude
#				- year (for merging with labelled data)
#				- DEM 
#				- surface water
#
#		Also combines it all
#

library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(lubridate)
library(raster)
library(ncdf4)
library(rgdal)

print("num cores:")
print(detectCores())

source("featureFunctions.R")

coef_files = list.files("../../data/training_coefs",
												full.names=T)

registerDoParallel(detectCores())

coefsdt = rbindlist(mclapply(coef_files, fread, mc.cores=detectCores()),fill=T,use.names=T)

# fill missing asp and slp
coefsdt[is.na(asp), asp:=median(asp, na.rm=T), by = set]
coefsdt[is.na(slp), slp:=median(slp, na.rm=T), by = set]

coefsdt = na.omit(coefsdt)
coefsdt = unique(coefsdt)

print(Sys.time())
print('read in coefficient data!')

# get spacetime info
coefsdt[,pixel:=paste(tile,px,py,sep="-")]
coefsdt[,yr_start:=year(days(start)+ymd("0001-1-1"))]
coefsdt[,yr_end:=year(days(end)+ymd("0001-1-1"))]
coefsdt[,yr_br:=year(days(br)+ymd("0001-1-1"))]

# some change-based information
coefsdt[,nbreaks:=(.N-1L),by=pixel]
coefsdt[,seglength:=(end-start)]

alllat = mclapply(1:nrow(coefsdt),llApplier,T,mc.cores=detectCores()) 
alllon = mclapply(1:nrow(coefsdt),llApplier,F,mc.cores=detectCores()) 

coefsdt[,lat:=unlist(alllat)]
coefsdt[,lon:=unlist(alllon)]

print(Sys.time())
print("calculated lat and lon!")

ecozones = mclapply(1:nrow(coefsdt),getEcozone, mc.cores=detectCores())

coefsdt[, ecozone := unlist(ecozones)]

print(Sys.time())
print("assigned ecozones!")

blocks = seq(1, nrow(coefsdt), by = 5000)
for(b in blocks){
	#msynth = mclapply(1:nrow(coefsdt), getMonthlySynthetics,mc.cores=detectCores())

	bstart = b
	bend = b + 4999

	if(bend > nrow(coefsdt))bend = nrow(coefsdt)

#	msynth = mclapply(5001:5010, getMonthlySynthetics)
	msynth = mclapply(bstart:bend, getMonthlySynthetics,mc.cores=detectCores())
	#

	msynthclass = unlist(lapply(lapply(msynth,class),"[[",1))

	retrylist1 = (bstart:bend)[which(msynthclass == "try-error")]
	retry1 = mclapply(retrylist1, getMonthlySynthetics, mc.cores=detectCores())
	retryclass1 = unlist(lapply(lapply(retry1,class),"[[",1))

	retrylist2 = retrylist1[which(retryclass1 == "try-error")]
	retry2 = mclapply(retrylist2, getMonthlySynthetics, mc.cores=detectCores())
	retryclass2 = unlist(lapply(lapply(retry2,class),"[[",1))

	retrylist3 = retrylist1[which(retryclass2 == "try-error")]
	retry3 = mclapply(retrylist3, getMonthlySynthetics, mc.cores=detectCores())
	retryclass3 = unlist(lapply(lapply(retry3,class),"[[",1))

	msynth = msynth[which(msynthclass != "try-error")]
	retry1 = retry1[which(retryclass1 != "try-error")]
	retry2 = retry2[which(retryclass2 != "try-error")]
	retry3 = retry3[which(retryclass3 != "try-error")]

	msynthdt = rbindlist(c(msynth,retry1,retry2,retry3), use.names=T,fill=T)
	
	print(Sys.time())
	print("generated synthetic metrics!")
	
	setkey(msynthdt, px, py, yr_start)
	setkey(coefsdt, px, py, yr_start)
	
	coefsdt2 = coefsdt[msynthdt]
	
	#write.csv(coefsdt2, "../../data/features/featureSet_20171102_hasFilterAndDiffAndEcoz.csv", row.names=F)
	#write.csv(coefsdt2, paste0("../../data/features/featureSet_20171102_hasFilterAndEcoz_bl",(b%/%5000+1),".csv"), row.names=F)
	#write.csv(coefsdt2, paste0("../../data/features/featureSet_20171201_hasFilterAndEcozAndDEMAndSW_bl",(b%/%5000+1),".csv"), row.names=F)
	#write.csv(coefsdt2, paste0("../../data/features/featureSet_20180111_hasFilterAndEcozAndDEMAndSWAndChrom_bl",(b%/%5000+1),".csv"), row.names=F)
	write.csv(coefsdt2, paste0("../../data/features/featureSet_20180118_coreSample_bl",(b%/%5000+1),".csv"), row.names=F)
	
	print(Sys.time())
	print("wrote data! block:")
	print(b)
}


# which i is failing at monthifying?
geti = function(x){

	if(class(x) == "data.table"){
		thei = x[1,i]
	}else{
		return(NULL)
	}
	return(thei)

}

ms_i = lapply(msynth, geti)
goodi = unlist(ms_i)
alli = 1:5000
badi = alli[!alli %in% goodi]

bms_i = lapply(msynth, geti)
bgoodi = unlist(bms_i)
alli = 1:5000
bbadi = alli[!alli %in% bgoodi]
