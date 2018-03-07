######
###
#
#		After building out features, why are some files missing?
#		Which ones are missing?
#

library(feather)
library(data.table)
library(parallel)

ti = commandArgs(TRUE)[1]

print("checking for missing and filling in gaps for:")
print(ti)

# get list of ccdc_feathers
coefs_files = list.files(paste0("../../data/ccdc_feathers/",ti),
												 full.names=T)

coefs_rows = unlist(lapply(strsplit(coefs_files, "_c|\\.npz"),"[[",2))

# get list of feature feathers
feat_files = list.files(paste0("../../data/features/",ti),
												full.names=T)

feat_rows = unlist(lapply(strsplit(feat_files, "_c|_features"),"[[",2))

# missing feature feathers
missing_rows = coefs_rows[!coefs_rows %in% feat_rows]
missing_row_search = paste0("_c",missing_rows,"\\.npz")

source('featureFunctions.R')

# see the missing featherfiles
featherfiles = list.files(paste0("../../data/ccdc_feathers/",ti),
													full.names=T)
featherfiles = featherfiles[grepl(paste(missing_row_search,collapse="|"), featherfiles)]

# loop through and featurify the feathers
for(fb in 1:length(featherfiles)){

	f = featherfiles[fb]

	fout = strsplit(basename(f),"\\.npz\\.ft")[[1]]
	fout = paste0(fout, "_features.ft")

	if(file.exists(fout))next

	print(Sys.time())
	print(basename(f))
	print(paste0('ffile: ', fb))

	coefsdt = as.data.table(read_feather(f))

	# drop unnecessary vars
	dropnames = names(coefsdt)[grepl("^a0|^a1|^b1|^c1|^a2|^b2|^a3|^b3|magnitude", names(coefsdt))]
	coefsdt[,(dropnames):=NULL]

	print(Sys.time())
	print('read in coefficient data!')

	# get spacetime info
	coefsdt[,pixel:=paste(px,py,sep="-")]
	coefsdt[,yr_start:=year(days(start)+ymd("0001-1-1"))]
	coefsdt[,yr_end:=year(days(end)+ymd("0001-1-1"))]

	# filter by year - only interested in 2010 for now
	coefsdt = coefsdt[yr_start <= 2010 & yr_end >= 2010,]
	print(paste0('there are ', nrow(coefsdt), ' rows!'))

	coefsdt[,tile:=ti]

	# some change-based information
	coefsdt[,nbreaks:=(.N-1L),by=pixel]

	coefsdt[,seglength:=(end-start)]

	# get spatial data
	alllat = mclapply(1:nrow(coefsdt),llApplier,T,mc.cores=detectCores())
	alllon = mclapply(1:nrow(coefsdt),llApplier,F,mc.cores=detectCores())

	coefsdt[,lat:=unlist(alllat)]
	coefsdt[,lon:=unlist(alllon)]

	# how many breaks does each pixel experience?
	#coefsdt[,avgbreaks := mean(nbreaks), by = pixel]
	#subdt = unique(coefsdt[,.(pixel, avgbreaks)])

	print(Sys.time())
	print("calculated lat and lon!")

	# get ecozone data
	ecozones = mclapply(1:nrow(coefsdt),getEcozone, mc.cores=detectCores())
	coefsdt[, ecozone := unlist(ecozones)]

	print(Sys.time())
	print("assigned ecozones!")

	msynth = mclapply(1:nrow(coefsdt), getMonthlySynthetics, mc.cores=detectCores())

	msynthdt = rbindlist(msynth, use.names=T, fill=T)

	print(Sys.time())
	print("calculated monthly features!")

	msynthdt[,pixel:=paste0(px,'-',py)]

	setkey(msynthdt,pixel,yr_start)
	setkey(coefsdt,pixel,yr_start)

	coefsdt2 = coefsdt[msynthdt]

	dir.create(paste0("../../data/features/",ti), showWarnings=F)
	write_feather(coefsdt2, paste0('../../data/features/',ti,'/',fout))
	print("feather written!")

	rm(coefsdt2)
	rm(msynthdt)
	rm(msynth)

}
