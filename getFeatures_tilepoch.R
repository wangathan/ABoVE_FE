##############
#
#	 To build out the features of a tile 
#	 	To do that, we need feature functions, ccdc outputs
#
#

library(data.table)
library(feather)
library(parallel)

ti = commandArgs(TRUE)[1]
b = commandArgs(TRUE)[2]
# Jan 27 2018 - after some optimizations seems like it might be fast enough to skip the epoch stuff? 
# just do each file fully
# thepoch = as.numeric(commandArgs(TRUE)[3])

if(suppressWarnings(!is.na(as.numeric(b))))
	b = as.numeric(b)

# get feature functions
source('featureFunctions.R')

# get coefs
featherfiles = list.files(paste0("../../data/ccdc_feathers/",ti),
													full.names=T)

print(ti)

# FYI
print("cores: ")
print(detectCores())


# filling last gaps
allFiles = 1:3000 -1
#if(thepoch == 2010)
	#doneFiles = list.files(paste0("../../data/features/",ti))
#}else{
#	doneFiles = list.files(paste0("../../data/features_",thepoch,"/",ti))
#}
#fileID = as.numeric(unlist(lapply(strsplit(doneFiles, "_c|_f"),"[[",2)))
#notDone = allFiles[!allFiles %in% fileID]
#toDoID = paste0("_c",notDone,".npz")


if(class(b) == "numeric"){
	bstart = (b-1)*100+1
	bend = bstart + 100
	if(bend > length(featherfiles))bend=length(featherfiles)

	blockdo = bstart:bend
}
#if(class(b) == "character"){
#	if(b == "leftovers"){
#		blockdo = unlist(lapply(toDoID,grep,featherfiles))
#	}
#}

# for those days when you need to be VERY specific
#filedo = paste0("c",715:725,"\\.")  
#filedo = paste(filedo, collapse="|")
#blockdo = grep(filedo, featherfiles)

# for testing 
#featherfiles = featherfiles[sample(1:length(featherfiles), 100)]
print(Sys.time())
for(fb in blockdo){#bstart:bend){

	f = featherfiles[fb]

	fout = strsplit(basename(f),"\\.npz\\.ft")[[1]]
	#fid = strsplit(fout, "_c")[[1]][2] # for getting the checklist
	fout = paste0(fout, "_features.ft")

#	if(thepoch == 2010){
		fout= paste0("../../data/features/",ti,"/",fout)
#	}else{
#		fout= paste0("../../data/features_",thepoch,"/",ti,"/",fout)
#	}
	if(file.exists(fout))next

	# to determine which segments need to be featurized
	# checklist = fread(paste0("../../data/checklist/",ti,"/check_c",fid,".csv"))
	# checklist = checklist[epoch == thepoch,]

	print(Sys.time())
	print(basename(f))
	print(paste0('--- ffile: ', fb))

	coefsdt = as.data.table(read_feather(f))

	# drop unnecessary vars
	dropnames = names(coefsdt)[grepl("^a0|^a1|^b1|^c1|^a2|^b2|^a3|^b3|magnitude", names(coefsdt))]
	coefsdt[,(dropnames):=NULL]

	#print(Sys.time())
	print('read in coefficient data!')

	# get spacetime info
	coefsdt[,yr_start:=year(days(start)+ymd("0001-1-1"))]
	coefsdt[,yr_end:=year(days(end)+ymd("0001-1-1"))]
	
	# filter by epoch
#	coefsdt[,epoch := assignEpoch(yr_start, yr_end), by=1:nrow(coefsdt)]
	# coefsdt = coefsdt[epoch == thepoch,]
	if(nrow(coefsdt)==0){
#		# ok so if there's nothing in this feather then skip... make sure to check for existence on featurelist
		next
	}
#	coefsdt = coefsdt[yr_start <= 2010 & yr_end >= 2010,]
	print(paste0('there are ', nrow(coefsdt), ' rows!'))

	# some spatial info
	coefsdt[,pixel:=paste(px,py,sep="-")]
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

	#print(Sys.time())
	print("calculated lat and lon!")

	# get ecozone data
	ecozones = mclapply(1:nrow(coefsdt),getEcozone, mc.cores=detectCores())
	coefsdt[, ecozone := unlist(ecozones)]

	#print(Sys.time())
	print("assigned ecozones!")

	msynth = mclapply(1:nrow(coefsdt), getMonthlySynthetics, mc.cores=detectCores())
#	system.time(
#		msynth <- mclapply(1:3000, getMonthlySynthetics, mc.cores=detectCores())
#	)
#
#	# some retries to make sure it's complete
#
	msynthclass = unlist(lapply(lapply(msynth,class),"[[",1))
	retrylist1 = (1:nrow(coefsdt))[which(msynthclass == "try-error")]
	retryCounter = 0

	while(length(retrylist1) > 0){

		print(retryCounter)

		# keep track of how much the list reduces in size each time
		startingList = length(retrylist1)

		retry1 = mclapply(retrylist1, getMonthlySynthetics, mc.cores=detectCores())
		retryclass1 = unlist(lapply(lapply(retry1,class),"[[",1))
		retrylist1 = retrylist1[which(retryclass1 == "try-error")]

		msynth = c(msynth, retry1)


		# see how much the retry list has fallen in size
		endingList = length(retrylist1)
		if(endingList < startingList)
			retryCounter = 0
		if(endingList >= startingList)
			retryCounter = retryCounter + 1

		if(retryCounter > 12)
			break

	}

	msynthclass = unlist(lapply(lapply(msynth,class),"[[",1))
	msynth = msynth[which(msynthclass != "try-error")]
	msynthdt = rbindlist(msynth, use.names=T, fill=T)

	print(Sys.time())
	print("calculated monthly features!")

	msynthdt[,pixel:=paste0(px,'-',py)]

	rm(msynth)

	setkey(msynthdt,pixel,yr_start)
	setkey(coefsdt,pixel,yr_start)

	coefsdt2 = coefsdt[msynthdt]
	rm(coefsdt)
	rm(msynthdt)

#	if(thepoch == 2010){
		dir.create(paste0("../../data/features/",ti), showWarnings=F, recursive = T)
#	}else{
 # 	dir.create(paste0("../../data/features_",thepoch,"/",ti), showWarnings=F, recursive = T)
#	}
	write_feather(coefsdt2, fout)
	print("feather written!")

	rm(coefsdt2)
	rm(msynthdt)
	rm(msynth)

}
