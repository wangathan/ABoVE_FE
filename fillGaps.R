##################################
#
#		Produce gap feature files
#
#		Basically for classifying snow
#
#

library(data.table)
library(feather)
library(parallel)

# for a tile get the files
# for each file determine the missing pixels
# check the NC data for each missing pixel 
# build gap data 

source("featureFunctions.R")

ti = commandArgs(TRUE)[1]

# build coefsdt 
# get coefs                                                      
featherfiles = list.files(paste0("../../data/ccdc_feathers/",ti),
													                          full.names=T) 

for(fb in 1:length(featherfiles)){
	print(fb)

	f = featherfiles[fb]

	fout = strsplit(basename(f),"\\.npz\\.ft")[[1]]
	#fid = strsplit(fout, "_c")[[1]][2] # for getting the checklist
	fout = paste0(fout, "_features_gap.ft")

	fout= paste0("../../data/features_gap/",ti,"/",fout)
	if(file.exists(fout))next

	# to determine which segments need to be featurized
	# checklist = fread(paste0("../../data/checklist/",ti,"/check_c",fid,".csv"))
	# checklist = checklist[epoch == thepoch,]

	print(Sys.time())
	print(basename(f))
	print(paste0('--- ffile: ', fb))

	coefsdt = as.data.table(read_feather(f))

	# some spatial info
	coefsdt[,pixel:=paste(px,py,sep="-")]
	coefsdt[,tile:=ti]

	missingPx = getMissingPx(coefsdt)
	print(paste0("tile ", ti, " is missing ", length(missingPx), " pixels!"))

	if(length(missingPx)==0){
		#save a dummy file
		missingdt = data.table(px = -1, py = -1)
		dir.create(paste0("../../data/features_gap/",ti), showWarnings=F, recursive = T)
		write_feather(missingdt, fout)
		next
	}

	system.time(
		missingdt_l <- mclapply(missingPx, fillGap, tile = ti, mc.cores=detectCores())
	)
	missingdt = rbindlist(missingdt_l)
	missingdt[is.na(snowmelt), snowmelt:=200]
	missingdt[is.na(snowfall), snowmelt:=200]

	missingdt[, tile := ti]

	coefsdt = missingdt

	# get spatial data                                                                                       
	alllat = mclapply(1:nrow(missingdt),llApplier,T,mc.cores=detectCores())                                    
	alllon = mclapply(1:nrow(missingdt),llApplier,F,mc.cores=detectCores())

	missingdt[,lat := unlist(alllat)]
	missingdt[,lon := unlist(alllon)]

	# get ecozone data                                                     
	ecozones = mclapply(1:nrow(missingdt),getEcozone, mc.cores=detectCores())
	missingdt[, ecozone := unlist(ecozones)]                                 

	dir.create(paste0("../../data/features_gap/",ti), showWarnings=F, recursive = T)
	write_feather(missingdt, fout)
}
