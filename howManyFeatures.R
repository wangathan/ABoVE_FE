############################
#
#		Just to look up how many features there are TOTAL
#
#		Loop through feathers, count each
#
#

library(feather)
library(data.table)
library(parallel)

tiledt = data.table(tile = fread("../analyzeRegion/coreTiles")$x)

rowCounter = function(ti){
	print(ti)

	# get all feathers
	featherlist = list.files(paste0("../../data/ccdc_feathers/",ti),
													 full.names=T)

#	featherlist = featherlist[1:100]

	# loop through them
	tilerows = mclapply(featherlist, function(x)feather_metadata(x)$dim[1], mc.cores=detectCores())
	totalrows = sum(unlist(tilerows), na.rm=T)
	return(totalrows)

}

counts = lapply(tiledt$tile, rowCounter)
tiledt[,rows := unlist(counts)]
fwrite(tiledt, "coreTilesFeaturesCount.csv")

tiledt[,rows_seconds := rows*0.0077]
tiledt[,rows_hours := rows_seconds / 60 / 60 / 30]
# 

