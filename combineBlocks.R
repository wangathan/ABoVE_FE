####
#
#	To combine the training set features
#
#

library(data.table)

featureFiles = list.files("../../data/features/",
													pattern = "20180411.*bl.*",
													full.names=T)

featuresdt = lapply(featureFiles, fread)
featuresdt = rbindlist(featuresdt, use.names=T, fill=T)
#featuresdt[,V1:=NULL]

# include toolik lake and chasmer fens
tl = fread("../../data/training_coefs/training_TL_coefs.csv")
lc = fread("../../data/training_coefs/training_LC_coefs.csv")
# drop the coefficients
toDrop = "a1|b1|c1|a0|a2|b2|a3|b3|magnitude"
dropNames = names(lc)[grep(toDrop, names(lc))]
tl[, (dropNames) := NULL]
lc[, (dropNames) := NULL]

# drop the land cover classes (can go look them up later)
tl[, Toolik_SurfaceClass_Mosaic_ABoVE:= NULL]
lc[, Scotty_classification_ABoVE := NULL]

featuresdt = rbindlist(list(featuresdt, tl, lc), fill=T)
featuresdt[,c("x","y","cell","i.yr_end", "V1", "i.px", "i.py"):=NULL]

# try imputing here, by ecozone
for(n in names(featuresdt)[!grepl("tile|set|pixel|ecozone", names(featuresdt))]){

  if(class(unlist(featuresdt[,n,with=F]))=="integer")
    featuresdt[,(n) := as.numeric(get(n))]
  featuresdt[, (n) := ifelse(is.na(get(n)), median(get(n),na.rm=T),get(n)), by=ecozone]

}

getMode = function(x){

  tb = table(x)
  index = which.max(tb)
  return(names(tb)[index])

}

# impute ecozone using the tiles
featuresdt[, ecozone := ifelse(is.na(ecozone), getMode(ecozone), ecozone),by = tile]

#featuresdt[, names(featuresdt)[grepl("NaN",names(featuresdt))]:=NULL]

write.csv(featuresdt,"../../data/features/featureSet_20180411_coreSample.csv")


