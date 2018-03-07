####
#
#	To combine the training set features
#
#

library(data.table)

featureFiles = list.files("../../data/features/",
													pattern = "20180118.*bl.*",
													full.names=T)

featuresdt = lapply(featureFiles, fread)
featuresdt = rbindlist(featuresdt, use.names=T, fill=T)
#featuresdt[,V1:=NULL]

#featuresdt[, names(featuresdt)[grepl("NaN",names(featuresdt))]:=NULL]

write.csv(featuresdt,"../../data/features/featureSet_20180119_coreSample.csv")


