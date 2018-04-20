######
#
#   Use file size to kill bad feathers
#
#     requires first find . -size -2M > smallFeatherFile
#
#
#

# keep track of which tiles are cleared
# only kill featuers

library(data.table)

smalls = fread("../../data/features/smallFilesList", header=F)

smalls[, whatever := tstrsplit(V1,"/",keep=1L)]
smalls[, tile := tstrsplit(V1,"/",keep=2L)]
smalls[, feat := tstrsplit(V1,"/",keep=3L)]

## only consider features
smalls = smalls[!is.na(feat),]
smalls = smalls[tile!="zipArchive",]

## construct
smalls[, filepath := paste0("../../data/features/",tile,"/",feat)]

## which tiles are erroneous
badTiles = unique(smalls$tile)

## remove them
file.remove(smalls$filepath)

write.table(badTiles, "tilesToRedo", quote=F, row.names=F, col.names=F)
