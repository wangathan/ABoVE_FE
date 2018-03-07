## Jan 30 2018

## to answer these questions about the netcdfs:
## 1: are the NC file lists the same for each tile?
## 2: do they cover all py values? i.e. are there any py orphans?

poop1 = list.files("../../../NC_FILES/Bh04v06/output")
poop2 = list.files("../../../NC_FILES/Bh09v05/output")
poop3 = list.files("../../../NC_FILES/Bh12v12/output")
poop4 = list.files("../../../NC_FILES/Bh06v06/output")
poop5 = list.files("../../../NC_FILES/Bh02v03/output")
poop6 = list.files("../../../NC_FILES/Bh11v11/output")

getpy = function(x){

	pylist = as.numeric(strsplit(x, "-|\\.")[[1]][c(1,2)])
	if(pylist[2]-pylist[1] > 1)
		pylist = c(pylist, pylist[2]-1)
	return(pylist)

}

plist1 = unique(unlist(lapply(poop1, getpy)))
plist6 = unlist(lapply(poop6, getpy))
# I think we're good
