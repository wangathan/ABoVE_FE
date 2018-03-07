##########
#
#		To calculate slope, aspect, etc for ASTER data
#
#
#
 
library(insol) # for calculating sun position
library(parallel)
library(raster)
library(rgdal)
#library(RSAGA)
#library(dynatopmodel) topographic wetness maybe? maybe later?
library(googleway) # for lazy calculating time zone 

#ti = commandArgs(TRUE)[1]

allTi = unlist(lapply(strsplit(list.files("../../../ABoVE_Xtra/ASTER_DEM"),
								 			"\\."),
											"[[",
											2))

calcAster = function(ti){
  # Read in ASTER
  elvout = paste0("../../data/aster/",ti,"/AST_DEM_",ti,"_elv.tif")
  aspout = paste0("../../data/aster/",ti,"/AST_DEM_",ti,"_asp.tif")
  slpout = paste0("../../data/aster/",ti,"/AST_DEM_",ti,"_slp.tif")
 
	if(file.exists(elvout) & file.exists(aspout) & file.exists(slpout))return(NULL)

  asterin = paste0("../../../ABoVE_Xtra/ASTER_DEM/AST_DEM.",ti,".tif")
  asterraster = raster(asterin)
  #asterlon_aea = mean(asterraster@extent@xmin,asterraster@extent@xmax)
  #asterlat_aea = mean(asterraster@extent@ymin,asterraster@extent@ymax)
  
  # convert to ll
  #d = data.frame(lon=asterlon_aea,lat=asterlat_aea)
  #coordinates(d) = c("lon", "lat")
  #proj4string(d) = CRS(projection(asterraster))
  #asterll = spTransform(d, CRS("+init=epsg:4326"))
  
  #sgridfi = paste0("../../data/aster/SAGA/AST_DEM.",ti)
  #writeRaster(asterraster,
  #						sgridfi,
  #						format='SAGA') # a few seconds
  
  # Calc things
  #terrainthings = rsaga.slope.asp.curv(in.dem = paste0(sgridfi,".sgrd"),
  #																		 out.we
  																		 
  
  																		 
  # clean up
  asterraster[asterraster %in% c(-32768, 32767)] = NA																		 
	asterraster[is.na(asterraster)] = median(values(asterraster), na.rm=T)
  # calc things and save
  asterasp <- terrain(asterraster, 
  										 opt = "aspect",
  										 unit = "degrees")
  asterslope <- terrain(asterraster, 
  										 opt = "slope",
  										 unit = "degrees")
  
 
  dir.create(paste0("../../data/aster/",ti), showWarnings=F)	
  writeRaster(asterraster, elvout, overwrite=T)
  writeRaster(asterasp, aspout, overwrite=T)
  writeRaster(asterslope, slpout, overwrite=T)

	print(paste0("Done with ti: ", ti))
}

mclapply(allTi, calcAster, mc.cores=detectCores())

###############################
#
#		Some testing to see if hillshade was viable. Seems like no, comes out as white noise :\
#

## time is 10:30 local? need to determine time zone from lon
## date is doy 120-300, calc all of them and take average
##
#toRadians = function(x)
#{
#	  return(x*pi/180)
#}
#
#toDegrees = function(x)
#{
#	  return(x*180/pi)
#}
#
#decl = function(d)
#{
#	  #make sure d is in Day of Year!
#	  x = 2*pi/365 * (d-172)
#  #x_radians = toRadians(x)
#  decl = 23.5 * cos(x) 
#	  return (decl) #in degrees!
#}
#hrAngle = function(t,long)
#{
#	  #t is in hours since solar noon 
#		#since the hour angle describes the east-west relative displacement of the sun,
#		#seems like there should be a correction for local longitude (though it's not in the notes)
#		
#		#And I've shifted the 0 point of t because I prefer my daily insolation graphs to have midnight be 0, noon be 12. 
#	  ha = 2*pi*(t-12)/24 +  toRadians(long)
#  return(ha)
#}

#zenith = function(lat,decl,ha)
#{
#	  latR = toRadians(lat)
#  declR = toRadians(decl)
#	  cosPsi = sin(latR)*sin(declR) + cos(latR)*cos(declR)*cos(ha)
#	  Zen = acos(cosPsi)
#		  return(Zen)
#}

#azimuth = function(lat,decl,ha)
#{
#	  Z = zenith(lat,decl,ha)
#  declR = toRadians(decl)
#	  sinA = -cos(declR)*sin(ha)/sin(Z)
#	  Azi = asin(sinA)
#		  return(Azi#)
#}

#tz = google_timezone(c(asterll@bbox[2,1], asterll@bbox[1,1]), key=NULL)
#sp = sunpos(sunvector((120:270)+(10.5/24), #julian days + fractional days/time 
#											asterll@bbox[2,1],
#											asterll@bbox[1,1], #
#											tz$rawOffset/3600)) #tz 
#
#allDecl=decl(120:300)
#theHA = hrAngle(-1.5+tz$rawOffset/3600, 0)
								#,asterll@bbox[1,1])

#theHA = mean(hourangle(120:300 + 0.4,asterll@bbox[1,1],tz$rawOffset/3600))
#allEle = toDegrees(pi/2 - zenith(asterll@bbox[2,1],allDecl,theHA))
#allAzi = toDegrees(azimuth(asterll@bbox[2,1],allDecl,theHA))

#astershade <- hillShade(asterslope, 
#												asterasp)#, 
#												#angle = mean(allEle),
#												#direction = mean(allAzi))

# calc hillshade??? Need sun elevation and azimuth. Can calculate average of those... each ASTER file as an average latlon that I'll use



# 

