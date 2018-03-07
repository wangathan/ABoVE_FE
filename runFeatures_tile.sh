#!/bin/bash -l
#$ -V
#$ -l mem_total=40G
#$ -l h_rt=24:00:00
#$ -pe omp 28

#module load mrsid_dsdk/9.1.0.4045
#module load sqlite3/3.17.0
#module load proj4/4.8.0
#module load libgta/1.0.5
#module load hdf/4.2.11_no_netcdf
#module load hdf5/1.8.16
#module load netcdf/4.4.0
#module load jasper/1.900.1
#module load ecw/3.3.20060906
#module load mrsid_dsdk/9.1.0.404
#module load xerces-c/3.1.1
#module load libkml/1.3.0
#module load libdap/3.12.0
#module load geos/3.6.1
#module load freexl/1.0.2
#module load libspatialite/4.3.0a
#module load epsilon/0.9.2
#module load webp/0.4.2
#module load postgresql/9.4.4
#module load gdal/2.1.3
#module load rgdal

module purge
source ~/.bashrc

Rscript getFeatures_tilepoch.R $1 $2
#Rscript getFeatures_tile.R $1 $2 
