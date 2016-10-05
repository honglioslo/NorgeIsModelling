#!/bin/sh
### step 1: prepare MET data
### extract pre and tem from SeNorge datasets
### input: 
/hdata/fou/personlig/guru/DEW/Dew/Source/pre_3.01b/stationMask control_mask.txt
mkdir -p /hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Met1km/Pre
mkdir -p /hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Met1km/Tem
export  METDATA=/hdata/grid/metdata
export  METMASK=/hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Met1km
/hdata/fou/personlig/guru/DEW/Dew/Source/pre_3.01b/ExtractSeNorge control_grid.txt

mkdir -p /hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Met100m/Pre
mkdir -p /hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Met100m/Tem
export  METDATA=/hdata/grid/metdata
export  METMASK=/hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Met100m
/hdata/fou/personlig/guru/DEW/Dew/Source/pre_3.01b/gridMaskVariable control_grid100.txt

### step 2: make discharge data
mkdir /hdata/fou/NorgeIsModelling/HardangenjokulenTest/500013Bjoreio/Disc
####################################################################
ts=timeseries

for file in `more $ts`
do
    reg=`awk -v f=$file 'BEGIN{split(f,y,".");print y[1]}'`
    mno=`awk -v f=$file 'BEGIN{split(f,y,".");print y[2]}'`
    
   echo $reg  $mno  
lescon_var -s 196101010000 -e 201501010000 -f timevalue 5 $reg $mno 0 1001 0 > $reg.$mno.0.1001.0.var
timeseriegraph $reg.$mno.0.1001.0.var
done

### step 3: make pest files for calibration 
####################################################################
# Preparing Qobs for pest-file and for instruction-files
####################################################################
# cd /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Aalfotbreen_Calibration_01/Disc
# Run
/hdata/fou/personlig/guru/DEW/Dew/Source/utilities/calibration_data_five timeseries  2000 2014 1 1
cp *_ins ..
cp obs_streamflow.txt ..
cp stations_nve.txt ..
cp pestfile_observed.txt ..

cp dew.pst ..


# output: 
#         obs_streamflow.txt
#         pestfile_observed.txt
#         stations_nve.txt
##
####################################################################
# Preparing pest-file
####################################################################
# 1. insert pestfile_observed.txt after "* observation datain" in
#     file dew.pst
# 2. insert *_ins file name after " * model input/output"
# 3. Change NPAR, NOBS and NINSFLE in pest control file
# output: dew.pst
####
####################################################################
# Run stationMaskone more time if calibration on small area
####################################################################
# 1. preparing control_mask.txt
# File with sub-catchment identifiers                    : watershed.txt
# File with grid cell sub-catchment identifiers          : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/watershed_area.asc
# Output file name                                       : stations.txt
# 
# 2. Run stationMask
### step 4: run the HBV model. There are four sub-routines
### step 4.1: mask
/hdata/fou/personlig/guru/DEW/Dew/Source/pre_3.01b/stationMask control_mask.txt
# output : 
#         stations.txt
#stationMask writes the data in one column.
# ncols 5
# nrows 8 ...
# xllcorner 0 
# yllcorner 0 
# cellsize 500
# NODATA_value -9999 
# 80001 
# 80001 
# -9999 
# -9999 
# -9999 
# 80001
# 8001
# 8001
####

####################################################################
### step 4: run the HBV model. There are four sub-routines
### step 4.2: flow

### step 4: run the HBV model. There are four sub-routines
### step 4.3: predew

# Run predew
####################################################################
# 1. Preparing control_pre.txt
# Model structure, HBV (0) or KinematicWave (1)                                  : 0
# Landscape elements hierarchy, flow direction network(N) or nested catchments(C): C
# Output file name                                                               : pre_out.txt
# File with meteorological stations                                              : met_stations.txt
# File with common parameters                                                    : dew_common_parameters.txt
# File with geographical analysis area                                           : stations.txt
# File with grid cell elevations                                                 : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/lokaldtm100.asc
# File with slope lengths                                                        : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/lokaldtm100.asc
# File with slope angles                                                         : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/lokaldtm100.asc
# File with slope aspects                                                        : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/lokaldtm100.asc
# File with lake percentage                                                      : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/innsjopro.asc
# File with forest percentage                                                    : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/skogarpro.asc
# File with bog percentage                                                       : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/myrarepro.asc
# File with glacier percentage                                                   : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/brearepro.asc
# File with glacier surface elevations                                           : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/lokaldtm100.asc
# File with glacier ice thickness                                                : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/istykkelse.asc
# File with tree levels                                                          : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/tregrense.asc
# File with sub-catchment hierarchy                                              : watershed.txt
# File with flow direction grid for landscape elements                           : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/flowdir.asc
# File with watercourse/sub-catchment identifiers                                : /hdata/fou/personlig/guru/DEW/Dew/Aalfotbreen/Geodata/watershed_area.asc
#------------------------------------------------------------------------------------------
# 2. copy a parameters file dew_common_parameters.txt
# 3. Create file met_stations.txt

cp ../../Results/Jotunheimen/Jotunheimen_Calibration_00/dew_common_parameters.txt .
cp ../../Results/Jotunheimen/Jotunheimen_Calibration_00/met_stations.txt .
cp ../../Results/Jotunheimen/Jotunheimen_Calibration_00/*.tpl .
cp ../../Results/Jotunheimen/Jotunheimen_Calibration_00/*_parameters.txt .
cp ../../Results/Jotunheimen/Jotunheimen_Calibration_00/*_elements.txt .
cp ../../Results/Jotunheimen/Jotunheimen_Calibration_00/catchment_correction.txt .
cp ../../Results/Aafotbreen_Calibration_1km_00/start_dew.bat . 

head -25 watershed_00.txt > subcatchment_output_elements.txt
# 4. Run predew
/hdata/fou/personlig/guru/DEW/Dew/Source/pre_3.01b/predew control_pre.txt
# output:
#        dew_landscape.txt
#        dew_grid_index.txt
#        dew_waterland.txt
#        pre_out.txt
###
####################################################################
# Copy files
####################################################################
# template-files
# kiwa_soil_parameters.txt
# kiwa_elements.txt 
# glacier_retreat_parameters.txt 
# hbv_elements.txt 
# catchment_correction.txt
# create subcatchment_output_elements.txt from watershed.txt
###
####################################################################
### step 4: run the HBV model. There are four sub-routines
### step 4.4: run dew

# Set environment variable
####################################################################
export  METDATA=/hdata/grid/metdata
export  METMASK=/hdata/fou/personlig/guru/DEW/Dew/"Aalfotbreen"/Met

