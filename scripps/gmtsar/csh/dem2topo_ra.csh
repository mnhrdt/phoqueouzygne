#!/bin/csh -f
# Matt WEI Feb 1 2010
# modified by Xiaopeng Feb 9 2010
# modifed by E. Feilding, DST, XT to add TSX data Jan 10 2014
#=======================================================================
#  script to make topography for interferograms 
#  The USGS elevations are height above WGS84 so this is OK.
#
alias rm 'rm -f'
unset noclobber
#
if ($#argv < 2) then
 echo " "
 echo "Usage: dem2topo_ra.csh master.PRM dem.grd" 
 echo " "
 echo "        Note: Works for TSX,ALOS,ERS,ENVISAT"
 echo " "
 exit 1
endif 

#
# local variables 
#
  set scale = -JX8i
#
#========================Mosaic topo data===============================

#-----------------------------------------------------------------------
# 

#------------------------Get bounds in radar coordinates----------------
set XMAX = `grep num_rng_bins $1 | awk '{print $3}'`
set yvalid = `grep num_valid_az $1 | awk '{print $3}'`
set num_patch = `grep num_patches $1 | awk '{print $3}'`
set YMAX = `echo "$yvalid $num_patch" | awk '{print $1*$2}'`
set SC = `grep SC_identity $1 | awk '{print $3}'`
#
# look for range sampling rate
#
  set rng_samp_rate = `grep rng_samp_rate $1 | awk 'NR == 1 {printf("%d", $3)}'`

# set the range spacing of simulation in units of image range pixel size

if($rng_samp_rate > 0 && $rng_samp_rate < 25000000) then
  set rng = 1
else if($rng_samp_rate >= 25000000 && $rng_samp_rate < 72000000 || $SC == 7 ) then
  set rng = 2
else if($rng_samp_rate >= 72000000) then
  set rng = 4
else
   echo "range sampling rate out of bounds"
   exit 0
endif
echo " range decimation is: " $rng

#   use special ALOS_llt2rat if this is ALOS otherwise use SAT_llt2rat

if ($SC == 5) then
  echo " processing for ALOS data"
  grd2xyz --D_FORMAT=%lf $2 -S | ALOS_llt2rat $1 -bod  > trans.dat
else 
  echo " processing generic data"
  grd2xyz --D_FORMAT=%lf $2 -S | SAT_llt2rat $1 -bod  > trans.dat
endif
  gmtconvert trans.dat -F0,1,2 -bid5 -bod3 | blockmedian -R0/$XMAX/0/$YMAX -I$rng/4 -bid -bod -V > temp.rat 
  surface temp.rat -R0/$XMAX/0/$YMAX -I$rng/4 -bid -T.50 -N1000 -Gnode.grd -V
# 
# resample and flip top to bottom for both ascending and descending passes
#  
  grdsample node.grd -T -Gtopo_rat.grd
  grdmath topo_rat.grd FLIPUD = topo_ra.grd
# 
# plotting
# 
  grd2cpt topo_ra.grd -Cgray -V -Z > topo_ra.cpt 
  grdimage topo_ra.grd $scale -X0.2i -P -Ctopo_ra.cpt -V > topo_ra.ps
#
#  clean up
#
rm node.grd temp.rat dem.xyz topo_rat.grd 
rm topo_ra.cpt
