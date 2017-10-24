#!/bin/csh -f
#
#  April 21, 1999 - David T. Sandwell
#
#  Script to calculate a table of parameters from master and
#  slave PRM files.
#
unset noclobber
#
#  Modified by M.Wei to add ALOS function, 9/27/06
#
 if ($#argv < 2) then
  echo " "
  echo "Usage: baseline_table master.PRM slave.PRM [GMT] "
  echo "               [GMT] creates table for pstext"
  echo " "
  echo " Output:"
  echo " sat_orb slave_time slave_days(1992ERS,2006ALOS) Bpl Bperp xshift yshift"
  echo " "
  exit 1
 endif
#
#  get the time information from the master
#
 set MT0 = `grep SC_clock_start $1 | awk '{print $3}'`
 set MTF = `grep SC_clock_stop $1 | awk '{print $3}'`
 set MSC = `grep SC_identity $1 | awk '{print $3}'`
#
#  get the time information from the slaver
#
 set ST0 = `grep SC_clock_start $2 | awk '{print $3}'`
 set STF = `grep SC_clock_stop $2 | awk '{print $3}'`
 set SSC = `grep SC_identity $2 | awk '{print $3}'`
# 
# convert the start time to days since 1992
#
 @ T0  = `grep SC_clock_start $2 | awk '{print $3}' | awk -F"." '{print $1}'`
 @ DAY = $T0 % 1000
#
if ($SSC == 1 || $SSC == 2) then
 ERS_baseline $1 $2 > temp
 @ YR = $T0  / 1000 - 1992
 @ YDAY = $YR * 365 + $DAY
else if ($SSC == 4) then
 ENVI_baseline $1 $2 > temp
 @ YR = $T0  / 1000 - 1992
 @ YDAY = $YR * 365 + $DAY
else 
 ALOS_baseline $1 $2 > temp
 @ YR = $T0  / 1000 - 2006
 @ YDAY = $YR * 365 + $DAY
endif
#
#  get the needed parameters from temp
#
 set BPL = `grep B_parallel temp | awk '{print $3}'`
 set BPR = `grep B_perpendicular temp | awk '{print $3}'`
 set XS  = `grep xshift temp | awk '{print $3}'`
 set YS  = `grep yshift temp | awk '{print $3}'`
 set NM  = `grep SC_identity $2 | awk '{print $3}'`
if ($SSC == 5) then
 set ORB = `grep input_file $2 | awk '{print $3}' | awk '{print substr($1,14,5)}'` 
else if ($SSC == 4) then
 set ORB = `grep input_file $2 | awk '{print $3}' | cut -c17-21`
else if ($SSC == 1 || $SSC == 2) then
 set ORB = `grep input_file $2 | awk '{print $3}' | cut -c1-8`
endif
#
 if ($#argv < 3) then
    echo $ORB $ST0 $YDAY $BPL $BPR $XS $YS 
 else 
    echo $YDAY $BPR '8' '0.' '0' '5' $ORB
 endif
#
# clean up
#
  rm temp
