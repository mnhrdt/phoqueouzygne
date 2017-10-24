#!/bin/csh -f 
#  D. Sandwell 1/12/07
#
alias rm 'rm -f'
unset noclobber
#
#
#  project an ASCII file from range/azimuth coordinates into lon/lat coordinates
#  this version only works with GMT V4.0 and higher
#
#  Input:
#  trans.dat    - file generated by llt_grid2rat  (r a topo lon lat)
#  phase_ra.txt - an ASCII file of phase or anything
#
#  Output:
#  phase_ll.txt - an ASCII file of phase in longitude/latitude coordinates
#
# check for number of arguments
#
 if ($#argv < 3) then
  echo " "
  echo "Usage: proj_ra2ll_ascii.csh trans.dat phase.txt phase_ll.txt" 
  echo "        trans.dat    - file generated by SAT_llt2rat  (r a topo lon lat)"
  echo "        phase_ra.txt - an ASCII file of phase or anything" 
  echo "        phase_ll.txt - output file in lon/lat-coordinates" 
  echo " "
  exit 1
 endif 
#
#   make grids of longitude and latitude versus range and azimuth unless they already exist
#
if (! -f raln.grd || ! -f ralt.grd ) then
  gmtconvert $1 -F0,1,3 -bid5 -bos3 > raln
  gmtconvert $1 -F0,1,4 -bid5 -bos3 > ralt
# awk '{ if ($4 > 180.) printf("%f %f %e \n",$1,$2,$4-360.); else printf("%f %f %e \n",$1,$2,$4) }' < $1 | gmtconvert -bos3 > raln
# awk '{ printf("%f %f %e \n",$1,$2,$5) }' < $1 | gmtconvert -bos3 > ralt
#
surface raln `minmax $2 -I16/32` -bis3 -I8/32 -T.50 -Graln.grd -V
surface ralt `minmax $2 -I16/32` -bis3 -I8/32 -T.50 -Gralt.grd -V
endif
#
grdtrack $2 -Q -Graln.grd > rapln
grdtrack rapln -Q -Gralt.grd > raplnlt
#
# get the lon, lat, phase columns and grid
#
awk '{print $4,$5,$3}' < raplnlt > $3
#
# clean
#
rm -f rap* llp llpb