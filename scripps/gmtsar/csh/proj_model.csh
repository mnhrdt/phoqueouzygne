#!/bin/csh -f

# Creator: Xiaopeng Tong, David Sandwell 
# Date:    July 6th 2013

alias rm 'rm -f'

if ($#argv != 7) then
  echo ""
  echo "Usage: proj_model.csh SAT ve.grd vn.grd vu.grd master.PRM dem.grd los.grd"
  echo ""
  echo "  project a simulated crust motion model into radar look directions "
  echo "  the radar look direction is spatial variable"
  echo ""
  echo "  SAT             --  can be ALOS or ERS or ENVI or generic SAT"
  echo "  ve, vn, vu.grd  --  input: east, north, vertical motion simulated from numerical model "
  echo "  master.PRM      --  PRM file of the radar image"
  echo "  dem.grd         --  DEM grid (wider than the coverage of the radar image)"
  echo "  los.grd         --  output: file name of the LOS grid "
  echo ""
  echo "  note that the grid of the grd files must be consistent" 
  echo "  note the program need master.PRM and its LED file"
  echo ""
  exit 1
endif

if ((! -e $2) || (! -e $3) || (! -e $4)) then
  echo ""
  echo "no input model file found: $2, $3, $4"
  echo ""
  exit 1
endif

if (! -e $5) then
  echo ""
  echo "no input PRM file found: $5"
  echo ""
  exit 1
endif

if (! -e $6) then
  echo ""
  echo "no input DEM grid found: $6"
  echo ""
  exit 1
endif

set SAT = $1
if ($SAT != ENVI && $SAT != ERS && $SAT != ALOS && $SAT != SAT) then
    echo ""
    echo " SAT can be ALOS ENVISAT or ERS or generic SAT"
    echo ""
    exit 1
endif

set ve = $2
set vn = $3
set vu = $4
set prm = $5
set dem = $6
set out = $7

set int = `grdinfo -I $ve`

grd2xyz $dem | $SAT"_look" $prm > look.xyz
awk '{print $1,$2,$4}' look.xyz > look.lle    
awk '{print $1,$2,$5}' look.xyz > look.lln    
awk '{print $1,$2,$6}' look.xyz > look.llu    

blockmedian look.lle `grdinfo -I- $dem` $int -bo | surface -bi -Gtmp.grd $int `grdinfo -I- $dem` -T0.5 
grdsample tmp.grd -T -Glle.grd
blockmedian look.lln `grdinfo -I- $dem` $int -bo | surface -bi -Gtmp.grd $int `grdinfo -I- $dem` -T0.5
grdsample tmp.grd -T -Glln.grd
blockmedian look.llu `grdinfo -I- $dem` $int -bo | surface -bi -Gtmp.grd $int `grdinfo -I- $dem` -T0.5
grdsample tmp.grd -T -Gllu.grd

set region = `grdinfo -I- llu.grd`
#echo $region
grdsample $ve $region $int -Gtmpve.grd 
grdsample $vn $region $int -Gtmpvn.grd
grdsample $vu $region $int -Gtmpvu.grd
grdmath tmpve.grd "lle.grd" MUL tmpvn.grd "lln.grd" MUL ADD tmpvu.grd "llu.grd" MUL ADD = $out 


#clean up
rm look.xyz look.lle look.lln look.llu tmp.grd
rm lle.grd llu.grd lln.grd tmpve.grd tmpvn.grd tmpvu.grd 

