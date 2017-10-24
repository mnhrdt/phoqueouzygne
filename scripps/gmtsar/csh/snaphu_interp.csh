#!/bin/csh -f
#
#
alias rm 'rm -f'
unset noclobber
#
  if ($#argv < 2) then
errormessage:
    echo ""
    echo "snaphu_interp.csh [GMT5SAR] - Unwrap the phase, after masking and nearest-neighbor interpolating over decorrelated areas"
    echo " "
    echo "This script is an alternative to snaphu.csh, and should run significantly faster and produce fewer unwrapping errors."
    echo " "
    echo "Usage: snaphu_interp.csh correlation_threshold maximum_discontinuity [<rng0>/<rngf>/<azi0>/<azif>]"
    echo ""
    echo "       correlation is reset to zero when < threshold"
    echo "       maximum_discontinuity enables phase jumps for earthquake ruptures, etc."
    echo "       set maximum_discontinuity = 0 for continuous phase such as interseismic "
    echo ""
    echo "Example: snaphu_interp.csh .12 40 1000/3000/24000/27000"
    echo ""
    echo "References:"
    echo "Unwrapping:"
    echo "Chen C. W. and H. A. Zebker, Network approaches to two-dimensional phase unwrapping: intractability and two new algorithms, Journal of the Optical Society of America A, vol. 17, pp. 401-414 (2000)."
    echo "Interpolation:" 
    echo "Agram, P. S. and H. A. Zebker, Sparse Two-Dimensional Phase Unwrapping Using Regular-Grid Methods, IEEE GEOSCIENCE AND REMOTE SENSING LETTERS, VOL. 6, NO. 2, APRIL 2009."
    exit 1
  endif
#
# prepare the files adding the correlation mask
#
if ($#argv == 3 ) then
   grdcut mask.grd -R$3 -Gmask_patch.grd
   grdcut corr.grd -R$3 -Gcorr_patch.grd
   grdcut phasefilt.grd -R$3 -Gphase_patch.grd
else
   ln -s mask.grd mask_patch.grd
   ln -s corr.grd corr_patch.grd
   ln -s phasefilt.grd phase_patch.grd
endif
#
# create landmask
#
if (-e landmask_ra.grd) then
  if ($#argv == 3 ) then 
    grdsample landmask_ra.grd -R$3 `grdinfo -I phase_patch.grd` -Glandmask_ra_patch.grd
  else 
    grdsample landmask_ra.grd `grdinfo -I phase_patch.grd` -Glandmask_ra_patch.grd
  endif
  grdmath phase_patch.grd landmask_ra_patch.grd MUL = phase_patch.grd -V
endif
#
# user defined mask 
#
if (-e mask_def.grd) then
  if ($#argv == 3 ) then
    grdcut mask_def.grd -R$3 -Gmask_def_patch.grd
  else
    cp mask_def.grd mask_def_patch.grd
  endif
  grdmath corr_patch.grd mask_def_patch.grd MUL = corr_patch.grd -V
endif
#
#convert files for input to snaphu
#
grdmath corr_patch.grd $1 GE 0 NAN mask_patch.grd MUL = mask2_patch.grd
grdmath corr_patch.grd 0. XOR 1. MIN  = corr_patch.grd
grdmath mask2_patch.grd corr_patch.grd MUL = corr_tmp.grd 
#use the interpolated phase for snaphu below (old: grd2xyz phase_patch.grd -ZTLf -N0 > phase.in)
grd2xyz corr_tmp.grd -ZTLf  -N0 > corr.in

#
# nearest neighbor interpolation with GDAL
#

#mask the phase by correlation
grdmath mask2_patch.grd phase_patch.grd MUL = phase_patch_mask.grd

# basenames of files
set in = 'phase_patch_mask'
set out = 'phase_patch_interp'

# get x,y bounds
set minx = `grdinfo -C $in.grd |cut -f 2`
set maxx = `grdinfo -C $in.grd |cut -f 3`
set nx = `grdinfo -C $in.grd |cut -f 10`
set boundsx = "$minx $maxx"
set miny = `grdinfo -C $in.grd |cut -f 4`
set maxy = `grdinfo -C $in.grd |cut -f 5`
set ny = `grdinfo -C $in.grd |cut -f 11`
# for some reason we have to reverse these two
set boundsy = "$maxy $miny"

# first convert to ascii
grd2xyz $in.grd -S -V > $in.gmt

# run gdal, then convert back to grd
gdal_grid -of GTiff -txe $boundsx -tye $boundsy -outsize $nx $ny -l $in -a nearest $in.gmt $out.tiff
gdal_translate -of GMT -ot Float32 $out.tiff $out.grd

# fix the grd header metadata
grdedit $out.grd -T
grdedit $out.grd -R$minx/$maxx/$miny/$maxy

# create the snaphu input file
grd2xyz phase_patch_interp.grd -ZTLf -N0 > phase.in

# we could do all of the above with GMT in one line (but it runs incredibly slow, possibly forever):
#grd2xyz $in.grd -bo -S | nearneighbor -bi -G${out}_gmt.grd -R$in.grd -S50 -N1


#
# run snaphu
#

echo "unwrapping phase with snaphu - higher threshold for faster unwrapping "

if ($2 == 0) then
  snaphu phase.in `grdinfo -C phase_patch.grd | cut -f 10` -f $GMTSAR/csh/snaphu.conf.brief -c corr.in -o unwrap.out -v -s
else
  sed "s/.*DEFOMAX_CYCLE.*/DEFOMAX_CYCLE  $2/g" $GMTSAR/csh/snaphu.conf.brief > snaphu.conf.brief
  snaphu phase.in `grdinfo -C phase_patch.grd | cut -f 10` -f snaphu.conf.brief -c corr.in -o unwrap.out -v -d
endif
#
#snaphu phase.in -s `grdinfo -C phase_patch.grd | cut -f 10` -f $GMTSAR/csh/snaphu.conf.brief -c corr.in -o unwrap.out -v -d 
#snaphu phase.in -s `grdinfo -C phase_patch.grd | cut -f 10` -f $GMTSAR/csh/snaphu.conf.brief -c corr.in -o unwrap.out -v -d --nproc 4 --tile 2 2 500 1000
#
# convert to grd
#
xyz2grd unwrap.out -ZTLf -F `grdinfo -I- phase_patch.grd` `grdinfo -I phase_patch.grd` -Gtmp.grd
grdmath tmp.grd mask2_patch.grd MUL = tmp.grd
#
# detrend the unwrapped phase. currently commented out
#
#grdtrend tmp.grd -N3r -Dunwrap.grd
mv tmp.grd unwrap.grd
#
# landmask
if (-e landmask_ra.grd) then
  grdmath unwrap.grd landmask_ra_patch.grd MUL = tmp.grd -V
  mv tmp.grd unwrap.grd
endif
#
# user defined mask
#
if (-e mask_def.grd) then
  grdmath unwrap.grd mask_def_patch.grd MUL = tmp.grd -V
  mv tmp.grd unwrap.grd
endif
#
#  plot the unwrapped phase
#
grdgradient unwrap.grd -Nt.9 -A0. -Gunwrap_grad.grd
set tmp = `grdinfo -C -L2 unwrap.grd`
set limitU = `echo $tmp | awk '{printf("%5.1f", $12+$13*2)}'`
set limitL = `echo $tmp | awk '{printf("%5.1f", $12-$13*2)}'`
set std = `echo $tmp | awk '{printf("%5.1f", $13)}'`
makecpt -Cseis -I -Z -T"$limitL"/"$limitU"/1 -D > unwrap.cpt
set boundR = `grdinfo unwrap.grd -C | awk '{print ($3-$2)/4}'`
set boundA = `grdinfo unwrap.grd -C | awk '{print ($5-$4)/4}'`
grdimage unwrap.grd -Iunwrap_grad.grd -Cunwrap.cpt -JX6.5i -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > unwrap.ps
psscale -D3.3/-1.5/5/0.2h -Cunwrap.cpt -B"$std":"unwrapped phase, rad": -O -E >> unwrap.ps
#
# clean up
#
rm tmp.grd corr_tmp.grd unwrap.out tmp2.grd unwrap_grad.grd 
rm phase.in corr.in 
#
#   cleanup more
#
rm wrap.grd corr_patch.grd phase_patch.grd mask_patch.grd mask2_patch.grd
rm phase_patch_interp.tiff phase_patch_mask.gmt phase_patch_mask.grd phase_patch_interp.grd

#rm wrap.grd corr_patch.grd phase_patch.grd mask_patch.grd mask2_patch.grd mask3.grd mask3.out
#

