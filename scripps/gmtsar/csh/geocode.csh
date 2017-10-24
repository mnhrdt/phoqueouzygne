#!/bin/csh -f
#
#  D. Sandwell FEB 10 20100
#
#
alias rm 'rm -f'
unset noclobber
#
  if ($#argv < 1) then
errormessage:
    echo ""
    echo "Usage: geocode.csh correlation_threshold"
    echo ""
    echo " phase is masked when correlation is less than correlation_threshold"
    echo ""
    echo "Example: geocode.csh .12"
    echo ""
    exit 1
  endif
#
#   first mask the phase and phase gradient using the correlation
#
grdmath corr.grd $1 GE 0 NAN mask.grd MUL = mask2.grd -V
grdmath phase.grd mask2.grd MUL = phase_mask.grd
if (-e xphase.grd) then
  grdmath xphase.grd mask2.grd MUL = xphase_mask.grd
  grdmath yphase.grd mask2.grd MUL = yphase_mask.grd
endif
if (-e unwrap.grd) then
  grdcut mask2.grd `grdinfo unwrap.grd -I-` -Gmask3.grd
  grdmath unwrap.grd mask3.grd MUL = unwrap_mask.grd
endif
if (-e phasefilt.grd) then
  grdmath phasefilt.grd mask2.grd MUL = phasefilt_mask.grd
endif
#
#   look at the masked phase
#
set boundR = `grdinfo display_amp.grd -C | awk '{print ($3-$2)/4}'`
set boundA = `grdinfo display_amp.grd -C | awk '{print ($5-$4)/4}'`
grdimage phase_mask.grd -JX6.5i -Cphase.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > phase_mask.ps
psscale -D3.3/-1.5/5/0.2h -Cphase.cpt -B1.57:"phase, rad": -O >> phase_mask.ps
if (-e xphase_mask.grd) then
  grdimage xphase_mask.grd -JX8i -Cphase_grad.cpt -X.2i -Y.5i -P > xphase_mask.ps
  grdimage yphase_mask.grd -JX8i -Cphase_grad.cpt -X.2i -Y.5i -P > yphase_mask.ps
endif
if (-e unwrap_mask.grd) then
  grdimage unwrap_mask.grd -JX6.5i -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cunwrap.cpt -X1.3i -Y3i -P -K > unwrap_mask.ps
  set std = `grdinfo -C -L2 unwrap_mask.grd | awk '{printf("%5.1f", $13)}'`
  psscale -D3.3/-1.5/5/0.2h -Cunwrap.cpt -B"$std":"unwrapped phase, rad": -O -E >> unwrap_mask.ps
endif
if (-e phasefilt_mask.grd) then
  grdimage phasefilt_mask.grd -JX6.5i -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cphase.cpt -X1.3i -Y3i -P -K > phasefilt_mask.ps
  psscale -D3.3/-1.5/5/0.2h -Cphase.cpt -B1.57:"phase, rad": -O >> phasefilt_mask.ps
endif
# line-of-sight displacement
if (-e unwrap_mask.grd) then
  set wavel = `grep wavelength *.PRM | awk '{print($3)}' | head -1 `
  grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
  grdgradient los.grd -Nt.9 -A0. -Glos_grad.grd
  set tmp = `grdinfo -C -L2 los.grd`
  set limitU = `echo $tmp | awk '{printf("%5.1f", $12+$13*2)}'`
  set limitL = `echo $tmp | awk '{printf("%5.1f", $12-$13*2)}'`
  set std = `echo $tmp | awk '{printf("%5.1f", $13)}'`
  makecpt -Cpolar -Z -T"$limitL"/"$limitU"/1 -D > los.cpt
  grdimage los.grd -Ilos_grad.grd -Clos.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -JX6.5i -X1.3i -Y3i -P -K > los.ps
  psscale -D3.3/-1.5/4/0.2h -Clos.cpt -B"$std":"LOS displacement, mm":/:"range decrease": -O -E >> los.ps
endif
#
#  now reproject the phase to lon/lat space
#
echo "geocode.csh"
echo "project correlation, phase, unwrapped and amplitude back to lon lat coordinates"
proj_ra2ll.csh trans.dat corr.grd corr_ll.grd
proj_ra2ll.csh trans.dat phase_mask.grd phase_mask_ll.grd
proj_ra2ll.csh trans.dat display_amp.grd display_amp_ll.grd
if (-e xphase_mask.grd) then
  proj_ra2ll.csh trans.dat xphase_mask.grd xphase_mask_ll.grd
  proj_ra2ll.csh trans.dat yphase_mask.grd yphase_mask_ll.grd
endif
if (-e unwrap_mask.grd) then
  proj_ra2ll.csh trans.dat unwrap_mask.grd unwrap_mask_ll.grd
endif
if (-e phasefilt_mask.grd) then
  proj_ra2ll.csh trans.dat phasefilt_mask.grd phasefilt_mask_ll.grd
endif
#
#   now image for google earth
#
echo "geocode.csh"
echo "make the KML files for Google Earth"
grd2kml.csh display_amp_ll display_amp.cpt
grd2kml.csh corr_ll corr.cpt
grd2kml.csh phase_mask_ll phase.cpt
#ln -s phasefilt_mask_ll.grd phase_mask_ll_bw.grd
#grd2kml.csh phase_mask_ll_bw phase_bw.cpt
#rm phase_mask_ll_bw.grd
if (-e xphase_mask_ll.grd) then
  grd2kml.csh xphase_mask_ll phase_grad.cpt
  grd2kml.csh yphase_mask_ll phase_grad.cpt
endif
if (-e unwrap_mask_ll.grd) then
  grd2kml.csh unwrap_mask_ll unwrap.cpt
endif
if (-e phasefilt_mask_ll.grd) then
  grd2kml.csh phasefilt_mask_ll phase.cpt
endif
if (-e unwrap_mask_ll.grd) then
  grdmath unwrap_mask_ll.grd $wavel MUL -79.58 MUL = los_ll.grd
  grd2kml.csh los_ll los.cpt
endif
