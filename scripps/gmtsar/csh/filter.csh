#!/bin/csh -f
#
#  Xiaopeng Tong and David Sandwell 
#  FEB 4 2010
#  Matt Wei May 4 2010, ENVISAT
#  DTS - May 26, 2010, added phase gadient
#  EF, DTS, XT - Jan 10 2014, TSX
# Convolve the real.grd and imag.grd with gauss filters. 
# Form amplitude, phase, phase gradient, and correlation images. 
#
#
alias rm 'rm -f'
#
#
# set grdimage options
#
  set scale = "-JX6.5i"
  set thresh = "5.e-21"
  gmtset COLOR_MODEL = hsv
  gmtset MEASURE_UNIT = inch

  if ($#argv != 4) then
errormessage:
    echo ""
    echo "Usage: filter.csh master.PRM slave.PRM filter decimation"
    echo ""
    echo " Apply gaussian filter to amplitude and phase images."
    echo " The gaussian filter and decimation for the amplitude and phase images are the same."
    echo " filter is the  name of the filter."
    echo " decimation control the size of the amplitude and phase images. It is either 1 or 2."
    echo " Set the decimation to be 1 if you want higher resolution images."
    echo " Set the decimation to be 2 if you want images with smaller file size."
    echo " "
    echo "Example: filter.csh IMG-HH-ALPSRP055750660-H1.0__A.PRM IMG-HH-ALPSRP049040660-H1.0__A.PRM gauss_alos_300m 2"
    echo ""
    exit 1
  endif
#
# define filter and decimation variables
#
  set filter2 = $GMTSAR/filters/$3
  set filter3 = $GMTSAR/filters/fill.3x3
  set filter4 = $GMTSAR/filters/xdir
  set filter5 = $GMTSAR/filters/ydir
  set dec  = $4
  set az_lks = 4 
  set PRF = `grep PRF *.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  if( $PRF < 1000 ) then
     set az_lks = 1
  endif
  echo "az_lks = " $az_lks
#
# look for range sampling rate
#
  set rng_samp_rate = `grep rng_samp_rate $1 | awk 'NR == 1 {printf("%d", $3)}'`

# set the range spacing in units of image range pixel size
  if ($?rng_samp_rate) then
    if ($rng_samp_rate > 110000000) then 
      set dec_rng = 4
      set filter1 = $GMTSAR/filters/gauss15x5
    else if ($rng_samp_rate < 110000000 && $rng_samp_rate > 20000000) then
      set dec_rng = 2
      set filter1 = $GMTSAR/filters/gauss15x5
    else  
      set dec_rng = 1
      set filter1 = $GMTSAR/filters/gauss15x3
    endif
  else
    echo "Undefined rng_samp_rate in the master PRM file"
    exit 1
  endif
#
# filter the two amplitude images
#
  echo "filter.csh"
  echo "making amplitudes..."
  conv $az_lks $dec_rng $filter1 $1 amp1_tmp.grd
  conv $dec $dec $filter2 amp1_tmp.grd amp1.grd
  rm amp1_tmp.grd
  conv $az_lks $dec_rng $filter1 $2 amp2_tmp.grd
  conv $dec $dec $filter2 amp2_tmp.grd amp2.grd
  rm amp2_tmp.grd
#
# filter the real and imaginary parts of the interferogram
# also compute gradients
#
  echo "filtering interferogram..."
  conv $az_lks $dec_rng $filter1 real.grd real_tmp.grd
  conv $dec $dec $filter2 real_tmp.grd realfilt.grd
#  conv $dec $dec $filter4 real_tmp.grd xreal.grd
#  conv $dec $dec $filter5 real_tmp.grd yreal.grd
  rm real_tmp.grd real.grd
  conv $az_lks $dec_rng $filter1 imag.grd imag_tmp.grd
  conv $dec $dec $filter2 imag_tmp.grd imagfilt.grd
#  conv $dec $dec $filter4 imag_tmp.grd ximag.grd
#  conv $dec $dec $filter5 imag_tmp.grd yimag.grd
  rm imag_tmp.grd imag.grd
#
# form amplitude image
#
  echo "making amplitude..."
  grdmath realfilt.grd imagfilt.grd HYPOT  = amp.grd 
  grdmath amp.grd 0.5 POW FLIPUD = display_amp.grd 
  set AMAX = `grdinfo -L2 display_amp.grd | grep stdev | awk '{ print 4*$5 }'`
  set boundR = `grdinfo display_amp.grd -C | awk '{print ($3-$2)/4}'`
  set boundA = `grdinfo display_amp.grd -C | awk '{print ($5-$4)/4}'`
  grd2cpt display_amp.grd -Z -D -L0/$AMAX -Cgray > display_amp.cpt
  echo "N  255   255   254" >> display_amp.cpt
  grdimage display_amp.grd -Cdisplay_amp.cpt $scale -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3 -Y3 -P > display_amp.ps
#
# form the correlation
#
  echo "making correlation..."
  grdmath amp1.grd amp2.grd MUL = tmp.grd
  grdmath tmp.grd $thresh GE 0 NAN = mask.grd
  grdmath amp.grd tmp.grd SQRT DIV mask.grd MUL FLIPUD = tmp2.grd=bf
  conv 1 1 $filter3 tmp2.grd corr.grd
  makecpt -T0./.5/0.1 -Cgray -Z -N > corr.cpt
  echo "N  255   255   254" >> corr.cpt
  grdimage corr.grd $scale -Ccorr.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3 -Y3 -P -K > corr.ps
  psscale -D3.3/-1.5/5/0.2h -Ccorr.cpt -B0.2:correlation: -O -E >> corr.ps
#
# form the phase 
#
  echo "making phase..."
  grdmath imagfilt.grd realfilt.grd ATAN2 mask.grd MUL FLIPUD = phase.grd
  makecpt -Crainbow -T-3.15/3.15/0.1 -Z -N > phase.cpt
# makecpt -Cgray -T-3.14/3.14/0.1 -Z -N > phase_bw.cpt
# echo "N  255   255   254" >> phase_bw.cpt
  grdimage phase.grd $scale -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cphase.cpt -X1.3 -Y3 -P -K > phase.ps
  psscale -D3.3/-1.5/5/0.2h -Cphase.cpt -B1.57:"phase, rad": -O >> phase.ps

#
#  make the Werner/Goldstein filtered phase
#
  echo "filtering phase..."
# phasefilt -imag imagfilt.grd -real realfilt.grd -amp1 amp1.grd -amp2 amp2.grd -psize 16 -complex_out
  phasefilt -imag imagfilt.grd -real realfilt.grd -amp1 amp1.grd -amp2 amp2.grd -psize 16 
  grdedit filtphase.grd `grdinfo mask.grd -I- --D_FORMAT=%.12lg` 
  grdmath filtphase.grd mask.grd MUL FLIPUD = phasefilt.grd
  rm filtphase.grd
  grdimage phasefilt.grd $scale -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cphase.cpt -X1.3 -Y3 -P -K > phasefilt.ps
  psscale -D3.3/-1.5/5/0.2h -Cphase.cpt -B1.57:"phase, rad": -O >> phasefilt.ps
# grdimage phasefilt.grd $scale -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cphase_bw.cpt -X1.3 -Y3 -P -K > phase_bw.ps
# psscale -D3.3/-1.5/5/0.2h -Cphase_bw.cpt -B1.57:"phase, rad": -O >> phase_bw.ps

# 
#  form the phase gradients
#
#  echo "making phase gradient..."
#  grdmath amp.grd 2. POW = amp_pow.grd
#  grdmath realfilt.grd ximag.grd MUL imagfilt.grd xreal.grd MUL SUB amp_pow.grd DIV mask.grd MUL FLIPUD = xphase.grd
#  grdmath realfilt.grd yimag.grd MUL imagfilt.grd yreal.grd MUL SUB amp_pow.grd DIV mask.grd MUL FLIPUD = yphase.grd 
#  makecpt -Cgray -T-0.7/0.7/0.1 -Z -N > phase_grad.cpt
#  echo "N  255   255   254" >> phase_grad.cpt
#  grdimage xphase.grd $scale -Cphase_grad.cpt -X.2 -Y.5 -P > xphase.ps
#  grdimage yphase.grd $scale -Cphase_grad.cpt -X.2 -Y.5 -P > yphase.ps
#
  mv mask.grd tmp.grd 
  grdmath tmp.grd FLIPUD = mask.grd
#
# delete files
 rm tmp.grd tmp2.grd ximag.grd yimag.grd xreal.grd yreal.grd 
