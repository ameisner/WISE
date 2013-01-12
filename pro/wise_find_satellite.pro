;+
; NAME:
;   wise_find_satellite
;
; PURPOSE:
;   identify satellite-like objects in an L1b exposure
;
; CALLING SEQUENCE:
;   sat = wise_find_satellite(raw, msk)
;
; INPUTS:
;   raw - raw 1016x1016 (W3) L1b image
;   msk - corresponding L1b image mask
;   
; OUTPUTS:
;   sat - bitmask marking satellite locations (0=no satellite,1=satellite)
;
; EXAMPLES:
;   see routines sat_test and t203 below
;
; COMMENTS:
;   I can imagine many approaches to finding satellites..here I'm
;   trying to identify satelites based on morphology without reference
;   to time-series information, whereas a different approach would be
;   to look for anomalies in the ISSA tile artifact masks and then map 
;   back to L1b exposures
;
;   one L1b exposure takes ~0.3 seconds to process
;
; REVISION HISTORY:
;   2012-Oct-14 - Aaron Meisner
;----------------------------------------------------------------------
function wise_find_satellite, raw, msk

  binfac = 4
  medsz = 13
  IMPIX = 1016 ; W3 L1b
  bthresh = 2000 ; DN
  vthresh = 1000
  ;vthresh = 5
  trim = (medsz/2)*binfac+binfac


  xbox = lindgen(IMPIX,IMPIX) MOD IMPIX
  ybox = lindgen(IMPIX,IMPIX) / IMPIX

  im = wise_l1b_maskinterp(raw, msk)
  reb = rebin(im, IMPIX/binfac, IMPIX/binfac)
  med = median(reb, medsz) ; outer medsz/2 pixels left unchanged
  med = rebin(med, IMPIX, IMPIX)
; --- zero out bad regions near edges
  med *= (xbox GE trim) AND (xbox LE (IMPIX-trim)) AND $ 
         (ybox GE trim) AND (ybox LE (IMPIX-trim))
  im  *= (xbox GE trim) AND (xbox LE (IMPIX-trim)) AND $ 
         (ybox GE trim) AND (ybox LE (IMPIX-trim))
  brt = (im-med) GT bthresh

  wbrt = where(brt, nbrt)

  if (nbrt EQ 0) then return, -1
  sat = bytarr(IMPIX, IMPIX)


  count = 0L
  while (total(brt) GT 0) do begin 
    print, 'count = ', count
    xbrt = wbrt MOD IMPIX
    ybrt = wbrt / IMPIX
    xstart = xbrt[0]
    ystart = ybrt[0]
    brt_sv = brt
    flood_fill, xstart, ystart, brt
    contiguous = where(brt NE brt_sv, nc)
    if (nc LE 2) then GOTO, ENDLOOP
; ----- this approach will have trouble if all x values same (vertical line)
    xc = contiguous MOD IMPIX
    yc = contiguous / IMPIX
    A = fltarr(2, nc)
    A[0,*] = float(xc)
    A[1,*] = 1
    YY = float(yc)
    coeff = invert(transpose(A) ## A) ## (transpose(A) ## YY)
    theta = atan(coeff[0]) ; rad
 ;   rot = [[cos(theta),-sin(theta)],[sin(theta), cos(theta)]]
    u1 = cos(theta)*xc+sin(theta)*yc
    u2 = -sin(theta)*xc+cos(theta)*yc
    var1 = variance(u1)
    var2 = variance(u2)
    print, var1, var2, var1/var2, n_elements(xc), $ 
        (var1 GT vthresh*var2), mean(xc), mean(yc)
    if (var1 GT vthresh*var2) then sat[contiguous] = 1
    ENDLOOP: 
    wbrt = where(brt, nbrt)
    count++
  endwhile
  ; do linear fit to x vs. y
  ; project onto linear fit line, and perpendicular to it
  ; if ratio of these projections greater than THRESH then it is satellite
  
  return, sat
end

pro sat_test, satmsk

; ----- this image has one obvious satellite
  fname = '/n/home08/dfink/wisedata/downloads/4band_p1bm_frm/7a/01557a/101/01557a101-w3-int-1b.fits'

; ----- two obvious satellites and a bright star
  fname = '/n/home08/dfink/wisedata/downloads/4band_p1bm_frm/2b/01542b/124/01542b124-w3-int-1b.fits'
  raw = readfits(fname)
  mname = repstr(fname, 'int', 'msk') + '.gz'
  msk = readfits(mname)
  t0 = systime(1)
  satmsk = wise_find_satellite(raw, msk)
  print, 'dt = ', systime(1)-t0

end

pro t203, indstart, nproc, outpath=outpath

  if ~keyword_set(outpath) then outpath = '/n/panlfs/ameisner/satellites'
; test all ~7000 images in tile , should be 6891 images

  tnum = 203
  tstr = mrdfits('$WISE_DATA/wisetile-index-allsky.fits', 1)
  euler, tstr[tnum-1].ra, tstr[tnum-1].dec, lgal, bgal, 1
  lb = [lgal, bgal]
  ind = wise_index_metadata(lb, /allsky, angle=9.)
  nfile = n_elements(ind)

  indend = (long(indstart)+nproc-1) < (nfile-1)

  for i=long(indstart), indend do begin
      fname = ind[i].fname
      print, i, '   ', fname
      raw = readfits(fname, /silent, h)
      mname = repstr(fname, 'int', 'msk') + '.gz'
      msk = readfits(mname, /silent)
      t0 = systime(1)
      satmsk = wise_find_satellite(raw, msk)
      print, 'dt = ', systime(1)-t0
      outname = concat_dir(outpath, fileandpath(fname))
      writefits, outname, satmsk, h
  endfor

end
