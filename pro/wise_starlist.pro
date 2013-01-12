;+
; NAME:
;   wise_starlist
;
; PURPOSE:
;   Get a list of the stars  in (and just off the edge of) an L1b image
;   that will be PSF subtracted
;
; CALLING SEQUENCE:
;   wise_starlist, h, xlist, ylist, maglist, catra, catdec, catmag, catflag
;
; INPUTS:
;   h      - FITS header containing WCS astrometry
;   catra  - W3 catalog RA values for w3 < 11.5
;   catdec - W3 catalog DEC values for w3 < 11.5
;   catmag - W3 catalog magnitude for w3 < 11.5, from field w3mpro
;   catflag - W3 catalog CC flag, from field w3cc_map
;
; OUTPUTS:
;   xlist   - list of x values of stars to be PSF subtracted
;   ylist   - list of y values of stars to be PSF subtracted
;   maglist - list of magnitudes of stars to be PSF subtracted
;   
; EXAMPLES:
;   cat = mrdfits('$WISE_DATA/w3_catalog.fits')
;   wise_starlist, h, xlist, ylist, maglist, cat.ra, cat.dec, cat.mag, cat.flag
;
; COMMENTS:
;   list stars returned is dependent on PSF parameters, which are hardcoded
;
;   binary_search function should eventually become its own separate routine
;   so that future routines that need it don't have to copy/paste
;
; REVISION HISTORY:
;   2011-Feb-16 - Written by Aaron Meisner
;
;----------------------------------------------------------------------

function binary_search, list, val
  
  n = n_elements(list) 
  step = n/2
  
  ind = n/2

  for bar=1, (alog(n)/alog(2))+1 do begin 
     step = (step/2) > 1
     ind += (list[ind] LT val ? step : -step)
     ind = (ind < (n-1)) > 0
  endfor

  return, ind
end

pro wise_starlist, h, xlist, ylist, maglist, allsky=allsky, $ 
                   mjdlim=mjdlim, wm=m, w4=w4

  par = psf_par_struc(allsky=allsky, w4=w4, /everything)
  IMPIX  = par.impix
  NPIX_BRIGHT = par.psfpix
  NPIX_FAINT = par.pfaint

  init_source_catalog, w4=w4, allsky=allsky
  COMMON CATALOG, catra, catdec, catmag, mjdmin, mjdmax, catm

; ----- eventually use psf_par_struc.pro to define pscl
  pscl = keyword_set(w4) ? 5.53 : 2.75 ; asec/pixel
  maxsep = sqrt(2)*(IMPIX/2 + NPIX_BRIGHT/2)*pscl/3600. ;deg
  
  expra  = sxpar(h, 'CRVAL1')
  expdec = sxpar(h, 'CRVAL2')
  indmin = binary_search(catdec, expdec - maxsep)
  indmax = binary_search(catdec, expdec + maxsep)
  if (indmax EQ indmin) then return
  wdec   = lindgen(indmax - indmin + 1)+indmin

  dangle = djs_diff_angle(catra[wdec], catdec[wdec], expra, expdec)
  
  wclose = where(dangle LE maxsep, nclose)
  if (nclose EQ 0) then return
  
  ra     = catra[wdec[wclose]]
  dec    = catdec[wdec[wclose]]
  mag    = catmag[wdec[wclose]]
  mjdl   = mjdmin[wdec[wclose]]
  mjdu   = mjdmax[wdec[wclose]]
  m      = catm[wdec[wclose]]
  
  bright = keyword_set(w4) ?  $ 
      (mag LT ((-0.2+2.5*alog10(sqrt(m > 1)) > 1.2) < 2)) : $
      (mag LT ((par.bpad+2.6+2.5*alog10(sqrt(m>1)) > 4) < 6))
  
  extast, h, astr
  ad2xy, ra, dec, astr, x, y
  ix     = round(x)
  iy     = round(y)
  
  wfaint = where((~bright) AND (ix GE -NPIX_FAINT/2) $ 
    AND (ix LT IMPIX+NPIX_FAINT/2) AND (iy GE -NPIX_FAINT/2) $ 
    AND (iy LT IMPIX+NPIX_FAINT/2),nfaint)

  wbright = where((bright) AND (ix GE -NPIX_BRIGHT/2) $ 
    AND (ix LT IMPIX+NPIX_BRIGHT/2) AND (iy GE -NPIX_BRIGHT/2) $ 
    AND (iy LT IMPIX+NPIX_BRIGHT/2), nbright)

   if ((nbright+nfaint) EQ 0) then return

   wsub = (nbright GT 0) ? ((nfaint GT 0) ? [wbright, wfaint] : wbright) $ 
     : wfaint
   xlist = x[wsub]
   ylist = y[wsub]
   maglist = mag[wsub]
   mjdlim = [[mjdl[wsub]],[mjdu[wsub]]]
   m = m[wsub]

end
