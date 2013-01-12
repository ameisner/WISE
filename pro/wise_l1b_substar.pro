;+
; NAME:
;  wise_l1b_substar
;   
; PURPOSE:
;  PSF subtract stars from L1b image, interpolate over the core residuals
;
; CALLING SEQUENCE:  
;  clean = wise_l1b_substar(intim, h, starmask=)
;
; INPUTS:
;  intim  - image with stars
;  h      - FITS header containing WCS astrometry
;
; OUTPUTS:
;   clean - image with stars PSF subtracted
;
; OPTIONAL OUTPUTS:
;   starmask - bit mask identifying saturated, bright, 
;              ghost-contaminated pixels and those pixels which
;              this routine interpolates over, see header of
;              wise_l1b_clean for bit details
;
; REVISION HISTORY:
;   2012-Feb-16 - Written by Aaron Meisner
;----------------------------------------------------------------------
function wise_l1b_substar, intim, h, msk, starmask=starmask, $ 
                           allsky=allsky, w4=w4, srcstr=srcstr

  clean=intim ;don't want to modify intim

  par = psf_par_struc(allsky=allsky, w4=w4, /everything)
  IMPIX = par.impix
  NPIX_BRIGHT = par.psfpix
  NPIX_FAINT = par.pfaint
  INPIX = par.szwings ; central non-zero portion of bright star PSF
  INHALF = INPIX/2
  GOFFS = -par.ygoffs ; sign convention semantics...
  GX = par.xgpix
  GY = par.ygpix
; ----- faint limit for L1b core interpolation
  interp_max = keyword_set(w4) ? 4.3  : 9.
  INTERP_SIZE = [23, 19, 15, 15, 11, 11, 11, 9, 7, 7, 5]
  MAGZP = par.magzp

  if ~keyword_set(srcstr) then $ 
      wise_starlist, h, xstar, ystar, maglist, allsky=allsky, w4=w4, $ 
          mjdlim=mjdlim,  wm=m $ 
  else begin
      xstar = srcstr.xlist
      ystar = srcstr.ylist
      maglist = srcstr.maglist
      mjdlim = srcstr.mjdlim
      m = srcstr.m
  endelse

  starmask = bytarr(IMPIX, IMPIX)
  if n_elements(mjdlim) EQ 0 then return, clean
  mjdrange = mjdlim[*, 1] - mjdlim[*, 0]
  multi_epoch = mjdrange GT 60. ; ~2 months
  bright = keyword_set(w4) ?  $ 
      (maglist LT ((-0.2+2.5*alog10(sqrt(m > 1)) > 1.2) < 2)) : $ 
      (maglist LT ((par.bpad+2.6+2.5*alog10(sqrt(m > 1)) > 4) < 6))
  flux = 10^((MAGZP-maglist)/2.5) ; 1 DN = MAGZP w?mag
  NSTAR=n_elements(xstar)

  t0=systime(1)
  for i=0L, NSTAR-1 do begin
    PSFPIX = bright[i] ? NPIX_BRIGHT : NPIX_FAINT
    PSFHALF=PSFPIX/2

    thismag = maglist[i]
    thisx = xstar[i]
    thisy = ystar[i]
    
    ix = round(thisx)
    iy = round(thisy)
    xstar[i] = ix
    ystar[i] = iy
    
    psf = wise_psf_cutout(thisx, thisy, psf_coeff=psf_coeff, bright=bright[i],$
                          allsky=allsky, w4=w4, flux=flux[i])

    if (bright[i] AND ~keyword_set(w4)) then begin  
    psf[(PSFHALF-INHALF):(PSFHALF+INHALF),(PSFHALF-INHALF):(PSFHALF+INHALF)]= $
      sshift2d(psf[(PSFHALF-INHALF):(PSFHALF+INHALF), $ 
                   (PSFHALF-INHALF):(PSFHALF+INHALF)],[thisx-ix,thisy-iy])
    endif else begin 
  if (thismag LT (par.interplim+(2.5*alog10(sqrt((m[i]> 1)/16.)) > 0))) then $ 
      psf = sshift2d(psf, [thisx-ix,thisy-iy])
    endelse

    clean[((ix-PSFHALF) > 0):((ix+PSFHALF) < (IMPIX-1)), $ 
    ((iy-PSFHALF) > 0):((iy+PSFHALF) < (IMPIX-1))] -= $ 
    psf[((PSFHALF-ix) > 0):((IMPIX-ix+PSFHALF-1) < (PSFPIX-1)), $ 
     ((PSFHALF-iy) > 0):((IMPIX-iy+PSFHALF-1) < (PSFPIX-1))]

    ; make interpolation mask for subtracted cores
    if (thismag LT interp_max) then begin
      masksize = keyword_set(w4) ? INTERP_SIZE[floor(thismag+6.7) > 0] : $ 
                                   INTERP_SIZE[(floor(thismag)+2) > 0]
      if (ix GE -masksize/2) AND (ix LT IMPIX+masksize/2) AND $ 
        (iy GE -masksize/2) AND (iy LT IMPIX+masksize/2) then begin
        starmask[((ix-masksize/2) > 0):((ix+masksize/2) < (IMPIX-1)), $ 
          ((iy-masksize/2) > 0):((iy+masksize/2) < (IMPIX-1))] OR= 16
      endif
    endif
 
    ; mask more things for bright stars
    if (bright[i] OR (keyword_set(w4) AND (thismag LT 1.6)))then begin
      thismask = wise_brightmask(psf, multi_epoch=multi_epoch[i], w4=w4, $ 
                                 faint=(~bright[i]))
      starmask[((ix-PSFHALF) > 0):((ix+PSFHALF) < (IMPIX-1)), $ 
        ((iy-PSFHALF) > 0):((iy+PSFHALF) < (IMPIX-1))] OR= $ 
        thismask[((PSFHALF-ix) > 0):((IMPIX-ix+PSFHALF-1) < (PSFPIX-1)), $ 
        ((PSFHALF-iy) > 0):((IMPIX-iy+PSFHALF-1) < (PSFPIX-1))]
    endif
   endfor
  
  intx = djs_maskinterp(clean,((starmask AND 16) NE 0) OR (msk),iaxis=0,/const)
  inty = djs_maskinterp(clean,((starmask AND 16) NE 0) OR (msk),iaxis=1,/const)
  clean = (intx+inty)/2

  print, 'star subtraction time: ', systime(1)-t0, ' s'
  
  return, clean
end
