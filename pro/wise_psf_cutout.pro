;+
; NAME:
;  wise_psf_cutout
;   
; PURPOSE:
;  generate a WISE PSF cutout for specified (x, y) L1b coordinates   
;
; CALLING SEQUENCE:
;  psf = wise_psf_cutout(x, y, BRIGHT=, psf_coeff=, allsky=, w4=,
;                        msk=, flux=)   
;
; INPUTS:
;  x      - x coordinate for PSF centroid
;  y      - y coordinate for PSF centroid
;
; KEYWORDS:  
;  BRIGHT    - set when requesting full PSF, including ghost, for a
;              bright star
;  psf_coeff - pass in pre-read PSF image, rather than re-reading the
;              file containing the PSF for every PSF generated
;  allsky    - set for allsky release (default preliminary release)
;  w4        - set for W4 (default W3)
;  msk       - L1b mask image, don't worry about this keyword as it's rarely
;              useful
;  flux      - source flux in DN, if particular PSF scaling is desired
;  
;
; OUTPUTS:
;   cutout - PSF cutout appropriate to supplied (x, y) L1b
;            coordinates, with L1b pixel scale
;
; EXAMPLES:
;   faint_cutout  = wise_psf_cutout(x, y, /allsky)
;   bright_cutout = wise_psf_cutout(x, y, /BRIGHT, /allsky)
;
; DEPENDENCIES:
;   psf_par_struc.pro
;   read_psf_coeff.pro
;   taper_cutout.pro
;   taper_weight.pro
;   wise_translate_ghost.pro
;
; COMMENTS:
;   This routine does not perform fractional pixel sinc shifting.
;   This routine tapers the PSF model edges gradually and
;   symmetrically to zero by default, and there is currently no keyword
;   to disable this feature.
;
; REVISION HISTORY:
;   2012-Feb-16 - Written by Aaron Meisner
;----------------------------------------------------------------------
function wise_psf_cutout, x, y, BRIGHT=BRIGHT, allsky=allsky, w4=w4, $ 
                          psf_coeff=psf_coeff, msk=msk, flux=flux

  par = psf_par_struc(w4=w4, allsky=allsky, /everything)
; ----- minimum allowed separation b/w PSF star centroid, image edge
  NPAD = keyword_set(w4) ? -0.5 : 114.5
  coord_min_faint = par.psfpix/2-par.pfaint/2 
  coord_max_faint = par.psfpix/2+par.pfaint/2

  if ~keyword_set(psf_coeff) then begin 
    psf_coeff = read_psf_coeff(allsky=allsky, w4=w4)
  endif

  if ~keyword_set(BRIGHT) then begin
; ----- faint star case
    cutout = psf_coeff[coord_min_faint:coord_max_faint, $ 
                       coord_min_faint:coord_max_faint, 0]
    cutout = taper_cutout(cutout, feat='wings', allsky=allsky, w4=w4)
  endif else begin
;----- bright star case
;----- don't extrapolate PSF beyond region over which it was fit
    xpsf   = (x > NPAD) < (par.impix-1-NPAD)
    ypsf   = (y > NPAD) < (par.impix-1-NPAD)

    dx     = xpsf - par.crpix
    dy     = ypsf - par.crpix
    cutout = psf_coeff[*, *, 0] + $ 
             psf_coeff[*, *, 1]*dx + $
             psf_coeff[*, *, 2]*dy + $
             psf_coeff[*, *, 3]*dx*dy + $
             psf_coeff[*, *, 4]*dx^2 + $ 
             psf_coeff[*, *, 5]*dy^2 + $ 
             psf_coeff[*, *, 6]*(dx^2)*dy + $
             psf_coeff[*, *, 7]*(dy^2)*dx + $
             psf_coeff[*, *, 8]*(dx^3) + $ 
             psf_coeff[*, *, 9]*(dy^3)
; ----- substitute in new ghost translation model
    if ~keyword_set(w4) then begin
      ghost_cutout = wise_translate_ghost(x, y, intshift=ishift, $ 
                                          allsky=allsky, w4=w4)
  cutout[(par.psfpix/2-par.xgpix/2+ishift):(par.psfpix/2+par.xgpix/2+ishift), $
         0:(par.ygpix-1)] = ghost_cutout
    endif
  endelse

; ----- if flux keyword set, then scale PSF to requested total DN
  if keyword_set(flux) then begin
      PSFHALF = keyword_set(BRIGHT) ? par.psfpix/2 : par.pfaint/2
      psf_ap_flux = $
    (djs_phot(PSFHALF+0., PSFHALF+0., 10, [PSFHALF/2,2+PSFHALF/2], cutout))[0]
      cutout = cutout*flux/psf_ap_flux ; rescale PSF
  endif

; ---- if msk passed, then interpolate over PSF to match interpolation
;      of image from which the PSF will be subtracted
;      this interpolation makes sense in limit that bad pixels are sparse, but
;      in cases where nearly all pixels are bad, this is probably not the
;      right thing to do
  if keyword_set(msk) then begin
      ix = round(x)
      iy = round(y)
      psfpix = keyword_set(BRIGHT) ? par.psfpix : par.pfaint
      mcut = $ 
        byte(wise_l1b_cutout(msk, round(x), round(y), psfpix, psfpix, w4=w4))
      mcut = (mcut AND (cutout NE 0))
      intx = djs_maskinterp(cutout, mcut, iaxis=0, /const)
      inty = djs_maskinterp(cutout, mcut, iaxis=1, /const)
      cutout = (intx+inty)/2
  endif
  return, cutout

end
