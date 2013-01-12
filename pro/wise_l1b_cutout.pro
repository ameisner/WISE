;+
; NAME:
;   wise_l1b_cutout
;
; PURPOSE:
;   retrieve cutout from L1b image
;
; CALLING SEQUENCE:
;   cutout = wise_l1b_cutout(im, ix, iy, xpix, ypix, wt=, w4=, bool=)
;
; INPUTS:
;   im     - L1b image from which to extract cutout
;   ix     - integer central x coordinate of cutout region
;   iy     - integer central y coordinate of cutout region
;   xpix   - x sidelength of desired cutout (pixels)
;   ypix   - y sidelenght of desired cutout (pixels)
;
; KEYWORDS:
;   w4     - set if considering w4 L1b image
;   bool   - set if you only want a bitflag indicating whether desired
;             cutout has nonzero overlap with (1=some overlap, 0=no overlap)
;
; OUTPUTS:
;   cutout - cutout with specified center and sidelengths, zero where
;            data is not available (outside of image boundaries)
;
; OPTIONAL OUTPUTS:
;   wt     - xpix by ypix coverage bitmask (1=in image, 0=outside of image)
;
; EXAMPLES:
;   see wise_model_psf.pro
;
; COMMENTS:
;   cutout need not be entirely contained within image, and
;   the center of the cutout need not be within the image boundaries
;
; REVISION HISTORY:
;   2012-Oct-15 - Aaron Meisner
;----------------------------------------------------------------------
function wise_l1b_cutout, im, ix, iy, xpix, ypix, wt=wt, w4=w4, bool=bool, $
                          full=full, silent=silent

; ----- eventually generalize to even sidelengths
  if ((xpix MOD 2) EQ 0) OR ((ypix MOD 2) EQ 0) then begin
      if ~keyword_set(silent) then $
          print, 'this routine assumes odd sidelengths'
      return, 0
  endif

  par = psf_par_struc(w4=w4)
  impix = par.impix
  xhalf = xpix/2
  yhalf = ypix/2

  cutout = fltarr(xpix, ypix)
  wt = bytarr(xpix, ypix)

  xmax = ix+xhalf
  xmin = ix-xhalf
  ymax = iy+yhalf
  ymin = iy-yhalf

  if (keyword_set(full) AND ((xmax GT (IMPIX-1)) OR (xmin LT 0) OR $ 
      (ymax GT (IMPIX-1)) OR (ymin LT 0))) then begin
      if ~keyword_set(silent) then $
          print, 'cutout is not fully contained within the image'
      return, 0
  endif

  if ((xmax LT 0) OR (xmin GT (impix-1)) OR $ 
      (ymax LT 0) OR (ymin GT (impix-1))) then begin
      if ~keyword_set(silent) then $ 
          print, 'desired cutout has no overlap with image footprint'
      return, 0
  endif
  if keyword_set(bool) then return, 1

  cutout[((xhalf-ix) > 0):((impix-ix+xhalf-1) < (xpix-1)), $
         ((yhalf-iy) > 0):((impix-iy+yhalf-1) < (ypix-1))] = $ 
      im[((ix-xhalf) > 0):((ix+xhalf) < (impix-1)), $ 
         ((iy-yhalf) > 0):((iy+yhalf) < (impix-1))] 
  wt[((xhalf-ix) > 0):((impix-ix+xhalf-1) < (xpix-1)), $
     ((yhalf-iy) > 0):((impix-iy+yhalf-1) < (ypix-1))] = 1

  return, cutout

end
