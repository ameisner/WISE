;+
; NAME:
;   wise_translate_ghost
;
; PURPOSE:
;   
; CALLING SEQUENCE:
;   
; INPUTS:
;   x - x coordinate of parent source
;   y - y coordinate of parent source
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   cutout - ghost cutout, with appropriate fractional pixel sinc 
;            interpolation applied
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   only intended for W3 at this point
;
; REVISION HISTORY:
;   2012-Nov-29 - Written by Aaron Meisner
;----------------------------------------------------------------------
function wise_translate_ghost, x, y, intshift=intshift, allsky=allsky, w4=w4

  par = psf_par_struc(allsky=allsky, w4=w4, /everything)

  COMMON GHOST, ghost_cutout, scalefac, xshift, yshift
  if n_elements(ghost_cutout) EQ 0 then begin
      ghost_cutout = readfits(par.fghost)
      scalefac = $ 
          readfits('/n/panlfs/ameisner/psf/results/w3_ghost-scalefac.fits')
      xshift = readfits('/n/panlfs/ameisner/psf/results/w3_ghost-shift.fits')
      yshift = $ 
          readfits('/n/panlfs/ameisner/psf/results/w3_ghost-shift.fits', ex=1)
  endif

; ----- bound coordinates within region over which ghost shift has been fit
  xx = (x < (par.impix-0.5)) > (-0.5)
  yy = (y < (par.impix+par.ygoffs-0.5)) > (par.ygoffs-0.5)

  sfac = interpolate(scalefac, (xx-7.5)/10., (yy-(par.ygoffs+7.5))/10.)
  dx = interpolate(xshift,(xx-7.5)/10., (yy-(par.ygoffs+7.5))/10.)
  dy = interpolate(yshift,(xx-7.5)/10., (yy-(par.ygoffs+7.5))/10)
;;  xcoeff = [-0.12814304, 0.010751015, 0.00011044494]
;;  ycoeff = [0.056743526, 0.0011828195, 0.0014948531]
    
;;  dx = xcoeff[0]+xcoeff[1]*(xx-par.crpix)+xcoeff[2]*(yy-par.crpix)
;;  dy = ycoeff[0]+ycoeff[1]*(xx-par.crpix)+ycoeff[2]*(yy-par.crpix)

  xghost = x + dx
;  yghost = y + dy

  ix = round(xghost)
;  iy = round(yghost)

;;  intshift = round((round(x)-x)+dx)
  intshift = round(xghost)-round(x)
; ----- y shift should be implemented eventually, but don't do it yet
  xyshift = [xghost-ix, y-round(y)+dy]

  cutout = sshift2d(ghost_cutout, xyshift)
; ----- taper after sinc shift
  cutout = taper_cutout(cutout, feat='ghost', allsky=allsky, w4=w4)
  cutout *= sfac
  return, cutout

end

pro test_translate, xshift, yshift

  xbox = (findgen(1200,1200) MOD 1200) - 100
  ybox = (findgen(1200,1200) / 1200) + 100

  xshift = fltarr(1200,1200)
  yshift = fltarr(1200,1200)
  for i = 0, 1199 do begin
      print, i
      for j = 0, 1199 do begin
          xyshift = wise_translate_ghost(xbox[i,j], ybox[i,j], /allsky)
          xshift[i,j] = xyshift[0]
          yshift[i,j] = xyshift[1]
      endfor
  endfor

end
