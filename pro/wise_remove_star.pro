;+
; NAME:
;   wise_remove_star
;
; PURPOSE:
;   Remove stars from a WISE image
;
; CALLING SEQUENCE:
;   wise_remove_star, image, x, y, flux, bthresh=bthresh, mask=mask
;
; INPUTS:
;   image      - image to remove stars from 
;   x, y       - coords of stars to remove
;   flux       - flux of stars (so we can start with the brightest)
;
; OPTIONAL INPUTS:
;   bthresh    - threshold for calling a star "bright"
;   mask       - mask [1=pixel modified] (created if not passed)
;
; OPTIONAL OUTPUTS:
;   mask       - mask [1=pixel modified]
;
; OUTPUTS:
;   image      - continuum-corrected image - overwritten
;
; COMMENTS:
;   replace point sources with median of annulus around them, one by one.
;   should always go brightest to faintest!!!
; 
;   Bright stars (flux > bthresh) get a bigger disc. 
;
; PROCEDURES CALLED:
;   linfit
;   reverse
;
; REVISION HISTORY:
;   2001-Nov-13  Written as halpha_remove_star by Douglas Finkbeiner, Princeton
;   2001-Dec-01  Only use unmasked pixels in gradient fit - DPF
;   2002-Sep-11  Can use non-square arrays - DPF
;   2011-Apr-16  Based on halpha_remove_star.pro by DPF
;-
;------------------------------------------------------------------------------
pro wise_remove_star, image, xin, yin, fluxin, mask=mask

  sx = (size(image, /dim))[0]
  sy = (size(image, /dim))[1]
  bthresh = 500.

; re-order the stars, brightest first
  sind = sort(-fluxin)
  x = xin[sind]
  y = yin[sind]
  flux = fluxin[sind]

; define postage stamp size
  nstar = n_elements(x)
  sz = 19
  xbox = (lindgen(sz, sz) MOD sz) - (sz/2)
  ybox = (lindgen(sz, sz) / sz) - (sz/2)
  
; create mask if it was not passed
  if n_elements(mask) EQ 0 then mask = bytarr(sx, sy)

  for i=0, nstar-1 do begin 
     ix = round(x[i])
     iy = round(y[i])
     if ((ix < iy) GT sz/2) AND (ix LT sx-1-sz/2) AND (iy LT sy-1-sz/2) then begin 
        stamp = image[ix-sz/2:ix+sz/2, iy-sz/2:iy+sz/2]
        smask = mask[ix-sz/2:ix+sz/2, iy-sz/2:iy+sz/2]
        d = sqrt((xbox-(x[i]-ix))^2+(ybox-(y[i]-iy) )^2)
        rad = flux[i] GT bthresh ? 6.0: 3.5
        if flux[i] GT 5000 then rad = 8.


        ; get index lists
        in = where(d lt rad)
        an = where((d ge rad) and (d lt (rad+1.5)))
        ; determine sky with annulus median (good enough for this)
        good = where(finite(stamp[an]) AND (smask[an] EQ 0), ngood)

        if ngood GE 10 then begin ; demand at least 10 annulus pixels

           sky = median(stamp[an[good]], /even)
           peak = max(abs(stamp[in[good]]-sky))    ; highest value
           if peak GT (.5*sky) then begin 
              ; determine local gradient by fitting pixels in annulus
              ndeg = 1
              nsig = 5
              poly_iter, xbox[an[good]], stamp[an[good]], ndeg, nsig, coeff=xcoeff
              poly_iter, ybox[an[good]], stamp[an[good]], ndeg, nsig, coeff=ycoeff

;              print, xcoeff, ycoeff, peak

              xslope = xcoeff[1]
              yslope = ycoeff[1]

              ; skip sources that are on too steep a slope. 
              if peak GT (sqrt(xslope^2 + yslope^2)*6) then begin 
                 stamp[in] = 0.5*(xcoeff[0]+ycoeff[0])+ $
                                  xbox[in]*xslope+ybox[in]*yslope
                 smask[in] = 1B
                 
                 image[ix-sz/2:ix+sz/2, iy-sz/2:iy+sz/2] = stamp
                 mask[ix-sz/2:ix+sz/2, iy-sz/2:iy+sz/2] = smask
              endif 
           endif 
        endif 
     endif 
  endfor

  return
end
