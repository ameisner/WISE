;+
; NAME:
;   akari_sourcemask
;
; PURPOSE:
;   create a mask marking the locations of Akari sources
;
; CALLING SEQUENCE:
;   akari_sourcemask, h, mask, long=, cat=, dilate=
;
; INPUTS:
;   h      - header of fits file for which a source mask is desired
;
; KEYWORDS:
;   long   - use Akari WideL catalog (if not set default to using
;            WideS catalog)
;   cat    - variable into which to store a structure containing all
;            fields present in Akari catalog, but only for the relevant
;            source subset
;   dilate - get a bitmask that masks a ~circular region about each source
;
; OUTPUTS:
;   mask   - mask marking the location of each source with Akari
;            catalog flux in Jy. If /DILATE set, a bit mask flagging
;            pixels nearby each source
;
; EXAMPLES:
;   h =  headfits('/raid14/data/akari/FITS/Release1.0/WideL/093_WideL.fits')
;   akari_sourcemask, h, mask, /LONG, cat=cat, /DILATE
;   help, cat, /ST
;
; COMMENTS:
;  some more work could be done on properly tuning amount of mask
;  dilation as a function of source flux. current dilation is very
;  generous.
;
; REVISION HISTORY:
;   2012-Jun-28 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro akari_sourcemask, h, mask, long=long, cat=cat, dilate=dilate

  nx = sxpar(h, 'NAXIS1')
  ny = sxpar(h, 'NAXIS2')

  akari_starlist, h, xlist, ylist, long=long, cat=cat

  if n_elements(xlist) EQ 0 then return

  ix = round(xlist)
  iy = round(ylist)

  if ~keyword_set(dilate) then begin
    mask = fltarr(nx, ny)
    wim = where((ix GE 0) AND (ix LT nx) AND (iy GE 0) AND (iy LT ny), nwim)
    if nwim EQ 0 then return
; ----- this has a problem if two distinct sources land on the same
;       pixel but in general it might be useful to know the brightness
;       of the source at each location
    mask[ix[wim], iy[wim]] = cat[wim].flux90 ; set to flux in Jy
  endif else begin
    kpix = 21
    mask = bytarr(nx, ny)
    kern = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)
    for i=0L, n_elements(cat)-1 do begin
      mask[((ix[i]-kpix/2) > 0):((ix[i]+kpix/2) < (nx-1)), $ 
               ((iy[i]-kpix/2) > 0):((iy[i]+kpix/2) < (ny-1))] OR= $ 
      kern[((kpix/2-ix[i]) > 0):((nx-ix[i]+kpix/2-1) < (kpix-1)), $ 
        ((kpix/2-iy[i]) > 0):((ny-iy[i]+kpix/2-1) < (kpix-1))]
    endfor
  endelse

end

; test if anything can be broken
pro test_sourcemask, long=long, dilate=dilate

  tilepath = '/raid14/data/akari/FITS/Release1.0/WideS'
  ntile = 430
  tnum = string(lindgen(ntile) + 1, format='(I03)')
  fname = tnum + '_WideS.fits'
  tname = concat_dir(tilepath, fname)
  if keyword_set(long) then tname = repstr(tname, 'WideS', 'WideL')

  for i=0, ntile-1 do begin
    print, tname[i]
    h = headfits(tname[i])
    akari_sourcemask, h, mask, long=long, dilate=dilate
  endfor
end
