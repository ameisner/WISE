;+
; NAME:
;   atlas_build_radeclist
;
; PURPOSE:
;   make an index structure for Doug's old preliminary release
;   Atlas downloads, with full file name and (ra, dec) coords
;
; CALLING SEQUENCE:
;   indstr = atlas_build_radeclist(dir=)
;
; KEYWORDS:
;   dir    - base directory into which L3a Atlas files were downloaded
;
; OUTPUTS:
;   indstr - index structure listing full file
;
; EXAMPLES:
;   indstr = atlas_build_radeclist()
;   mwrfits, indstr, '$WISE_DATA/index-prelim-L3a.fits'
;
; COMMENTS:
;   intended for odyssey cluster, not tested for other computers e.g. nebel
;
;   is there a useful way to sort the index structure?
;
; REVISION HISTORY:
;   2012-Sep-13 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function atlas_build_radeclist, dir=dir

  if ~keyword_set(dir) then dir='$WISE_DATA/L3a'
  flist = file_search(dir+'/*/*w3-int-3.fits', count=nfile)

  ra = fltarr(nfile)
  dec = fltarr(nfile)

  for i=0L, nfile-1 do begin
      if (i MOD 20) EQ 0 then print, i
      h = headfits(flist[i])
      ra[i] = sxpar(h, 'CRVAL1')
      dec[i] = sxpar(h, 'CRVAL2')
  endfor

  indstr = replicate({fname:'', ra:0., dec:0.}, nfile)

  indstr.fname = flist
  indstr.ra = ra
  indstr.dec = dec

  return, indstr

end
