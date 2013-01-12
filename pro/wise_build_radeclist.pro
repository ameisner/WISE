;+
; NAME:
;   wise_build_radeclist
;
; PURPOSE:
;   Build ra, dec index list
;
; CALLING SEQUENCE:
;   wise_build_radeclist
;
; INPUTS:
;   
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2011-Nov-26 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_build_radeclist

  t0 = systime(1)

  dir = '$WISE_DATA'
  flist = file_search(dir+'/L1b/p1bm_frm/*/*/*/*int-1b.fits', count=nfile)
  ra  = fltarr(nfile)
  dec = fltarr(nfile)
  print, nfile, ' files'
  for i=0L, nfile-1 do begin 
     hdr = headfits(flist[i])
     ra[i]  = sxpar(hdr, 'CRVAL1')
     dec[i] = sxpar(hdr, 'CRVAL2')

     if (i mod 5000) eq 0 then begin
        print, i
        euler, ra, dec, l, b, 1
        plot, l, b, ps=3, yr=[-90, 90], xr=[360, 0], /xst, /yst
     endif
  endfor

  str = replicate({fname:'', ra:0., dec:0.}, nfile)
  str.fname = flist
  str.ra = ra
  str.dec = dec
  mwrfits, str, concat_dir(dir, 'index-L1b.fits'), /create

  print, systime(1)-t0, ' seconds'

  return
end
