;+
; NAME:
;   wisetile_build_radeclist
;
; PURPOSE:
;   Build ra, dec index list
;
; CALLING SEQUENCE:
;   wisetile_build_radeclist
;
; INPUTS:
;   <files>
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
;   2012-Feb-05 - Written (as wise_build_radeclist by Douglas Finkbeiner, CfA
;   2012-Mar-13 - modified to do ISSA tiles of WISE data
;
;----------------------------------------------------------------------
pro wisetile_build_radeclist

  t0 = systime(1)

;  dir = '$WISETILE_DATA'
  dir = '/n/home08/dfink/runwise/tile'
  flist = file_search(dir+'/wise_???.fits', count=nfile)
  if nfile EQ 0 then message, 'No files found'

  ra  = fltarr(nfile)
  dec = fltarr(nfile)
  print, nfile, ' files'
  for i=0L, nfile-1 do begin 
     hdr = headfits(flist[i])
     ra[i]  = sxpar(hdr, 'CRVAL1')
     dec[i] = sxpar(hdr, 'CRVAL2')

     if (i mod 10) eq 0 then print, '.', format='($,A)'

  endfor
  print

  str = replicate({fname:'', ra:0., dec:0.}, nfile)
  str.fname = fname
  str.ra  = ra
  str.dec = dec
  mwrfits, str, concat_dir(dir, 'wisetile-index.fits'), /create

  print, systime(1)-t0, ' seconds'

  return
end
