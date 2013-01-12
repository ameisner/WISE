;+
; NAME:
;   wise_write_wcs_file
;
; PURPOSE:
;   Gather WCS header info for all WISE fields into one file
;
; CALLING SEQUENCE:
;   wise_write_wcs_file, datadir=
;
; INPUTS:
;   datadir  - wise data root
;
; OUTPUTS:
;   <file> wise-wcs.fits
;
; COMMENTS:
;   Run this once to gather all the header info.
;   
; REVISION HISTORY:
;   2011-May-15 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_write_wcs_file, datadir=datadir

  if file_test('$WISE_DIR')  eq 0 then message, 'Please set $WISE_DIR'
  if ~ keyword_set(datadir) then begin 
     if file_test('$WISE_DATA') eq 0 then message, 'Please set $WISE_DATA'
     wise_root = '$WISE_DATA/L3a'
  endif else wise_root = datadir

  print, 'Getting file list...'
  flist = file_search(concat_dir(wise_root, '*/*-w3-int-3.fits'), count=nfile)
  print, nfile, ' files found'

; -------- get an example astrometry structure, insert nonsense values
  h = headfits(flist[0])
  extast, h, astr
  astr.ctype = ['Undefined', ' ']
  astr.cdelt = [0.d, 0.d]
  astr.crval = [-999d, -999]
  astr2 = struct_addtags(astr, {idstr:'None'})
  arr = replicate(astr2, nfile)

  for i=0, nfile-1 do begin 
     filename = flist[i]
     idstr = strmid(fileandpath(filename), 0, 8)
     h = headfits(filename)
     
     extast, h, astr
     astr2 = struct_addtags(astr, {idstr:idstr})
     arr[i] = astr2
     if (i mod 100) eq 0 then print, i
  endfor

  mwrfits, arr, '$WISE_DIR/etc/wise-wcs.fits', /create

  return
end
