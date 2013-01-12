;+
; NAME:
;   wise_read_imcube
;
; PURPOSE:
;   Read many wise images into a data cube for convenience
;
; CALLING SEQUENCE:
;   imcube = wise_read_imcube(lb, nimage=nimage)
;
; INPUTS:
;   lb      - coordinates
;   
; KEYWORDS:
;   nimage  - number of images to read
;   exten   - FITS extension (default to primary HDU)
;   cleanpath - input directory
;
; OUTPUTS:
;   imcube  - image cube
;
; EXAMPLES:
;   
; COMMENTS:
;   Should not hardwire index file!!!!
;
; REVISION HISTORY:
;   2011-Dec-08 - Written by Douglas Finkbeiner, CfA
;   2012-Jan-25 - Added cleanpath keyword to wise_name call - Aaron Meisner
;
;----------------------------------------------------------------------
function wise_read_imcube, lb, nimage=nimage, cleanpath=cleanpath

  t0 = systime(1)

; -------- set defaults
  if ~ keyword_set(lb) then $
    message, 'imcube = wise_read_imcube(lb, nimage=, cleanpath=)'
  if ~ keyword_set(nimage) then nimage = 1900
  if ~ keyword_set(cleanpath) then cleanpath = '$WISE_DATA/clean'

; -------- read metadata table
  indstr = wise_index_metadata(lb, nimage=nimage)
  flist  = indstr.fname

; -------- read header of first file to determine image size
  hname = wise_name('clean', flist[0], cleanpath=cleanpath)

  if file_test(hname) EQ 0 then message, ' Cannot find file '+hname
  hdr = headfits(hname)
  nx = sxpar(hdr, 'NAXIS1')
  ny = sxpar(hdr, 'NAXIS2')

; -------- allocate image cube array
  imcube = fltarr(nx, ny, nimage)

; -------- loop over images and read
  for i=0L, nimage-1 do begin 

     fname = wise_name('clean', flist[i], cleanpath=cleanpath)
     imcube[*, *, i] = $
       readfits(fname, h1, /silent)
     if (i mod 50) eq 0 then print, i
  endfor

  splog, 'Read time:', systime(1)-t0, ' sec'

  return, imcube
end
