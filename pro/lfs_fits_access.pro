;+
; NAME:
;   lfs_fits_access
;
; PURPOSE:
;   try accessing a fits file once, wait a few seconds if file not
;   found, and then try a second and final time to read the file
;
; CALLING SEQUENCE:
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
;   if headeronly keyword not set, this routine wraps readfits
;   if headeronly keyword is set, this routine wraps headfits
;   tries at most TWICE to read target fits file
;
; REVISION HISTORY:
;   2012-May-18 - Written by Aaron Meisner
;----------------------------------------------------------------------
function lfs_fits_access, filename, h, headeronly=headeronly, silent=silent, $ 
                          exten_no=exten_no, pause=pause
  firstattempt = 1

  RETRY:
  if keyword_set(headeronly) then begin
    out = headfits(filename, silent=silent, exten=exten_no)
  endif else begin
    if arg_present(h) then $
      out = readfits(filename, h, silent=silent, exten_no=exten_no) $ 
    else $ 
      out = readfits(filename, silent=silent, exten_no=exten_no)
  endelse

;----- on "Unable to locate file" errors readfits, headfits return -1
  if (n_elements(out) EQ 1) AND firstattempt then begin
    print, 'LFS READ FAILED ', filename, ' ', systime(/UTC)
    wait, keyword_set(pause) ? pause : 10
    firstattempt = 0
    GOTO, RETRY
  endif

  return, out
end
