;+
; NAME:
;   wise_name
;
; PURPOSE:
;   Return the name for SDSS data file including path information.
;
; CALLING SEQUENCE:
;   fullname = wise_name( ftype, /no_path ] )
;
; INPUTS:
;   ftype      - File type; supported types are:
;                  clean
;
;   no_path    - If set, then do not set the default path for this file name.
;
; OUTPUTS:
;   fullname   - Full file name (which may not actually exist on disk)
;
; COMMENTS:
;   modeled on sdss_name.pro by D. J. Schlegel.
;
; REVISION HISTORY:
;   16-Dec-2011  Written by Douglas Finkbeiner, CfA
;   25-Jan-2011  Added cleanpath keyword - Aaron Meisner
;-
;------------------------------------------------------------------------------
FUNCTION wise_name, ftype, fname, no_path=no_path, cleanpath=cleanpath

; -------- 
  if ftype EQ 'clean' then begin 
     if ~ keyword_set(cleanpath) then cleanpath = '$WISE_DATA/clean'
     thisfile = fileandpath(fname)
     thisdir  = strmid(thisfile, 0, 4)
     thispath = concat_dir(cleanpath, thisdir)
     fullname = concat_dir(thispath, thisfile)
     
  endif else begin
     message, 'Unknown FTYPE=' + ftype
  endelse
  
  return, fullname
end
