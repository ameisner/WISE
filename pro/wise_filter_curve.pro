;+
; NAME:
;   wise_filter_curve
;
; PURPOSE:
;  read in WISE filter curve structure (also a place to save routine for
;  building filter curve structures, wise_build_filter)
;   
; CALLING SEQUENCE:
;   trans = wise_filter_curve(band=)
;
; KEYWORDS:
;   band  - WISE band of interest, integer 1-4
;
; OUTPUTS:
;   trans - structure containing filter curve information for desired WISE band
;
; COMMENTS:
;   .txt transmission curves obtained from WISE all-sky explanatory 
;   supplement, section IV.4.h
;    
;   http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/figures/RSR-W?.txt
;
; REVISION HISTORY:
;   2012-Mar-20 - Written by Aaron Meisner
;----------------------------------------------------------------------
function wise_filter_curve, band=band

  if ~keyword_set(band) then band=3

  filtfile = concat_dir(getenv('WISE_DATA'), 'wise-filters.fits')
  trans = mrdfits(filtfile, band)
  return, trans
end

pro wise_build_filter

  wisedata = getenv('WISE_DATA')
  
  for band = 1, 4 do begin
    sband  = string(band, format='(I1)')
    filterfile = 'RSR-W' + sband + '.txt' 
 
    readcol, concat_dir(wisedata, filterfile), lambda_um, response, unc, $
      F = 'F, F, I', COMMENT='#'
  
    npt    = n_elements(lambda_um)
    filtstr = replicate({lambda_um:0., response: 0., unc:0}, npt)
  
    filtstr.lambda_um = lambda_um
    filtstr.response = response
    filtstr.unc = unc ; uncertainty provided as integer parts per thousand

    mwrfits, filtstr, concat_dir(wisedata, 'wise-filters.fits')
   endfor
  
end
