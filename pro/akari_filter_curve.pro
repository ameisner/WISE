;+
; NAME:
;   akari_filter_curve
;
; PURPOSE:
;  read in AKARI filter curve structure (also a place to save routine for
;  building filter curve structures, akari_build_filter)
;   
; CALLING SEQUENCE:
;   trans = akari_filter_curve(band)
;
; INPUTS:
;   band  - AKARI band of interest, one of ['WideS', 'WideL', 'N60', 'N160']
;
; OUTPUTS:
;   trans - structure containing filter curve info for desired AKARI band
;
; COMMENTS:
;   ascii .dat transmission curves obtained from:
;    
;   http://svo.cab.inta-csic.es/theory/fps3/index.php?mode=browse
;
;   note that filter curves are "relative" and each band has its own
;   arbitrary normalization
; REVISION HISTORY:
;   2012-Mar-20 - Written by Aaron Meisner
;----------------------------------------------------------------------
function akari_filter_curve, band

  bands  = ['WideS', 'WideL', 'N60', 'N160']
  w      = where(bands EQ band)
  filterfile = concat_dir(getenv('WISE_DATA'), 'akari-filters.fits')
  trans  = (w[0] GE 0) ? $ 
    mrdfits(filterfile, w[0] + 1, /silent) : 0

  return, trans
end

pro akari_build_filter

  bands = ['WideS', 'WideL', 'N60', 'N160']

  for i = 0, n_elements(bands) - 1 do begin
    filterfile = 'AKARI-FIS.' + bands[i] + '.dat' 
    wisedata = getenv('WISE_DATA')
  
    readcol, concat_dir(wisedata, filterfile), lambda_aa, response, $
     F = 'F, F'
  
    npt    = n_elements(lambda_aa)
    lambda_um = lambda_aa/(1e4)
    filtstr = replicate({lambda_um:0., response: 0.}, npt)
  
    filtstr.lambda_um = lambda_um
    filtstr.response = response
  
    mwrfits, filtstr, concat_dir(wisedata, 'akari-filters.fits')
  endfor

end
