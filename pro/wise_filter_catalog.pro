;+
; NAME:
;   wise_filter_catalog
;
; PURPOSE:
;   apply cuts to sources in WISE source catalog
;
; CALLING SEQUENCE:
;   wgood = wise_filter_catalog(cat, ngood=, w4=)
;
; INPUTS:
;   cat    - catalog structure containing superset of desired sources
;
; KEYWORDS:
;   w4     - set to make cuts specific to W4 (rather than default W3),
;            not yet implemented
;
; OUTPUTS:
;   wgood  - indices of cat structure containing sources that pass cuts
;
; OPTIONAL OUTPUTS:
;   ngood  - number of sources that pass cuts
;
; EXAMPLES:
;   see wise_starlist
;
; COMMENTS:
;   right now this routine doesn't do much but as source cuts become
;   more elaborate it will be good to have a dedicated routine to ensure
;   same cuts are always being made
;
; REVISION HISTORY:
;   2012-Sep-13 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function wise_filter_catalog, cat, ngood=ngood, w4=w4

; ----- check curve of growth, mask: 1=good, 0=bad
  mask = wise_filter_cog(cat, faint=faint, w4=w4)

; ----- no spurious sources of any origin
  cc_map = keyword_set(w4) ? cat.w4cc_map : cat.w3cc_map
  wgood = where((mask) AND $ 
               ~(((cc_map AND 256) EQ 256) AND (~faint)), ngood)
  return, wgood

end
