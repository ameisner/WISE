;+
; NAME:
;   match_astr_params
;
; PURPOSE:
;   gather parameters of header necessary to form astrometry structure
;   with extast
;
; CALLING SEQUENCE:
;   match_astr_params, fname
;
; INPUTS:
;   fname  - file name for which header parameters are desired
;   
; OUTPUTS:
;   hastr  - 78 element array of 32-character strings, performing
;            extast on this array is identical to performing
;            extast on the full header
;
; EXAMPLES:
;   see wise_build_hdr
;
; COMMENTS:
;   it would be good to experiment with further reducing size of
;   strings in hastr output, since right now an index structure with
;   header information appended is dominated by the header astrometry
;   in terms of storage space...
;   see also .header field of $WISE_DATA/index-meta-astr-L1b.fits
;   
; REVISION HISTORY:
;   2011-Mar-17 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function match_astr_params, fname
 
;------- subset of header keywords  necessary for extast

  COMMON KASTR, kastr

  if (n_elements(kastr) EQ 0) then $ 
  kastr = ['NAXIS1', 'NAXIS2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'WCDELT1', $
          'WCDELT2', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CTYPE1', $ 
          'CTYPE2', 'A_0_0', 'A_0_1', 'A_0_2', 'A_0_3', 'A_0_4', 'A_1_0', $ 
          'A_1_1', 'A_1_2', 'A_1_3', 'A_2_0', 'A_2_1', 'A_2_2', 'A_3_0', $ 
          'A_3_1', 'A_4_0', 'A_ORDER', 'AP_0_0', 'AP_0_1', 'AP_0_2', $ 
          'AP_0_3', 'AP_0_4', 'AP_1_0', 'AP_1_1', 'AP_1_2', 'AP_1_3', $ 
          'AP_2_0', 'AP_2_1', 'AP_2_2', 'AP_3_0', 'AP_3_1', 'AP_4_0', $ 
          'AP_ORDER', 'B_0_0', 'B_0_1', 'B_0_2', 'B_0_3', 'B_0_4', 'B_1_0', $ 
          'B_1_1', 'B_1_2', 'B_1_3', 'B_2_0', 'B_2_1', 'B_2_2', 'B_3_0', $ 
          'B_3_1', 'B_4_0', 'B_ORDER', 'BP_0_0', 'BP_0_1', 'BP_0_2', $ 
          'BP_0_3', 'BP_0_4', 'BP_1_0', 'BP_1_1', 'BP_1_2', 'BP_1_3', $ 
          'BP_2_0', 'BP_2_1', 'BP_2_2', 'BP_3_0', 'BP_3_1', 'BP_4_0', $ 
          'BP_ORDER']

  h     = headfits(fname)
  kall  = strtrim(strmid(h,0,8))

  matchlist, kastr, kall, m1, m2
  if (n_elements(m1) NE n_elements(kastr)) then print, 'MISSING KEYWORD'

;------- only leading 32 characters to avoid storing same comments 1M times
  hastr = strmid(h[m2], 0, 32)
  
  return, hastr
end
