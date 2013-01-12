;+
; NAME:
;   read_psf_coeff
;
; PURPOSE:
;   handle reading and caching of PSF model
;
; CALLING SEQUENCE:
;   psf_coeff = read_psf_coeff(allsky=)
;
; KEYWORDS:
;   allsky    - set if All-sky release PSF model is desired
;
; OUTPUTS:
;   psf_coeff - PSF model, suitable to be input to wise_psf_cutout
;
; EXAMPLES:
;   see wise_l1b_substar.pro
;
; COMMENTS:
;   should generalize to w4 eventually
;
; REVISION HISTORY:
;   2012-Sep-22 - Written by Aaron Meisner
;----------------------------------------------------------------------
function read_psf_coeff, allsky=allsky, w4=w4

  par = psf_par_struc(allsky=allsky, w4=w4, /everything)
  coeff_file = par.fpsf

  COMMON PSFIMAGE, psf_coeff, coeff_file_sav
  if (n_elements(psf_coeff) EQ 0) || (coeff_file NE coeff_file_sav) then begin 
      psf_coeff = readfits(coeff_file)
      ncoeff = (size(psf_coeff,/DIM))[2]
      for i=0, ncoeff-1 do $
          psf_coeff[*,*,i] = $ 
              taper_cutout(psf_coeff[*,*,i], feat='wings', /bright, $ 
                           allsky=allsky, w4=w4)
      coeff_file_sav = coeff_file
  endif

  return, psf_coeff

end
