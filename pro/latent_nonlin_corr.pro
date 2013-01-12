;+
; NAME:
;   latent_nonlin_corr
;
; PURPOSE:
;   return latent nonlinearity correction factor
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
;   see wise_correct_latent.pro
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Oct-21 - Aaron Meisner
;----------------------------------------------------------------------
function latent_nonlin_corr, mag, allsky=allsky, w4=w4

  par = psf_par_struc(w4=w4, allsky=allsky, /everything)
  corr_ulim_allsky = 1.804

;----- correction = nlin_coeff[0]+w?mag*nlin_coeff[1] + ...
  nlin_coeff = keyword_set(w4) ? [1., -0.0118, -0.0267] : $ 
(keyword_set(allsky) ? [1., 0.2207946, -0.02749962, -0.01057454] : [1, 0.228])

  if keyword_set(w4) then begin
      nlin_corr = (nlin_coeff[0] + $ 
                  ((mag > par.latmin) < par.latmax)*nlin_coeff[1] + $ 
                  (((mag > par.latmin) < par.latmax)^2)*nlin_coeff[2])
  endif else begin
      if keyword_set(allsky) then begin
          fitmin = -3.
          fitmax = 3.1
          nlin_corr = nlin_coeff[0] + $ 
                      ((mag > fitmin) < fitmax)*nlin_coeff[1] + $ 
                      (((mag > fitmin) < fitmax)^2)*nlin_coeff[2] + $ 
                      (((mag > fitmin) < fitmax)^3)*nlin_coeff[3]
      endif else begin
          nlin_corr = nlin_coeff[0] + $ 
                      ((mag > par.latmin) < par.latmax)*nlin_coeff[1]
      endelse
  endelse

  return, nlin_corr

end
