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
function latent_nonlin_corr, mag, allsky=allsky, w4=w4, band=band

  par = psf_par_struc(w4=w4, allsky=allsky, /everything, band=band)

  case band of
      1 : nlin_coeff = [2.7763787, -2.4806535, 0.81217395]
      2 : nlin_coeff = [1.2812997, -0.97785225, 1.0474848]
      3 : nlin_coeff = keyword_set(allsky) ? $ 
                       [1., 0.2207946, -0.02749962, -0.01057454] : [1, 0.228]
      4 : nlin_coeff = [1., -0.0118, -0.0267]
  endcase

  mags_u = [3.30, 2.50, keyword_set(allsky) ? 3.1 : par.latmax, par.latmax]
  mags_l = [1.95, 0.65, keyword_set(allsky) ? -3. : par.latmin, par.latmin]

  mag_u = mags_u[band-1]
  mag_l = mags_l[band-1]

  if (band NE 3) then begin
      nlin_corr = (nlin_coeff[0] + $ 
                  ((mag > mag_l) < mag_u)*nlin_coeff[1] + $ 
                  (((mag > mag_l) < mag_u)^2)*nlin_coeff[2])
  endif else begin
; ----- W3 case
      if keyword_set(allsky) then begin
          nlin_corr = nlin_coeff[0] + $ 
                      ((mag > mag_l) < mag_u)*nlin_coeff[1] + $ 
                      (((mag > mag_l) < mag_u)^2)*nlin_coeff[2] + $ 
                      (((mag > mag_l) < mag_u)^3)*nlin_coeff[3]
      endif else begin
          nlin_corr = nlin_coeff[0] + $ 
                      ((mag > mag_l) < mag_u)*nlin_coeff[1]
      endelse
  endelse

  return, nlin_corr

end
