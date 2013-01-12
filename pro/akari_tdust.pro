;+
; NAME:
;   akari_tdust
;   
; REVISION HISTORY:
;   2012-Jun-25 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function integrate_response, const=const, alpha=alpha, T=T, long=long

  band = keyword_set(long) ? 'WideL': 'WideS'
  if ~keyword_set(alpha) then alpha = 2.
  if ~keyword_set(T) then T = 18. ; K
  filt = akari_filter_curve(band)
  nsamp = n_elements(filt)
  filt = filt[reverse(lindgen(nsamp))] ; want to order by increasing frequency
  freq_THz = (1e-12)*(3e8)/((1e-6)*filt.lambda_um)
  
  indsamp = lindgen(nsamp-1)+0.5
  indlower = lindgen(nsamp-1) 
  indupper = lindgen(nsamp-1)+1
  
  dfreq = freq_THz[indupper]-freq_THz[indlower]
  trans = interpolate(filt.response, indsamp)

; ----- sample at center of frequency bin
  freq_GHz = 1000*interpolate(freq_THz, indsamp)

  if keyword_set(const) then return, total(dfreq*trans)

  b_nu = djs_planck(T, freq_GHz, units='GHz', /MJy)

; ----- the 1e12 makes it so that there are no dimensions of Hz
;       remaining, not sure what the order of magnitude of the result (~1e17)
;       really means though
  return, total((b_nu*1e12)*((freq_GHz/1000.)^alpha)*dfreq)
  ;stop

end

; generate a lookup table of f(T) where WideS = f*WideL
; eventually decide on the range/resolution of T values worth computing
pro akari_gen_rat, Tvals, f, alpha=alpha

  if ~keyword_set(alpha) then alpha=2.

; ----- 16 K -> 25 K in steps of 0.01 K
  Tvals = 0.01*lindgen(901) + 16.

  nT = n_elements(Tvals)
; ----- WideS = f*WideL
  f = fltarr(nT)

  norm_short = integrate_response(alpha=alpha, /const) 
  norm_long = integrate_response(/const, /long, alpha=alpha)
  for i=0, nT-1 do begin
    integral_short = integrate_response(T=Tvals[i], alpha=alpha)
    integral_long = integrate_response(T=Tvals[i], /long, alpha=alpha)
    f[i] = integral_short*norm_long/(integral_long*norm_short)
  endfor

end

; use lookup table built by akari_gen_rat to interpolate
pro akari_tdust, alpha=alpha, tvals, ratvals

  akari_gen_rat, tvals, ratvals, alpha=alpha
; ----- take ISSA tile 116 as an example for now
  ratmap, ratshort, ratlong, corrshort, corrlong
; ----- only try to compute a temperature where both WideS and WideL
;       plausibly show a correlation with WISE 12um
  wgood = where((corrshort GT 0.3) AND (corrlong GT 0.3))
; ---- convention of ratio = WideS/WideL to be consistent with
;      akari_gen_rat
  imrat = ratshort[wgood]/ratlong[wgood]
; ----- don't extrapolate or do anything with NaNs
  wrat = $ 
   where(finite(imrat) AND (imrat LT max(ratvals)) AND (imrat GT min(ratvals)))
; ----- now identify the T value in T vs. ratio curve corresponding to each
;       good ratio, which is annoying because ratio values
;       aren't sampled uniformly
  coeff = poly_fit(ratvals, tvals, 2)
  t_out = coeff[0] + imrat[wrat]*coeff[1] + (imrat[wrat]^2)*coeff[2]

  imtemp = fltarr(n_elements(ratshort)) ; 1500^2
  imtemp[wgood[wrat]] = t_out
  sz = size(ratshort, /DIM)
  imtemp = reform(imtemp, sz[0], sz[1])
end
