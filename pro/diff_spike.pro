;+
; NAME:
;   diff_spike
;
; PURPOSE:
;   for given L1b exposure, determine where diffraction spikes are and
;   return a corresponding mask
;   
; CALLING SEQUENCE:
;   spkmsk = diff_spike(h, allsky=allsky, w4=w4)
;
; INPUTS:
;   h  - header for L1b image of interest
;
; KEYWORDS:
;   allsky - set for all-sky release
;   w4     - set for W4 rather than W3 (not yet tested)
;
; OUTPUTS:
;   spkmsk - bitmask marking diffraction spike locations (1=spike, 0
;
; EXAMPLES:s
;   see wise_l1b_clean.pro
;
; COMMENTS:
;   right now only marking optical spikes, have not gotten around to 
;   calibrating extent of electronic bleeding spikes as a function of source
;   brightness
;
; REVISION HISTORY:
;   2012-Nov-19 - Aaron Meisner
;----------------------------------------------------------------------
function rot_corr, beta, w4=w4

; make simple-minded correction that accounts for the fact that at
; high ecliptic latitudes, observations at different epochs have 
; different orientations, thus spreading out and reducing apparent extent of
; diffractions spikes upon stacking

  par = psf_par_struc(w4=w4)
  spar = diff_spike_par(w4=w4)
; ----- fiducial angular extent of single-epoch diffraction spikes, based
;       on examination of PSF model
  dtheta0 = spar.dtheta0 ; radians

; ----- sidelength of L1b expsosure in degrees
  sz = par.pscl*par.impix/3600. ; deg

; ----- assume beta specified in degrees
  dtheta = 2.*!pi*(sz/(360.*cos(beta/!radeg))) < !pi/2 ; radians
  
; ----- factor by which diff spike amplitude reduced by changing orientation
  sbfac = dtheta/dtheta0
  corr = 2.5*alog10(sbfac > 1)

  return, corr
end

function effective_mag, mag, m, w4=w4, i100rms=i100rms, bound=bound, $ 
                        radec=radec
  par = diff_spike_par(w4=w4)
  meff = mag - (keyword_set(i100rms) ? ~((i100rms GT 5.) AND $ 
(m GT par.covmed)) : 1.)*2.5*alog10(sqrt(float((m < par.mmax) > 1)/par.covmed))
  if keyword_set(i100rms) then meff += par.i100fac*alog10((i100rms-2.) > 1)
  if keyword_set(radec) then begin
      euler, radec[*,0], radec[*,1], lambda, beta, 3
      rcorr = rot_corr(beta, w4=w4)
      meff += rcorr
  endif
  if keyword_set(bound) then meff = (meff < par.meff_max) > par.meff_min

  return, meff

end

function spike_half_length, mag, m, w4=w4, i100rms=i100rms, radec=radec

  par = diff_spike_par(w4=w4)

  meff = effective_mag(mag, m, w4=w4, i100rms=i100rms, radec=radec, /bound)
; ----- units of arcminutes
  coeff = par.coeff
  len = (coeff[0]+coeff[1]*meff+coeff[2]*(meff^2)+coeff[3]*(meff^3))
  
  return, len
end

function diff_spike, h, allsky=allsky, w4=w4

  init_bright, allsky=allsky, w4=w4
  COMMON BRIGHT, ra_brt, dec_brt, mag_brt, m_brt, len_deg, i100rms

  par = psf_par_struc(w4=w4, allsky=allsky, /everything)
  spar = diff_spike_par(w4=w4)
  spmsk = bytarr(par.impix, par.impix)

  racen = sxpar(h,'CRVAL1')
  deccen = sxpar(h,'CRVAL2')
  maxsep = (par.impix/sqrt(2))*(par.pscl/3600.)+len_deg
  binary_search, dec_brt, deccen+max(maxsep), indmax
  binary_search, dec_brt, deccen-max(maxsep), indmin
  if (indmax EQ -1) then indmax = n_elements(mag_brt)-1
  indmin >= 0
  if (indmax EQ indmin) then return, spmsk
  wdec   = lindgen(indmax - indmin + 1)+indmin

  dangle = djs_diff_angle(racen, deccen, ra_brt[wdec], dec_brt[wdec])
  wspk = where(dangle LT maxsep[wdec], nspk)

  if nspk EQ 0 then return, spmsk

  extast, h, astr
  ad2xy, ra_brt[wdec[wspk]], dec_brt[wdec[wspk]], astr, x, y
  ix = round(x)
  iy = round(y)
  for i=0L, nspk-1 do begin
; ---- what is right/optimal amount of dilation ????
      thismsk = diff_spike_mask(mag_brt[wdec[wspk[i]]], m_brt[wdec[wspk[i]]], $
          w4=w4, sz=sz, dil=spar.dil, i100rms=i100rms[wdec[wspk[i]]], $ 
          radec=[[ra_brt[wdec[wspk[i]]]], [dec_brt[wdec[wspk[i]]]]])
      incl = wise_l1b_cutout(_, ix[i], iy[i], sz, sz, w4=w4, /bool, /silent)
      if ~incl then continue
      spmsk[((ix[i]-sz/2) > 0):((ix[i]+sz/2) < (par.impix-1)), $ 
            ((iy[i]-sz/2) > 0):((iy[i]+sz/2) < (par.impix-1))] OR= $ 
      thismsk[((sz/2-ix[i]) > 0):((par.impix-ix[i]+sz/2-1) < (sz-1)), $ 
              ((sz/2-iy[i]) > 0):((par.impix-iy[i]+sz/2-1) < (sz-1))]
  endfor

  return, spmsk

end
