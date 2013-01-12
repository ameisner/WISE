;+
; NAME:
;   init_bright
;
; PURPOSE:
;   handle initialization/caching of WISE bright star catalog used for 
;   diffraction spike masking
;
; CALLING SEQUENCE:
;   init_bright, allsky=, w4=
;
; KEYWORDS:
;   allsky - set for all-sky release
;   w4     - set for W4 (default W3)
;
; EXAMPLES:
;   see diff_spike.pro
;
; COMMENTS:
;   dependent upon Eddie's heal_interp.pro
;
; REVISION HISTORY:
;   2012-Nov-22 - Aaron Meisner
;----------------------------------------------------------------------
pro init_bright, allsky=allsky, w4=w4

; ----- cache list of sources bright enough to require diffraction
;       spike masking
  COMMON BRIGHT, ra_brt, dec_brt, mag_brt, m_brt, len_deg, i100rms
  if n_elements(ra_brt) EQ 0 then begin
      init_source_catalog, allsky=allsky, w4=w4
      COMMON CATALOG, catra, catdec, catmag, mjdmin, mjdmax, catm
      spar = diff_spike_par(w4=w4)
      rms_heal = sfd_i100_rms()
      euler, catra, catdec, lgal, bgal, 1
      i100rms = heal_interp(rms_heal, lgal, bgal)
      meff = effective_mag(catmag, catm, w4=w4, i100rms=i100rms, $ 
                           radec=[[catra],[catdec]])
      wbrt = where(meff LT spar.meff_max)
      ra_brt = catra[wbrt]
      dec_brt = catdec[wbrt]
      mag_brt = catmag[wbrt]
      m_brt = catm[wbrt]
      i100rms = i100rms[wbrt]
      len_deg = $ 
          spike_half_length(mag_brt, m_brt, w4=w4, i100rms=i100rms, $ 
                            radec=[[ra_brt],[dec_brt]])/60. ; deg
  endif

end
