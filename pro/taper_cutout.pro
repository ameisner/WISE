;+
; NAME:
;   taper_cutout
;
; PURPOSE:
;   taper the edges of an image cutout so that its profile transitions smoothly
;   to zero
;
; CALLING SEQUENCE:
;   tap = taper_cutout(cutout, feat=, bright=, allsky=, w4=)
;
; INPUTS:
;   cutout - image cutout to be tapired
;
; KEYWORDS:
;   feat   - PSF feature being tapired, must be one of 'wings',
;          'latent', 'latent2', 'ghost' ... here 'wings' includes
;          extended wings+core
;   bright - set for bright PSF wings (default faint)
;   allsky - set for all-sky release
;   w4     - set for W4 (default W3)
;
; OUTPUTS:
;   tap    - cutout tapered smoothly/symmetrically to zero at edges
;
; EXAMPLES:
;   
; COMMENTS:
;   not yet meant to be general enough to handle W4
;
; REVISION HISTORY:
;   2011-Nov-25 - Aaron Meisner
;----------------------------------------------------------------------
function taper_cutout, cutout, feat=feat, bright=bright, allsky=allsky, $ 
                       w4=w4, band=band

  w4 = keyword_set(w4)
  par = psf_par_struc(w4=w4, allsky=allsky, /everything, band=band)
; feat should either be 'wings' (here meant to include
; wings+core), 'latent', 'latent2', or 'ghost'
  if ~keyword_set(feat) then feat = 'wings'
  case feat of
      'wings'   : begin
                      if keyword_set(bright) then begin
                          r0 = w4 ? 142 : 162 ; L1b pix
                          r1 = w4 ? 132 : 147
                          psfhalf = par.psfpix/2
                          whalf = par.szwings/2
                          weight = taper_weight(par.szwings, par.szwings, $ 
                                                r0, r1)
                          tap = cutout
                          tap[(psfhalf-whalf):(psfhalf+whalf), $ 
                              (psfhalf-whalf):(psfhalf+whalf)] *= weight
                      endif else begin
                          r0 = 57 ; same value now applies to both W3, W4 !!!
                          r1 = 42 ; same value now applies to both W3, W4 !!! 
                          weight = taper_weight(par.pfaint, par.pfaint, $ 
                                                r0, r1)
                          tap = weight*cutout
                      endelse
                  end
      'ghost'   : begin
                      r0 = 43
                      r1 = 38
                      weight = taper_weight(par.xgpix, par.ygpix, r0, r1)
                      tap = weight*cutout
                  end
      'latent'  : begin
; ----- eventually absorb these into psf_par_struc.pro and generalize
;       to w4
                      r0s = [50, -1, 162, 95] ; L1b pix
                      r1s = [40, -1, 147, 85]
                      r0 = r0s[band-1]
                      r1 = r1s[band-1]
                      weight = taper_weight(par.szlat, par.szlat, r0, r1)
                      tap = weight*cutout
                  end 
      'latent2' : begin
                      r0 = w4 ? 94 : 135
                      r1 = w4 ? 84 : 120
                      weight = taper_weight(par.szlat2, par.szlat2, r0, r1)
                      tap = weight*cutout
                  end
  endcase

  return, tap
end
