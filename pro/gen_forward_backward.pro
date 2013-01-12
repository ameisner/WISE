;+
; NAME:
;   gen_forward_backward
;
; PURPOSE:
;   Generate a flag (to be appended to index structure) that indicates
;   whether a given pointing is in the "forward" or "backward" direction,
;   by comparing ecliptic longitude of Earth as viewed from
;   Sun to the ecliptic longitude of the pointing. Forward is 
;   defined as being in the direction of earth's orbit.
;
; CALLING SEQUENCE:
;   gen_forward_backward, forward, allsky=
;
; KEYWORDS:
;   allsky - if set, compute forward/backward flag for all pointings
;            in all-sky release (default is just prelim release)
;
; OUTPUTS:
;   forward - flag indicating forward (1) backward (0) pointing
;
; REVISION HISTORY:
;   2012-Jun-12 - Written by Aaron Meisner
;----------------------------------------------------------------------
pro gen_forward_backward, forward, allsky=allsky

  wisedata = getenv('WISE_DATA')
  index = keyword_set(allsky) ? 'index-allsky-L1b.fits' : $ 
    'index-metadata-L1b.fits'
  indstr = mrdfits(concat_dir(wisedata, index), 1)
  lambda_earth = $ 
    wise_ephemeris_model(indstr.mjd, /INTERP, beta_earth=beta_earth)

  euler, indstr.ra, indstr.dec, lambda_pointing, beta_pointing, 3
  
;----- this is the unit vector FROM sun TO earth
  uv_earth = ll2uv([[lambda_earth], [beta_earth]])
  uv_pointing = ll2uv([[lambda_pointing], [beta_pointing]])

;----- take cross product of sun->earth direction and pointing direction
;      z component sign tells you forward vs. backward
  cross = crossprod(uv_earth, uv_pointing)
  forward = cross[*, 2] GT 0

end
