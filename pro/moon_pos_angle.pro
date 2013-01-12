;+
; NAME:
;   moon_pos_angle
;
; PURPOSE:
;   given a list of pointing (lon, lat) and corresponding list of
;   moon (lon, lat), compute the position angle east of north of the great
;   circle from each pointing direction to each corresponding moon direction
;
; CALLING SEQUENCE:
;   pos_angle = moon_pos_angle(lon, lat, moonra, moondec)
;
; INPUTS:
;   lon     - list of pointing longitudes (deg)  
;   lat     - list of pointing latitudes (deg)
;   moonlon  - list of moon longitudes (deg)
;   moonlat - list of moon latitudes (deg)
;  
; OUTPUTS:
;   moon position angle east of north in degrees, in range [0, 360)
;
; EXAMPLES:
;   see pa_test and wise_moon_radec wrapper routines below
;
; COMMENTS:
;   Also included routine moon_ra_dec, which I used to calculate the
;   moon (ra, dec) as viewed from WISE for each WISE pointing, thus
;   generating the ra_moon_wise and dec_moon_wise fields in metadata
;   index structures. These values are not the same as simply using
;   moonpos to get geocentric (ra, dec) of Moon, and can differ by
;   up to ~1 deg
;
; REVISION HISTORY:
;   2011-Jun-12 - written by Aaron Meisner
;----------------------------------------------------------------------
function moon_pos_angle, lon, lat, moonlon, moonlat

  n = n_elements(lon)
  ind = lindgen(n,3)

  u_pointing = ll2uv([[lon], [lat]])
  u_moon = ll2uv([[moonlon], [moonlat]])

;----- u_pointing cross u_moon gives normal to the plane containing pointing
;      unit vector and moon pointing vector, call this moon_normal
  moon_normal = crossprod(u_pointing, u_moon) 
  moon_normal_norm = sqrt(total(moon_normal^2, 2))
  moon_normal[ind] = moon_normal[ind]/moon_normal_norm[ind MOD n]

;----- u_pointing cross (0,0,1) gives normal to the plane containing pointing
;      unit vector and +z axis (north pole), call this pole_normal
  north_pole = [[0.], [0.], [1.]]
  pole_normal = crossprod(u_pointing, north_pole)
  pole_normal_norm = sqrt(total(pole_normal^2, 2))
  pole_normal[ind] = pole_normal[ind]/pole_normal_norm[ind MOD n]
;----- u_pointing cross moon_normal = line given by intersection of 
;      moon/pointing plane and plane perpendicular to pointing = moon_dir
  moon_dir = crossprod(moon_normal, u_pointing) ;guaranteed unit norm
;----- u_pointing cross pole_normal = line given by intersection of 
;      pole/pointing plane and plane perpendicular to pointing = pole_dir
  pole_dir = crossprod(pole_normal, u_pointing) ;guaranteed unit norm

  cos_pa = total(moon_dir*pole_dir, 2)
  sin_pa = total((crossprod(moon_dir, pole_dir))*(u_pointing), 2)
  pos_angle = (atan(sin_pa, cos_pa)*!radeg + 360.) MOD 360.

  return, pos_angle
end

;----- Routine to infer moon (ra, dec) as viewed from WISE spacecraft
;      given the information available in the L1b fits headers, the 
;      angular separation of the moon and the position angle east of
;      north of the great circle along which moon lies.
;----- This routine is not vectorized, but I don't care since I only
;      want to compute moon ra/dec once for each pointing, then store result.
function moon_ra_dec, ra_pointing, dec_pointing, theta, moonsep

;----- theta is moon direction east of celestial north, in degrees
  ll = [[ra_pointing], [dec_pointing]]

;----- u is the unit vector from origin to direction of pointing
  u = ll2uv(ll, /DOUBLE)

  theta = -theta/!radeg  ; conventions
  moonsep = moonsep/!radeg
;----- u_pole = rotation of north pole by theta about u direction
  u_pole = transpose([u[0]*u[2]*(1-cos(theta))+u[1]*sin(theta), $ 
                    u[1]*u[2]*(1-cos(theta))-u[0]*sin(theta), $ 
                    cos(theta) + (u[2]^2)*(1-cos(theta))])

;----- u_moon = direction of unit vector to follow in plane
;      perpendicular to u, in order to arrive at the moon position
  u_moon = $
    crossprod(u, crossprod(u_pole, u))
;----- ensure that u_moon is a unit vector, possible division by zero??
  u_moon = u_moon/sqrt(total(u_moon^2))
;----- now compute unit vector perpendicular to plane containing u and u_moon
  ax = crossprod(u, u_moon)
;----- then rotate u about ax by moon_sep degrees
;      uv_moon = vector at the end of which the moon lies on unit sphere
;      with lon, lat <-> ra, dec
  rot = fltarr(3, 3)
  rot[0, 0] = cos(moonsep) + (ax[0]^2)*(1-cos(moonsep))
  rot[0, 1] = ax[1]*ax[0]*(1-cos(moonsep)) + ax[2]*sin(moonsep)
  rot[0, 2] = ax[2]*ax[0]*(1-cos(moonsep)) - ax[1]*sin(moonsep)
  rot[1, 0] = ax[0]*ax[1]*(1-cos(moonsep)) - ax[2]*sin(moonsep)
  rot[1, 1] = cos(moonsep) + (ax[1]^2)*(1-cos(moonsep))
  rot[1, 2] = ax[2]*ax[1]*(1-cos(moonsep)) + ax[0]*sin(moonsep)
  rot[2, 0] = ax[0]*ax[2]*(1-cos(moonsep)) + ax[1]*sin(moonsep)
  rot[2, 1] = ax[1]*ax[2]*(1-cos(moonsep)) - ax[0]*sin(moonsep)
  rot[2, 2] = cos(moonsep) + (ax[2]^2)*(1-cos(moonsep))

  uv_moon = rot ## u
  return, uv2ll(uv_moon)

end

pro wise_moon_radec, moonra, moondec, allsky=allsky

  if ~keyword_set(allsky) then $ 
    indstr = mrdfits('$WISE_DATA/index-metadata-L1b.fits', 1) $ 
  else $
    indstr = mrdfits('$WISE_DATA/index-allsky-L1b.fits', 1)

  n = n_elements(indstr)
  moonra = fltarr(n)
  moondec = fltarr(n)
  for i=0L, n-1 do begin
    ll = moon_ra_dec(indstr[i].ra, indstr[i].dec, indstr[i].moonpa, $ 
                      indstr[i].moon_sep)
    if (i MOD 10000) EQ 0 then print, i
    moonra[i] = ll[0]
    moondec[i] = ll[1]
  endfor
end

pro pa_test

  indstr = mrdfits('$WISE_DATA/index-metadata-L1b.fits', 1)
  pos_angle = moon_pos_angle(indstr.ra, indstr.dec, indstr.ra_moon_wise, $ 
                             indstr.dec_moon_wise)

;----- indstr.moonpa is the parameter specified in header as moon position
;      angle east of north, want to compare with my calculation
  diff = indstr.moonpa - pos_angle
  print, 'Worst agreement = ', max(abs(diff)), ' deg'
;----- apparently there are no issues at high absolute dec, or even when
;      the moon separation gets very close to 180 deg. (which is problematic
;      since pointing and moon directions no longer define a plane, and
;      there isn't a well defined pos. angle of the moon anyway, as
;      following any great circle for 180 deg works)
end
