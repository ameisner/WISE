;+
; NAME:
;   wise_forced_photometry
;
; PURPOSE:
;   Read WISE image and call djs_phot for a list of positions
;
; CALLING SEQUENCE:
;   str = wise_forced_photometry(ra, dec, int_name, unc_name, aper=, skyrad=)
;
; INPUTS:
;   int_name  - intensity map name  (FITS file)
;   ra, dec   - coordinates of points to photometer
;    
; OPTIONAL INPUTS:
;   unc_name  - uncertainty map name
;
; KEYWORDS:
;   aper      - array of aperture radii [asec]
;   skyrad    - 2 elemnt array of sky annulus radii [asec]
;
; OUTPUTS:
;   str       - structure array of flux, sky, etc. 
;
; EXAMPLES:
;   see wisephot.pro
;
; COMMENTS:
;   We do not recenter on the (ra, dec) -- this is forced aperture
;                                          photometry 
;   All output flux quantities are nMgy (but not AB system!)
;
;   Only call this with positions that are actually inside the image, or
;     you will get zeros returned!
;
; REVISION HISTORY:
;   2011-May-14 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function wise_forced_photometry, ra, dec, int_name, unc_name, $
                                 aper=aper, skyrad=skyrad

; -------- read WISE image
  if file_test(int_name) EQ 0 then message, 'Cannot find file '+int_name
  im  = readfits(int_name, hdr, /silent)

; -------- If NaNs, just bail for now.  Should handle this more gracefully
  if total(finite(im) EQ 0) NE 0 then begin
     message, 'NaN alert!', /info
     return, 0
  endif

; -------- define inverse variance map for flux uncertainties
  if keyword_set(unc_name) then begin 
     if file_test(unc_name) EQ 0 then begin 
        print, 'Skipping uncertainty file' 
     endif else begin
        unc = readfits(unc_name, /silent)
        if min(unc LE 0) then message, 'uncertainty less LE zero???'
        ivar = 1./unc^2
     endelse
  endif
  sz = size(im, /dimen)

; -------- extra WCS astrometry from header
  extast, hdr, astr
  asecperpix = abs(astr.cdelt[0]*3600)
  magzp   = sxpar(hdr, 'MAGZP')   ; magnitude zero point (Vega system)
  coaddid = sxpar(hdr, 'COADDID') ; WISE coadd id

; -------- transform RA,dec to x,y
  ad2xy, ra, dec, astr, x, y
  if min(x) le -1 then print, 'WARNING: out of bounds'
  if min(y) le -1 then print, 'WARNING: out of bounds'
  if max(x) GE sz[0] then print, 'WARNING out of bounds'
  if max(y) GE sz[1] then print, 'WARNING out of bounds'

; --------  do forced photometry at each (x,y)
  flux = djs_phot(x, y, aper/asecperpix, skyrad/asecperpix, $
                  im, ivar, $
                  calg='none', flerr=flerr, skyval=skyval, skyrms=skyrms)

  f = fltarr(n_elements(aper)) 
  str0 = {flux: f, $
          flerr:f, $
          skyval: 0., $
          skyrms: 0., $
          coaddid: ''}
  
  str = replicate(str0, n_elements(skyval))

; -------- convert all fluxes to nMgy, like SDSS fluxes
;           Note that these are NOT AB magnitudes. 
  nmgyfactor = 10.^((22.5-magzp)/2.5)

  str.flux    = transpose(flux)*nmgyfactor
  str.flerr   = transpose(flerr)*nmgyfactor
  str.skyval  = skyval*nmgyfactor
  str.skyrms  = skyrms*nmgyfactor
  str.coaddid = coaddid

  return, str
end


pro callit


  cat = query_irsa_cat([41.0,-0.30],radius=400, cat=catname)

  wise_photometry, cat.ra, cat.dec, 'w1', str

  return
end


pro gaussphot

  big = readfits(fname, h)
  wise_rebin, big, h, im, hout
  
  fwhm = 8.5/2.75
  sm = fastconv(im, fwhm, 1, /nodisp)


  catname = 'wise_prelim_p3as_psd'
  cat = query_irsa_cat([150.8,69.7],radius=3000, cat=catname)


  return
end
