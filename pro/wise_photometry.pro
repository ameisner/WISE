;+
; NAME:
;   wise_photometry
;
; PURPOSE:
;   Do forced aperture photometry on WISE for a list of (RA, dec)
;
; CALLING SEQUENCE:
;   wise_photometry, ra, dec, bandname, str
;
; INPUTS:
;   ra, dec    - locations to make measurements
;   bandname   - ['1', '2', '3', '4'] for 3.5, 4.6, 12, 22 micron
;
; OUTPUTS:
;   str        - output structure containing:
;     wise_flux      [naper, nband]: aperture flux, nMgy
;     wise_flux_ivar [naper, nband]: inverse variance of flux
;     wise_skyflux   [nband]       : flux in sky annulus [nMgy/pixel]
;     wise_skyrms    [nband]       : sky flux RMS in annulus [nMgy/pixel]
;     wise_coaddid                 : WISE Coadd ID, e.g. '2178m016_aa11'
;
; EXAMPLES:
;   See wisephot.pro 
;
; COMMENTS:
;   No re-centering is done. 
;
; REVISION HISTORY:
;   2011-May-15 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_photometry, ra, dec, bandname, str

  dW = [2.683, 3.319, 5.242, 6.604]
  ABfac = 10.^(-dW/2.5)  ; multiple vegamaggies by this number to get AB

  nband = n_elements(bandname)

; FWHM values measured off images are (8.5,9.5,12,19 asec)
; -------- Aperture radii (these are 88, 95, 98% containment at 3.5 mu)
  aper = [7., 10., 12.5,  15., 20.] ; radius [arcsec ]
  skyrad = [20.0, 28.0] ; arcsec

  sigasec = [8.5, 9.5, 12.0, 19.0]
  sigpix0 = sigasec/1.375/2 ; factor of 2 is for 2048 images
  sigpix = sqrt(sigasec^2-6.0^2)/1.375/2
  errfac = sqrt(4*!pi*sigpix^2)
  t0 = systime(1)

; -------- define structure for WISE photometry
  f0 = fltarr(n_elements(aper), nband) 
  s0 = fltarr(nband)
  str0 = {wise_flux: f0, $
          wise_flux_ivar:f0, $
          wise_skyflux: s0, $
          wise_skyrms: s0, $
          wise_coaddid: ''}

  str = replicate(str0, n_elements(ra))

; -------- Read WISE WCS information
  wcs = mrdfits('$WISE_DIR/etc/wise-wcs.fits', 1)

; -------- Match WISE pointings with SDSS objects
  wise_findimage, ra, dec, wcs, ind1, ind2

; -------- Loop over unique list of WISE plates
  plateind = ind2[uniq(ind2, sort(ind2))]

  for i=0L, n_elements(plateind)-1 do begin 

; -------- WISE Atlas directory is hardwired to 2048 images here
     idstr = wcs[plateind[i]].idstr
     dir = '$WISE_DATA/L3a_2048/'+idstr+'_aa11'
     int_name = concat_dir(dir, idstr+'_aa11-w'+bandname[0]+'-int-3.fits')

; -------- read WCS header
     hdr = headfits(int_name)

; -------- get WCS astrometry from header
     extast, hdr, astr
     sz = [sxpar(hdr, 'NAXIS1'), sxpar(hdr, 'NAXIS2')]

; -------- transform RA,dec to x,y to see if objects are away from edge
     objind = ind1[where(ind2 EQ plateind[i])] 

     ad2xy, ra[objind], dec[objind], astr, x, y

     print, 'Considering ', n_elements(x), ' objects'
     pad = 24  ; pixels
     x0 = pad
     x1 = sz[0]-pad-1
     
     w = where((round(x) GE x0) AND (round(x) LE x1) AND $
               (round(y) GE x0) AND (round(y) LE x1), nw)
     print, 'idstr: ', idstr, nw, ' objects'
     
; -------- if there are objects more than pad pixels from edge, call djs_phot
     if nw GT 0 then begin 
        objw = objind[w]

        for iband=0L, nband-1 do begin 

           bandstr = bandname[iband]
           int_name = concat_dir(dir, idstr+'_aa11-w'+bandstr+'-int-3.fits')
           unc_name = concat_dir(dir, idstr+'_aa11-w'+bandstr+'-unc-3.fits')

           if file_test(int_name) EQ 0 then message, 'Cannot find file: '+int_name
           print, int_name

           phot = wise_forced_photometry(ra[objw], dec[objw], int_name, unc_name, $
                                         aper=aper, skyrad=skyrad)

; -------- if we get objects back, then put them in output structure
           if keyword_set(phot) then begin 
              str[objw].wise_flux[*, iband]      = phot.flux*ABfac[iband]
              str[objw].wise_flux_ivar[*, iband] = 1./(phot.flerr*errfac[iband]*ABfac[iband])^2
              str[objw].wise_skyflux[iband]      = phot.skyval*ABfac[iband]
              str[objw].wise_skyrms[iband]       = phot.skyrms*ABfac[iband]
              str[objw].wise_coaddid             = phot.coaddid
           endif 
        endfor
     endif 
  endfor

; -------- declare victory
  print, n_elements(ra), ' objects processed on', n_elements(plateind), ' images.'
  print, 'Elapsed time: ', systime(1)-t0

  return
end
