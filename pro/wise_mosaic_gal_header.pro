;+
; NAME:
;   wise_mosaic_gal_header
;
; PURPOSE:
;   Make FITS header with WCS info
;
; CALLING SEQUENCE:
;   hdr = wise_mosaic_gal_header(flux, cd, lgal=, bgal=)
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
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2011-Nov-27 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function wise_mosaic_gal_header, flux, cd, lgal=lgal, bgal=bgal

  if n_elements(lgal) NE 1 then message, 'must set lgal, bgal'

  sz = (size(flux, /dimen))[0]
; -------- define header for mosaic
  mkhdr, hdr, flux

  sxaddpar, hdr, 'CTYPE1', 'GLON-ZEA'
  sxaddpar, hdr, 'CTYPE2', 'GLAT-ZEA'
  sxaddpar, hdr, 'CD1_1', -cd ; / Degrees / Pixel
  sxaddpar, hdr, 'CD2_1', 0.0
  sxaddpar, hdr, 'CD1_2', 0.0
  sxaddpar, hdr, 'CD2_2', cd  ;  / Degrees / Pixel
  sxaddpar, hdr, 'CRPIX1', sz/2 ; Reference Pixel in X
  sxaddpar, hdr, 'CRPIX2', sz/2 ; Reference Pixel in Y
  sxaddpar, hdr, 'CRVAL1', double(lgal) ;  / R.A. (degrees) of reference pixel
  sxaddpar, hdr, 'CRVAL2', double(bgal); / Declination of reference pixel
  sxaddpar, hdr, 'LONPOLE', 180.000000000 ; / Native longitude of Celestial pole
  sxaddpar, hdr, 'LATPOLE', 90.0000000000 ; / Celestial latitude of native pole

  return, hdr
end
