;+
; NAME:
;   wise_l1b_clean
;
; PURPOSE:
;   Clean a single level 1b image
;
; CALLING SEQUENCE:
;   image = wise_l1b_clean(raw, msk, h, hout=, bitmask=, dirty=,
;     binfac=, starmask=, interp1b=)
;
; INPUTS:
;   raw    - raw image from int file
;   msk    - bit mask from msk file
;   h      - FITS header containing WCS astrometry
;
; KEYWORDS:
;   binfac - bin down by this factor (should be power of 2)
;   interp1b - interpolate over star artifact mask before smoothing/rebinning
;
; OUTPUTS:
;   image   - image with stars and CRs removed, bad regions interpolated
;             over, smoothed to 15 arcsec FWHM, and rebinned by binfac.
;   hout    - output header, takes account of removal of npad pixels at
;             edge and rebinning
;
; OPTIONAL OUTPUTS:
;   bitmask - deprecated for now, but should keep track of missing regions
;   dirty   - image without point sources cleaned up, but rebinned and
;             smoothed 
;   starmask - bit mask for masking related to stars
;              2^0 = saturated stellar core
;              2^1 = bright star ghost
;              2^2 = bright region of stellar profile, right now with 
;                    threshold of 250 DN above background
;              2^3 = central pixels of latent
;              2^4 = pixels near core of PSF subtraction which have
;                    already been interpolated over by
;                    wise_l1b_substar
;              2^5 = pixel is bright part of a ghost, where "bright" 
;                    threshold is lower than that for 2^2 bit
;              2^6 = pixels at location of known SSO have been interpolated
;                    over
;              2^7 = marks location of an artifact which is not
;                    present at all epochs of observation
;              2^8 = location of second latent
;              2^9 = location of third latent
;              2^10 = location of fourth latent
;              2^11 = bright SSO ghost
;              2^12 = bright SSO latent
;              2^13 = bright compact source optical diffraction spike
;              2^14 = good pixels (not in flagged in static mask) that are
;                     saturated for any reason. Intended to flag
;                     extremely bright, saturated nebulosity. Should
;                     use "AND" mask for this bit.
;                     
; EXAMPLES:
;   see calls in wise_write_clean1b
;
; COMMENTS:
;   This is a first attempt to fill NaNs, remove CRs, and remove faint
;     stars.  Bright stars still need attention. 
;  
;   Interpretation of the msk array is entirely empirical.  Read the
;     docs some time!
; 
;   dirty part still not tested !!!
;
; REVISION HISTORY:
;   2011-Nov-24 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function wise_l1b_clean, raw, msk, h, hout=hout, bitmask=bitmask, $
                         dirty=dirty, binfac=binfac, starmask=starmask, $
                         interp1b=interp1b, allsky=allsky, w4=w4

  ; ------- interpolate over NaNs, CRs, and other bad pixels
  intim = wise_l1b_maskinterp(raw, msk, badmask=badmask, satmask=satmask)
  ; ------- subtract stars
  clean = wise_l1b_substar(intim, h, badmask, starmask=starmask, $ 
      allsky=allsky, w4=w4)
  ; ------- correct latents
  wise_correct_latent, clean, h, latentmask=latentmask, allsky=allsky, w4=w4
  ; ------- create diffraction spike mask
  spkmsk = diff_spike(h, allsky=allsky, w4=w4)
  ; ------- remove moving objects
  wise_remove_sso, clean, h, smask=ssomask, w4=w4, allsky=allsky
  starmask = fix(starmask) ; latent can have some 2^[>7] bits set
  starmask = starmask OR latentmask OR 64*(ssomask AND 1) OR $ 
             2048*((ssomask AND 2) NE 0) OR 4096*((ssomask AND 4) NE 0) OR $ 
             8192*(spkmsk NE 0)
  ; ------- if requested, perform star mask interpolation at L1b level
  if keyword_set(interp1b) then begin
    kpix = 23 ;  best choice for this value?
    kern = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)
    brt  = dilate((starmask AND 44) NE 0, kern)
    intx = djs_maskinterp(clean, brt, iaxis=0, /const)
    inty = djs_maskinterp(clean, brt, iaxis=1, /const)
    clean = (intx+inty)/2
  endif
  ; ------- "smooth" masks
  kpix   = 5
  kern   = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)
  starmask = dilate(starmask AND 1, kern) + $ 
    2*dilate(starmask AND 2, kern) + 4*dilate(starmask AND 4, kern) + $ 
    8*dilate(starmask AND 8, kern) + 16*dilate(starmask AND 16, kern) + $ 
    32*dilate(starmask AND 32, kern) + 64*dilate(starmask AND 64, kern) + $
    128*dilate(starmask AND 128, kern) + $ 
    256*dilate((starmask AND 256) NE 0, kern) + $
    512*dilate((starmask AND 512) NE 0, kern) + $ 
    1024*dilate((starmask AND 1024) NE 0, kern) + $
    2048*dilate((starmask AND 2048) NE 0, kern) + $
    4096*dilate((starmask AND 4096) NE 0, kern) + $ 
    8192*dilate((starmask AND 8192) NE 0, kern) + $ 
    16384*satmask
; -------- 
; Because of the funny PSF, the images smooth as if they have FWHM 3.6
; pix.  sqrt((15./2.754)^2-3.6^2) = 4.08
  fwhm_pix = keyword_set(w4) ? 4.18 : 4.08 ; smoothing kernel FWHM in pixels
  smth = fastconv(clean, fwhm_pix, 1, /nodisp, /silent)

; -------- if dirty image requested, smooth it also
  if arg_present(dirty) then $ 
      dirty = fastconv(intim, fwhm_pix, 1, /nodisp, /silent)

; -------- the smoothing wrecks the outer 8 pixels of the image
  npad = 8
  sz = size(smth, /dimen)
  outsz = sz-npad*2
  out   = smth[npad:sz[0]-1-npad, npad:sz[1]-1-npad]
  starmask = starmask[npad:sz[0]-1-npad, npad:sz[1]-1-npad] ; trim this as well
  if arg_present(dirty) then dirty = dirty[npad:sz[0]-1-npad, npad:sz[1]-1-npad]
; -------- modify header for removal of outer npad pixels
  hout = h
  sxaddpar, hout, 'NAXIS1', outsz[0]
  sxaddpar, hout, 'NAXIS2', outsz[1]
  sxaddpar, hout, 'CRPIX1', outsz[0]/2.0d + 0.5d
  sxaddpar, hout, 'CRPIX2', outsz[1]/2.0d + 0.5d

; -------- now rebin by binfac, updating header accordingly
  if keyword_set(binfac) then begin 
     if binfac NE 2 then message, 'Warning: Only tested for binfac = 2'
     out1 = temporary(out)
     hdr = hout
; note: this seems to be good only at the .02 pixel level - not sure why.
     hrebin, out1, hdr, out, hout, outsize=outsz/binfac
     if keyword_set(dirty) then $
       dirty = rebin(dirty, outsz[0]/binfac, outsz[1]/binfac)
  ; ------- rebin the star mask
     starmask = $ 
 (rebin(float(starmask AND 1),outsz[0]/binfac,outsz[1]/binfac) GE 0.5) $ 
+(rebin(float(starmask AND 2),outsz[0]/binfac,outsz[1]/binfac) GE 1)*2 $ 
+(rebin(float(starmask AND 4),outsz[0]/binfac,outsz[1]/binfac) GE 2)*4 $ 
+(rebin(float(starmask AND 8),outsz[0]/binfac,outsz[1]/binfac) GE 4)*8 $ 
+(rebin(float(starmask AND 16),outsz[0]/binfac,outsz[1]/binfac) GE 8)*16 $ 
+(rebin(float(starmask AND 32),outsz[0]/binfac,outsz[1]/binfac) GE 16)*32 $
+(rebin(float(starmask AND 64),outsz[0]/binfac,outsz[1]/binfac) GE 32)*64 $
+(rebin(float(starmask AND 128),outsz[0]/binfac,outsz[1]/binfac) GE 64)*128 $
+(rebin(float(starmask AND 256),outsz[0]/binfac,outsz[1]/binfac) GE 128)*256 $
+(rebin(float(starmask AND 512),outsz[0]/binfac,outsz[1]/binfac) GE 256)*512 $
+(rebin(float(starmask AND 1024),outsz[0]/binfac,outsz[1]/binfac) GE 512)*1024+(rebin(float(starmask AND 2048),outsz[0]/binfac,outsz[1]/binfac) GE 1024)*2048+(rebin(float(starmask AND 4096),outsz[0]/binfac,outsz[1]/binfac) GE 2048)*4096+(rebin(float(starmask AND 8192),outsz[0]/binfac,outsz[1]/binfac) GE 4096)*8192+(rebin(float(starmask AND 16384),outsz[0]/binfac,outsz[1]/binfac) GE 8192)*16384
     starmask = fix(starmask)
  endif 

  return, out
end
