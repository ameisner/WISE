;+
; NAME:
;   wise_filter_cog
;
; PURPOSE:
;   filter WISE source catalog based on reported aperture magnitudes
;
; CALLING SEQUENCE:
;   mask = wise_filter_cog(cat, blimit=, thresh=)
;
; INPUTS:
;   cat  - structure contatining a  subset of rows from WISE source 
;          catalog, and having fields w3mpro, w3mag_[1-6], w3sigm_1 
;
; KEYWORDS:
;   blimit - w3mpro value for brightest source which will be subject
;            to filtering based on its curve of growth
;
; OUTPUTS:
;   mask - array with same number of elements as there are rows in 
;          input cat structure, where 1 indicates that a source is
;          "good" and 0 indicates that a source is "bad"
;
; EXAMPLES:
;   mask = wise_filter_cog(cat)
;   wgood = where(mask)
;   cat = cat[wgood]
;
; COMMENTS:
;  this method of filtering WISE sources derives from the basic ideas of: 
;
;  http://wise2.ipac.caltech.edu/staff/jarrett/wise/nebula/fragmented_xsrc.html
;
;  An implementation of those exact cuts can now be found in the subroutine
;  jarrett_filter_cog
;
;  updated wise_filter_cog rejects sources that are better fit
;  by a constant surface brightness C.O.G. than a PSF C.0.G., instead
;  of simply cutting on number of points in the C.0.G. brighter than 
;  expected for a PSF by thresh mags
; 
;  model fits are assessed in terms of sum of absolute deviations
;  rather than a formal chi^2 or likelihood since so many of the aperture 
;  mag error bars are undefined or suspect
;
; REVISION HISTORY:
;   2011-Sep-16 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function aper_par_struc, w4=w4

; ----- WISE source catalog aperture radii in arcsec
;         w3mag_1 w3mag_2 w3mag_3 w3mag_4 w3mag_5 w3mag_6 w3mag_7 w3mag_8
  aper = [5.5,    8.25,   11.0,   13.75,  16.5,   19.25,  22.0,   24.75] 
  if keyword_set(w4) then aper *= 2
  par_str = {aper: aper}
  return, par_str

end

pro model_const_sb, mag, w4=w4

 ; simplistic model of C.O.G resulting from measurements in
 ; WISE apertures over a region of constant surface brightness

  par = aper_par_struc(w4=w4)
  aper = par.aper

  flux = aper^2
  mag = -2.5*alog10(flux)
  mag -= mag[0]

end

function wise_filter_cog, cat, blimit=blimit, faint=faint, w4=w4

  w4 = keyword_set(w4)
  if ~keyword_set(blimit) then blimit = w4 ? 2.5 : 7.2
  faint = w4 ? (cat.w4mpro GT blimit) : (cat.w3mpro GT blimit)
  uncbad = w4 ? (cat.w4sigm_1 EQ -9999) : (cat.w3sigm_1 EQ -9999)

; ----- preliminary release
; dm = [-0.586, -0.897, -1.057, -1.139, -1.182]
; ----- allsky release
  dm_real = w4 ? [-0.504, -0.792, -0.940, -1.003, -1.031] : $ 
                 [-0.585, -0.894, -1.053, -1.135, -1.177]

  model_const_sb, dm_const, w4=w4

  dm_const = dm_const[1:5]

  mag_1 = w4 ? cat.w4mag_1 : cat.w3mag_1
  dm_21 = w4 ? (cat.w4mag_2-mag_1) : (cat.w3mag_2-mag_1)
  dm_31 = w4 ? (cat.w4mag_3-mag_1) : (cat.w3mag_3-mag_1)
  dm_41 = w4 ? (cat.w4mag_4-mag_1) : (cat.w3mag_4-mag_1)
  dm_51 = w4 ? (cat.w4mag_5-mag_1) : (cat.w3mag_5-mag_1)
  dm_61 = w4 ? (cat.w4mag_6-mag_1) : (cat.w3mag_6-mag_1)

  totdev_real = abs(dm_21-dm_real[0])+abs(dm_31-dm_real[1])+$
                abs(dm_41-dm_real[2])+abs(dm_51-dm_real[3])+$ 
                abs(dm_61-dm_real[4])
  totdev_const = abs(dm_21-dm_const[0])+abs(dm_31-dm_const[1])+$
                 abs(dm_41-dm_const[2])+abs(dm_51-dm_const[3])+$ 
                 abs(dm_61-dm_const[4])

  bad = ((uncbad) OR (totdev_real GT totdev_const)) AND (faint)

  return, ~bad
end

function jarrett_filter_cog, cat, blimit=blimit, thresh=thresh

;  preliminary release
;  dm = [-0.586, -0.897, -1.057, -1.139, -1.182]
; allsky release
  dm = [-0.585, -0.894, -1.053, -1.135, -1.177]
  
  if ~keyword_set(blimit) then blimit = 7.2
  if ~keyword_set(thresh) then thresh=0.5

  faint = (cat.w3mpro GT blimit)
  badtot = ((cat.w3mag_2-cat.w3mag_1) LT (dm[0]-thresh)) + $ 
           ((cat.w3mag_3-cat.w3mag_1) LT (dm[1]-thresh)) + $ 
           ((cat.w3mag_4-cat.w3mag_1) LT (dm[2]-thresh)) + $
           ((cat.w3mag_5-cat.w3mag_1) LT (dm[3]-thresh)) + $ 
           ((cat.w3mag_6-cat.w3mag_1) LT (dm[4]-thresh))

  bad = ((badtot GE 2) OR (cat.w3sigm_1 EQ -9999)) AND faint
  
  return, ~bad

end
