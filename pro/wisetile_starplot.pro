;+
; NAME:
;   wisetile_starplot
;
; PURPOSE:
;   overplot WISE catalog source positions on an ISSA tile
;
; CALLING SEQUENCE:
;   wisetile_starplot, tnum, tpath_prelim= , tpath_allsky= , allsky= , 
;       maglim= , noplot= , xp= , yp= , xa= , ya= , dirty= 
;
; INPUTS:
;   tnum         - ISSA tile number in interval [1, 430]
;
; KEYWORDS:
;   tpath_prelim - directory containing preliminary release tiles
;   tpath_allsky - directory containing allsky release tiles
;   allsky       - if plotting, set to view allsky rather than prelim
;                  release tile
;   maglim       - 2-element array giving lower magnitude limit
;                  (element 0), and upper magnitude limit (element 1)
;   noplot       - set to suppress image display
;   dirty        - display tile with stars instead of clean version
;
; OPTIONAL OUTPUTS:
;   xp           - retrieve preliminary release source x positions within tile
;   yp           - retrieve preliminary release source y positions within tile
;   xa           - retrieve allsky release source x positions within tile
;   ya           - retrieve allsky release source y positions within tile
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Aug-04 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro wisetile_starplot, tnum, tpath_prelim=tpath_prelim, $ 
                       tpath_allsky=tpath_allsky, allsky=allsky, $
                       maglim=maglim, noplot=noplot, xp=xp, yp=yp, $ 
                       xa=xa, ya=ya, dirty=dirty, _extra=extra

  if ~keyword_set(tpath_prelim) then tpath_prelim='/n/panlfs/ameisner/tilef'
  if ~keyword_set(tpath_allsky) then $ 
      tpath_allsky='/n/panlfs/ameisner/tile-allsky-flat'
  ext_no = keyword_set(dirty)

  fprelim = $ 
      concat_dir(tpath_prelim, 'wise_'+string(tnum, format='(I03)')+'.fits')
  fallsky = $
      concat_dir(tpath_allsky, 'wise_'+string(tnum, format='(I03)')+'.fits')

  inprelim = file_test(fprelim)
  if ~inprelim then $ 
      print, 'tile '+string(tnum, format='(I03)')+' not in preliminary release'

  COMMON catalogs, catp, cata 
  if (n_elements(catp) EQ 0) then begin
      catp = mrdfits('$WISE_DATA/w3_catalog.fits', 1)
      cata = mrdfits('$WISE_DATA/w3_catalog-allsky.fits', 1)
  endif

  ima = readfits(fallsky, ha, ex=ext_no)

  nx = sxpar(ha, 'NAXIS1')
  ny = sxpar(ha, 'NAXIS2')
  racen = sxpar(ha, 'CRVAL1')
  deccen = sxpar(ha, 'CRVAL2')
  extast, ha, astr
  
  magl = keyword_set(maglim) ? maglim[0]  : -1e6
  magu = keyword_set(maglim) ? maglim[1]  : 1e6
  maxsep = sqrt((nx*sxpar(ha, 'CD1_1')/2)^2+(ny*sxpar(ha, 'CD2_2')/2)^2)

  if inprelim then begin
      dangp = djs_diff_angle(catp.ra, catp.dec, racen, deccen)
      w = where(dangp LT maxsep)
      ad2xy, catp[w].ra, catp[w].dec, astr, xp, yp
      wgood = where((xp LT nx-0.5) AND (xp GT -0.5) AND $  
                    (yp LT ny-0.5) AND (yp GT -0.5) AND $ 
                    ((catp[w].flag AND 256) EQ 0) AND $ 
                    (catp[w].mag LT magu) AND (catp[w].mag GT magl))
      xp = xp[wgood] & yp = yp[wgood]
  endif

  danga = djs_diff_angle(cata.ra, cata.dec, racen, deccen)
; ----- don't bother with any binary searching...
  w = where(danga LT maxsep)
  ad2xy, cata[w].ra, cata[w].dec, astr, xa, ya
  wgood = where((xa LT nx-0.5) AND (xa GT -0.5) AND $  
                (ya LT ny-0.5) AND (ya GT -0.5) AND $ 
                ((cata[w].flag AND 256) EQ 0) AND $ 
                (cata[w].mag LT magu) AND (cata[w].mag GT magl))
  xa = xa[wgood] & ya = ya[wgood]

  if ~keyword_set(noplot) then begin
      if (~keyword_set(allsky) AND inprelim) then begin
          imp = readfits(fprelim, hp, ex=ext_no)
          atv, imp, header=hp, _extra=extra
          atvplot, xa, ya, psym=1
          atvplot, xp, yp, psym=1, color=djs_icolor('green')
      endif else begin
          atv, ima, header=ha, _extra=extra
          atvplot, xa, ya, psym=1
          if inprelim then $ 
              atvplot, xp, yp, psym=1, color=djs_icolor('green')
      endelse 
  end

end
