;+
; NAME:
;   init_source_catalog
;
; PURPOSE:
;   read in/assemble/filter WISE compact source catalog
;
; CALLING SEQUENCE:
;   init_source_catalog, w4=, allsky=
;
; KEYWORDS:
;   w4   - set for W4 catalog (default is W3)
;   
; EXAMPLES:
;   see wise_starlist
;
; COMMENTS:
;   creates a common block with necessary catalog parameters
;   
; REVISION HISTORY:
;   2012-Nov-8 - Aaron Meisner
;----------------------------------------------------------------------
pro init_source_catalog, w4=w4, allsky=allsky

  COMMON CATALOG, catra, catdec, catmag, mjdmin, mjdmax, catm
  if n_elements(catra) EQ 0 then begin
      par = psf_par_struc(w4=w4, allsky=allsky)
      cat  = mrdfits(par.catfile, 1)
      wgood = wise_filter_catalog(cat, w4=w4)
      cat      = cat[wgood]
      catra    = cat.ra
      catdec   = cat.dec
      catmag   = keyword_set(w4) ? cat.w4mpro : cat.w3mpro
      wise_substitute_bright, catmag, catra, catdec, fluxflag, w4=w4
      wkeep    = where(~fluxflag)
      catmag   = catmag[wkeep]
      catra    = catra[wkeep]
      catdec   = catdec[wkeep]
      mjdmin   = keyword_set(w4) ? cat[wkeep].w4mjdmin : cat[wkeep].w3mjdmin
      mjdmax   = keyword_set(w4) ? cat[wkeep].w4mjdmax : cat[wkeep].w3mjdmax
      catm     = keyword_set(w4) ? cat[wkeep].w4m : cat[wkeep].w3m
      sind     = sort(catdec)
      catra    = catra[sind]
      catdec   = catdec[sind]
      catmag   = catmag[sind]
      mjdmin   = mjdmin[sind]
      mjdmax   = mjdmax[sind]
      catm     = catm[sind]
      delvarx, cat ; free this memory ASAP
  endif

end
