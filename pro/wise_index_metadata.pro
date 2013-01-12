;+
; NAME:
;   wise_index_metadata
;
; PURPOSE:
;   Read WISE L1b image metadata table
;
; CALLING SEQUENCE:
;   indstr = wise_index_metadata(lb, nimage=nimage)
;
; INPUTS:
;   lb      - (l,b) of center of region of interest
;
; OPTIONAL INPUTS:
;   nimage  - return data for nimage images nearest to (l,b) 
;   NOSORT  - keyword indicating sorted index file is not desired
;
; OUTPUTS:
;   indstr  - structure with metadata
;
; EXAMPLES:
;   See wise_write_clean1b.pro
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Feb-22 - Written by Douglas Finkbeiner, CfA
;   2012-Mar-08 - added angle keyword - DPF
;
;----------------------------------------------------------------------
function wise_index_metadata, lb, nimage=nimage, NOSORT=NOSORT, angle=angle, $ 
                              allsky=allsky, w4=w4

  common wise_index_meta_cache, indexfile_cache, ind_sav
  par = psf_par_struc(w4=w4, allsky=allsky)
  indexfile = par.indexfile

; -------- read index file  
  if keyword_set(indexfile_cache) && (indexfile EQ indexfile_cache) then begin 
     splog, 'Restoring meta index from cache'
     indstr = ind_sav
  endif else begin 
     splog, 'Reading index'
     meta = mrdfits(indexfile, 1)
     indexfile_cache = indexfile

     euler, meta.ra, meta.dec, lambda, beta, 3
     msknumsat = keyword_set(w4) ? meta.w4msknumsat : meta.w3msknumsat
; -------- make some cuts, pretend these don't exist
     wgood = where((meta.moon_sep GT 12) AND (meta.saa_sep GT -5) $ 
       AND (msknumsat LT par.nsatmax) AND ~((beta LT -76) AND $ 
       (lambda GT 180) AND (meta.dtanneal LT 1000)) $ 
       AND ~((beta GT 76) AND (lambda LT 180) AND (meta.dtanneal LT 1000)), nw)

     print, nw, ' of', n_elements(meta), ' frames pass cuts'
     indstr = meta[wgood]
     ind_sav = indstr
  endelse

; -------- sort list, trim if nside is set
  if (~ keyword_set(NOSORT)) OR keyword_set(angle) then begin 
     euler, lb[0], lb[1], racen, deccen, 2
     dangle = djs_diff_angle(indstr.ra, indstr.dec, racen, deccen)
     angsind = sort(dangle)

     if keyword_set(nimage) then angsind = angsind[0:nimage-1]
     indstr = indstr[angsind]
  endif

; -------- trip with angle cut if angle is set
  if keyword_set(angle) then begin
     wangle = where(dangle[angsind] LE angle, nangle)
     if nangle GT 0 then indstr = indstr[wangle] else indstr = 0
  endif 

  return, indstr
end
