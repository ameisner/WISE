;+
; NAME:
;   wise_build_metadata
;
; PURPOSE:
;   build index file with metadata information from WISE metadata table
;   
; CALLING SEQUENCE:
;   wise_build_metadata, indstr, w4=w4
;
; OUTPUTS:
;   index structure with metadata information for all frames in
;   all-sky released. This structure also
;   containts full file names for files downloaded from preliminary release
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;   2011-Mar-17 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro wise_build_metadata, indstr, w4=w4, missing=missing

  w4 = keyword_set(w4)
  metafile = $ 
      '/n/panlfs/ameisner/wise_allsky.wise_allsky_4band_p1bs_frm4293.tbl'
  dldir = w4 ? '/n/wise/dfink/wise/band4/4band_p1bm_frm' : $ 
               '/n/home08/dfink/wisedata/downloads/4band_p1bm_frm'
  F = w4 ? 'A,I,F,F,D,D,X,A,D,X,A,F' : 'A,I,F,F,D,D,A,X,D,A,X,F'
  readcol, metafile, scan_id, frame_num, moon_sep, saa_sep, ra, dec, bandrun, $
      mjd, msknumsat, dtanneal, F=F, skipline=43

  ntot = n_elements(scan_id)
  wgood = where(bandrun NE 'null', ngood)

  indstr = replicate({scan_id:'', frame_num:0, ra:0.d, dec:0.d,  mjd: 0.d, $ 
                      fname:'',  dtanneal:0., moon_sep:0., saa_sep:0.}, ntot) 

  indstr.scan_id = scan_id
  indstr.frame_num = fix(frame_num)
  indstr.moon_sep = moon_sep
  indstr.saa_sep = saa_sep
  indstr.ra = ra
  indstr.dec = dec
  indstr.mjd = mjd
  indstr.dtanneal = dtanneal
  if total(msknumsat EQ 'null') NE 0 then $ 
      msknumsat[where(msknumsat EQ 'null')] = '-9999'
  if w4 then begin 
      addstr = replicate({w4msknumsat:0L}, ntot)
      addstr.w4msknumsat = long(msknumsat)
  endif else begin
      addstr = replicate({w3msknumsat:0L}, ntot)
      addstr.w3msknumsat = long(msknumsat)
  endelse
  indstr = struct_addtags(indstr, addstr)
; ----- full file location
  fname = dldir+'/'+strmid(scan_id,4,2)+'/'+scan_id+'/'+$
      string(frame_num,format='(I03)')+'/'+scan_id+$
      string(frame_num,format='(I03)')+'-'+(w4 ? 'w4': 'w3')+'-int-1b.fits'
  indstr.fname = fname
  indstr = indstr[wgood]
; ----- depending on file system, this may be a very bad way to 
;       identify downloaded intensity images, presumably could also be
;       more clever with find, xargs or do brute force file_test
  spawn, 'find '+dldir+" -name '*int-1b.fits'", fname_dl
  
  matchlist, fileandpath(indstr.fname), fileandpath(fname_dl), mall, mdl

  nmissing = ngood-n_elements(mall)
  print, '# of L1b exposures missing from download: ', nmissing

  if (nmissing NE 0) then begin
      missing = bytarr(ngood)+1
      missing[mall] = 0
      missing = fileandpath(indstr[where(missing)].fname)
  endif

  indstr = indstr[mall]
  sind = sort(indstr.mjd)
  indstr = indstr[sind]

end


