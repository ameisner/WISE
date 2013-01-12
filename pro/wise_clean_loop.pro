;+
; NAME:
;   wise_clean_loop
;
; PURPOSE:
;   given an arbitrary L1b index structure, loop over its file names
;   and write out cleaned files
;
; CALLING SEQUENCE:
;   wise_clean_loop, indstr, outpath=, allsky=
;
; INPUTS:
;   indstr  - index structure with field fname giving list of files to
;             be processed
;
; KEYWORDS:
;   outpath - directory into which to write cleaned files
;
; EXAMPLES:
;   see wise_clean_batch and wise_write_clean1b
;
; COMMENTS:
;   i made this code a separate routine for the sake of modularity,
;   since it may not always be the case that we want to clean
;   files by a range of MJD (as in wise_clean_batch) or by proximity
;   to a some (l,b) pointing (wise_write_clean1b)
;
;   also, it seems better to have a common place for the file writing
;   code rather than two nearly identical versions in wise_clean_batch
;   and wise_l1b_clean
;
; REVISION HISTORY:
;   2012-Sep-16 - Code written predominantly by Doug in
;                 wise_write_clean1b and wise_clean_batch moved here for
;                 the sake of modularity by Aaron
;
;----------------------------------------------------------------------
pro wise_clean_loop, indstr, outpath=outpath, allsky=allsky, w4=w4, $ 
                     tmpdir=tmpdir

  nimage = n_elements(indstr) 
  splog, 'Working on ', nimage, ' images...'
  print
  flist = indstr.fname
  ct = 0L

  t0 = systime(1)
; -------- loop over files
  for i=0L, nimage-1 do begin 

; -------- read image
     thisfile = fileandpath(flist[i])
     print, i, '  ', thisfile
     thisdir = strmid(thisfile, 0, 4)
     thispath = concat_dir(outpath, thisdir)
     file_mkdir, thispath  ; 10 musec if dir already exists
     outname = concat_dir(thispath, thisfile)

; -------- see if we need to work on file
     if file_test(outname) EQ 0 then begin 
; -------- touch the file (takes ~ 1 msec).  This allows multiple
;          threads to play nice together. 
        spawn, ['touch', outname], /noshell
        
        ct++
        raw = lfs_fits_access(flist[i], h, /silent)

        mname = repstr(flist[i], 'int', 'msk')
        mname = repstr(mname, '.fits', '.fits.gz')
        msk = lfs_fits_access(mname, /silent)
; -------- if intensity image or mask is missing, give up on this file
        if (n_elements(raw) EQ 1) OR (n_elements(msk) EQ 1) then begin 
            print, 'missing mask or intensity image !!!'
; -------- in case something strange is going on with the file system
;          such that files cannot be read but can be created, wait to ensure
;          you aren't touching thousands of files in rapid-fire succession
            wait, 5
            continue
        endif
        
; -------- process image
        clean = wise_l1b_clean(raw, msk, h, hout=hout, binfac=2, $ 
                               dirty=dirty, starmask=starmask, allsky=allsky, $
                               w4=w4)
        
; -------- write image
        randstr = string(randomu(iseed, /double)*1e12, format='(I12.12)')
        tempext = '.temp.'+randstr
        tempname = outname+tempext

; -------- could do this to build file on local disk and then move...
; tempname = '/tmp/'+fileandpath(outname)+tempext
        if keyword_set(tmpdir) then $ 
            tempname = concat_dir(tmpdir, fileandpath(tempname))

        t2 = systime(1)
        writefits, tempname, clean, hout, /checksum

; -------- fix up header for extension HDUs
        sxaddpar, hout, 'XTENSION', 'IMAGE', ' IMAGE extension',before='SIMPLE'
        sxdelpar, hout, 'SIMPLE'
        writefits, tempname, dirty, hout, /checksum, /append
        sxaddpar, hout, 'BITPIX', 16
        sxaddpar, hout, 'BZERO', 32768, after='NAXIS2'
        writefits, tempname, uint(starmask), hout, /checksum, /append
        file_move, tempname, outname, /overwrite
        print, 'Write time:', systime(1)-t2
     endif

  endfor

  print, ct, ' images processed in', systime(1)-t0, ' seconds'

end
