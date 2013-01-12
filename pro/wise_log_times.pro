;+
; NAME:
;   wise_log_times, str
;
; PURPOSE:
;   grep image processing times from the wise log files
;
; CALLING SEQUENCE:
;   wise_log_times, str
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
;   2012-Mar-01 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_log_times, str

  logdir = '~ameisner/wise/pro'
  flist = file_search(concat_dir(logdir, 'b*.log'), count=nfile)

  print, 'Reading', nfile, ' log files...'

  str = replicate({fname:'', time:fltarr(100000), ntime:0L}, nfile)

  for i=0L, nfile-1 do begin 

     cmd = 'cat '+flist[i]+' | grep -v latent | grep -v HREBIN | grep -v fits'
     print, cmd
     spawn, cmd, result
     w = where(strlen(result) EQ 16, ntime)
     time = float(result[w])
    
;     plot,time,yr=[0.1,60],/ylog
     str[i].fname = fileandpath(flist[i])
     str[i].ntime = ntime
     str[i].time[0:ntime-1] = time
  endfor

  return
end


pro doplot, str

  !p.multi = [0, 1, 9]
  nfile = n_elements(str) 
  yr=[0.5,120]
  xr = [0, max(str.ntime)]
  yl = 60./2^findgen(4)

  colors = ['', 'red', 'green', 'cyan']

  k = 0L
  while 1 do begin 
     plot, str[k].time[0:str[k].ntime-1], xr=xr, yr=yr, /ylog, /yst, $
       chars=1.5, ymargin=[1.5, 1], /xst
     djs_xyouts, xr[1], yl[0], str[k].fname, align=1.2
     k++
     while 1 do begin 
        k4 = k mod 4
        if (k EQ nfile) or (k4 EQ 0) then break
        djs_oplot, str[k].time[0:str[k].ntime-1], color=colors[k4]
        djs_xyouts, xr[1], yl[k4], str[k].fname, align=1.2, color=colors[k4]
        k++
     endwhile
     
     if k EQ nfile then break
  endwhile

  totfiles = round(total(str.ntime))
  print, 'Total files processed so far:', totfiles
  w = where(str.time NE 0)
  print, 'Total processing time [s]:', total((str.time)[w])
  print, 'Total processing time [h]:', total((str.time)[w])/3600
  print, 'Max time per image:   ', max((str.time)[w])
  print, 'Mean time per image:  ', mean((str.time)[w])
  print, 'Median time per image:', median((str.time)[w])
  return
end

