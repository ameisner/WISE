;+
; NAME:
;   wise_pairwise_reject
;
; PURPOSE:
;   Use pairwise overlap variance to reject bad images
;   
; CALLING SEQUENCE:
;   wise_pairwise_reject, pair, imcube=imcube
;
; INPUTS:
;   pair    - structure made by wise_pariwise_compare
;   imcube  - image cube from wise_read_imcube()
;
; OUTPUTS:
;   pair    - updated - ndiff=0 means bad
;
; COMMENTS:
;   See wise_make_mosaic for usage
;   May need to increase maxndiff some day
;   thresholds are arbitrary
; 
; REVISION HISTORY:
;   2011-Nov-27 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_pairwise_reject, pair, imcube=imcube

;  sz = size(imcube, /dimen)
;  small=rebin(imcube,125,125,sz[2])

  varcut = 150
  maxndiff = 20

  imtot = n_elements(pair)

  while 1 do begin   ; repeat until we break out

; -------- find worst match

;     worst = max(pair.medvar, maxind)
     worst = max(pair.medvar, maxind)
     fworst = fileandpath(pair[maxind].fname)
     print, 'Worst: ', worst, '    maxind:', maxind, '    name: ', fworst
     if worst LT varcut then break   ; know when to stop

     w = where(pair.ind EQ maxind, nw)
     

; -------- remove those matches from list
     vind = w mod maxndiff
     pind = w / maxndiff
     for i=0L, nw-1 do begin 
        pi = pind[i]
        vi = vind[i]
        pair[pi].var[vi] = -1
        pair[pi].ind[vi] = -1
        thisind = where(pair[pi].ind GE 0, nind)
        if nind GE 3 then begin 
           if pair[pi].ndiff GT 0 then begin 
              pair[pi].ndiff --
              pair[pi].totvar = total(pair[pi].var[thisind])/pair[pi].ndiff
              pair[pi].medvar = median(pair[pi].var[thisind])
           endif 
        endif else begin 
           pair[pi].ndiff = 0
           pair[pi].medvar = -1
           pair[pi].totvar = -1
           message, 'not enough comparisons left!', /info
        endelse
     endfor
     
; -------- remove bad one
     pair[maxind].ndiff = 0
     pair[maxind].medvar = -1
     pair[maxind].totvar = -1

     wgood = where(pair.ndiff GT 0, ngood)

;     if ((i mod 10) EQ 0) or done then begin 
;        sind = sort(-pair[wgood].totvar)
;        atv, grid(small[*,*,wgood[sind]], /med),/al, /s
;        wait, .1
;     endif 
  endwhile

  splog, 'ngood =', round(total(pair.ndiff NE 0)), ' of ', imtot, $ 
    ' L1b exposures'
        
  return
end
