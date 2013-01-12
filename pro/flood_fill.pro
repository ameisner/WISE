;+
; NAME:
;   flood_fill
;
; PURPOSE:
;   determine region connected to some starting pixel within boolean array   
;
; CALLING SEQUENCE:
;   flood_fill, x, y, bitmap
;
; INPUTS:
;   x      - integer x value of starting node
;   y      - integer y value of starting node
;   bitmap - boolean image in which you want to know the contiguous
;            region connected to (x, y)
;   
; OUTPUTS:
;   bitmap - modified version of input in which the contiguous region
;            desired has been negated
;   
; EXAMPLES:
;   see routine flood_fill_example below
;
; COMMENTS:
;   this is "four-way" flood fill, where pixels only touching at the
;   one corner do not count as connected, as opposed to "eight-way"
;   flood fill
;
; REVISION HISTORY:
;   2011-Sep-23 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro flood_fill, x, y, bitmap

  forward_function flood_fill
  sz = size(bitmap, /DIM)

  if (x LT 0) OR (x GE sz[0]) OR (y LT 0) OR (y GE sz[1]) then return
  if (bitmap[x, y] EQ 0) then return

  bitmap[x, y] = 0
  flood_fill, x+1, y, bitmap
  flood_fill, x-1, y, bitmap
  flood_fill, x, y+1, bitmap
  flood_fill, x, y-1, bitmap

end

pro flood_fill_example

  sz = 11
  bitmap = bytarr(sz, sz)
  bitmap[(sz/2-1):(sz/2+1), (sz/2-1):(sz/2+1)] = 1
  bitmap[0,*] = 1
  bitmap[sz-1,*] = 1
  bitmap[*,0] = 1
  bitmap[*,sz-1] = 1
  
;---- bitmap
;   1   1   1   1   1   1   1   1   1   1   1
;   1   0   0   0   0   0   0   0   0   0   1
;   1   0   0   0   0   0   0   0   0   0   1
;   1   0   0   0   0   0   0   0   0   0   1
;   1   0   0   0   1   1   1   0   0   0   1
;   1   0   0   0   1   1   1   0   0   0   1
;   1   0   0   0   1   1   1   0   0   0   1
;   1   0   0   0   0   0   0   0   0   0   1
;   1   0   0   0   0   0   0   0   0   0   1
;   1   0   0   0   0   0   0   0   0   0   1
;   1   1   1   1   1   1   1   1   1   1   1
;----
  orig = bitmap
  flood_fill, sz/2, sz/2, bitmap
  contiguous = bitmap NE orig  

;---- contiguous
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   1   1   1   0   0   0   0
;   0   0   0   0   1   1   1   0   0   0   0
;   0   0   0   0   1   1   1   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;   0   0   0   0   0   0   0   0   0   0   0
;----
end
