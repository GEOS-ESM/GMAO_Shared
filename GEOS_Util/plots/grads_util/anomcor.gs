function anomcor (args)

f      = subwrd(args,1)
a      = subwrd(args,2)
c      = subwrd(args,3)
lonbeg = subwrd(args,4)
lonend = subwrd(args,5)
latbeg = subwrd(args,6)
latend = subwrd(args,7)
level  = subwrd(args,8)
file   = subwrd(args,9)
tag    = subwrd(args,10)


say 'f = 'f
say 'a = 'a
say 'c = 'c
say 'lonbeg = 'lonbeg
say 'lonend = 'lonend
say 'latbeg = 'latbeg
say 'latend = 'latend
say 'level  = 'level
say 'file   = 'file

'set dfile 'file
'sett'
'set lev   'level

* Compute Region Boundaries
* -------------------------
'set lon 'lonbeg
'run getinfo 'xpos
              xbeg = result
'set lon 'lonend
'run getinfo 'xpos
              xend = result

'set lat 'latbeg
'run getinfo 'ypos
              ybeg = result
'set lat 'latend
'run getinfo 'ypos
              yend = result

* Define Anomaly Correlation
* --------------------------
'set x 1'
'set y 1'
'define fmcbar = aave('f'-'c',x='xbeg',x='xend',y='ybeg',y='yend')'
'define amcbar = aave('a'-'c',x='xbeg',x='xend',y='ybeg',y='yend')'

'define varfmc  = aave( pow( ('f'-'c')-fmcbar,2),x='xbeg',x='xend',y='ybeg',y='yend')'
'define varamc  = aave( pow( ('a'-'c')-amcbar,2),x='xbeg',x='xend',y='ybeg',y='yend')'
'define covfac  = aave(    ( ('f'-'c')-fmcbar )*( ('a'-'c')-amcbar ),x='xbeg',x='xend',y='ybeg',y='yend')'
'define stdfmc  = sqrt( varfmc )'
'define stdamc  = sqrt( varamc )'

'define acc'tag' = covfac/(stdfmc*stdamc) '

return
