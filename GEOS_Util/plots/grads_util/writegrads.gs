function writegrads(args)

'numargs  'args
 numargs = result

'q gxout'
   gxout = sublin(result,4)
   gxout = subwrd(gxout,6)

'run getenv ZFLIP'
            zflip = result

*******************************************
* Note: To Write Data in Reverse Order
*       run setenv ZFLIP ON
*       args:  list of variables to write
*             (Default, all variables)  
*******************************************

'setx'
'sety'
'setz'

'getinfo xdim'
         xdim = result
'getinfo ydim'
         ydim = result
'getinfo zdim'
         zdim = result
'getinfo edim'
         edim = result
'getinfo tdim'
         tdim  = result

'getinfo dlon'
         dlon = result
'getinfo dlat'
         dlat = result

if( xdim = 1 ) ; dlon = 1 ; endif

'set x 1'
'set y 1'
'getinfo lon'
         lonbeg = result
'getinfo lat'
         latbeg = result

'setx'
'sety'

z = 1
levels = ''
while( z<=zdim )
  'set z 'z
  'getinfo level'
           level = result
  levels = levels' 'level
z = z + 1
endwhile

'set t 1'
'getinfo date'
      begdate = result
'getinfo tinc'
         tinc = result
'getinfo tunit'
         tunit = result
'getinfo label'
         label = result

'getinfo nvars'
         nvars = result

if( tunit = 'year'  ) ; tunit = yr ; endif
if( tunit = 'month' ) ; tunit = mo ; endif
if( tunit = 'day'   ) ; tunit = dy ; endif
if( tunit = 'hour'  ) ; tunit = hr ; endif

* Write grads.ctl
* ---------------
'!remove grads.ctl'
'!touch  grads.ctl'
'!echo dset ^grads.data                           >> grads.ctl'
'!echo title 'label'                              >> grads.ctl'
'!echo undef 1e15                                 >> grads.ctl'
'!echo xdef 'xdim' linear 'lonbeg' 'dlon'         >> grads.ctl'
'!echo ydef 'ydim' linear 'latbeg' 'dlat'         >> grads.ctl'
'!echo zdef 'zdim' levels 'levels'                >> grads.ctl'
'!echo tdef 'tdim' linear 'begdate' 'tinc''tunit' >> grads.ctl'

if( numargs = 0 )
  '!echo vars 'nvars'   >> grads.ctl'
else
  '!echo vars 'numargs' >> grads.ctl'
endif

* --------------------------------------------------------------------

'set gxout fwrite'
'set fwrite grads.data'

'set undef 1e15'

t=1
while(t<=tdim)
  'set t 't

n=1
while(n<=nvars)

'run getvarz 'n' noprint'
    var  = result
   name  = subwrd(result,1)
   zdim  = subwrd(result,2)

   found = true

if( numargs != 0 )
    found = false
    k=1
    while( k<=numargs ) 
      if( subwrd(args,k) = name ) ; found = true ; endif
    k=k+1
    endwhile
endif

 if( found = true )

  'q file'
   loc = 6 + n 
   varinfo = sublin(result,loc)
  '!remove DESC.txt'
  '!echo "'varinfo'" > DESC.txt'
  'run getenv "DESC"'
    desc = result
  '!remove DESC.txt'
   say 'Writing: 'desc

   if( name = slp ) ; zdim = 0 ; endif

   if( zdim = 0 )
      '!echo "'desc'" >> grads.ctl'
*      say 'Writing Variable: 'name
      'set z 1'
       e = 1
       while( e<=edim )
      'set e 'e
      'd 'name
       e = e + 1
       endwhile
   else
      '!echo "'desc'" >> grads.ctl'

      if(zflip != 'ON' )
         z=1
         while(z<=zdim)
            'set z 'z
            'getinfo level'
             lev = result
*            say 'Writing Variable: 'name' for Level: 'lev
             e = 1
             while( e<=edim )
            'set e 'e
            'd 'name
             e = e + 1
             endwhile
         z=z+1
         endwhile
      else
         z=zdim
         while(z>=1)
            'set z 'z
            'getinfo level'
             lev = result
*            say 'Writing Variable: 'name' for Level: 'lev
             e = 1
             while( e<=edim )
            'set e 'e
            'd 'name
             e = e + 1
             endwhile
         z=z-1
         endwhile
      endif

   endif
*      say '  '
 endif

n=n+1
endwhile
t=t+1
endwhile

'!echo endvars >> grads.ctl'
'disable fwrite'
'set gxout 'gxout
return
