function formula (args)

'numargs  'args
 numargs = result

* OBS EXPORTS & GCS
* -----------------
if( subwrd(args,1) = '-OBS' )
    return "MSDWLWRFCS:ERA5"
endif

* Formula:  lwgdwnc = lwgdwnc
* ---------------------------
       n  = 1
while( n <= numargs )
       oname.n = subwrd(args,n)
       n  = n + 1
endwhile

say 'define qobs = 'oname.1
    'define qobs = 'oname.1

return
