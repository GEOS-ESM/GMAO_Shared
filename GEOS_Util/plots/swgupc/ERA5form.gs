function formula (args)

'numargs  'args
 numargs = result

* OBS EXPORTS & GCS
* -----------------
if( subwrd(args,1) = '-OBS' )
    return "SLRSFC:SOLAR RSCS:SOLAR"
endif

* Formula:  swgupc = msdwswrf - msnswrf
* ---------------------------------
       n  = 1
while( n <= numargs )
       oname.n = subwrd(args,n)
       n  = n + 1
endwhile

say 'define qobs = 'oname.1' - 'oname.2
    'define qobs = 'oname.1' - 'oname.2

return
