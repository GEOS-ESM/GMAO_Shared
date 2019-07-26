function formula (args)

'numargs  'args
 numargs = result

* OBS EXPORTS & GCS
* -----------------
if( subwrd(args,1) = '-OBS' )
    return "MSNLWRF:ERA5"
endif

* Formula:  lwgnet = lwgnet
* ---------------------------
       n  = 1
while( n <= numargs )
       oname.n = subwrd(args,n)
       n  = n + 1
endwhile

say 'define qobs = 'oname.1
    'define qobs = 'oname.1

return
