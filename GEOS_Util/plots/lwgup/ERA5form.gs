function formula (args)

'numargs  'args
 numargs = result

* OBS EXPORTS & GCS
* -----------------
if( subwrd(args,1) = '-OBS' )
    return " MSDWLWRF:ERA5 MSNLWRF:ERA5"
endif

* Formula:  lwgup = msdwlwrf - (-1 * msnlwrf)
* Import of msnlwrf is scaled by -1 in VERIFICATION.ERA5.rc for comparison with lwgnet
* ------------------------------------
       n  = 1
while( n <= numargs )
       oname.n = subwrd(args,n)
       n  = n + 1
endwhile

say 'define qobs = 'oname.1' + 'oname.2
    'define qobs = 'oname.1' + 'oname.2

return
