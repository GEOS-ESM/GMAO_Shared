function formula (args)

'numargs  'args
 numargs = result

* OBS EXPORTS & GCS
* -----------------
if( subwrd(args,1) = '-OBS' )
    return "MSDWLWRFCS:ERA5 MSNLWRFCS:ERA5"
endif

* Formula:  lwgup = msdwlwrfcs - (-1 * msnlwrfcs)
* Import of msnlwrfcs is scaled by -1 in VERIFICATION.ERA5.rc for comparison with lwgnetc
* ------------------------------------
       n  = 1
while( n <= numargs )
       oname.n = subwrd(args,n)
       n  = n + 1
endwhile

say 'define qobs = 'oname.1' + 'oname.2
    'define qobs = 'oname.1' + 'oname.2

return
