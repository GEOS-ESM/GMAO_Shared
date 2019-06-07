function makplot (args)

'numargs  'args
 numargs = result

* Initialize INPUT Parameters
* ---------------------------
        num = 0
while ( num < numargs )
        num = num + 1

if( subwrd(args,num) = '-MVAR'     ) ; mvar     = subwrd(args,num+1) ; say 'mvar  = 'mvar  ; endif
if( subwrd(args,num) = '-MNAME'    ) ; mname    = subwrd(args,num+1) ; say 'mname = 'mname ; endif
if( subwrd(args,num) = '-MFILE'    ) ; mfile    = subwrd(args,num+1) ; say 'mfile = 'mfile ; endif
if( subwrd(args,num) = '-MDESC'    ) ; mdesc    = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-MBEGDATE' ) ; bdate    = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-MENDDATE' ) ; edate    = subwrd(args,num+1) ; endif

if( subwrd(args,num) = '-OVAR'     ) ; ovar     = subwrd(args,num+1) ; say 'ovar  = 'ovar  ; endif
if( subwrd(args,num) = '-ONAME'    ) ; oname    = subwrd(args,num+1) ; say 'oname = 'oname ; endif
if( subwrd(args,num) = '-OFILE'    ) ; ofile    = subwrd(args,num+1) ; say 'ofile = 'ofile ; endif
if( subwrd(args,num) = '-ODESC'    ) ; odesc    = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-OBEGDATE' ) ; bdateo   = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-OENDDATE' ) ; edateo   = subwrd(args,num+1) ; endif

if( subwrd(args,num) = '-EXPID'    ) ; expid    = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-PREFIX'   ) ; prefix   = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-SEASON'   ) ; season   = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-OUTPUT'   ) ; output   = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-CLIMATE'  ) ; climate  = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-GC'       ) ; gridcomp = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-MATH'     ) ; math     = subwrd(args,num+1) ; endif

endwhile

* ---------------------------
if(     math = NULL ) ;     math = '' ; endif
if(   season = NULL ) ;   season = '' ; endif
if( gridcomp = NULL ) ; gridcomp = '' ; endif

'set t 1'
'run getenv "GEOSUTIL"'
             geosutil = result

'run getenv "LEVTYPE"'
             LEVTYPE = result
         if( LEVTYPE = 'NULL' ) ; LEVTYPE = DLEVS ; endif

if( prefix != NULL )
    PFX = prefix'_'
else
    PFX = ''
endif
say ''

title = 'NULL'
clevs = 'NULL'
ccols = 'NULL'
dlevs = 'NULL'
dcols = 'NULL'
ccint = 'NULL'

* Check for Existance of NAME Specific plot.rc
* --------------------------------------------
PLOTRC = geosutil'/plots/'mname'/plot.rc'

'!remove   grads.txt'
'!listfile 'PLOTRC' > grads.txt'
checkrc = sublin ( read(grads.txt),2 )
             rc = close(grads.txt)

* Check Variable Attributes from NAME Specific PLOTRC
* ---------------------------------------------------
if( checkrc = PLOTRC )
                        'getresource 'PLOTRC' 'PFX'TITLE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      TITLE' ; endif
                                                   title  = result
                        'getresource 'PLOTRC' 'PFX'CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CINT' ; endif
                                                   ccint  = result
                        'getresource 'PLOTRC' 'PFX'CLEVS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CLEVS' ; endif
                                                   clevs  = result
                        'getresource 'PLOTRC' 'PFX'CCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CCOLS' ; endif
                                                   ccols  = result
                        'getresource 'PLOTRC' 'PFX''LEVTYPE
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      'LEVTYPE ; endif
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PFX'DLEVS'   ; endif
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      DLEVS'   ; endif
                                                   dlevs  = result
                        'getresource 'PLOTRC' 'PFX'DCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      DCOLS' ; endif
                                                   dcols  = result
factor = 1

else

* Check Variable Attributes from Generic PLOTRC
* ---------------------------------------------
PLOTRC = geosutil'/plots/grads_util/plot.rc'

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CBSCALE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CBSCALE' ; endif
                                                                cbscale = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_FACTOR'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_FACTOR' ; endif
                                                                factor = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_TITLE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_TITLE' ; endif
                                                                title = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CINT' ; endif
                                                                ccint = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CCOLS' ; endif
                                                                ccols = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CLEVS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CLEVS' ; endif
                                                                clevs = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_REGRID'
                                                                method = result
endif

'run getenv "CINTDIFF"'
             CINTDIFF  = result
         if( CINTDIFF != 'NULL' ) ; dcols = 'NULL' ; endif

say ''
if( factor = 'NULL' ) ; factor = 1 ; endif
if( title  = 'NULL' )
   'getdesc 'mname
             desc = result
    title = mname':'gridcomp'  'desc
   "rmstring '"title"' '[column]'"
    title = result
   "rmstring '"title"' '__ENSEMBLE__'"
    title = result
endif


* Perform Mathematics if necessary
* --------------------------------
'define qmod = 'mvar''season'*'factor
'define qobs = 'ovar''season'*'factor

'd log10('factor')'
        m = subwrd(result,4)
if( ccols = NULL )
   'set gxout stat'
   'd qmod'
   qminmax = sublin(result,8)
   qmin    = subwrd(qminmax,4)
   qmax    = subwrd(qminmax,5)
   say 'Original  QMIN: 'qmin
   say 'Original  QMAX: 'qmax
   say 'Original STATS: 'result
   'set gxout shaded'
   'd abs('qmin')'
           qmin = subwrd(result,4)
   'd abs('qmax')'
           qmax = subwrd(result,4)
   if( qmin > qmax ) ; qmax = qmin ; endif
   if( qmax > 0 )
      'd log10('qmax')'
       m = subwrd(result,4)
   else
       m = 0
   endif
   say '    Log Factor: 'm
   if( m<0 ) ; m = m-2 ; endif
   'getint 'm
            m = result
   if( m>0 )
       if( m<=2 )
           m = 0
       else
           m = m-2
       endif
   endif
   say 'Field Scaling Factor: 'm
   minv = -m
   'define qmod = qmod * 1e'minv
   'define qobs = qobs * 1e'minv
   'set gxout stat'
   'd qmod'
   qminmax = sublin(result,8)
   qmin    = subwrd(qminmax,4)
   qmax    = subwrd(qminmax,5)
   say 'Final  QMIN: 'qmin
   say 'Final  QMAX: 'qmax
   say 'Final STATS: 'result
   'set gxout shaded'
endif

if( math = LOG )
   'define qmod = log(qmod+0.00001)'
   'define qobs = log(qobs+0.00001)'
endif

'set vpage off'
'set parea off'
'set grid  off'
'set mproj scaled'
'set frame on'
'set xlopts 1 3 .11'
'set ylopts 1 3 .11'

* Count Seasons
* -------------
'set dfile 'ofile
'count "'season'" 'bdateo' 'edateo
 nobs = result

'set dfile 'mfile
'count "'season'" 'bdate' 'edate
 nmod = result

'set t 1'

* Top Panel
* ---------
'set vpage 0 8.5 0.0 11'
'set parea 1.5 7.0 7.70 10.50'
'set grads off'
if( ccols != NULL )
   'set clevs 'clevs
   'set ccols 'ccols
else
   if( ccint != NULL )
      'shades qmod 0 -cint 'ccint
   else
      'shades qmod 0'
   endif
endif
'd qmod'
'set parea 0 8.5 7.0 11'
'cbarn -vert'
'set parea off'

* Middle Panel
* ------------
'set vpage 0 8.5 0.0 11'
'set parea 1.5 7.0 4.30 7.10'
'set grads off'
if( ccols != NULL )
   'set clevs 'clevs
   'set ccols 'ccols
else
   if( ccint != NULL )
      'shades qmod 0 -cint 'ccint
   else
      'shades qmod 0'
   endif
endif
'd qobs'
'set parea off'

* Bottom Panel
* ------------
'set vpage 0 8.5 0.0 11'
'set parea 1.5 7.0 0.90 3.70'
'set grads off'
'rgbset'
'getinfo lon'
         lon = result
*'define obsg = regrid2( qobs,1,1,bs_p1,'lon',-90)'
*'define modg = regrid2( qmod,1,1,bs_p1,'lon',-90)'
'define obsg = qobs'
'define modg = qmod'
'define difg = maskout( modg-obsg,abs(obsg) )'

n = 0
if( dcols != NULL )
   'set clevs 'dlevs
   'set ccols 'dcols
   'd difg'
else
   'stats difg'
     avgdif = subwrd(result,1)
     stddif = subwrd(result,2)
       qmax = stddif/3
   if( qmax > 0 )
      'd log10('qmax')'
       n = subwrd(result,4)
   else
       n = 0
   endif
   say '    Log Factor: 'n
   if( n<0 ) ; n = n-2 ; endif
   'getint 'n
            n = result
   if( n>0 )
       if( n<=2 )
           n = 0
        else
           n = n+2
        endif
   endif
   say 'Diff Scaling Factor: 'n
      'd 'qmax'/1e'n
       cint = subwrd(result,4)
      'shades 'cint
      'define difg = difg/1e'n
      'd difg'
endif
'cbarn -snum 0.55 -xmid 4.25 -ymid 0.4'

'stats maskout(modg,abs(obsg))'
 avgmod = subwrd(result,1)
 stdmod = subwrd(result,2)
'stats maskout(obsg,abs(obsg))'
 avgobs = subwrd(result,1)
 stdobs = subwrd(result,2)
'stats difg'
 avgdif = subwrd(result,1)
 stddif = subwrd(result,2)

'set vpage off'
'set string 1 l 4'
'set strsiz .065'
'draw string 0.05 0.08 ( EXPID:  'expid' )'

'set string 1 c 6'
'set strsiz .14'
'draw string 4.25 10.85 'math' 'title
'set strsiz .10'

if( m != 0 )
   if( m>0 )
      'draw string 4.25 10.62 'mdesc'  'season' ('nmod')  (x 10** -'m')'
   else
      'draw string 4.25 10.62 'mdesc'  'season' ('nmod')  (x 10**'m')'
   endif
else
   'draw string 4.25 10.62 'mdesc'  'season' ('nmod')'
endif
   'draw string 4.25 7.22 'odesc'  'season' ('nobs')  ('climate')'

if( n != 0 )
   'draw string 4.25 3.80 Difference (Top-Middle)  (x 10**'n')'
else
   'draw string 4.25 3.80 Difference (Top-Middle)'
endif

* Print Beginning and Ending Dates
* --------------------------------
                date = getdate (bdate)
bmnthm = subwrd(date,1)
byearm = subwrd(date,2)
                date = getdate (edate)
emnthm = subwrd(date,1)
eyearm = subwrd(date,2)
                date = getdate (bdateo)
bmntho = subwrd(date,1)
byearo = subwrd(date,2)
                date = getdate (edateo)
emntho = subwrd(date,1)
eyearo = subwrd(date,2)

* --------------------------------

'set string 1 l 4'
'set strsiz .08'
'draw string 0.050 10.50 Beg: 'bmnthm' 'byearm
'draw string 0.050 10.35 End: 'emnthm' 'eyearm
'draw string 0.050 7.10 Beg: 'bmntho' 'byearo
'draw string 0.050 6.95 End: 'emntho' 'eyearo

'draw string 0.050 9.85 Mean: 'avgmod
'draw string 0.050 9.70  Std: 'stdmod
'draw string 0.050 6.45 Mean: 'avgobs
'draw string 0.050 6.30  Std: 'stdobs
'draw string 0.050 3.05 Mean: 'avgdif
'draw string 0.050 2.90  Std: 'stddif

if( output != 'NULL' )
if( math = LOG )
    mathparm = '_LOG'
else
    mathparm = ''
endif
'myprint -name 'output'/'mname''mathparm'_'PFX''oname'.'season
endif

'set mproj latlon'
return

function getdate (date,month,year)
       num = 1
       bit = substr(date,num,1)
while( bit != '' )
       num = num+1
       bit = substr(date,num,1)
endwhile
       loc = num-7
     month = substr(date,loc  ,3)
      year = substr(date,loc+3,4)
return month' 'year
