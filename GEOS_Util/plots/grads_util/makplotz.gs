function makplotz (args)

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
if( subwrd(args,num) = '-RGFILE'   ) ; rfile    = subwrd(args,num+1) ; endif

endwhile
* ---------------------------

if( math = NULL ) ; math = '' ; endif

'set t 1'
'run getenv "GEOSUTIL"'
             geosutil = result

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
axlim = 'NULL'
ylab  = 'NULL'
grid  = 'NULL'

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
                        'getresource 'PLOTRC' 'PFX'CLEVS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CLEVS' ; endif
                                                   clevs  = result
                        'getresource 'PLOTRC' 'PFX'DLEVS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      DLEVS' ; endif
                                                   dlevs  = result
                        'getresource 'PLOTRC' 'PFX'CCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CCOLS' ; endif
                                                   ccols  = result
                        'getresource 'PLOTRC' 'PFX'DCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      DCOLS' ; endif
                                                   dcols  = result
                        'getresource 'PLOTRC' 'PFX'AXLIM'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      AXLIM' ; endif
                                                   axlim  = result
                        'getresource 'PLOTRC' 'PFX'YLAB' 
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      YLAB'  ; endif
                                                   ylab   = result
                        'getresource 'PLOTRC' 'PFX'GRID' 
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      GRID'  ; endif
                                                   grid   = result
else

* Check Variable Attributes from Generic PLOTRC
* ---------------------------------------------
PLOTRC = geosutil'/plots/grads_util/plot.rc'

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CBSCALE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CBSCALE' ; endif
                                                                cbscale = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_TITLE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_TITLE' ; endif
                                                                title = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CCOLS' ; endif
                                                                ccols = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CLEVS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CLEVS' ; endif
                                                                clevs = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_REGRID'
                                                                method = result

if( axlim = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_AXLIM'math ; endif ; axlim =  result

endif

say ''

if( axlim != NULL )
axmin = subwrd(axlim,1)
axmax = subwrd(axlim,2)
endif


* Perform Mathematics if necessary
* --------------------------------
'define qmod = 'mvar''season
'define qobs = 'ovar''season

m = 0
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
'set mproj scaled'
'set grid  on'
'set frame on'
'set xlopts 1 3 .11'
'set ylopts 1 3 .11'

* Count Seasons
* -------------
'set dfile 'mfile
'count "'season'" 'bdate' 'edate
 nmod = result

if( title  = 'NULL' )
   'getdesc 'mname
             desc = result
    title = mname':'gridcomp'  'desc
   "rmstring '"title"' '[column]'"
    title = result
   "rmstring '"title"' '__ENSEMBLE__'"
    title = result
endif

'set dfile 'ofile
'count "'season'" 'bdateo' 'edateo
 nobs = result

* Set DFILE to REGRID File to enable zonal averages with proper dimensions
* ------------------------------------------------------------------------
'set dfile 'rfile
'setlons'
'setlats'
'define qobs = maskout(qobs,abs(qmod))'
'define qmod = maskout(qmod,abs(qobs))'

'makez qobs z'
'makez qmod z'

'setlons'
'setlats'

'set dfile 'mfile
'set t   1'
'set lon 0'

'set vpage 0 8.5 0 11'
'set parea .8 8 4 9'
if( ylab  != NULL ) ; 'set ylab  'ylab  ; endif
if( grid  != NULL ) ; 'set grid  'grid  ; endif

if( axlim = NULL ) 
'set dfile 'ofile
'set t   1'
'set lon 0'
   'run minmax qobsz'
           qmax = subwrd(result,1)
           qmin = subwrd(result,2)
'set dfile 'rfile
'set t   1'
'set lon 0'
   'run minmax qmodz'
           dmax = subwrd(result,1)
           dmin = subwrd(result,2)
       if( dmax > qmax ) ; qmax = dmax ; endif
       if( dmin < qmin ) ; qmin = dmin ; endif
       qave = (qmax + qmin)/2
       qdel = (qmax - qave)*1.10
      axmax =  qave + qdel
      axmin =  qave - qdel
      axlim = ''axmin' 'axmax
endif

'set axlim 'axlim

'set grads off'
'set cmark  0'
'set cstyle 1'
'set ccolor 1'
'd qobsz'
'set cmark  0'
'set cstyle 1'
'set ccolor 4'
'd qmodz'

* Draw Zero Line
* --------------
if( axlim != NULL )
if( axmin*axmax < 0 )
'set cmark  0'
'set cstyle 1'
'set cthick 1'
'set ccolor 2'
'd lon-lon'
endif
endif

'setlons'

'set vpage off'
'set string 4 c 6'
'set strsiz .11'
*'xlabel 1 4.25 10.5'
'draw string 4.25  10.5 EXPID: 'expid'  'mdesc
'draw string 4.25  9.95 'math'  'title' 'season' ('nmod')  (blue)'
'set string 1 c 6'
'draw string 4.25  9.70 vs'
'draw string 4.25  9.45 'odesc'  'season' ('nobs')  ('climate')  (black)'

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

'set string 4 l 4'
'set strsiz .08'
'draw string 2.50  3.37 Beg: 'bmnthm' 'byearm
'draw string 2.50  3.24 End: 'emnthm' 'eyearm
'set string 1 l 4'
'draw string 2.50  3.50 Mod Dates:'
'draw string 5.00  3.50 Obs Dates:'
'draw string 5.00  3.37 Beg: 'bmntho' 'byearo
'draw string 5.00  3.24 End: 'emntho' 'eyearo
'set string 1 c 6'
* --------------------------------

if( math = LOG )
    mathparm = 'LOG_'
else
    mathparm = ''
endif
'myprint -name 'output'/'mname'_z_'mathparm''PFX''oname'.'season

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
