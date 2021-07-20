function closeness (args)

* Initialize INPUT Parameters
* ---------------------------
    cvar = getarg (args,CVAR)
    mvar = getarg (args,MVAR)
    ovar = getarg (args,OVAR)

   cname = getarg (args,CNAME)
   mname = getarg (args,MNAME)
   oname = getarg (args,ONAME)

   cdesc = getarg (args,CDESC)
   mdesc = getarg (args,MDESC)
   odesc = getarg (args,ODESC)

   mfile = getarg (args,MFILE)
   bdate = getarg (args,MBEGDATE)
   edate = getarg (args,MENDDATE)

   ofile = getarg (args,OFILE)
  bdateo = getarg (args,OBEGDATE)
  edateo = getarg (args,OENDDATE)

   expid = getarg (args,EXPID)
gridcomp = getarg (args,GC)
  prefix = getarg (args,PREFIX)
  season = getarg (args,SEASON)
  output = getarg (args,OUTPUT)
 climate = getarg (args,CLIMATE)
    math = getarg (args,MATH)
   level = getarg (args,LEVEL)

* ---------------------------
if(     math = NULL ) ;     math = '' ; endif
if(   season = NULL ) ;   season = '' ; endif
if( gridcomp = NULL ) ; gridcomp = '' ; endif

'set t 1'
'run getenv "GEOSUTIL"'
             geosutil = result

'run getenv "LEVTYPE"'
             LEVTYPE = result
say 'GETENV LEVTYPE = 'LEVTYPE
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

                        'getresource 'PLOTRC' 'PFX'FIXED_PLOT_FACTOR'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      FIXED_PLOT_FACTOR' ; endif
                                                   fixpltfact = result

                        'getresource 'PLOTRC' 'PFX'CLOSE_PLOT_FACTOR'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CLOSE_PLOT_FACTOR' ; endif
                                                   clspltfact = result

                        'getresource 'PLOTRC' 'PFX'FIXED_PLOT_CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      FIXED_PLOT_CINT' ; endif
                                                   fixpltcint = result

                        'getresource 'PLOTRC' 'PFX'CLOSE_PLOT_CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC'      CLOSE_PLOT_CINT' ; endif
                                                   clspltcint = result

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

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_FIXED_PLOT_FACTOR'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_FIXED_PLOT_FACTOR' ; endif
                                                                 fixpltfact = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CLOSE_PLOT_FACTOR'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CLOSE_PLOT_FACTOR' ; endif
                                                                 clspltfact = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_TITLE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_TITLE' ; endif
                                                                title = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CINT' ; endif
                                                                ccint = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_FIXED_PLOT_CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_FIXED_PLOT_CINT' ; endif
                                                                 fixpltcint = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CLOSE_PLOT_CINT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CLOSE_PLOT_CINT' ; endif
                                                                 clspltcint = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CCOLS' ; endif
                                                                ccols = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_CLEVS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_CLEVS' ; endif
                                                                clevs = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_REGRID'
                                                                method = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_DCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_DCOLS' ; endif
                                                                 dcols = result

                        'getresource 'PLOTRC' 'mname'_'gridcomp'_'level'_'LEVTYPE
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_'level' 'DLEVS ; endif
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'gridcomp'_'DLEVS         ; endif
                                                                  dlevs = result
endif

'run getenv "CINTDIFF"'
             CINTDIFF  = result
         if( CINTDIFF != 'NULL' ) ; dcols = 'NULL' ; dlevs = 'NULL' ; endif

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

if( dcols = 'NULL' )
     dcols = '55  49  47  45  44  36  34  33  32  0  21  22  23  24  25  26  27  28 69'
endif

* Perform Mathematics if necessary
* --------------------------------
'define qmod = 'mvar''season'*'factor
'define cmod = 'cvar''season'*'factor
'define qobs = 'ovar''season'*'factor

if( math = LOG )
   'define qmod = log(qmod+0.00001)'
   'define cmod = log(cmod+0.00001)'
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

'set vpage 0 8.5 0.0 11'
'set grads off'
'rgbset'

'getinfo lon'
         lon = result
'define obsg = qobs'
'define modg = qmod'
'define codg = cmod'
'define difg = maskout( modg-obsg,abs(obsg) )'
'define cifg = maskout( codg-obsg,abs(obsg) )'

* Top Panel
* ---------
'set parea 1.5 7.0 7.70 10.50'

   'set gxout stat'
   'd modg'
   qmodminmax = sublin(result,8)
   qmodmin    = subwrd(qmodminmax,4)
   qmodmax    = subwrd(qmodminmax,5)
   say 'QMOD_Max Value: 'qmodmax
   say 'QMOD_Min Value: 'qmodmin
   'set gxout shaded'
   'd abs('qmodmin')'
          aqmodmin = subwrd(result,4)
   'd abs('qmodmax')'
          aqmodmax = subwrd(result,4)
   if( aqmodmin > aqmodmax ) ; aqmodmax = aqmodmin ; endif
   say 'Absolute QMOD_MAX: '   aqmodmax

ntop = 0
    say 'TOP CLEVS: 'dlevs
    say 'TOP CCOLS: 'dcols
if( dcols != NULL & dlevs != NULL )
   'set clevs 'dlevs
   'set ccols 'dcols
   'd difg'
else
   'stats difg'
     avgdif = subwrd(result,1)
     stddif = subwrd(result,2)

say 'avgdif = 'avgdif
say 'stddif = 'stddif
     dqmax  =  stddif/3

     dqrel = dqmax / aqmodmax  * 100

say 'dqrel = 'dqrel
    'getint 'dqrel*100
             dqrel = result/100
say 'dqrel = 'dqrel

     dpct  = 0.1
     say 'Absolute DQMAX: 'dqmax'  QMOD_MAX: 'aqmodmax
     say 'Relative Percent Difference: 'dqrel' (100*DQMAX/QMOD_MAX)'
     say ' Default Percent Difference for  Plots: 'dpct

     if( dqrel < dpct )
         dqrel = dpct
     endif
         dqmax = dqrel * aqmodmax / 100
         say 'Setting Diff CINT using Relative Percent Difference: 'dqrel'%  dqmax = 'dqmax

   if( dqmax > 0 )
      'd log10('dqmax')'
       ntop = subwrd(result,4)
   else
       ntop = 0
   endif
   say '    Log Factor: 'ntop
   if( ntop<0 ) ; ntop = ntop-2 ; endif
   'getint 'ntop
            ntop = result
   if( ntop>0 )
       if( ntop<=2 )
           ntop = 0
        else
           ntop = ntop+2
        endif
   endif

   if( fixpltfact != NULL )
    'd 'fixpltfact
     ntop = subwrd(result,4)
   endif
   say 'Diff Scaling Factor: 'ntop

   if( fixpltcint != NULL )
    'd 'fixpltcint
        fixpltcint = subwrd(result,4)
        cint = fixpltcint
   else
      if( ntop < 0 )
          ztop = -1 * ntop
         'd 'dqmax'*1e'ztop
      else
      '   d 'dqmax'/1e'ntop
      endif
        cint = subwrd(result,4)
   endif

       say 'TOP CINT: 'cint
      'shades 'cint
      if( ntop < 0 )
          ztop = -1 * ntop
         'define difg = difg*1e'ztop
      else
         'define difg = difg/1e'ntop
      endif
      'd difg'
endif

   'set gxout stat'
   'd difg'
   mdifminmax = sublin(result,8)
   mdifmin    = subwrd(mdifminmax,4)
   mdifmax    = subwrd(mdifminmax,5)
   say 'MDIF_Max Value: 'mdifmax
   say 'MDIF_Min Value: 'mdifmin
   'set gxout shaded'
   'd abs('mdifmin')'
          amdifmin = subwrd(result,4)
   'd abs('mdifmax')'
          amdifmax = subwrd(result,4)
   if( amdifmin > amdifmax ) ; amdifmax = amdifmin ; endif
   say 'Absolute MDIF_MAX: '   amdifmax

'set parea 0 8.5 7.0 11'
'cbarn -vert'

* Middle Panel
* ------------
'set parea 1.5 7.0 4.30 7.10'

nmid = 0
if( dcols != NULL & dlevs != NULL )
   'set clevs 'dlevs
   'set ccols 'dcols
   'd cifg'
else
       nmid = ntop
      'shades 'cint
      if( ntop < 0 )
          ztop = -1 * ntop
         'define cifg = cifg*1e'ztop
      else
         'define cifg = cifg/1e'ntop
      endif
      'd cifg'
endif

   'set gxout stat'
   'd cifg'
   qobsminmax = sublin(result,8)
   qobsmin    = subwrd(qobsminmax,4)
   qobsmax    = subwrd(qobsminmax,5)
   say 'QOBS_Max Value: 'qobsmax
   say 'QOBS_Min Value: 'qobsmin
   'set gxout shaded'
   'd abs('qobsmin')'
          aqobsmin = subwrd(result,4)
   'd abs('qobsmax')'
          aqobsmax = subwrd(result,4)
   if( aqobsmin > aqobsmax ) ; aqobsmax = aqobsmin ; endif
   say 'Absolute QOBS_MAX: '   aqobsmax

* Bottom Panel
* ------------
'set parea 1.5 7.0 0.90 3.70'

'define closeness = abs(difg)-abs(cifg)'

   'set gxout stat'
   'd closeness'
   closminmax = sublin(result,8)
   closmin    = subwrd(closminmax,4)
   closmax    = subwrd(closminmax,5)
   say 'CLOSE_Max Value: 'closmax
   say 'CLOSE_Min Value: 'closmin
   'set gxout shaded'
   'd abs('closmin')'
          aclosmin = subwrd(result,4)
   'd abs('closmax')'
          aclosmax = subwrd(result,4)
   if( aclosmin > aclosmax ) ; aclosmax = aclosmin ; endif
   say 'Absolute CLOSE_MAX: '  aclosmax

   'stats closeness'
     avgdif = subwrd(result,1)
     stddif = subwrd(result,2)
      dqmax = stddif/3

     dqrel = dqmax / amdifmax  * 100
say 'raw dqrel = 'dqrel
    'getint 'dqrel*100
             dqrel = result/100
say 'int dqrel = 'dqrel

     dpct  = 0.1
     say 'DQMAX: 'dqmax'  QMOD_MAX: 'amdifmax
     say 'Relative Percent Difference: 'dqrel' (100*DQMAX/QMOD_MAX)'
     say ' Default Percent Difference for  Plots: 'dpct

     if( dqrel < dpct )
         dqrel = dpct
     endif
         dqmax = dqrel * amdifmax / 100
         say 'Setting Diff CINT using Relative Percent Difference: 'dqrel'%  dqmax = 'dqmax

   if( dqmax > 0 )
      'd log10('dqmax')'
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

   if( clspltfact != NULL )
    'd 'clspltfact
     n = subwrd(result,4)
   endif
   say 'Diff Scaling Factor: 'n

   if( clspltcint != NULL )
    'd 'clspltcint
        clspltcint = subwrd(result,4)
        cint = clspltcint
   else
      if( n < 0 )
          ztop = -1 * n
         'd 'dqmax'*1e'ztop
      else
         'd 'dqmax'/1e'n
      endif
        cint = subwrd(result,4)
   endif

      'shades 'cint

      if( n < 0 )
          ztop = -1 * n
         'define closeness = closeness*1e'ztop
      else
         'define closeness = closeness/1e'n
      endif
      'd closeness'

'cbarn -snum 0.55 -xmid 4.25 -ymid 0.4'

'stats difg'
 avgmod = subwrd(result,1)
 stdmod = subwrd(result,2)
'stats cifg'
 avgobs = subwrd(result,1)
 stdobs = subwrd(result,2)
'stats closeness'
 avgdif = subwrd(result,1)
 stddif = subwrd(result,2)

'set vpage off'
'set string 1 l 4'
'set strsiz .065'
'draw string 0.05 0.08 ( EXPID:  'expid' )'

'set string 1 c 6'
'set strsiz .125'
if( level = '' | level = 0 )
   'draw string 4.25 10.85 'math' 'title
else
   'draw string 4.25 10.85 'level'-mb 'math' 'title
endif
'set strsiz .10'

if( ntop != 0 )
   if( ntop>0 )
      'draw string 4.25 10.62 'mdesc' - 'oname'  'season' ('nmod')  (x 10** -'ntop')'
   else
      'draw string 4.25 10.62 'mdesc' - 'oname'  'season' ('nmod')  (x 10**'ntop')'
   endif
else
   'draw string 4.25 10.62 'mdesc' - 'oname'  'season' ('nmod')'
endif

if( nmid != 0 )
   if( nmid>0 )
      'draw string 4.25 7.22 'cdesc' - 'oname'  'season' ('nobs')  (x 10** -'nmid')  ('climate')'
   else
      'draw string 4.25 7.22 'cdesc' - 'oname'  'season' ('nobs')  (x 10**'nmid')  ('climate')'
   endif
else
   'draw string 4.25 7.22 'cdesc' - 'oname'  'season' ('nobs')  ('climate')'
endif

*if( n != 0 )
if( ntop != 0 )
   if( ntop>0 )
      'draw string 4.25 3.80 Closeness to 'oname':  ABS(Top)-ABS(Middle)  (x 10**-'ntop')'
   else
      'draw string 4.25 3.80 Closeness to 'oname':  ABS(Top)-ABS(Middle)  (x 10**'ntop')'
   endif
else
   'draw string 4.25 3.80 Closeness to 'oname':  ABS(Top)-ABS(Middle)'
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
'draw string 0.050 10.30 Beg: 'bmnthm' 'byearm
'draw string 0.050 10.15 End: 'emnthm' 'eyearm
'draw string 0.050 9.85  Max: 'mdifmax
'draw string 0.050 9.70  Min: 'mdifmin
'draw string 0.050 9.40 Mean: 'avgmod
'draw string 0.050 9.25  Std: 'stdmod

'draw string 0.050 6.90 Beg: 'bmntho' 'byearo
'draw string 0.050 6.75 End: 'emntho' 'eyearo
'draw string 0.050 6.45  Max: 'qobsmax
'draw string 0.050 6.30  Min: 'qobsmin
'draw string 0.050 6.00 Mean: 'avgobs
'draw string 0.050 5.85  Std: 'stdobs

'draw string 0.050 3.50 Beg: 'bmnthm' 'byearm
'draw string 0.050 3.35 End: 'emnthm' 'eyearm
'draw string 0.050 3.05  Max: 'closmax
'draw string 0.050 2.90  Min: 'closmin
'draw string 0.050 2.60 Mean: 'avgdif
'draw string 0.050 2.45  Std: 'stddif

*if( CINTDIFF != 'NULL' )
   'set strsiz .07'
   'draw string 0.050 1.77 Plot represents'
   'draw string 0.050 1.62 values > 'dqrel' %'
   'draw string 0.050 1.47 Relative Difference'
   'draw string 0.050 1.32 ( DQ/QMax )'
*endif


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

function getarg (args,name)
'numargs  'args
 numargs = result
        num = 0
while ( num < numargs )
        num = num + 1
if( subwrd(args,num) = '-'name ) 
    arg = subwrd(args,num+1)
    say name' = 'arg
    return arg
endif
endwhile
return

