function stdiff (args)

********************************************************************************
****                                                                        ****
**** Note:  Forecast           =>  F                                        ****
****        Analysis           =>  A                                        ****
****                                                                        ****
****        Mean Square Error  =>  MSE  =   1/N * SUM[ (F-A)**2 ]           ****
****        Mean Error         =>  BIAS =   1/N * SUM[ (F-A) ]              ****
****        Mean Error Squared =>  MES  =   BIAS**2                         ****
****        Root Mean  Square  =>  RMS  = SQRT[ MSE ]                       ****
****        Variance           =>  VAR  =   1/N * SUM[ (F-A-BIAS)**2 ]      ****
****        Standard Deviation =>  STD  = SQRT[ VAR ]                       ****
****                                                                        ****
****        F Mean             =>  FBAR =   1/N * SUM[  F  ]                ****
****        A Mean             =>  ABAR =   1/N * SUM[  A  ]                ****
****        F Variance         =>  FVAR =   1/N * SUM[ (F-FBAR)**2 ]        ****
****        A Variance         =>  AVAR =   1/N * SUM[ (A-ABAR)**2 ]        ****
****        CoVariance         =>  COV  =   1/N * SUM[ (F-FBAR)*(A-ABAR) ]  ****
****                                                                        ****
****        BIAS      Error    =>  BIA  = [ FBAR - ABAR ]**2                ****
****        Amplitude Error    =>  AMP  = [ FSTD - ASTD ]**2                ****
****        Phase     Error    =>  PHZ  = 2*[ FSTD*ASTD - COV ]             ****
****                                                                        ****
****        Mean Square Error  =   BIAS- + Amplitude- + Phase-Error         ****
****                                                                        ****
****                      MSE  =   BIA + AMP + PHZ                          ****
****                                                                        ****
********************************************************************************

                           k = 1
field = subwrd  (args,k) ; k = k + 1
type  = subwrd  (args,k) ; k = k + 1
exp1  = subwrd  (args,k) ; k = k + 1
exp2  = subwrd  (args,k) ; k = k + 1
numf  = subwrd  (args,k) ; k = k + 1
title = subwrd  (args,k) ; k = k + 1

xloc  = subwrd  (args,k) ; k = k + 1
yloc  = subwrd  (args,k) ; k = k + 1
xmax  = subwrd  (args,k) ; k = k + 1
ymax  = subwrd  (args,k) ; k = k + 1
 
'fixname 'exp1
          tag1 = result
'fixname 'exp2
          tag2 = result
          tag2 = EXP

say 'TYPE: 'type
say 'exp1: 'exp1' tag1: 'tag1
say 'exp2: 'exp2' tag2: 'tag2

********************************************************************************
* if( type = std ) ; string = "Std.Dev.  (Forecast-Analysis)"           ; endif
* if( type = var ) ; string = "Variance  (Forecast-Analysis)"           ; endif
* if( type = rms ) ; string = "RMS  (Forecast-Analysis)"                ; endif
* if( type = fma ) ; string = "MEAN  (Forecast-Analysis)"               ; endif
* if( type = fmc ) ; string = "MEAN  (Forecast-Climatology)"            ; endif
* if( type = mse ) ; string = "Mean_Square_Error  (Forecast-Analysis)"  ; endif
* if( type = mes ) ; string = "Mean_Error_Squared  (Forecast-Analysis)" ; endif
* if( type = rmes) ; string = "Root_Mean_Error_Squared  (Forecast-Analysis)" ; endif
********************************************************************************

'run getinfo level'
             level = result
'run getinfo time'
             time  = result
'run getinfo tinc'
             tinc = result

if( tinc = "NULL" ) ; tinc = 6 ; endif
                      hour = (time-1)*tinc
if( hour < 10  )    ; hour = 0hour ; endif
if( hour < 100 )    ; hour = 0hour ; endif

'run getenv "SYSCMP_TDIM"'
                    tdim  = result

'parea 'xloc' 'yloc' 'xmax' 'ymax
 xmid = subwrd(result,1)
 ybot = subwrd(result,2)
 ytop = subwrd(result,3)

if( xloc = 1 & yloc = 1 )
    x1 = 1.014 + ( 4.12 / 4 )
    x2 = 5.135 + ( 4.12 / 4 )
    y1 = 4.58
    y2 = 7.62
   'set parea 'x1' 'x2' 'y1' 'y2
    xmid = ( x1 + x2 )/2
endif

'set datawarn off'
'set grid  off'
'set grads off'
'set gxout shaded'
'set ccols 92 0'
'set clevs 0.5'
'd mask'

********************************************************************************
*         TVAL = SQRT(N) x NAME_DIFF_MEAN / NAME_DIFF_STD     *
*                                                             *
*       astudt (N-1) 0.10  [90% Confidence]                   *
*       astudt (N-1) 0.05  [95% Confidence]                   *
*       astudt (N-1) 0.04  [96% Confidence]                   *
*       astudt (N-1) 0.02  [98% Confidence]                   *
*       astudt (N-1) 0.01  [99% Confidence]                   *
*                                                             *
*       q defval astudtout 1 1                                *
*       critval=subwrd(result,3)                              *
*                                                             *
*       and then CONTOUR TVAL using CLEVS = critval          
********************************************************************************

    tipe = type
if( type = Dres ) 
    type = Dmse
endif

'set gxout shaded'

* ---------------------------------------------------------------

'set t 'tdim
'define delDmse = 'field'Dmse'tag2
'define delDmes = 'field'Dmes'tag2
'define delDamp = 'field'Damp'tag2
'define delDphz = 'field'Dphz'tag2

if( tipe = Dres )
   'define  delDmse = delDmse - ( delDmes + delDamp + delDphz)'
   'define dumm = regrid2(  delDmse,.25, .25, bs_p1, 0, -90 )'
else
*   To Unify Contour Interval and Scaling, use: delDmse , otherwise use: del'type'
*   ------------------------------------------------------------------------------
   'define dumm = regrid2(  delDmse,  .25, .25, bs_p1, 0, -90 )'
*  'define dumm = regrid2(  del'type',.25, .25, bs_p1, 0, -90 )'
endif

       dummy  = getstuff( 'dumm' )
   diff_cint  = subwrd(dummy,1)
   diff_scale = subwrd(dummy,2)
        diffm = subwrd(dummy,3)

* ---------------------------------------------------------------

'set t 'time
'define delDmse = 'field'Dmse'tag2
'define delDmes = 'field'Dmes'tag2
'define delDamp = 'field'Damp'tag2
'define delDphz = 'field'Dphz'tag2

if( tipe = Dres )
   'define  delDmse = delDmse - ( delDmes + delDamp + delDphz)'
endif

'define diff = regrid2(  del'type',.25, .25, bs_p1, 0, -90 )'

* Compute Confidence Interval
* ---------------------------
    numfm1 = numf - 1

if( type = Dmse ) ; 'define tval = sqrt('numfm1') * 'field'Dmse'tag2' / sqrt('field'DDmse'tag2')' ; endif
if( type = Dmes ) ; 'define tval = sqrt('numfm1') * 'field'Dmes'tag2' / sqrt('field'DDmse'tag2')' ; endif
if( type = Damp ) ; 'define tval = sqrt('numfm1') * 'field'Damp'tag2' / sqrt('field'DDmse'tag2')' ; endif
if( type = Dphz ) ; 'define tval = sqrt('numfm1') * 'field'Dphz'tag2' / sqrt('field'DDmse'tag2')' ; endif

   'define tval = regrid2( tval,0.25,0.25,bs_p1,0,-90 )'

         ttest = 0.10
    confidence = 100 * (1-ttest)
   'astudt 'numfm1' 'ttest
   'q defval astudtout 1 1'
      critval=subwrd(result,3)
      say 'critval: 'critval

if( tipe = Dres )
   'shades 'diff_cint
else
   'shades 'diff_cint' -quad'
            diff_cint = result
endif

'define  diff0 = maskout( maskout( diff*'diff_scale',abs(diff*'diff_scale')-'diff_cint' ), abs(tval)-'critval')'
'define  diff  =          maskout( diff*'diff_scale',abs(diff*'diff_scale')-'diff_cint' )'

'd diff'

* Create New File containing Ratio: DIFF/DIFF0 for contouring region of significance
* ----------------------------------------------------------------------------------
'getinfo file'
      curfile = result

if( tipe != Dres )
'define diffr = diff/diff0'
'getinfo undef'
         undef = result
'set undef 0.0'
'set sdfwrite -5d diffr.nc4'
'sdfwrite diffr'
'undefine diffr'

'sdfopen  diffr.nc4'
'getinfo  numfiles'
          diffile = result
'set dfile 'diffile
'q ctlinfo'
   ctlinfo = result

          fname = 'diffr.ctl'
'!remove 'fname
             n = 1
             line = sublin(ctlinfo,n)
             write(fname,line)
      while( line != 'endvars' )
             n = n + 1
             line = sublin(ctlinfo,n)
             word = subwrd(line,1)
         if( word = 'undef' )
             line = 'undef 1e15'
         endif
             write(fname,line,append)
      endwhile
      close = close(fname)

'close 'diffile
'open  'fname
'getinfo numfiles'
         diffile = result
'set dfile 'diffile

'set gxout contour'
'set clab  off'
'set ccolor 1'
'set cthick 1'
'set clevs 0.5'
'd diffr'
'close 'diffile
endif

* Create New File containing DIFF with zeroes rather than UNDEF for Global Mean Metrics
* -------------------------------------------------------------------------------------
'define diff2 = diff0/'diff_scale
'set sdfwrite -5d diff'type'.nc4'
'set undef 0.0'
'sdfwrite diff2'
'sdfopen  diff'type'.nc4'
'getinfo  numfiles'
          diffile = result
'set dfile 'diffile
'q ctlinfo'
   ctlinfo = result
   
          fname = 'diff'type'.ctl'
'!remove 'fname
             n = 1
             line = sublin(ctlinfo,n)
             write(fname,line)
      while( line != 'endvars' )
             n = n + 1
             line = sublin(ctlinfo,n)
             word = subwrd(line,1)
         if( word = 'undef' )
             line = 'undef 1e15'
         endif
             write(fname,line,append)
      endwhile
      close = close(fname)

'close 'diffile
'open  'fname
'getinfo numfiles'
         diffile = result
'set dfile 'diffile
'define NHEM  = aave( diff2.'diffile',lon=0,lon=360,lat= 20,lat= 80 )'
'define SHEM  = aave( diff2.'diffile',lon=0,lon=360,lat=-80,lat=-20 )'
'define GLOB  = aave( diff2.'diffile',lon=0,lon=360,lat=-90,lat= 90 )'
'define TROP  = aave( diff2.'diffile',lon=0,lon=360,lat=-20,lat= 20 )'
'close 'diffile

'set dfile 'curfile
'set undef 'undef
* -------------------------------------------------------------


'run getenv MONTHLAB'
            month = result
say 'MONTH_LABEL: 'month

'getinfo year'
         year  = result

'set string 1 c 5'

if( tipe = Dmse )
   'set strsiz 0.10'
   'draw string 'xmid' 'ytop' MSE: Mean Square Error  (x10**'diffm')' ; 
   'cbarn -sbar  0.5   -snum 0.35 -xmid 'xmid' -ymid 'ybot' -scaley 0.5 -scalex 1.0'
endif

if( tipe = Dres )
   'set strsiz 0.08'
   'draw string 'xmid' 'ytop' MSE - [ BIAS + AMPLITUDE + PHASE] (x10**'diffm')' ; 
   'cbarn -sbar  0.4   -snum 0.35 -xmid 'xmid' -ymid 'ybot' -scaley 0.5 -scalex 0.9'
endif
if( type = Dmes )
   'set strsiz 0.08'
   'draw string 'xmid' 'ytop' BIAS Error Difference  (x10**'diffm')' ; 
   'cbarn -sbar  0.4   -snum 0.35 -xmid 'xmid' -ymid 'ybot' -scaley 0.5 -scalex 0.9'
endif
if( type = Damp )
   'set strsiz 0.08'
   'draw string 'xmid' 'ytop' AMPL Error Difference  (x10**'diffm')' ; 
   'cbarn -sbar  0.4   -snum 0.35 -xmid 'xmid' -ymid 'ybot' -scaley 0.5 -scalex 0.9'
endif
if( type = Dphz )
   'set strsiz 0.08'
   'draw string 'xmid' 'ytop' PHASE Error Difference  (x10**'diffm')' ; 
   'cbarn -sbar  0.4   -snum 0.35 -xmid 'xmid' -ymid 'ybot' -scaley 0.5 -scalex 0.9'
endif


if( tipe = Dres )
'set vpage off'
'set parea off'
'set string 1 c 6'

'set strsiz 0.16'
'run uppercase 'field
                UFIELD = result
'draw string 5.50 8.40 'level'-mb 'UFIELD'   'month' 'year'   Forecast Hour: 'hour
'set strsiz 0.12'
'draw string 5.50 8.18 'exp2' - 'exp1'  'numf'-member Ensemble   (Contour > 'confidence'% Confidence)'
'set strsiz .10'
'draw string 5.50 7.95 'title

endif

return

function getstuff( q )

'minmax.simple 'q
   qmax = subwrd(result,1)
   qmin = subwrd(result,2)
   cint = (qmax-qmin)/18

'stats 'q
  cint = subwrd(result,2)

say 'Inside getstuff for 'q', cint: 'cint
if( cint = 0 )
    fact = 1
   icint = 0
   scale = 0
    cmax = 0
    cmin = 0
else

'd log10('cint')'
   log10cint = subwrd(result,4)
'getint 'log10cint
         scale = result

if( scale <= 0 )
   'd pow(10,abs('scale'))'
    fact = subwrd(result,4)
else
   'd pow(10,-'scale')'
    fact = subwrd(result,4)
endif

     'getint 'cint*fact
      icint = result

say ' scale: 'scale'  fact: 'fact'  icint: 'icint
while( icint < 1  )
   if( scale <= 0 )
       fact = fact*10
      scale = scale - 1
   else
       fact = fact/10
      scale = scale + 1
   endif

   'getint 'cint*fact
     icint = result
say ' scale: 'scale'  fact: 'fact'  icint: 'icint
endwhile

'minmax.simple 'q
   qmax = subwrd(result,1)
   qmin = subwrd(result,2)

'getint 'qmax*fact/icint
        dqmax = result
'getint 'qmin*fact/icint
        dqmin = result
cmax = icint*dqmax
cmin = icint*dqmin

say 'qmax: 'qmax'  cmax: 'cmax
say 'qmin: 'qmin'  cmin: 'cmin

endif
return icint' 'fact' 'scale' 'cmax' 'cmin
