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
****        Mean Square Error  =   Variance + BIAS**2                       ****
****                      MSE  =   VAR + MES                                ****
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

'run getinfo time'
             time  = result
'run getinfo tinc'
             tinc = result

if( tinc = "NULL" ) ; tinc = 6 ; endif
'run getinfo time'
             time = result
                      hour = (time-1)*tinc
if( hour < 10  )    ; hour = 0hour ; endif
if( hour < 100 )    ; hour = 0hour ; endif

'run getenv "SYSCMP_TDIM"'
                    tdim  = result
say 'TDIM: 'tdim

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

'getinfo zmin'
         zmin = result
'getinfo zmax'
         zmax = result

'set z 'zmin
'getinfo level'
         levmin = result
'set z 'zmax
'getinfo level'
         levmax = result

 if( levmax < levmin )
     levmin = levmax
 endif
'set z 'zmin' 'zmax

 if( levmin < 100 & field != q )
    'set zlog on'
    'setlevs'
 else
    'set zlog off'
 endif


*******************************************************************************
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
*******************************************************************************

    say 'TYPE: 'type
    say 'exp1: 'exp1' tag1: 'tag1
    say 'exp2: 'exp2' tag2: 'tag2

    tipe = type
if( type = Dres ) 
    type = Dmse
endif

* ----------------------------------------------------------------------------------

'set t 'tdim

'define delDmse = 'field'Dmse'tag2'z'
'define delDmes = 'field'Dmes'tag2'z'
'define delDamp = 'field'Damp'tag2'z'
'define delDphz = 'field'Dphz'tag2'z'

* To Unify Contour Interval and Scaling, use: delDmse , otherwise use: del'type'
* ------------------------------------------------------------------------------
*      dummy  = getstuff( 'delDmse' )
       dummy  = getstuff(  del''type'' )
   diff_cint  = subwrd(dummy,1)
   diff_scale = subwrd(dummy,2)
        diffm = subwrd(dummy,3)

if( tipe = Dres )
       'define  delDmse = delDmse - ( delDmes + delDamp + delDphz)'

        dummy2  = getstuff( 'delDmse' )
        diff_cint2  = subwrd(dummy2,1)
        diff_scale2 = subwrd(dummy2,2)
        diffm2 = subwrd(dummy2,3)
        if( diffm2 != diffm )
            diff_cint  = diff_cint2
            diff_scale = diff_scale2
            diffm = diffm2
        endif
endif

'set t 'time

'define delDmse = 'field'Dmse'tag2'z'
'define delDmes = 'field'Dmes'tag2'z'
'define delDamp = 'field'Damp'tag2'z'
'define delDphz = 'field'Dphz'tag2'z'

if( tipe = Dres )
       'define  delDmse = delDmse - ( delDmes + delDamp + delDphz)'
endif

'define diff  = del'type

* Compute Confidence Interval
* ---------------------------
    numfm1 = numf - 1

if( type = Dmes ) ; 'define tval = sqrt('numfm1') * 'field'Dmes'tag2'z / sqrt('field'DDmse'tag2'z)' ; endif
if( type = Damp ) ; 'define tval = sqrt('numfm1') * 'field'Damp'tag2'z / sqrt('field'DDmse'tag2'z)' ; endif
if( type = Dphz ) ; 'define tval = sqrt('numfm1') * 'field'Dphz'tag2'z / sqrt('field'DDmse'tag2'z)' ; endif
if( type = Dmse ) ; 'define tval = sqrt('numfm1') * 'field'Dmse'tag2'z / sqrt('field'DDmse'tag2'z)' ; endif

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
'define  diff1 =          maskout( diff*'diff_scale',abs(diff*'diff_scale')-'diff_cint' )'

'd diff*'diff_scale

* Create New File containing Ratio: DIFF/DIFF0 for contouring region of significance
* ----------------------------------------------------------------------------------
'getinfo file'
      curfile = result

if( tipe != Dres )
'define diffr = diff1/diff0'
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
'set csmooth off'
'set clab  off'
'set ccolor 1'
'set cthick 1'
'set clevs 0.999999'
'd diffr'
'close 'diffile

'set dfile 'curfile
'set undef 'undef
endif

* -------------------------------------------------------------


'run getenv MONTHLAB'
            month = result
say 'MONTH_LABEL: 'month

'getinfo year'
         year  = result

'run uppercase 'field
                UFIELD = result

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
'set string 1 c 6'

'set strsiz 0.16'
'draw string 5.50 8.40 Zonal Mean 'UFIELD'   'month' 'year'   Forecast Hour: 'hour
'set strsiz 0.12'
'draw string 5.50 8.18 'exp2' - 'exp1'  'numf'-member Ensemble   (Contour > 'confidence'% Confidence)'
'set strsiz .10'
'draw string 5.50 7.95 'title

endif

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
