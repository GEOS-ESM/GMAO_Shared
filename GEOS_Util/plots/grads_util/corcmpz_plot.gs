function corcmpz (args)

'numargs  'args
 numargs = result

'run getenv SOURCE'
        SOURCE = result

field = h
desc  = ''
PLOT  = TRUE

       num = 0
while( num < numargs )
       num = num + 1
if( subwrd(args,num)='-field'  ) ; field  = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-numexp' ) ; numexp = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-desc'   ) ; desc   = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-debug'  ) ; debug  = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-NOPLOT' ) ; PLOT   = FALSE              ; endif
endwhile
                                   mexps  = numexp-1
       num = 0
while( num < numargs )
       num = num + 1
       m = 0
       while( m<=mexps )
       if( subwrd(args,num)='-desc'm  ) ; expdsc.m = subwrd(args,num+1) ; endif
       m = m + 1
       endwhile
endwhile

***********************************************************
*****                                                 *****
*****  Note:  numexp is the Total Number of           *****
*****         Experiments (including the Control)     *****
*****         being compared.                         *****
*****                                                 *****
***********************************************************

'getinfo numfiles'
files_per_month = result/numexp

expcmp.0 = 1
expcol.0 = 1
       m = 1
while( m <= mexps )
expcmp.m = m+1
expcol.m = m+1
       m = m+1
endwhile


say '   Exp0: 'expcmp.0'  Color: 'expcol.0
       m = 1
while( m <= mexps )
say '   Exp'm': 'expcmp.m'  Color: 'expcol.m
       m = m+1
endwhile

say 'MEXPs = 'mexps
say 'FPM: 'files_per_month

'rgbset'
'run setenv "ANTIALIAS" NULL'


* Initialize Plot Values
* ----------------------
if( field = "h" )
    name  = "Heights"
    unit  = "(m)"
    label = "hght"
    cint  = 0.02
  rbrange = "0.70 1.0"
endif
if( field = "u" )
    name  = "U-Wind"
    unit  = "(m/sec)"
    label = "uwnd"
    cint  = 0.02
  rbrange = "0.70 1.0"
endif
if( field = "v" )
    name  = "V-Wind"
    unit  = "(m/sec)"
    label = "vwnd"
    cint  = 0.02
  rbrange = "0.70 1.0"
endif
if( field = "t" )
    name  = "Temperature"
    unit  = "(K)"
    label = "tmpu"
    cint  = 0.02
  rbrange = "0.70 1.0"
endif
if( field = "q" )
    name  = "Specific Humidity"
    unit  = "(g/g)"
    label = "sphu"
    cint  = 0.02
  rbrange = "0.70 1.0"
endif


'getinfo xpos'
         xpos  = result
'getinfo ypos'
         ypos  = result

if( xpos =  1 ) ; region = "Global"                                     ;  reg = "GLO"  ; endif
if( xpos =  2 ) ; region = "N.Hem. ExtraTropics (Lats: 20,80)"          ;  reg = "NHE"  ; endif
if( xpos =  3 ) ; region = "Tropics (Lats: -20,20)"                     ;  reg = "TRO"  ; endif
if( xpos =  4 ) ; region = "S.Hem. ExtraTropics (Lats: -20,-80)"        ;  reg = "SHE"  ; endif
if( xpos =  5 ) ; region = "N.W. Quadrant (Lons:-180,0  Lats: 0, 90)"   ;  reg = "NWQ"  ; endif
if( xpos =  6 ) ; region = "N.E. Quadrant (Lons: 0,180  Lats: 0, 90)"   ;  reg = "NEQ"  ; endif
if( xpos =  7 ) ; region = "S.W. Quadrant (Lons:-180,0  Lats: 0,-90)"   ;  reg = "SWQ"  ; endif
if( xpos =  8 ) ; region = "S.E. Quadrant (Lons: 0,180  Lats: 0,-90)"   ;  reg = "SEQ"  ; endif
if( xpos =  9 ) ; region = "North America (Lons:-140,-60  Lats: 20,60)" ;  reg = "NAM"  ; endif
if( xpos = 10 ) ; region = "Europe (Lons:-10,30  Lats: 30,60)"          ;  reg = "EUR"  ; endif
if( xpos = 11 ) ; region = "N.Polar (Lats: 60,90)"                      ;  reg = "NPO"  ; endif
if( xpos = 12 ) ; region = "S.Polar (Lats: -60,-90)"                    ;  reg = "SPO"  ; endif
if( xpos = 13 ) ; region = "X.Polar (Lats: -60,60)"                     ;  reg = "XPO"  ; endif


filebeg  = 1
fileend  = filebeg + files_per_month - 1
numfiles = fileend-filebeg+1


* Compute Beginning Times Relative to Control (1st File) across ALL Experiments
* -----------------------------------------------------------------------------
       m = 0
while( m<=mexps )
'set dfile 1'
'set t 1'
            n   = filebeg
            n.m = n + m*numfiles
'set dfile 'n.m

'q dims'
   tline    = sublin(result,5)
   toff.m.n = subwrd(tline,9)
say 'toffset:  m = 'm'  n = 'n'  toff = 'toff.m.n
m = m+1
endwhile

m = 0
while( m<=mexps )
        n   = filebeg
while ( n  <= fileend )
        n.m = n + m*numfiles
'set dfile 1'
'set t 1'
'set dfile 'n.m
'q dims'
tline       = sublin(result,5)
toffset.n.m = subwrd(tline,9)
*say 'toffset:  m: 'm'  n: 'n'  file(n.m) = 'n.m'  toffset = 'toffset.n.m
n = n + 1
endwhile
m = m + 1
endwhile

* Determine Months from Ensembles
* -------------------------------
month  = ''
months = ''
        n   = filebeg
while ( n  <= fileend )
'set dfile 'n
'set t 1'
'getinfo date'
         date = result
        dummy = substr(date,6,3)
    if( dummy != month )
        month  = dummy
        if( months = '' )
            months = month
         else
            months = months'-'month
         endif
     endif
n = n + 1
endwhile
say 'Months Used in Forecasts: 'months

* Define NDAY and NDAYMAX across ALL Experiments
* ----------------------------------------------
 ndaymax = 999
       m = 0
while( m<=mexps )
            n   = filebeg
            n.m = n + m*numfiles
'set dfile 'n.m
'set lev 1000 100'

'run getinfo tinc'
             tinc = result
         if( tinc = "NULL" ) ; tinc = 6 ; endif

'run getinfo tdim'
             tdum = result - toff.m.n
             nday = tdum * tinc / 24
         if( nday < ndaymax ) ; ndaymax = nday ; endif
m = m+1
endwhile

'run getenv "NDAY"'
             nday = result
         if( nday = "NULL" ) ; nday = ndaymax ; endif
         if( nday > ndaymax) ; nday = ndaymax ; endif
say 'NDAY: ' nday

* Determine TDIM based on TINC from Experiment Files for Diff Calculations
* ------------------------------------------------------------------------
   m = 0
while( m<=mexps )
            n   = filebeg
            n.m = n + m*numfiles
'set dfile 'n.m
'set lev 1000 100'
'getinfo tdim'
         tdim = result
'getinfo tinc'
         tinc = result

         ndayloc = (tdim-toff.m.n) * tinc / 24
         tmax    = toff.m.n + ndayloc*(24/tinc)
         tbeg.m  = toff.m.n -(tmax-tdim)
         tdim.m = tbeg.m + nday*(24/tinc)

if( tdim.m < tdim.0 )
    tdif.m = tdim.m
    ddif.m = n.m
else
    tdif.m = tdim.0
    ddif.m = 1
endif
say 'tbeg.'m': 'tbeg.m'  tdim'm': 'tdim.m'  tinc: 'tinc'  tdif.'m': 'tdif.m

m = m+1
endwhile

* Find LCD for TDIM among all Experiment Files for Synoptic Average Stats
* -----------------------------------------------------------------------
 tmin = tdim.0
dfile = 1
    m = 1
while( m<=mexps )
            n   = filebeg
            n.m = n + m*numfiles
            d.m = 1 + m*numfiles
if( tdim.m < tmin )
    tmin   = tdim.m
    dfile  = d.m
endif
m = m+1
endwhile

'set dfile 1'
'getinfo undef'
         undef = result
'    set undef ' undef

* Initialize Variables to Zero
* ----------------------------
'set dfile '1
'set lev 1000 100'
'set t  'tbeg.0' 'tdim.0
'define  zero = 0.0'

say 'TBEG: 'tbeg.0'  TEND: 'tdim.0
say 'XPOS: 'xpos  '  YPOS: 'ypos
say ' '
say 'Initialize zave and zvar ...'
   m = 0
while( m<=mexps )
            n   = filebeg
            n.m = n + m*numfiles
            d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t  'tbeg.m' 'tdim.m
'define  zave'm' = 0.0'
'define  zvar'm' = 0.0'
m = m+1
endwhile


say 'Initialize zaved and zvard ...'
   m = 0
while( m<=mexps )
            n   = filebeg
            n.m = n + m*numfiles
'set dfile 'ddif.m
'set lev 1000 100'
'set t  'tbeg.m' 'tdif.m
'define zaved'm' = 0.0'
'define zvard'm' = 0.0'
m = m+1
endwhile

* Plot each individual forecast while computing mean and Fisher Transform (to force Gaussian Distribution)
* --------------------------------------------------------------------------------------------------------

*say 'Redfine field'
m = 0
while( m<=mexps )
        n   = filebeg
while ( n  <= fileend )
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t  'tbeg.m' 'tdim.m
 toffset = 1-toffset.n.m
*say 'dfile 'd.m'  tbeg: 'tbeg.m'  tend: 'tdim.m'  toffset: 'toffset'  'field'cor'n.m' = 'field'cor.'n.m'(t+'toffset')'

'define 'field'cor'n.m' = 'field'cor.'n.m'(t+'toffset')'
n = n + 1
endwhile
m = m + 1
endwhile

* Define New Fisher Transform Variable (to force Gaussian Distribution)
* ---------------------------------------------------------------------
say ' Define New Fisher Transform Variable znem ...'
m = 0
while( m<=mexps )
        n   = filebeg
while ( n  <= fileend )
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t  'tbeg.m' 'tdim.m
'define     q'm' = 'field'cor'n.m
'define z'n'e'm' = 0.5*log( (1+ q'm')/(1- q'm'+5.0e-6) )'
*say '   DFILE: 'n.m'   Defined q'm' and z'n'e'm
n = n + 1
endwhile
m = m + 1
endwhile

say ' Define New Fisher Transform Variable zdnem ...'
m = 0
while( m<=mexps )
        n   = filebeg
while ( n  <= fileend )
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles

'set dfile 'n
'set dfile 'ddif.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
'define ctl = 'field'cor'n
*say '   DFILE: 'n'  CTL defined for TBEG = 'tbeg.0' to TEND: 'tdif.0

'set dfile 'd.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
'define dum  = 'field'cor'n.m
*say   '   DFILE: 'n.m'  DUM defined for TBEG = 'tbeg.0' to TEND: 'tdif.0

'set dfile 'ddif.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
if( m = 0 )
*say 'Computing DumDiff for EXP: 'm'   File: 'ddif.m' for TBEG = 'tbeg.0' to TEND: 'tdif.0
   'define dumdiff = dum'
else
*say 'Computing DumDiff for EXP: 'm'   File: 'n.m' for TBEG = 'tbeg.0' to TEND: 'tdif.0
   'makezdif2 -q1 dum -file1 'd.m' -q2 zero    -file2   1 -ptop 100 -name dum'
endif

'set dfile 'n
'set dfile 'ddif.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
'define     dq'm' = 0.5*(ctl-dumdiff)'
'define zd'n'e'm' = 0.5*log( (1+dq'm')/(1-dq'm') )'
n = n + 1
endwhile
m = m + 1
endwhile

* Compute Mean
* ------------
say ' Compute mean zave ...'
m = 0
while( m<=mexps )
        n  = filebeg
while ( n <= fileend )
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdim.m
'define  zave'm' =  zave'm' +  z'n'e'm
n = n + 1
endwhile
'define  zave'm' =  zave'm'/'numfiles
m = m + 1
endwhile
*say 'Hit Enter to continue ...'
*pull flag

say ' Compute mean zdave ...'
m = 0
while( m<=mexps )
        n  = filebeg
while ( n <= fileend )
'set dfile 'ddif.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
'define zaved'm' = zaved'm' + zd'n'e'm
n = n + 1
endwhile
'define zaved'm' = zaved'm'/'numfiles
m = m + 1
endwhile
*say 'Hit Enter to continue ...'
*pull flag


* Compute Variance and Standard Deviations
* ----------------------------------------
say ' Compute variance zvar and zstd ...'
m = 0
while( m<=mexps )
        n  = filebeg
while ( n <= fileend )
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdim.m
'define  zvar'm' =  zvar'm' + pow(  z'n'e'm'- zave'm',2 )'
n = n + 1
endwhile
'define  zvar'm' =  zvar'm'/('numfiles'-1)'
'define  zstd'm' = sqrt(  zvar'm' )'
m = m + 1
endwhile
*pull flag

say ' Compute variance zdvar and zdstd ...'
m = 0
while( m<=mexps )
        n  = filebeg
while ( n <= fileend )
'set dfile 'ddif.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
'define zvard'm' = zvard'm' + pow( zd'n'e'm'-zaved'm',2 )'
n = n + 1
endwhile
'define zvard'm' = zvard'm'/('numfiles'-1)'
'define zstdd'm' = sqrt( zvard'm' )'
m = m + 1
endwhile
*pull flag


* Compute Fisher Mean
* -------------------
say ' Compute Fisher Mean rave ...'
m = 0
while( m<=mexps )
        n  = filebeg
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdim.m
'define rave'm' = (exp(2*zave'm')-1) / (exp(2*zave'm')+1)'
m = m + 1
endwhile
*pull flag


* Copmute Confidence Intervals for Two-Tailed Students T-Test Distribution
* ------------------------------------------------------------------------
 dof = numfiles-1    ;* Degrees of Freedom (dof)

'astudt 'dof' 0.0001'  ;* 99.99% Confidence, and Minimum Value used for Shading
'q defval astudtout 1 1'
critval9999=subwrd(result,3)

'astudt 'dof' 0.01'  ;* 99% Confidence
'q defval astudtout 1 1'
critval99=subwrd(result,3)

'astudt 'dof' 0.05'  ;* 95% Confidence
'q defval astudtout 1 1'
critval95=subwrd(result,3)

'astudt 'dof' 0.10'  ;* 90% Confidence
'q defval astudtout 1 1'
critval90=subwrd(result,3)

'astudt 'dof' 0.32'  ;* 68% Confidence
'q defval astudtout 1 1'
critval68=subwrd(result,3)

* Estimate Statistically Significant Range for Zero-Mean Hypothesis in a Paired t-Test)
* -------------------------------------------------------------------------------------
say ' Estimate T-Test Range rUpm rLpm ...'
m = 1
while( m<=mexps )
'set dfile 'ddif.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdif.m
'define se  = sqrt( zvard'm'/'numfiles' )'

'define dx = se*'critval9999
'define rUp9999'm' = 2*(exp( 2*dx)-1)/(exp( 2*dx)+1)'
'define rLp9999'm' = 2*(exp(-2*dx)-1)/(exp(-2*dx)+1)'

'define dx = se*'critval99
'define rUp99'm' = 2*(exp( 2*dx)-1)/(exp( 2*dx)+1)'
'define rLp99'm' = 2*(exp(-2*dx)-1)/(exp(-2*dx)+1)'

'define dx = se*'critval95
'define rUp95'm' = 2*(exp( 2*dx)-1)/(exp( 2*dx)+1)'
'define rLp95'm' = 2*(exp(-2*dx)-1)/(exp(-2*dx)+1)'

'define dx = se*'critval90
'define rUp90'm' = 2*(exp( 2*dx)-1)/(exp( 2*dx)+1)'
'define rLp90'm' = 2*(exp(-2*dx)-1)/(exp(-2*dx)+1)'

'define dx = se*'critval68
'define rUp68'm' = 2*(exp( 2*dx)-1)/(exp( 2*dx)+1)'
'define rLp68'm' = 2*(exp(-2*dx)-1)/(exp(-2*dx)+1)'

m = m + 1
endwhile


************************************************************************
if( PLOT = TRUE )
************************************************************************

* Plot Fisher Mean for Experiments
* --------------------------------
say '  Plot Fisher Mean for Experiments'
'getinfo year'
         year = result
'getinfo date'
         date = result

        m  = 0
while( m<=mexps )
        n  = filebeg
        n.m = n + m*numfiles
        d.m = 1 + m*numfiles
'set dfile 'd.m
'set lev 1000 100'
'set t 'tbeg.m' 'tdim.m

'set vpage off'
'set grads off'
'set gxout contour'
'set parea 2.25 9.75 1.0 7.5'
'set clab  on'
'set xaxis 0 'nday' .5'
'set cmark  0'
'set cthick 8'
'set cstyle 1'
'set ccolor rainbow'
'set cint 'cint
'set rbrange 'rbrange
'd rave'm

'draw ylab Pressure (hPa)'
'set  string 1 c 6 0'
'set  strsiz .17'
'draw string 6.0 8.15 'expdsc.m' ('numfiles')   Anomaly Correlation (CINT: 'cint')'
'draw string 6.0 7.8  'name'  'region
'set  strsiz .12'
'draw string 6.0 0.72 Forecast Day'

'set  string 1 c 6 90'
'set  strsiz .18'
'draw string 0.80 4.25 'months' 'year

say 'EXP'm'  Field: 'name'  Region: 'region

'!/bin/mkdir -p 'SOURCE'/corcmp'
if( nday = ndaymax )
   'myprint -name 'SOURCE'/corcmp/'expdsc.m'_'expdsc.m'_stats_'label'_corcmp_'reg'_z_'months' -rotate 90 -density 100x100'
else
   'myprint -name 'SOURCE'/corcmp/'expdsc.m'_'expdsc.m'_stats_'label'_corcmp_'reg'_z_'months'_'nday'DAY -rotate 90 -density 100x100'
endif
if( debug = "TRUE" )
    say "Hit ENTER for next plot"
    pull flag
endif
'c'

m = m + 1
endwhile

************************************************************************
endif
************************************************************************


* Plot Difference plus Significance
* ---------------------------------
m = 1
while( m<=mexps )
       n  = filebeg
       d.m = 1 + m*numfiles

mfile = 1 + m*files_per_month


* Compute RMS difference between ravem and rave0: ravediff
* --------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m
'makezdif3 -q1  rave'm' -file1 'mfile' -q2 rave0  -file2 1  -ptop 100 -name rave'
'getinfo numfiles'
         newfile = result
'close ' newfile

* Compute difference between 99.99% Confidence and zero: rUp9999diff
* ------------------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m
'makezdif3 -q1   rUp9999'm' -file1 'mfile' -q2  zero  -file2 1  -ptop 100 -name  rUp9999'
'getinfo numfiles'
         newfile = result
'close ' newfile

* Compute difference between 99% Confidence and zero: rUp99diff
* -------------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m
'makezdif3 -q1   rUp99'm' -file1 'mfile' -q2  zero  -file2 1  -ptop 100 -name  rUp99'
'getinfo numfiles'
         newfile = result
'close ' newfile

* Compute difference between 95% Confidence and zero: rUp95diff
* -------------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m
'makezdif3 -q1   rUp95'm' -file1 'mfile' -q2  zero  -file2 1  -ptop 100 -name  rUp95'
'getinfo numfiles'
         newfile = result
'close ' newfile

* Compute difference between 90% Confidence and zero: rUp90diff
* -------------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m
'makezdif3 -q1   rUp90'm' -file1 'mfile' -q2  zero  -file2 1  -ptop 100 -name  rUp90'
'getinfo numfiles'
         newfile = result
'close ' newfile

* Compute difference between 68% Confidence and zero: rUp68diff
* -------------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m
'makezdif3 -q1   rUp68'm' -file1 'mfile' -q2  zero  -file2 1  -ptop 100 -name  rUp68'
'getinfo numfiles'
         newfile = result
'close ' newfile

'set gxout shaded'
'rgbset'

* Define New Variables for Montage Plots
* -------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m

* For ravediff > 0
* ----------------
'define sigdiffp90 = 1000 * ( ravediff-rUp90diff )'
'define sigdiffp95 = 1000 * ( ravediff-rUp95diff )'
'define sigdiffp99 = 1000 * ( ravediff-rUp99diff )'

* For ravediff < 0
* ----------------
'define sigdiffm90 = 1000 * ( ravediff+rUp90diff )'
'define sigdiffm95 = 1000 * ( ravediff+rUp95diff )'
'define sigdiffm99 = 1000 * ( ravediff+rUp99diff )'

'define maskm90 = ( sigdiffm90 - abs(sigdiffm90) )/2'
'define maskm95 = ( sigdiffm95 - abs(sigdiffm95) )/2'
'define maskm99 = ( sigdiffm99 - abs(sigdiffm99) )/2'

'define maskp90 = ( sigdiffp90 + abs(sigdiffp90) )/2'
'define maskp95 = ( sigdiffp95 + abs(sigdiffp95) )/2'
'define maskp99 = ( sigdiffp99 + abs(sigdiffp99) )/2'

'define sigdiff90 = maskm90 + maskp90'
'define sigdiff95 = maskm95 + maskp95'
'define sigdiff99 = maskm99 + maskp99'

'define sigdiffp9999 = 1000 * ( rUp9999diff-rUp90diff )'

* Find maximum value of sigdiff90 across all levels and times
* -----------------------------------------------------------
'set dfile 'ddif.m
'set t 'tbeg.m' 'tdif.m

' minmax sigdiff90'

    qmax = subwrd(result,1)
    qmin = subwrd(result,2)

    xmax = subwrd(result,3)
    ymax = subwrd(result,4)
    zmax = subwrd(result,7)
    tmax = subwrd(result,9)

    xmin = subwrd(result,5)
    ymin = subwrd(result,6)
    zmin = subwrd(result,8)
    tmin = subwrd(result,10)

    qmax = math_abs(qmax)
    qmin = math_abs(qmin)

if( qmin > qmax )
    qmax = qmin
    xmax = xmin
    ymax = ymin
    zmax = zmin
    tmax = tmin
endif

* Compare maximum value of sigdiff90 to 99.99% Confidence Difference
* ------------------------------------------------------------------
'set t 'tmax
'set z 'zmax
'getinfo level'
         level = result
'd sigdiff90'
     qtmp90 = subwrd(result,4)
'd sigdiffp9999'
     qtmp = subwrd(result,4)
say 'initial qmax = 'qmax
say 'tmax = 'tmax
say 'zmax = 'zmax' level: 'level
say 'sigdiffp9999: 'qtmp
say 'sigdiff90: 'qtmp90

if( qtmp > qmax )
    qmax = qtmp
endif
say '  final qmax = 'qmax

'set t 'tbeg.m' 'tdif.m
'set lev 1000 100'

   dcint = qmax / 9

if( PLOT = FALSE )
  if( xpos=2 ) ; 'run setenv DCINT_'field'_rms'rms'.'xpos' 'dcint ; say 'DCINT_'field'_rms'rms'.'xpos' = 'dcint ; endif
  if( xpos=3 ) ; 'run setenv DCINT_'field'_rms'rms'.'xpos' 'dcint ; say 'DCINT_'field'_rms'rms'.'xpos' = 'dcint ; endif
  if( xpos=4 ) ; 'run setenv DCINT_'field'_rms'rms'.'xpos' 'dcint ; say 'DCINT_'field'_rms'rms'.'xpos' = 'dcint ; endif
endif

************************************************************************
if( PLOT = TRUE )
************************************************************************

'run getenv DCINT_'field'_rms'rms'.2'
                             DCINT.2 = result
'run getenv DCINT_'field'_rms'rms'.3'
                             DCINT.3 = result
'run getenv DCINT_'field'_rms'rms'.4'
                             DCINT.4 = result

dcint = DCINT.2
    if( DCINT.3 > dcint ) ; dcint = DCINT.3 ; endif
    if( DCINT.4 > dcint ) ; dcint = DCINT.4 ; endif

flag = ''
while( flag = '' )
'set vpage off'
'set grads off'
'set grid  off'
'set parea 2.25 9.75 1.0 7.5'
'set xaxis 0 'nday' .5'
'set string 1 c 6 0'

'set datawarn off'

* Shade where sigdiff > 90% confidence error bar (color shaded)
* -------------------------------------------------------------
 ' rgbset'
 ' set gxout grfill '
 ' set csmooth off'
 ' shades 'dcint
 ' run getenv SHADES_CLEVS'
              SHADES_CLEVS = result

  clevs = subwrd( SHADES_CLEVS,1 )
  ii = 2
  while( ii <= 9 )
  clev  = subwrd( SHADES_CLEVS,ii )
  clevs = clevs' 'clev
  ii = ii + 1
  endwhile
  clevs = clevs'  -0.0001 0.0001 '
  ii = 10
  while( ii <= 18 )
  clev  = subwrd( SHADES_CLEVS,ii )
  clevs = clevs' 'clev
  ii = ii + 1
  endwhile

 'set clevs 'clevs
 'set ccols 59 57 55 47 44 37 36 34 33 31 0 20 21 22 23 24 25 26 27 28 29'

*' set gxout shaded '
 ' d sigdiff90 '
 ' cbarn -xmid 6 -snum 0.70 -ndot 1'

* Contour sigdiff that is = 90, 95, & 99% confidence diffs (black lines without label)
* ------------------------------------------------------------------------------------
'set gxout contour'
'set csmooth on'
'set clab off'

'set cstyle 2'
'set cthick 8'
'set ccolor 1'
'set clevs  0'
'd sigdiffp99'
'set cstyle 2'
'set cthick 8'
'set ccolor 1'
'set clevs  0'
'd sigdiffm99'

'set cstyle 6'
'set cthick 4'
'set ccolor 1'
'set clevs  0'
'd sigdiffp95'
'set cstyle 6'
'set cthick 4'
'set ccolor 1'
'set clevs  0'
'd sigdiffm95'

'set cstyle 1'
'set cthick 1'
'set ccolor 1'
'set clevs  0'
'd sigdiffp90'
'set cstyle 1'
'set cthick 1'
'set ccolor 1'
'set clevs  0'
'd sigdiffm90'

* reset some background values
* ----------------------------
'set datawarn on'
'set clab on'

* Draw Labels
* -----------
'draw ylab Pressure (hPa)'
'set  string 1 c 6 0'

dcintx = dcint * 100
'getint 'dcintx
         dcintx = result / 100

'set  strsiz .132'
'draw string 6.0 8.15 'expdsc.m' - 'expdsc.0' ('numfiles')   'name'   'region
'draw string 6.0 7.90 Anomaly Correlation Difference  (x10`a-3`n)'
'set  strsiz .125'
'draw string 6.0 7.65 'desc
'set  strsiz .12'
'draw string 6.0 0.72 Forecast Day'

'run uppercase 'field
                FIELD = result

'set  string 1 l 8 0'
'set  strsiz .32'
'draw string 0.70 7.38 'reg
'draw string 0.70 7.01 'FIELD

'set  string 1 l 3 0'
'set  strsiz .087'
'draw string 0.23 1.50 Thin  Solid Line (>90%)'
'draw string 0.23 1.35 Dot -Dash Line (>95%)'
'draw string 0.23 1.20 Long-Dash Line (>99%)'

'set  string 1 c 6 90'
'set  strsiz .18'
'draw string 0.80 4.25 'months' 'year

say 'EXP'm'  Field: 'name'  Region: 'region

'!/bin/mkdir -p 'SOURCE'/corcmp'
if( nday = ndaymax )
   'myprint -name 'SOURCE'/corcmp/'expdsc.0'_'expdsc.m'_stats_'label'_corcmp_'reg'_z_'months' -rotate 90 -density 100x100'
else
   'myprint -name 'SOURCE'/corcmp/'expdsc.0'_'expdsc.m'_stats_'label'_corcmp_'reg'_z_'months'_'nday'DAY -rotate 90 -density 100x100'
endif
if( debug = "TRUE" )
    say "Hit ENTER for next plot"
    pull flag
else
    flag = 'c'
endif
'c'
endwhile

************************************************************************
endif
************************************************************************

m = m + 1
endwhile

return

function getlength (string)
tb = ""
i = 1
while (i<=80)
blank = substr(string,i,1)
if( blank = tb )
length = i-1
i = 81
else
i = i + 1
endif
endwhile
return length

