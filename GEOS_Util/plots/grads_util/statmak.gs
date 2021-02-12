function statmak (args)
field  = subwrd  (args,1)
tag    = subwrd  (args,2)
ctl    = subwrd  (args,3)

'numargs  'args
 numargs = result

say ' '
say 'Inside STATMAK, field: 'field
say '                  tag: 'tag
say '                  ctl: 'ctl
say ' '

'run getenv "GEOSUTIL"'
             geosutil = result

if( field = q )
    scale = 1000
else
    scale = 1
endif

* Get Total Number of Files for Subset Experiments
* ------------------------------------------------------------------------
'getinfo numfiles'
         numfiles = result

* Query Level Environment
* -----------------------
'getinfo zfreq'
         zfreq = result

if( zfreq = 'fixed' )
   'getinfo level'
            levmin = result
            levmax = result

   'open 'geosutil'/plots/grads_util/lwmask1440721.tabl'
   'getinfo numfiles'
            newfile = result

   'set dfile 'newfile
   'set t 1'
   'set z 1'
   'define mask = lwmask'
   'close 'newfile
   'set dfile 1'
endif

* Set Proper DFILE for Varying Levels
* -----------------------------------
if( zfreq = 'varying' )
  'getinfo zmin'
           zmin = result
  'getinfo zmax'
           zmax = result
  'set z ' zmin
  'getinfo level'
           levmin = result
  'set z ' zmax
  'getinfo level'
           levmax = result
endif

* Get Proper File for TINC or ZMIN
* --------------------------------
'run getenv "TINCFILE" '
             tincfile = result
'open '      tincfile
'getinfo     numfiles'
             newfile = result
'set dfile ' newfile
'set t 1'
'setlons'
'sety'
'set lev 'levmin' 'levmax

if( zfreq = 'varying' )
   'close 'newfile
   'run getenv "ZMINFILE" '
                zminfile = result
   'open '      zminfile
   'getinfo     numfiles'
                newfile = result
   'set dfile ' newfile
   'set z ' zmin
   'getinfo level'
            levmin = result
   'set z ' zmax
   'getinfo level'
            levmax = result
    if( levmax < levmin )
        levmin = levmax
    endif
   'set lev 1000 'levmin
endif


* Define Number of Forecast Days and Time Interval (hrs)
* ------------------------------------------------------

'run getenv "SYSCMP_TINC"'
                    tinc  = result
'run getenv "SYSCMP_TDIM"'
                    tdim  = result
                    tdum  = tdim - 1
                    ndaymax = tdum * tinc / 24
                       nmax = 1 + ndaymax*(24/tinc)

'run getenv "NDAY"'
             nday = result
         if( nday = "NULL" ) ; nday = ndaymax ; endif
         if( nday > ndaymax) ; nday = ndaymax ; endif

         tbeg = 1-(nmax-tdim)
'set t ' tbeg ' 'tdim

'q file'
say 'Default File: 'result
'q dims'
say 'Default File DIMS: 'result


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
****        Mean Square Error  =   BIAS      Error                          ****
****                           +   Amplitude Error                          ****
****                           +   Phase     Error                          ****
****                                                                        ****
****                      MSE  =   BIA + AMP + PHZ                          ****
****                                                                        ****
********************************************************************************

* Compute forecast statistics
* ---------------------------
*  fma: forecast minus analysis
*  mes: mean error squared
*  mse: mean square error
*  rms: root mean square error
*  std: standard deviation
* --------------------------------------------------

* Compute 2D Fields
* -----------------
'setlons'
'sety'

'define 'field'fm'tag'  = lat-lat + lon-lon + lev-lev'
'define 'field'am'tag'  = lat-lat + lon-lon + lev-lev'
'define 'field'fma'tag' = lat-lat + lon-lon + lev-lev'
'define 'field'mse'tag' = lat-lat + lon-lon + lev-lev'

n = 1
while ( n <= numfiles )
say 'Processing Field: 'field' for File: 'n'  Tag: 'tag

* Note: Add and Subtract uf & vf to force similar UNDEF locations
* ---------------------------------------------------------------
   if( field = chi )
       'define  uaa'n' = ua.'n'+uf.'n'-uf.'n
       'define  vaa'n' = va.'n'+vf.'n'-vf.'n
       'define chif'n' = fish_chi(uf.'n',vf.'n')'
       'define chia'n' = fish_chi(uaa'n',vaa'n')'
       'define chif'n' = chif'n'-aave(chif'n',g)'
       'define chia'n' = chia'n'-aave(chia'n',g)'
   endif
   if( field = psi )
       'define  uaa'n' = ua.'n'+uf.'n'-uf.'n
       'define  vaa'n' = va.'n'+vf.'n'-vf.'n
       'define psif'n' = fish_psi(uf.'n',vf.'n')'
       'define psia'n' = fish_psi(uaa'n',vaa'n')'
       'define psif'n' = psif'n'-aave(psif'n',g)'
       'define psia'n' = psia'n'-aave(psia'n',g)'
   endif

   if( field  = chi | field  = psi )
       'define 'field'fm'tag'  = 'field'fm'tag'  +     'field'f'n
       'define 'field'am'tag'  = 'field'am'tag'  +     'field'a'n
       'define 'field'fma'tag' = 'field'fma'tag' +     'field'f'n'-'field'a'n
       'define 'field'mse'tag' = 'field'mse'tag' + pow('field'f'n'-'field'a'n',2)'
       'define 'field'mse'tag''n' = pow('field'f'n'-'field'a'n',2)'
   else
       'define 'field'fs  = 'field'f.'n'*'scale
       'define 'field'as  = 'field'a.'n'*'scale
       'define 'field'cs  = 'field'c.'n'*'scale

       'define 'field'fm'tag'  = 'field'fm'tag'  +     'field'fs'
       'define 'field'am'tag'  = 'field'am'tag'  +     'field'as'
       'define 'field'fma'tag' = 'field'fma'tag' +     'field'fs-'field'as'
       'define 'field'mse'tag' = 'field'mse'tag' + pow('field'fs-'field'as,2)'
       'define 'field'mse'tag''n' = pow('field'fs-'field'as,2)'
   endif

n = n + 1
endwhile

'define 'field'fm'tag'  = 'field'fm'tag' /'numfiles
'define 'field'am'tag'  = 'field'am'tag' /'numfiles
'define 'field'fma'tag' = 'field'fma'tag'/'numfiles
'define 'field'mse'tag' = 'field'mse'tag'/'numfiles

'undefine 'field'fs'
'undefine 'field'as'



'define 'field'mes'tag'  = lat-lat + lon-lon + lev-lev'
'define 'field'varf'tag' = lat-lat + lon-lon + lev-lev'
'define 'field'vara'tag' = lat-lat + lon-lon + lev-lev'
'define 'field'cov'tag'  = lat-lat + lon-lon + lev-lev'
'define 'field'rnd'tag'  = lat-lat + lon-lon + lev-lev'
n = 1
while ( n <= numfiles )
'define 'field'mes'tag'  = 'field'mes'tag'  + pow(  'field'fm'tag'-'field'am'tag',2)'
'define 'field'varf'tag' = 'field'varf'tag' + pow(  'field'f.'n'*'scale'-'field'fm'tag',2)'
'define 'field'vara'tag' = 'field'vara'tag' + pow(  'field'a.'n'*'scale'-'field'am'tag',2)'
'define 'field'cov'tag'  = 'field'cov'tag'  +      ('field'f.'n'*'scale'-'field'fm'tag') * ('field'a.'n'*'scale'-'field'am'tag')'
'define 'field'rnd'tag'  = 'field'rnd'tag'  + pow( ('field'f.'n'*'scale'-'field'fm'tag') - ('field'a.'n'*'scale'-'field'am'tag') , 2)'
n = n + 1
endwhile
'define 'field'mes'tag'  = 'field'mes'tag' /'numfiles
'define 'field'varf'tag' = 'field'varf'tag'/'numfiles
'define 'field'vara'tag' = 'field'vara'tag'/'numfiles
'define 'field'cov'tag'  = 'field'cov'tag' /'numfiles
'define 'field'rnd'tag'  = 'field'rnd'tag' /'numfiles


'define 'field'stdf'tag' = sqrt('field'varf'tag')'
'define 'field'stda'tag' = sqrt('field'vara'tag')'

'define 'field'ampl'tag' = pow( 'field'stdf'tag'-'field'stda'tag',2 )'
'define 'field'phaz'tag' =  2*( 'field'stdf'tag'*'field'stda'tag' - 'field'cov'tag')'
'define 'field'ramp'tag' = sqrt( 'field'ampl'tag' )'
'define 'field'rphz'tag' = sqrt( 'field'phaz'tag' )'
'define 'field'rrnd'tag' = sqrt( 'field'rnd'tag' )'

'define 'field'rmes'tag' = sqrt('field'mes'tag')'
'define 'field'rms'tag'  = sqrt('field'mse'tag')'
'define 'field'var'tag'  =      'field'mse'tag'-'field'mes'tag
'define 'field'std'tag'  = sqrt('field'mse'tag'-'field'mes'tag')'

* -----------------------------------------------------------------


'define 'field'rmes'tag' = sqrt('field'mes'tag')'
'define 'field'ramp'tag' = sqrt('field'ampl'tag')'
'define 'field'rphz'tag' = sqrt('field'phaz'tag')'


* -----------------------------------------------------------------

if( zfreq = 'varying' )
   'makez  'field'fm'tag'   z'
   'makez  'field'am'tag'   z'
   'makez  'field'fma'tag'  z'
   'makez  'field'mes'tag'  z'
   'makez  'field'mse'tag'  z'
   'makez  'field'ampl'tag' z'
   'makez  'field'phaz'tag' z'
   'define 'field'rmes'tag'z = sqrt('field'mes'tag'z)'
   'define 'field'ramp'tag'z = sqrt('field'ampl'tag'z)'
   'define 'field'rphz'tag'z = sqrt('field'phaz'tag'z)'
   'define 'field'rms'tag'z  = sqrt('field'mse'tag'z)'
   'define 'field'var'tag'z  =      'field'mse'tag'z-'field'mes'tag'z'
   'define 'field'std'tag'z  = sqrt('field'mse'tag'z-'field'mes'tag'z)'
endif


* Compute TAG & CTL Variables: X & Y, and Diff Variable: Z = X-Y
* For: mean-error-squared(mes), mean-square-error(mse), variance(var)
* -------------------------------------------------------------------

if( tag != ctl )
    say 'Computing Difference Variables for Field: 'field' and TAGs: 'tag' and 'ctl
    say '----------------------------------------- '

* Define Difference-Variables
* ---------------------------
* RMS:    Dmse  =   mse_tag - mse_ctl   =  { 1/N * SUM[ (F-A)   **2      ] }_tag
*                                       -  { 1/N * SUM[ (F-A)   **2      ] }_ctl
* BIAS:   Dmes  =   mes_tag - mes_ctl   =  { 1/N * SUM[ (F-A)            ] }_tag **2
*                                       -  { 1/N * SUM[ (F-A)            ] }_ctl **2
* VAR:    Dvar  =  Xvar_tag - Xvar_ctl  =  { 1/N * SUM[ ( (F-A)-(FBAR-ABAR) )**2 ] }_tag
*                                       -  { 1/N * SUM[ ( (F-A)-(FBAR-ABAR) )**2 ] }_ctl
* AMP:    Damp  =  ampl_tag - ampl_ctl
* PHZ:    Dphz  =  phaz_tag - phaz_ctl
* ----------------------------------------------------------------------

          'define 'field'Dmse'tag' = 'field'mse'tag'  - 'field'mse'ctl
          'define 'field'Dmes'tag' = 'field'mes'tag'  - 'field'mes'ctl
          'define 'field'Damp'tag' = 'field'ampl'tag' - 'field'ampl'ctl
          'define 'field'Dphz'tag' = 'field'phaz'tag' - 'field'phaz'ctl

           n  = 1
   while ( n <= numfiles )
          'define 'field'Dmse'tag''n' = 'field'mse'tag''n' - 'field'mse'ctl''n
           n = n + 1
   endwhile

* Define Variances of MSE Differences
* -----------------------------------
          'define 'field'DDmse'tag' = lat-lat + lon-lon + lev-lev'
           n  = 1
   while ( n <= numfiles )
          'define 'field'DDmse'tag' = 'field'DDmse'tag' + pow( 'field'Dmse'tag''n'-'field'Dmse'tag',2 )'
           n = n + 1
   endwhile
          'define 'field'DDmse'tag' = 'field'DDmse'tag' / 'numfiles

* --------------------------------------------------------------
* --------------------------------------------------------------

if( zfreq = 'varying' )
   say 'Computing Zonal Mean of Difference Variables for Field: 'field' and TAGs: 'tag' and 'ctl
   say '------------------------------------------------------ '
   'makez 'field'Dmes'tag'  z'
   'makez 'field'Dmse'tag'  z'
   'makez 'field'Damp'tag'  z'
   'makez 'field'Dphz'tag'  z'

           n  = 1
   while ( n <= numfiles )
          'makez 'field'Dmse'tag''n' z'
           n = n + 1
   endwhile

          'set lon 0'
          'define 'field'DDmse'tag'z = lat-lat + lev-lev'
           n  = 1
   while ( n <= numfiles )
          'define 'field'DDmse'tag'z = 'field'DDmse'tag'z + pow( 'field'Dmse'tag''n'z-'field'Dmse'tag'z,2 )'
           n = n + 1
   endwhile
          'define 'field'DDmse'tag'z = 'field'DDmse'tag'z / 'numfiles
endif

* --------------------------------------------------------------
* --------------------------------------------------------------

* End CTL Test
* ------------
endif
'close 'newfile


'set dfile 1'
'set t 'tbeg' 'tdim
'setlons'
'sety'
'set lev 'levmin' 'levmax

return
