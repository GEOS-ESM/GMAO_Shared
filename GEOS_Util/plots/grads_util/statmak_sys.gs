function statmak (args)
field  = subwrd  (args,1)
tag    = subwrd  (args,2)

say ' '
say 'Inside STATMAK_SYS, field: 'field
say '                      tag: 'tag
say ' '

if( field = q )
    scale = 1000
else
    scale = 1
endif

'numargs  'args
 numargs = result

* Query Level Environment
* -----------------------
'getinfo zfreq'
         zfreq = result
     if( zfreq = 'varying' )
         'getinfo zmin'
                  zmin = result
         'getinfo zmax'
                  zmax = result
     endif
     if( zfreq = 'fixed' )
         'getinfo zpos'
          zmin = result
          zmax = result
     endif

'run getenv "GEOSUTIL"'
             geosutil = result

'getinfo numfiles'
         numfiles = result

'open 'geosutil'/plots/grads_util/lwmask1440721.tabl'
'getinfo numfiles'
         newfile = result

'set dfile 'newfile
'set t 1'
'set z 1'
'define mask = lwmask'
'close 'newfile

'set dfile 1'
'set t 1'
'set z 'zmin' 'zmax

* Define Number of Forecast Days and Time Interval (hrs)
* ------------------------------------------------------
'run getinfo tinc'
             tinc = result
         if( tinc = "NULL" ) ; tinc = 6 ; endif

'run getinfo tdim'
             tdim = result
             tdum = tdim - 1
             ndaymax = tdum * tinc / 24
                nmax = 1 + ndaymax*(24/tinc)

'run getenv "NDAY"'
             nday = result
         if( nday = "NULL" ) ; nday = ndaymax ; endif
         if( nday > ndaymax) ; nday = ndaymax ; endif

         tbeg = 1-(nmax-tdim)
'set t ' tbeg ' 'tdim

* Compute forecast statistics
* ---------------------------
*  fma: forecast minus analysis
*  fmc: forecast minus climatology
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
'define 'field'cm'tag'  = lat-lat + lon-lon + lev-lev'
'define 'field'fma'tag' = lat-lat + lon-lon + lev-lev'
'define 'field'fmc'tag' = lat-lat + lon-lon + lev-lev'
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
       'define chic'n' = fish_chi(uc.'n',vc.'n')'
       'define chif'n' = chif'n'-aave(chif'n',g)'
       'define chia'n' = chia'n'-aave(chia'n',g)'
       'define chic'n' = chic'n'-aave(chic'n',g)'
   endif
   if( field = psi )
       'define  uaa'n' = ua.'n'+uf.'n'-uf.'n
       'define  vaa'n' = va.'n'+vf.'n'-vf.'n
       'define psif'n' = fish_psi(uf.'n',vf.'n')'
       'define psia'n' = fish_psi(uaa'n',vaa'n')'
       'define psic'n' = fish_psi(uc.'n',vc.'n')'
       'define psif'n' = psif'n'-aave(psif'n',g)'
       'define psia'n' = psia'n'-aave(psia'n',g)'
       'define psic'n' = psic'n'-aave(psic'n',g)'
   endif

   if( field  = chi | field  = psi )
       'define 'field'fm'tag'  = 'field'fm'tag'  +     'field'f'n
       'define 'field'am'tag'  = 'field'am'tag'  +     'field'a'n
       'define 'field'cm'tag'  = 'field'cm'tag'  +     'field'c'n
       'define 'field'fma'tag' = 'field'fma'tag' +     'field'f'n'-'field'a'n
       'define 'field'fmc'tag' = 'field'fmc'tag' +     'field'f'n'-'field'c'n
       'define 'field'mse'tag' = 'field'mse'tag' + pow('field'f'n'-'field'a'n',2)'
   else
       'define 'field'fs  = 'field'f.'n'*'scale
       'define 'field'as  = 'field'a.'n'*'scale
       'define 'field'cs  = 'field'c.'n'*'scale

       'define 'field'fm'tag'  = 'field'fm'tag'  +     'field'fs'
       'define 'field'am'tag'  = 'field'am'tag'  +     'field'as'
       'define 'field'cm'tag'  = 'field'cm'tag'  +     'field'cs'
       'define 'field'fma'tag' = 'field'fma'tag' +     'field'fs-'field'as'
       'define 'field'fmc'tag' = 'field'fmc'tag' +     'field'fs-'field'cs'
       'define 'field'mse'tag' = 'field'mse'tag' + pow('field'fs-'field'as,2)'
   endif

n = n + 1
endwhile

'define 'field'fm'tag'  = 'field'fm'tag' /'numfiles
'define 'field'am'tag'  = 'field'am'tag' /'numfiles
'define 'field'cm'tag'  = 'field'cm'tag' /'numfiles
'define 'field'fma'tag' = 'field'fma'tag'/'numfiles
'define 'field'fmc'tag' = 'field'fmc'tag'/'numfiles
'define 'field'mse'tag' = 'field'mse'tag'/'numfiles

'undefine 'field'fs'
'undefine 'field'as'
'undefine 'field'cs'



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
   'makez  'field'cm'tag'   z'
   'makez  'field'fma'tag'  z'
   'makez  'field'fmc'tag'  z'
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

return

function region(args)

reg = subwrd(args,1)

       'set dfile 1'

   if( reg = 'GLO'   )
       'setx'
       'sety'
       'run getinfo xdim'
                    xend = result
                    xbeg = 1
       'run getinfo ydim'
                    yend = result
                    ybeg = 1
   endif
   if( reg = 'NHE'   )
       'setx'
       'set lat 20 90'
       'run getinfo xdim'
                    xend = result
                    xbeg = 1
       'run getinfo ymax'
                    yend = result
       'run getinfo ymin'
                    ybeg = result
   endif
   if( reg = 'TRO'   )
       'setx'
       'set lat -20 20'
       'run getinfo xdim'
                    xend = result
                    xbeg = 1
       'run getinfo ymax'
                    yend = result
       'run getinfo ymin'
                    ybeg = result
   endif
   if( reg = 'SHE'   )
       'setx'
       'set lat -90 -20'
       'run getinfo xdim'
                    xend = result
                    xbeg = 1
       'run getinfo ymax'
                    yend = result
       'run getinfo ymin'
                    ybeg = result
   endif

return xbeg' 'xend' 'ybeg' 'yend

function masks (name,latbeg,latend)
'set dfile 1'
'set lon -180 360'
'sety'
'setz'
'define masc = maskout( 1    ,lat-'latbeg' )'
'define masc = maskout( masc,'latend'-lat )'
'set sdfwrite -4d mask.nc4'
'set undef 0.0'
'sdfwrite masc'

'sdfopen mask.nc4'
'getinfo  numfiles'
          newfile = result
'set dfile 'newfile

'q ctlinfo'
   ctlinfo = result

          fname = mask.ctl
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

'close 'newfile
'open  'fname
'getinfo numfiles'
         newfile = result
'set dfile 'newfile
'define 'name' = masc.'newfile
'close 'newfile
return

