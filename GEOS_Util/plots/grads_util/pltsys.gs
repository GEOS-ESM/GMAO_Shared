function pltsys (args)

****************************************************************
*****                                                      *****
*****  This script is the driver for calling: pltsys       *****
*****  pltsys computes the ensemble mean systematic error  *****
*****  for forecasts defined in the stats.rc file.         *****
*****                                                      *****
*****  Optional Arguments: -exp -field -level -rc          *****
*****  Examples:                                           *****
*****         1) pltsys.gs                                 *****
*****         2) pltsys.gs -exp 0                          *****
*****         3) pltsys.gs -field v -level 500             *****
*****         4) pltsys.gs -exp 1 -field h                 *****
*****         5) pltsys.gs -exp 1 -field h -rc rcfile      *****
*****                                                      *****
*****  Notes:                                              *****
*****  --------------------------------------------------  *****
*****  Defaults: plot all fields for all experiments       *****
*****       exp: 0 1 2 3 ... n                             *****
*****     field: p u v t q h chi psi                       *****
*****     level: HORIZ   or   ZONAL  or  value             *****
*****    rcfile: stats.rc                                  *****
*****                                                      *****
****************************************************************

'run getenv GEOSUTIL'
        geosutil = result
'run getenv SOURCE'
        SOURCE   = result

if( SOURCE = "NULL" )
   'run getenv "PWD"'
    SOURCE = result
   'run setenv "SOURCE" 'SOURCE
endif

****************************************************************
****************************************************************

'numargs  'args
 numargs = result

rcfile = 'stats.rc'
field  = 'NULL'
level0 = 'NULL'
target = ''

       num = 0
while( num < numargs )
       num = num + 1
if( subwrd(args,num)='-exp'   ) ; target = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-field' ) ; field  = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-level' ) ; level0 = subwrd(args,num+1) ; endif
if( subwrd(args,num)='-rc'    ) ; rcfile = subwrd(args,num+1) ; endif
endwhile

*******************************************************
****             Read Stats Resource File          ****
*******************************************************

n = 0

'getresource 'rcfile' DESC' ;  title = result
'getresource 'rcfile' EXP'n ;  exp.n = result
'getresource 'rcfile' DSC'n ; desc.n = result
 if( exp.n != NULL | desc.n  != NULL )
     n = n+1
 else
     say 'You must supply a CONTROL and COMPARISON experiment list'
     say 'in the file:  stats.rc'
     return
 endif

while( n >= 0 )
'getresource 'rcfile' EXP'n ;  exp.n = result
'getresource 'rcfile' DSC'n ; desc.n = result
 if( exp.n != NULL | desc.n  != NULL )
     n = n+1
 else
     ntot = n-1
        n = -999
 endif
endwhile

if( target = '' )
    k    = 0
    ktot = ntot
else
    k    = target
    ktot = target
endif

while( k <= ktot )

*******************************************************
****                Open kth-Experiment            ****
*******************************************************

'run setenv "ANTIALIAS" TRUE'

'reinit'
'set display color white'
'rgbset'
'c'

'fixname 'desc.k
          DESC.k = result

* Open Files
* ----------
'!'geosutil'/plots/grads_util/make_globl_ctl1 'exp.k' 'DESC.k
'run getenv MONTHLAB'
            month = result
say 'MONTH_LABEL: 'month

'pltsys_open 'desc.k' 'DESC.k
        numfiles = result

* Remove possible BLANKS from Experiment Description
* --------------------------------------------------
DESC = ''
length = getlength(desc.k)
i = 1
while( i<=length )
  bit = substr(desc.k,i,1)
  if( bit != ' ' ) 
      if( DESC = '' ) 
          DESC = bit
      else
          DESC = DESC''bit
      endif
  endif
i = i+1
endwhile
'!/bin/mkdir -p 'SOURCE'/'DESC

*******************************************************
****            Determine Fields in File           ****
*******************************************************

'set dfile 1'
'q ctlinfo'
say 'ctlinfo 'result
'getinfo nvars'
         nvars = result
m=0
n=1
while(n<=nvars)

* Determine Number of Levels Defined for each Field
* -------------------------------------------------
'run getvarz 'n
   name  = subwrd(result,1)
   levs  = subwrd(result,2)
   len   = strlen(name)
   root  = substr(name,1,len-1)
 suffix  = substr(name,len,len)
 if( ( field  = 'NULL' )  | ( field = root ) | ( field = chi & root = u ) | ( field = psi & root = u ) )
     if( suffix = 'f' )
         if( field = chi )
             root  = chi
         endif
         if( field = psi )
             root  = psi
         endif
              m =  m + 1
        field.m = root
         levs.m = levs
     endif
 endif
n = n + 1
endwhile
numflds = m

*******************************************************
****         Begin Systematic Error Plots ...      ****
*******************************************************

'set datawarn off'
'set t 1'

n = 1
while ( n<=numflds )
              field = field.n
               name = field
'set dfile 1'
'setlons'
'sety'

  if( field = p   ) ; name = slp  ; endif
  if( field = h   ) ; name = hght ; endif
  if( field = u   ) ; name = uwnd ; endif
  if( field = v   ) ; name = vwnd ; endif
  if( field = t   ) ; name = tmpu ; endif
  if( field = q   ) ; name = sphu ; endif
  if( field = chi ) ; name = chi  ; endif
  if( field = psi ) ; name = psi  ; endif

*******************************************************
****         Begin Systematic Error Plots ...      ****
*******************************************************

say 'Processing 'field' for pltsys:'

if( level0 = 'NULL' | level0 = 'HORIZ' | level0 = 'ZONAL' )

        levmin  = 1000
        numlevs = levs.n
    if( numlevs > 1 )

*       Create Array of Levels associated with Field
*       --------------------------------------------
        'set z 1'
        'getinfo level'
         levels = result
         z = 2
         while( z<=numlevs )
           'set z 'z
           'getinfo level'
                    level  = result
                 if(level >= 1 )
                    if( level < levmin ) ; levmin = level ; endif
                    levels = levels' 'level
                 endif
         z = z + 1
         endwhile
        'set lev 1000 'levmin

    else
        numlevs = 1
        levels = '1000'
       'set z 1'
    endif

else
        levmin  = level0
        numlevs = 1
        levels = level0

*     Determine if Input Level is Valid
*     ---------------------------------
       'getlevs 'field'f'
          nlevs = result
          valid = false

       if(nlevs > 1)
          z = 1
          while( z<=nlevs )
            'set z 'z
            'getinfo level'
                     level = result
            say 'Checking 'level' against 'level0
            if( level = level0 )
                valid = true
            endif
            z = z + 1
          endwhile
       else
          if( level0 = 1000 )
              valid = true
          endif
       endif

       if( valid = true )
          'set lev 'level0
       else
           return
       endif
endif

say 'NUMLEVS = 'numlevs
say '   LEVS = 'levels

say 'Calling statmak_sys for 'field' 'DESC.k
'q dims'
say 'DIMS: 'result

'statmak_sys 'field'  DESC'


************** 

* Horizonal Plots
* ---------------
if( level0 != 'ZONAL' )
  z = 1
  while ( z<=numlevs )
              level = subwrd(levels,z)
           if(level >= levmin )
             'set dfile 1'
             'setlons'
             'sety'
             'set lev 'level
              say 'Calling Horizontal Movie Statplt for Field: 'field
              say 'Z: 'z
              say 'LEVEL: 'level
             'c'
             'movie statplt "'field' -desc 'DESC' -nfcst 'numfiles' -title 'title'" -print -rotate 90 -name 'SOURCE'/'DESC'/stats_'name'_all_GLO_'level'_'month
             'c'
             '!sleep 60 ; convert -loop 0 -delay 30 'SOURCE'/'DESC'/stats_'name'_all_GLO_'level'_'month'.*.gif 'SOURCE'/'DESC'/stats_'name'_all_GLO_'level'_'month'.gif &'

*              Write Systematic Error File
*              ---------------------------
              'set dfile 1'
              'getinfo undef'
                       undef = result
              'set     undef ' undef
              'setlons'
              'getinfo xmin'
                       xmin = result
              'getinfo xmax'
                       xmax = result - 1
                      level = subwrd(levels,1)
              'sett'
              'set lat -90 90'
              'set x 'xmin' 'xmax
              'set lev 'level

              if( field = p   ) ; name = slp  ; scale = 1    ; endif
              if( field = h   ) ; name = hght ; scale = 1    ; endif
              if( field = u   ) ; name = uwnd ; scale = 1    ; endif
              if( field = v   ) ; name = vwnd ; scale = 1    ; endif
              if( field = t   ) ; name = tmpu ; scale = 1    ; endif
              if( field = q   ) ; name = sphu ; scale = 1000 ; endif
              
              'define 'name' = 'field'fmaDESC/'scale
              'set sdfwrite -4d -flt -nc3 'SOURCE'/'DESC'/'DESC'.'field'fma_'level'.'month'.nc3'
                  'sdfwrite 'name
           endif
  z = z + 1
  endwhile
endif


* Zonal Mean Plots
* ----------------
if( level0 = 'NULL' | level0 = 'ZONAL' )
if( numlevs > 1 )
    say 'Calling Zonal Movie Statpltz for Field: 'field
   'set dfile 1'
   'set x 1'
   'sety'

    if( levmin <= 100 )
       'set lev 1000 100'
       'set zlog off'
       'c'
       'movie statpltz "'field' -desc 'DESC' -nfcst 'numfiles' -title 'title'" -print -rotate 90 -name 'SOURCE'/'DESC'/stats_'name'_all_GLO_z_'month
       'c'
       '!sleep 60 ; convert -loop 0 -delay 30 'SOURCE'/'DESC'/stats_'name'_all_GLO_z_'month'.*.gif 'SOURCE'/'DESC'/stats_'name'_all_GLO_z_'month'.gif &'
    endif

    if( levmin <= 10 & field != q )
       'set lev 1000 10'
       'set zlog off'
       'c'
       'movie statpltz "'field' -desc 'DESC' -nfcst 'numfiles' -title 'title'" -print -rotate 90 -name 'SOURCE'/'DESC'/stats_'name'_all_GLO_zlog10_'month
       'c'
       '!sleep 60 ; convert -loop 0 -delay 30 'SOURCE'/'DESC'/stats_'name'_all_GLO_zlog10_'month'.*.gif 'SOURCE'/'DESC'/stats_'name'_all_GLO_zlog10_'month'.gif &'
    endif

    if( levmin <= 1 & field != q )
       'set lev 1000 1'
       'set zlog off'
       'c'
       'movie statpltz "'field' -desc 'DESC' -nfcst 'numfiles' -title 'title'" -print -rotate 90 -name 'SOURCE'/'DESC'/stats_'name'_all_GLO_zlog1_'month
       'c'
       '!sleep 60 ; convert -loop 0 -delay 30 'SOURCE'/'DESC'/stats_'name'_all_GLO_zlog1_'month'.*.gif 'SOURCE'/'DESC'/stats_'name'_all_GLO_zlog1_'month'.gif &'
    endif

endif
endif

************** 

**************
** Loop over field
**************

n = n + 1
endwhile

**************
** Loop over EXP
**************

k = k+1
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
