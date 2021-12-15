function getvar_2G (args)

EXPORT = subwrd(args,1)
GC     = subwrd(args,2)
SOURCE = subwrd(args,3)

* Check for SUFFIX appended to EXPORT (as EXPORT:SUFFIX)
* ------------------------------------------------------
SUFFIX = ''
         m = 0
         n = 1
       bit = substr(EXPORT,n,1)
while( bit != '' )
         n = n + 1
       bit = substr(EXPORT,n,1)
   if( bit = ':' ) ; m = n ; endif
endwhile
if( m != 0 )
         n = m + 1
       bit = substr(EXPORT,n,1)
    SUFFIX  = ''
    while( bit != '' )
       if( n = m+1 )
           SUFFIX = bit
       else
           SUFFIX = SUFFIX''bit
       endif
         n = n + 1
       bit = substr(EXPORT,n,1)
    endwhile
endif

say ' '
say 'Inside GETVAR_2G,    INPUT_EXPORT = 'EXPORT'  GC = 'GC'  SUFFIX = 'SUFFIX


*'!remove EXPORT.txt'
*'! echo 'EXPORT' | cut -d: -f1 > EXPORT.txt'
*'run getenv "EXPORT"'
*             EXPORT = result
*say 'Inside GETVAR_2G, CHCKHIST_EXPORT = 'EXPORT
*say ' '


* Ensure UpperCase EXPORT and GC
* ------------------------------
*'run uppercase 'EXPORT
*                EXPORT = result
*'run uppercase 'GC
*                GC     = result

* Get Info from HISTORY.rc File
* -----------------------------
'run getenv "GEOSUTIL"'
             geosutil = result

'run getenv "PLOTS_DIR"'
             PLOTS_DIR = result

say 'Running chckhist for EXPORT: 'EXPORT'  GC: 'GC
'!./chckhist 'EXPORT' 'GC' 'SOURCE
say ' finish chckhist for EXPORT: 'EXPORT'  GC: 'GC

say 'CAT hist.txt:'
'!cat hist.txt'

expdsc = sublin( read(hist.txt),2 )
qname  = sublin( read(hist.txt),2 )
qfile  = sublin( read(hist.txt),2 )
scale  = sublin( read(hist.txt),2 )
format = sublin( read(hist.txt),2 )
base   = sublin( read(hist.txt),2 )
           rc = close(hist.txt)

if( qfile = "NULL" )
    say ' '
    say EXPORT' from 'GC' not found in HISTORY.rc!'
    say ' '
    return 'NULL NULL 1 NULL NULL NULL'
endif

'getfile 'qfile
           file = result

* Open New Experiment Dataset
* ---------------------------
say 'Opening: 'qfile
say '   Desc: 'expdsc
say '         '

if( file               = "NULL" )
if( substr(format,1,4) = "flat" ) ; '   open 'qfile ; endif
if( substr(format,1,4) = "CFIO" ) ; 'xdfopen 'qfile ; endif
if( substr(format,1,3) = "HDF"  ) ; 'sdfopen 'qfile ; endif

'getinfo numfiles'
         numfiles = result

* Set CASHE SIZE
* --------------
'getinfo file'
      curfile = result

'set dfile 'numfiles
'getinfo xdim'
         xdim = result
'getinfo ydim'
         ydim = result
         cash = xdim * ydim * 4

'run getenv CASHESZ'
            CASHESZ = result

if( CASHESZ = 'NULL' )
    say 'setting cachesf 'cash
   'set cachesf 'cash
   'run setenv "CASHESZ" 'cash
endif

if( CASHESZ != 'NULL' & cash > CASHESZ )
       say 'updating cachesf 'cash
      'set cachesf 'cash
      'run setenv "CASHESZ" 'cash
endif

* Check for Alias within File
* ---------------------------
  'lowercase  'qname
               qnamelc = result
*if( SUFFIX != 'NULL' )
*               qnamelc = qnamelc''SUFFIX
*                 qname =   qname''SUFFIX
*endif

say 'Looking for 'qnamelc
  'query file'
  numvar = sublin(result,6)
  numvar = subwrd(numvar,5)
  flag = false
  n = 1
  while ( n<numvar+1 )
  field = sublin(result,6+n)
  field = subwrd(field,1)
        if( qnamelc = field )
               flag = true
        endif
  n = n+1
  endwhile

if( flag = false ) ; say qnamelc' not found!' ; endif
'set dfile 'curfile

endif

say 'returning: 'qname' 'numfiles' 'scale' 'expdsc' 'base
return qname' 'numfiles' 'scale' 'expdsc' 'base
