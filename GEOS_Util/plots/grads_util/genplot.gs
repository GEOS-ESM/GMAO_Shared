function genplot (args)

*******************************************************
****                 INPUT Variables               ****
*******************************************************

'numargs  'args
 numargs = result

PREFIX = NULL
TAYLOR = FALSE
LAND   = FALSE
OCEAN  = FALSE
SCALE  = 1.0
RC     = NULL
LEVEL  = NULL
MATH   = NULL

          m = 0
          n = 0
        num = 0
while ( num < numargs )
        num = num + 1

if( subwrd(args,num) = '-DIR'    ) ; DIR    = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-NAME'   ) ; NAME   = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-MATH'   ) ; MATH   = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-EXPID'  ) ; EXPID  = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-OUTPUT' ) ; OUTPUT = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-DEBUG'  ) ; DEBUG  = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-SCALE'  ) ; SCALE  = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-PREFIX' ) ; PREFIX = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-RC'     ) ; RC     = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-LEVEL'  ) ; LEVEL  = subwrd(args,num+1) ; endif
if( subwrd(args,num) = '-LAND'   ) ; LAND   = TRUE               ; endif
if( subwrd(args,num) = '-OCEAN'  ) ; OCEAN  = TRUE               ; endif
if( subwrd(args,num) = '-TAYLOR' ) ; TAYLOR = TRUE               ; endif

* Read Model EXPORT:GC
* --------------------
if( subwrd(args,num) = '-EXPORT' )
              m = m + 1
       EXPORT.m = subwrd(args,num+m   )
           word = subwrd(args,num+m+1 )
            bit = checkbit(word)
       while( bit != '-' )
              m = m + 1
       EXPORT.m = subwrd(args,num+m   )
           word = subwrd(args,num+m+1 )
            bit = checkbit(word)
       endwhile
endif

* Read Verification OBS:OBSGC
* ---------------------------
if( subwrd(args,num) = '-OBS' )
              n = n + 1
          OBS.n = subwrd(args,num+n   )
           word = subwrd(args,num+n+1 )
            bit = checkbit(word)
       while( bit != '-' )
              n = n + 1
          OBS.n = subwrd(args,num+n   )
           word = subwrd(args,num+n+1 )
            bit = checkbit(word)
       endwhile
endif

* Read SEASONS
* -----------
if( subwrd(args,num) = '-SEASON' )
     seasons = ''
           k = 1
    while( k > 0 )
           L = num + k
        season  = subwrd(args,L)
    if( season  = '' )
        k = -1
    else
        bit = substr(season,1,1)
        if( bit = '-' )
              k = -1
        else
              seasons = seasons % ' ' % season
              k = k+1
        endif
    endif
    endwhile
endif

endwhile


* Construct Model GCs from Model EXPORTS
* --------------------------------------
numGCs = m
        m  = 0
        k  = 1
while ( k <= numGCs )
        EX = ''
         j = 1
       bit = substr(EXPORT.k,j,1)
       while(bit != ':' & bit != '')
        EX = EX''bit
         j = j + 1
       bit = substr(EXPORT.k,j,1)
       endwhile
       if( EX != EXPORT.k )
         m = m + 1
         j = j + 1
       GC.m = ''
       bit = substr(EXPORT.k,j,1)
       while(bit != '')
       GC.m = GC.m''bit
         j = j + 1
       bit = substr(EXPORT.k,j,1)
       endwhile
       EXPORT.k = EX
       endif
k = k + 1
endwhile


* Construct OBS GCs from OBS EXPORTS
* ----------------------------------
    numn = n
if( numn = 0 )
    numn = numGCs
           n  = 1
   while ( n <= numn )
       OBS.n = EXPORT.n
     OBSGC.n =     GC.n
           n = n + 1
   endwhile
else
        m  = 0
        k  = 1
while ( k <= numn )
        EX = ''
         j = 1
       bit = substr(OBS.k,j,1)
       while(bit != ':' & bit != '')
        EX = EX''bit
         j = j + 1
       bit = substr(OBS.k,j,1)
       endwhile
       if( EX != OBS.k )
         m = m + 1
         j = j + 1
       OBSGC.m = ''
       bit = substr(OBS.k,j,1)
       while(bit != '')
       OBSGC.m = OBSGC.m''bit
         j = j + 1
       bit = substr(OBS.k,j,1)
       endwhile
       OBS.k = EX
       endif
k = k + 1
endwhile
endif

**************************************************
**************************************************

* Initialize
* ----------
'reinit'
'set display color white'
'set clab off'
'c'

'run getenv "GEOSUTIL"'
             geosutil = result

'run getenv "VERIFICATION"'
             verification = result

'run uppercase 'seasons
                seasons = result

'run getenv "CMPEXP_ONLY"'
             cmpexp_only = result

**************************************************
****            Echo Calling Sequence         ****
**************************************************

say ' '
say 'NAME   = 'NAME
say 'EXPID  = 'EXPID
say 'OUTPUT = 'OUTPUT
say 'DEBUG  = 'DEBUG
say 'SCALE  = 'SCALE
say 'RC     = 'RC
say 'LAND   = 'LAND
say 'OCEAN  = 'OCEAN
say 'TAYLOR = 'TAYLOR
say 'SEASON = 'seasons
say ' '
m = 1
while( m<=numGCs )
say 'EXPORT.'m' = 'EXPORT.m
say '    GC.'m' = 'GC.m
m = m + 1
endwhile
say ' '
n = 1
while( n<=numn )
say '   OBS.'n' = 'OBS.n
say ' OBSGC.'n' = 'OBSGC.n
n = n + 1
endwhile
say ' '

**************************************************
**************************************************

* Get Model Variables
* -------------------
        m  = 1
while ( m <= numGCs )
'run getvar 'EXPORT.m' 'GC.m
        qname.m = subwrd(result,1)
        qfile.m = subwrd(result,2)
        qscal.m = subwrd(result,3)
        qdesc.m = subwrd(result,4)
         qtag.m = subwrd(result,5)
    if( qname.m = 'NULL' ) ; return ; endif
         m  = m + 1
endwhile


* Set proper ZDIM
* ---------------
 'set dfile 'qfile.1
 'getlevs   'qname.1
    nlevs = result
if( nlevs != 'NULL' )
   'run setenv "ZDIM" 'nlevs
endif


* Ensure NAMES have no underscores
* --------------------------------
        m=1
while ( m<numGCs+1 )
'fixname 'qname.m
          alias.m = result
say 'Alias #'m' = 'alias.m
        m = m+1
endwhile

* Set Geographic Environment from Model Dataset
* ---------------------------------------------
'set dfile 'qfile.1
'setlons'
'setlats'
'getinfo dlon'
         dlon = result
'getinfo dlat'
         dlat = result
'getinfo lonmin'
         lonmin = result
'getinfo latmin'
         latmin = result

* Check for Model Name Consistency
* --------------------------------
 m = 1
while( m <= numGCs )
'set dfile 'qfile.m
if( LEVEL = "NULL" )
   'set z 1'
else
   'set lev 'LEVEL
endif
'sett'
if( qname.m != alias.m ) ; 'rename 'qname.m ' 'alias.m ; endif
m = m + 1
endwhile

* Set BEGDATE and ENDDATE for Seasonal Calculations
* -------------------------------------------------
'setdates'

* Extract Beginning and Ending Dates for Plots
* --------------------------------------------
'run getenv "BEGDATE"'
             begdate  = result
'run getenv "ENDDATE"'
             enddate  = result
if( begdate = "NULL" )
   'set dfile 'qfile.1
   'set t    '1
   'getinfo date'
         begdate = result
endif
if( enddate = "NULL" )
   'set dfile 'qfile.1
   'getinfo tdim'
            tdim     = result
   'set t  'tdim
   'getinfo date'
         enddate = result
endif


* Land/Water Masks
* ----------------
if( LAND = 'TRUE' | OCEAN = 'TRUE' )
   'set t 1'
   'setmask     mod'
   'define lwmaskmod = regrid2( lwmaskmod,0.25,0.25,bs_p1,'lonmin','latmin')'
   'define  omaskmod = maskout( 1, lwmaskmod-0.5 )'
   'define  lmaskmod = maskout( 1, 0.5-lwmaskmod )'
endif


* Perform Model Formula Calculation
* ---------------------------------
'set dfile 'qfile.1
'sett'
if( numGCs = 1 )
   'define qmod = 'alias.1'.'qfile.1'*'qscal.1
else
    filename  = geosutil'/plots/'NAME'/modform.gs'
    ioflag    = sublin( read(filename),1 )
    if(ioflag = 0)
       close  = close(filename)
       mstring = ''
       m  = 1
       while ( m <= numGCs )
          if( qname.m != alias.m )
              mstring = mstring' 'alias.m'*'qscal.m
          else
              mstring = mstring' 'alias.m'.'qfile.m'*'qscal.m
          endif
              m  = m + 1
       endwhile
      'run 'geosutil'/plots/'NAME'/modform 'mstring
    else
       mstring = NAME
       m  = 1
       while ( m <= numGCs )
          if( qname.m != alias.m )
              mstring = mstring' 'alias.m'*'qscal.m
          else
              mstring = mstring' 'alias.m'.'qfile.m'*'qscal.m
          endif
              m  = m + 1
       endwhile
      'run 'geosutil'/plots/'DIR'/'NAME'.gs 'mstring
      'define qmod = 'NAME'.1'
    endif
endif
                          'define qmod = regrid2( qmod,0.25,0.25,bs_p1,'lonmin','latmin')'
if(    LAND  = 'TRUE'   |  OCEAN = 'TRUE' )
   if( LAND  = 'TRUE' ) ; 'define qmod = maskout( 'SCALE'*qmod,lmaskmod )' ; endif
   if( OCEAN = 'TRUE' ) ; 'define qmod = maskout( 'SCALE'*qmod,omaskmod )' ; endif
else
                          'define qmod =          'SCALE'*qmod'
endif

'seasonal qmod'

* Create Dummy File with REGRID Dimensions
* ----------------------------------------
  'define qmodr = qmod'
  'getinfo undef'
           undef = result
  'set sdfwrite -5d regrid.nc4'
  'set undef 'undef
  'sdfwrite qmodr'
  'sdfopen regrid.nc4'
  'getinfo    numfiles'
              rgfile = result


***********************************************************************************
*              Loop over Possible Experiment Datasets for Comparison
***********************************************************************************

'!/bin/mv HISTORY.T HISTORY.Tmp'
'run getenv "CMPEXP"'
         cmpexp = result
         numexp = 1

          dummy = get_cmpexp (cmpexp,numexp)
            exp = subwrd(dummy,1)
           type = subwrd(dummy,2)

while( exp != 'NULL' )
say ' '
say 'Comparing with: 'exp

* analysis = false  EXP=M CMP=M  => ALEVS
* analysis = false  EXP=M CMP=A  => DLEVS
* analysis = true   EXP=A CMP=A  => ALEVS
* analysis = true   EXP=A CMP=M  => DLEVS

if( analysis != "false" )
    if( type = A )
       'run setenv "LEVTYPE" 'ALEVS
    else
       'run setenv "LEVTYPE" 'DLEVS
    endif
else
    if( type = A )
       'run setenv "LEVTYPE" 'DLEVS
    else
       'run setenv "LEVTYPE" 'ALEVS
    endif
endif

'!chckfile 'exp'/.HOMDIR'
 'run getenv CHECKFILE'
         CHECKFILE  = result
     if( CHECKFILE != 'NULL' )
        '!/bin/cp `cat 'exp'/.HOMDIR`/HISTORY.rc .'
     else
        '!/bin/cp 'exp'/HISTORY.rc .'
     endif
'!remove CHECKFILE.txt'

'!cat HISTORY.rc | sed -e "s/,/ , /g" | sed -e "s/*/@/g" > HISTORY.T'

* Get EXP Comparison Variables
* ----------------------------
FOUND = TRUE
        m  = 1
while ( m <= numGCs & FOUND = TRUE )
'run getvar 'EXPORT.m' 'GC.m' 'exp
        cname.numexp.m = subwrd(result,1)
        cfile.numexp.m = subwrd(result,2)
        cscal.numexp.m = subwrd(result,3)
        cdesc.numexp.m = subwrd(result,4)
         ctag.numexp.m = subwrd(result,5)
say ''
    if( cname.numexp.m = 'NULL' ) ; FOUND = FALSE ; endif
         m  = m + 1
endwhile

if( FOUND = TRUE )
'setlons'
'setlats'

* Land/Water Masks
* ----------------
if( LAND = 'TRUE' | OCEAN = 'TRUE' )
   'set dfile 'cfile.numexp.1
if( LEVEL = "NULL" )
   'set z 1'
else
   'set lev 'LEVEL
endif
   'set t 1'
*  'setmask     obs'
*  'define lwmaskobs = regrid2( lwmaskobs,0.25,0.25,bs_p1,'lonmin','latmin')'
   'define lwmaskobs = lwmaskmod'
   'define  omaskobs = maskout( 1, lwmaskobs-0.5 )'
   'define  lmaskobs = maskout( 1, 0.5-lwmaskobs )'
endif

'set dfile 'cfile.numexp.1
if( LEVEL = "NULL" )
   'set z 1'
else
   'set lev 'LEVEL
endif
    'getdates'
     begdateo = subwrd(result,1)
     enddateo = subwrd(result,2)

* Perform OBS Formula Calculation
* -------------------------------
if( numGCs = 1 )
   'define cmod'numexp' = 'cname.numexp.1'.'cfile.numexp.1'*'cscal.numexp.1
else
    filename  = geosutil'/plots/'NAME'/modform.gs'
    ioflag    = sublin( read(filename),1 )
    if(ioflag = 0)
       close  = close(filename)
       cstring = ''
       n  = 1
       while ( n <= numGCs )
          cstring = cstring' 'cname.numexp.n'.'cfile.numexp.n'*'cscal.numexp.n
               n  = n + 1
       endwhile
      'run 'geosutil'/plots/'NAME'/modform 'cstring
      'define cmod'numexp' = qmod'
    else
       cstring = NAME
       n  = 1
       while ( n <= numGCs )
          ostring = ostring' 'cname.numexp.n'.'cfile.numexp.n'*'cscal.numexp.n
               n  = n + 1
       endwhile
      'run 'geosutil'/plots/'DIR'/'NAME'.gs 'ostring
      'define cmod'numexp' = 'NAME'.1'
    endif
endif

* Check for MERRA Resolution problem
* ----------------------------------
'getinfo xdim'
         xdim = result
if( xdim != 540 )

                              'define cmod'numexp' = regrid2( cmod'numexp',0.25,0.25,bs_p1,'lonmin','latmin')'
    if(    LAND  = 'TRUE'   |  OCEAN = 'TRUE' )
       if( LAND  = 'TRUE' ) ; 'define cmod'numexp' = maskout( 'SCALE'*cmod'numexp',lmaskobs )' ; endif
       if( OCEAN = 'TRUE' ) ; 'define cmod'numexp' = maskout( 'SCALE'*cmod'numexp',omaskobs )' ; endif
    else
                              'define cmod'numexp' =          'SCALE'*cmod'numexp
    endif

else

*                             'define cmod'numexp' = regrid2( cmod'numexp','dlon','dlat',bs_p1,'lonmin','latmin')'
                              'define cmod'numexp' = regrid2( cmod'numexp',0.25,0.25,bs_p1,'lonmin','latmin')'
    if(    LAND  = 'TRUE'   |  OCEAN = 'TRUE' )
       if( LAND  = 'TRUE' ) ; 'define cmod'numexp' = maskout( 'SCALE'*cmod'numexp',lmaskobs )' ; endif
       if( OCEAN = 'TRUE' ) ; 'define cmod'numexp' = maskout( 'SCALE'*cmod'numexp',omaskobs )' ; endif
    else
                              'define cmod'numexp' =          'SCALE'*cmod'numexp
    endif

endif

* Compute Seasonal Means
* ----------------------
'seasonal cmod'numexp


'run getenv "CLIMATE"'
             climate = result

* Loop over Seasons to Process
* ----------------------------
       m = 1
while( m > 0 )
    season = subwrd(seasons,m)
if( season = '' )
         m = -1
else
         m = m+1
         say 'Processing Season: 'season

'set dfile 'qfile.1
'set gxout shaded'
'rgbset'

* Horizontal Plot
* ---------------
       mathparm  =  MATH
while( mathparm != 'DONE' )
        flag = ""
while ( flag = "" )

'makplot -MVAR 'qmod' -MNAME 'NAME ' -MFILE 'qfile.1' -MDESC 'qdesc.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OVAR 'cmod''numexp' -ONAME 'ctag.numexp.1' -OFILE 'cfile.numexp.1' -ODESC 'cdesc.numexp.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm

 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile ;* END While_FLAG Loop
       if( mathparm != 'NULL' )
           mathparm  = 'NULL'
       else
           mathparm  = 'DONE'
       endif
endwhile ;* END While_MATH Loop

* Zonal Mean Plot
* ---------------
       mathparm  =  MATH
while( mathparm != 'DONE' )
        flag = ""
while ( flag = "" )

'makplotz -MVAR 'qmod' -MNAME 'NAME ' -MFILE 'qfile.1' -MDESC 'qdesc.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OVAR 'cmod''numexp' -ONAME 'ctag.numexp.1' -OFILE 'cfile.numexp.1' -ODESC 'cdesc.numexp.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm' -RGFILE 'rgfile

 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile ;* END While_FLAG Loop
       if( mathparm != 'NULL' )
           mathparm  = 'NULL'
       else
           mathparm  = 'DONE'
       endif
endwhile ;* END While_MATH Loop


* End Season Test
* ---------------
endif
* End Season Loop
* ---------------
endwhile
* End FOUND Test
* --------------
endif

* Check next Comparison Experiment Dataset
* ----------------------------------------
 numexp = numexp + 1
  dummy = get_cmpexp (cmpexp,numexp)
    exp = subwrd(dummy,1)
   type = subwrd(dummy,2)

endwhile
 numexp = numexp - 1

'!/bin/mv HISTORY.Tmp HISTORY.T'

* ---------------------------------------------------------
* Now that we have computed plots for each experiment,
* we can compute the Closeness plots to MERRA-2
* ---------------------------------------------------------

* Find MERRA2 experiment
* ----------------------
  MERRA2  = 0
       n  = 1
while( n <= numexp )
if( ctag.n.1 = "MERRA-2" )
    MERRA2 = n
endif
         n = n + 1
endwhile

if( MERRA2 != 0 )

* Loop over Seasons to Process
* ----------------------------
       m = 1
while( m > 0 )
    season = subwrd(seasons,m)
if( season = '' )
         m = -1
else
         m = m+1
         say 'Processing Season: 'season

'set dfile 'qfile.1
'set gxout shaded'
'rgbset'
'run setenv "LEVTYPE" 'DLEVS

* Horizontal Plot
* ---------------
       mathparm  =  MATH
while( mathparm != 'DONE' )

* Closeness Plot (Experiment_vs_Comparison to MERRA-2)
* ----------------------------------------------------
       n  = 1
while( n <= numexp )
if( ctag.n.1 != "NULL" & ctag.n.1 != "merra" & ctag.n.1 != "MERRA-2" )
say 'Closeness plot between  exp: 'qtag.1
say '                       cexp: 'ctag.n.1
say '                        obs: 'ctag.MERRA2.1
say ''
        flag = ""
while ( flag = "" )

'closeness -CVAR 'cmod''n' -MVAR 'qmod' -OVAR 'cmod''MERRA2' -CNAME 'ctag.n.1' -MNAME 'NAME' -ONAME 'ctag.MERRA2.1' -CDESC 'cdesc.n.1' -MDESC 'qdesc.1' -ODESC 'cdesc.MERRA2.1' -MFILE 'qfile.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OFILE 'cfile.MERRA2.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm

if( mathparm != NULL )
    MTH = '_'mathparm
else
    MTH = ''
endif
if( PREFIX != NULL )
    PFX = PREFIX'_'
else
    PFX = ''
endif
'myprint -name 'OUTPUT'/'NAME''MTH'_'ctag.n.1'_closeness_'PFX''ctag.MERRA2.1'.'season

 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile ;* END While_FLAG Loop 
endif
       n  = n + 1
endwhile ;* END While_N Loop 
       if( mathparm != 'NULL' )
           mathparm  = 'NULL'
       else
           mathparm  = 'DONE'
       endif
endwhile ;* END While_MATH Loop

* End Season Test
* ---------------
endif
* ---------------
endwhile ;* END While_m>0 Loop

endif ;* END MERRA-2 Test

if( cmpexp_only = TRUE ) ; return ; endif

***********************************************************************************
*                          Loop over Verification Datasets
***********************************************************************************

'getnumrc 'geosutil'/plots/'NAME
     rcinfo = result
     numrc  = subwrd( rcinfo,1 )
say 'Initial RCINFO: 'rcinfo

if( numrc = 0 )
   'getnumrc 'geosutil'/plots/'DIR
    rcinfo = result
    numrc  = subwrd( rcinfo,1 )
say '  Final RCINFO: 'rcinfo
endif

         k  = 1
while(   k <= numrc )
        loc = k + 1
     rcfile = subwrd( rcinfo,loc )
     RCFILE = rcfile
     if(RC != 'NULL') ; RCFILE = geosutil'/plots/'NAME'/VERIFICATION.'RC'.rc' ; endif
 if( RCFILE =  rcfile )

 say 'k = 'k'  RCFILE = 'RCFILE

* Check for VERIFICATION Formula File
* -----------------------------------
'!remove OBSFORM.txt'
'!echo `basename 'RCFILE' | cut -d. -f2 > OBSFORM.txt`'
'run getenv OBSFORM'
            OBSNAME = result
            OBSFORM = OBSNAME'form.gs'

say 'Checking for OBSFORM: 'geosutil'/plots/'NAME'/'OBSFORM

'!remove CHECKFILE.txt'
'!chckfile 'geosutil'/plots/'NAME'/'OBSFORM
 'run getenv CHECKFILE'
             CHECKFILE  = result

say 'CHECKFILE = 'CHECKFILE
if( CHECKFILE != 'NULL' )
        
         LOBSFORM = 'TRUE'
        'run 'geosutil'/plots/'NAME'/'OBSFORM' -OBS'
         OBS = result
         say 'OBSFORM OBS: 'OBS
        'numargs  'OBS
         numobs = result
         say '     NUMOBS: 'numobs

* Read Verification OBS:OBSGC
* ---------------------------
              n = 1
          OBSNAME.n = subwrd(OBS,n   )
           next = subwrd(OBS,n+1 )
            bit = checkbit(next)
       say 'OBSNAME.'n': 'OBSNAME.n
       while( bit != '' )
              n = n + 1
          OBSNAME.n = subwrd(OBS,n   )
           next = subwrd(OBS,n+1 )
            bit = checkbit(next)
       say 'OBSNAME.'n': 'OBSNAME.n
       endwhile

* Construct OBS GCs from OBS EXPORTS
* ----------------------------------
        m  = 0
        n  = 1
while ( n <= numobs )
        EX = ''
         j = 1
       bit = substr(OBSNAME.n,j,1)
       while(bit != ':' & bit != '')
        EX = EX''bit
         j = j + 1
       bit = substr(OBSNAME.n,j,1)
       endwhile
       if( EX != OBSNAME.n )
         m = m + 1
         j = j + 1
       OBSNAMEGC.m = ''
       bit = substr(OBSNAME.n,j,1)
       while(bit != '')
       OBSNAMEGC.m = OBSNAMEGC.m''bit
         j = j + 1
       bit = substr(OBSNAME.n,j,1)
       endwhile
       OBSNAME.n = EX
       endif
n = n + 1
endwhile

* Get Verification Variables
* --------------------------
FOUND = TRUE
        n  = 1
while ( n <= numobs & FOUND = TRUE )
'run getobs 'OBSNAME.n' 'OBSNAMEGC.n' 'rcfile
        oname.n = subwrd(result,1)
        ofile.n = subwrd(result,2)
        oscal.n = subwrd(result,3)
        odesc.n = subwrd(result,4)
         otag.n = subwrd(result,5)

say 'VERIFICATION_EXPORT_name: 'oname.n
say '             description: 'odesc.n
say '                     tag: ' otag.n
say '                    file: 'ofile.n
say '                 scaling: 'oscal.n
say ''
    if( oname.n = 'NULL' ) ; FOUND = FALSE ; endif
         n  = n + 1
endwhile

* Save Original numn for other RCs
* --------------------------------
  numnorig = numn
  numn     = numobs

else

  LOBSFORM = 'FALSE'
   OBSFORM = 'obsform.gs'

* Get Verification Variables
* --------------------------
FOUND = TRUE
        n  = 1
while ( n <= numn & FOUND = TRUE )
'run getobs 'OBS.n' 'OBSGC.n' 'rcfile
        oname.n = subwrd(result,1)
        ofile.n = subwrd(result,2)
        oscal.n = subwrd(result,3)
        odesc.n = subwrd(result,4)
         otag.n = subwrd(result,5)

say 'VERIFICATION_EXPORT_name: 'oname.n
say '             description: 'odesc.n
say '                     tag: ' otag.n
say '                    file: 'ofile.n
say '                 scaling: 'oscal.n
say ''
    if( oname.n = 'NULL' ) ; FOUND = FALSE ; endif
         n  = n + 1
endwhile

endif

if( FOUND = TRUE )
'setlons'
'setlats'

* Land/Water Masks
* ----------------
if( LAND = 'TRUE' | OCEAN = 'TRUE' )
   'set dfile 'ofile.1
if( LEVEL = "NULL" )
   'set z 1'
else
   'set lev 'LEVEL
endif
   'set t 1'
   'setmask     obs'
   'define lwmaskobs = regrid2( lwmaskobs,0.25,0.25,bs_p1,'lonmin','latmin')'
   'define  omaskobs = maskout( 1, lwmaskobs-0.5 )'
   'define  lmaskobs = maskout( 1, 0.5-lwmaskobs )'
endif


'set dfile 'ofile.1
if( LEVEL = "NULL" )
   'set z 1'
else
   'set lev 'LEVEL
endif
    'getdates'
     begdateo = subwrd(result,1)
     enddateo = subwrd(result,2)

* Perform OBS Formula Calculation
* -------------------------------
if( numn = 1 )
   'define qobs = 'oname.1'.'ofile.1'*'oscal.1
else
    filename  = geosutil'/plots/'NAME'/'OBSFORM
    ioflag    = sublin( read(filename),1 )
    if(ioflag = 0)
       close  = close(filename)
       ostring = ''
       n  = 1
       while ( n <= numn )
          ostring = ostring' 'oname.n'.'ofile.n'*'oscal.n
               n  = n + 1
       endwhile
      'run 'geosutil'/plots/'NAME'/'OBSFORM' 'ostring
    else
       ostring = NAME
       n  = 1
       while ( n <= numn )
          ostring = ostring' 'oname.n'.'ofile.n'*'oscal.n
               n  = n + 1
       endwhile
      'run 'geosutil'/plots/'DIR'/'NAME'.gs 'ostring
      'define qobs = 'NAME'.1'
    endif
endif

* Check for MERRA Resolution problem
* ----------------------------------
'getinfo xdim'
         xdim = result
if( xdim != 540 )

                              'define qobs = regrid2( qobs,0.25,0.25,bs_p1,'lonmin','latmin')'
    if(    LAND  = 'TRUE'   |  OCEAN = 'TRUE' )
       if( LAND  = 'TRUE' ) ; 'define qobs = maskout( qobs,lmaskobs )' ; endif
       if( OCEAN = 'TRUE' ) ; 'define qobs = maskout( qobs,omaskobs )' ; endif
    endif

else

*                             'define qobs = regrid2( qobs,'dlon','dlat',bs_p1,'lonmin','latmin')'
                              'define qobs = regrid2( qobs,0.25,0.25,bs_p1,'lonmin','latmin')'
    if(    LAND  = 'TRUE'   |  OCEAN = 'TRUE' )
       if( LAND  = 'TRUE' ) ; 'define qobs = maskout( qobs,lmaskmod )' ; endif
       if( OCEAN = 'TRUE' ) ; 'define qobs = maskout( qobs,omaskmod )' ; endif
    endif

endif

* Compute Seaonal Means
* ---------------------
'seasonal qobs'

'run getenv "CLIMATE"'
             climate = result


* Perform Taylor Plots
* --------------------
'set dfile 'qfile.1
'run getenv "TAYLOR"'
         taylor = result
if(      taylor = 'true' & TAYLOR = 'TRUE' )

'taylor qmodmdjf qobsdjf djf 'EXPID
'taylor qmodmjja qobsjja jja 'EXPID
'taylor qmodmson qobsson son 'EXPID
'taylor qmodmmam qobsmam mam 'EXPID
'taylor qmodmann qobsann ann 'EXPID

'taylor_write 'EXPID' 'NAME' 'OUTPUT
'taylor_read   GFDL   'NAME' 'verification
'taylor_read   CAM3   'NAME' 'verification
'taylor_read   e0203  'NAME' 'verification
                                                                                                   
"taylor_plt 4 CAM3 GFDL e0203 "EXPID" "OUTPUT" "NAME" '"EXPID" "NAME" vs "obsnam"' "DEBUG
endif


* Loop over Seasons to Process
* ----------------------------
       m = 1
while( m > 0 )
    season = subwrd(seasons,m)
if( season = '' )
         m = -1
else
         m = m+1
         say 'Processing Season: 'season

'set dfile 'qfile.1
'set gxout shaded'
'rgbset'
'run setenv "LEVTYPE" 'DLEVS

* Horizontal Plot
* ---------------
       mathparm  =  MATH
while( mathparm != 'DONE' )

* Standard Plot (Experiment_vs_Verification)
* ------------------------------------------
        flag = ""
while ( flag = "" )
'makplot -MVAR 'qmod' -MNAME 'NAME ' -MFILE 'qfile.1' -MDESC 'qdesc.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OVAR 'qobs' -ONAME 'otag.1' -OFILE 'ofile.1' -ODESC 'odesc.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm
 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile ;* END While_FLAG Loop 

* Closeness Plot (Experiment_vs_Comparison to Verification)
* ---------------------------------------------------------
       n  = 1
while( n <= numexp )
if( ctag.n.1 != "NULL" & ctag.n.1 != "merra" & ctag.n.1 != "MERRA-2" )
say 'Closeness plot between  exp: 'qtag.1
say '                       cexp: 'ctag.n.1
say '                        obs: 'otag.1
say ''
        flag = ""
while ( flag = "" )
'closeness -CVAR 'cmod''n' -MVAR 'qmod' -OVAR 'qobs' -CNAME 'ctag.n.1' -MNAME 'NAME' -ONAME 'otag.1' -CDESC 'cdesc.n.1' -MDESC 'qdesc.1' -ODESC 'odesc.1' -MFILE 'qfile.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OFILE 'ofile.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm

if( mathparm != NULL )
    MTH = '_'mathparm
else
    MTH = ''
endif
if( PREFIX != NULL )
    PFX = PREFIX'_'
else
    PFX = ''
endif
'myprint -name 'OUTPUT'/'NAME''MTH'_'ctag.n.1'_closeness_'PFX''otag.1'.'season

 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile ;* END While_FLAG Loop 

endif
       n  = n + 1
endwhile ;* END While_N Loop 
       if( mathparm != 'NULL' )
           mathparm  = 'NULL'
       else
           mathparm  = 'DONE'
       endif
endwhile ;* END While_MATH Loop

* Zonal Mean Plot
* ---------------
       mathparm  =  MATH
while( mathparm != 'DONE' )
        flag = ""
while ( flag = "" )

'makplotz -MVAR 'qmod' -MNAME 'NAME ' -MFILE 'qfile.1' -MDESC 'qdesc.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OVAR 'qobs' -ONAME 'otag.1' -OFILE 'ofile.1' -ODESC 'odesc.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm' -RGFILE 'rgfile

 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile ;* END While_FLAG Loop 
       if( mathparm != 'NULL' )
           mathparm  = 'NULL'
       else
           mathparm  = 'DONE'
       endif
endwhile ;* END While_MATH Loop


* End Season Test
* ---------------
endif
* End Season Loop
* ---------------
endwhile
* End FOUND Test
* --------------
endif


* End RC=NULL Test
* ----------------
endif
* Update Verification Loop Index
* ------------------------------
k = k + 1

* Restore Original numn for other RCs
* -----------------------------------
if( LOBSFORM = 'TRUE' )
        numn = numnorig
endif

* End Verification Loop
* ---------------------
endwhile



*******************************************************
****            No Verification Case               ****
*******************************************************

if( numrc = 0 )

* Loop over Seasons to Process
* ----------------------------
       m = 1
while( m > 0 )
    season = subwrd(seasons,m)
if( season = '' )
         m = -1
else
         m = m+1
         say 'Processing Season: 'season

'set dfile 'qfile.1
'set gxout shaded'
'rgbset'

* Horizontal Plot
* ---------------
        flag = ""
while ( flag = "" )

'uniplot 'NAME'  'EXPID' 'PREFIX' 'season' 'OUTPUT' 'qfile.1' 'qdesc.1' 'begdate' 'enddate' 'begdateo' 'enddateo' 'climate
 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile

* Zonal Mean Plot
* ---------------
        flag = ""
while ( flag = "" )
'uniplotz 'NAME'  'EXPID' 'PREFIX' 'season' 'OUTPUT' 'qfile.1' 'qdesc.1' 'begdate' 'enddate' 'begdateo' 'enddateo' 'climate
 if( DEBUG = "debug" )
     say "Hit ENTER to repeat plot, or NON-BLANK to continue"
     pull flag
 else
     flag = "next"
 endif
'c'
endwhile


* End Season Test
* ---------------
endif
* End Season Loop
* ---------------
endwhile


* End Test for NUMRC
* ------------------
endif
'quit'
return

* Get Next EXP from CMPEXP List
* -----------------------------
function get_cmpexp (cmpexp,num)
      exp  = subwrd(cmpexp,num)
      len = get_length (exp)
      bit = substr(exp,len-1,1)
      if( bit = ":" )
          type = substr(exp,len,1)
          exp  = substr(exp,1,len-2)
      else
          type = M
      endif
return exp' 'type

function get_length (string)
tb = ""
i = 1
while (i<=256)
blank = substr(string,i,1)
if( blank = tb )
length = i-1
i = 999
else
i = i + 1
endif
endwhile
return length

* To Prevent Problem with BIT: E
* ------------------------------
function checkbit (word)
      bit = substr(word,1,1)
      dum = bit'TEST'
      if( dum = "ETEST" ) ; bit = A ; endif
return bit
