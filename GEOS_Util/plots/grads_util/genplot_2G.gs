function genplot_2G (args)

*******************************************************
****                 INPUT Variables               ****
*******************************************************

*****************************************************************************************************************************************************
**** Note:  Updates for genplot                                                                                                                  ****
****        -------------------                                                                                                                  ****
****                                                                                                                                             ****
**** Original GOCART   required HISTORY.rc to be written as  'EXPORT001' , 'GC' , for each output BIN001.                                        ****
****                                                         'EXPORT002' , 'GC' , for each output BIN002.                                        ****
****                                                         'EXPORT003' , 'GC' , for each output BIN003.                                        ****
****                                                                                                                                             ****
****          GOCART_2G allowed HISTORY.rc to be written as  'EXPORT' , 'GC'  which would automatically yield EXPORT001, EXPORT002, EXPORT003    ****
****         or, alternatively, HISTORY.rc to be written as  'EXPORT' , 'GC' , 'ALIAS1;ALIAS2,ALIAS3' which would yield ALIAS1, ALIAS2, ALIAS3   ****
****                                                                                                                                             ****
**** QUICKPLOT needs to be able to compare experiments derived from:  GOCART vs GOCART, GOCART_2G vs GOCART_2G, and GOCART vs GOCART_2G          ****
**** To accomplish this, genplot was modified to:                                                                                                ****
****                                                                                                                                             ****
****                             1) First check for  EXPORT:GC         as the OLD EXPORT and GC names                                            ****
****                             2) Then  check for  EXPORT:GC:SUFFIX  as the NEW EXPORT and GC names  if OLD names are not found                ****
****                                                                                                                                             ****
****                                The SUFFIX is then used as an index to determine which ALIAS (if provided) to use, or ...                    ****
****                                be concatentated with EXPORT to create the actual BIN name written to the output file.                       ****
****                                                                                                                                             ****
*****************************************************************************************************************************************************


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

       mold = 0
       mnew = 0
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
              m = mold
if( subwrd(args,num) = '-EXPORT' )
              m = m + 1
   OLD_EXPORT.m = subwrd(args,num+m   )
   NEW_EXPORT.m = NULL
           word = subwrd(args,num+m+1 )
            bit = checkbit(word)
       while( bit != '-' )
              m = m + 1
   OLD_EXPORT.m = subwrd(args,num+m   )
   NEW_EXPORT.m = NULL
           word = subwrd(args,num+m+1 )
            bit = checkbit(word)
       endwhile

*    Construct Model GCs from Model EXPORTS
*    --------------------------------------
     numGCs = m
             m  = 0
             k  = 1
     while ( k <= numGCs )
             EX = ''
              j = 1
            bit = substr(OLD_EXPORT.k,j,1)
            while(bit != ':' & bit != '')
             EX = EX''bit
              j = j + 1
            bit = substr(OLD_EXPORT.k,j,1)
            endwhile
            if( EX != OLD_EXPORT.k )
                  m = m + 1
                  j = j + 1
                OLD_GC.m = ''
                bit = substr(OLD_EXPORT.k,j,1)
                while(bit != '')
                OLD_GC.m = OLD_GC.m''bit
                  j = j + 1
                bit = substr(OLD_EXPORT.k,j,1)
                endwhile
                OLD_EXPORT.k = EX
                    NEW_GC.k = NULL
            endif
     k = k + 1
     endwhile
endif

* Read Alternate Model EXPORT:GC[:SUFFIX]
* ---------------------------------------
              m = mnew
if( subwrd(args,num) = '-ALT_EXPORT' )
     say ' '
     say 'Reading ALT_EXPORTs'
     say ' '
              m = m + 1
   ALT_EXPORT.m = subwrd(args,num+m   )
           word = subwrd(args,num+m+1 )
            bit = checkbit(word)
     say '  ALT_EXPORT.'m'  'ALT_EXPORT.m

       while( bit != '-' )
              m = m + 1
   ALT_EXPORT.m = subwrd(args,num+m   )
           word = subwrd(args,num+m+1 )
            bit = checkbit(word)
     say '  ALT_EXPORT.'m'  'ALT_EXPORT.m
       endwhile
       numGCs = m

              k  = 1
              j  = 2
       while( k <= numGCs )

      '!remove NEW_EXPORT.txt'
      '! echo 'ALT_EXPORT.k' | cut -d: -f1 > NEW_EXPORT.txt'
      'run getenv NEW_EXPORT'
                  NEW_EXPORT.k = result

      '!remove NEW_GC.txt'
      '! echo 'ALT_EXPORT.k' | cut -d: -f2 > NEW_GC.txt'
      'run getenv NEW_GC'
                  NEW_GC.k = result

      '!remove SUFFIX.txt'
      '! echo 'ALT_EXPORT.k' | grep : | cut -d: -f3 > SUFFIX.txt'
      'run getenv SUFFIX'
                  SUFFIX.k = result
      
              say '  NEW_EXPORT.'k': 'NEW_EXPORT.k
              say '      NEW_GC.'k': 'NEW_GC.k
              say '      SUFFIX.'k': 'SUFFIX.k

              if( SUFFIX.k = NULL )
                  SUFFIX.k = ''
              say 'Reset SUFFIX.'k': 'SUFFIX.k
              endif

                    k  = k + 1
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


* Construct OBS GCs from OBS EXPORTS
* ----------------------------------
    numn = n
if( numn = 0 )

*   Set OBS and OBSGC based on Original (OLD) Model EXPORTS and GC
*   --------------------------------------------------------------
    numn = numGCs
           n  = 1
   while ( n <= numn )
       OBS.n = OLD_EXPORT.n
     OBSGC.n =     OLD_GC.n
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

**************************************************
**************************************************

* Get Model Variables
* -------------------
        m  = 1
while ( m <= numGCs )
   say 'run 1st call to getvar_2G for SOURCE, orig EXPORT: 'OLD_EXPORT.m'  orig GC: 'OLD_GC.m
   'run getvar_2G 'OLD_EXPORT.m' 'OLD_GC.m
        qname.m = subwrd(result,1)
        qfile.m = subwrd(result,2)
        qscal.m = subwrd(result,3)
        qdesc.m = subwrd(result,4)
         qtag.m = subwrd(result,5)
       EXPORT.m = OLD_EXPORT.m
           GC.m = OLD_GC.m
    say 'qname from first getvar call = 'qname.m

*   If OLD_EXPORT is not found, try searching for NEW_EXPORT
*   --------------------------------------------------------
    if( qname.m = 'NULL' )
       say ' '
      if( SUFFIX.m = '' )
       say 'run 2nd call to getvar_2G for SOURCE,  alt EXPORT: 'NEW_EXPORT.m'  alt_GC: 'NEW_GC.m
       'run getvar_2G 'NEW_EXPORT.m' 'NEW_GC.m
      else
       say 'run 2nd call to getvar_2G for SOURCE,  alt EXPORT:SUFFIX 'NEW_EXPORT.m':'SUFFIX.m'  alt_GC: 'NEW_GC.m
       'run getvar_2G 'NEW_EXPORT.m':'SUFFIX.m' 'NEW_GC.m
      endif
        qname.m = subwrd(result,1)
        qfile.m = subwrd(result,2)
        qscal.m = subwrd(result,3)
        qdesc.m = subwrd(result,4)
         qtag.m = subwrd(result,5)
    say 'qname from second getvar call = 'qname.m
      if( SUFFIX.m = '' )
       EXPORT.m = NEW_EXPORT.m
      else
       EXPORT.m = NEW_EXPORT.m':'SUFFIX.m
      endif
           GC.m = NEW_GC.m
         if( qname.m = 'NULL' )
             return 
         endif
    endif
m  = m + 1
endwhile


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
'q dims'
say 'Running regrid on qmod, DIMS; '
say result
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
     exp.numexp = exp
    type.numexp = type

while( exp != 'NULL' )
say ' '
say 'Comparing with: 'exp

* analysis = false  EXP=M CMP=M  => ALEVS
* analysis = false  EXP=M CMP=A  => DLEVS
* analysis = true   EXP=A CMP=A  => ALEVS
* analysis = true   EXP=A CMP=M  => DLEVS

* INPUT Experiment is an Analysis
*********************************
if( analysis != "false" )
    if( type = A | type = V )
*   CMP Experiment is an Analysis
       'run setenv "LEVTYPE" 'ALEVS
       'run setenv "DIFFTYPE" 'A

    else
*   CMP Experiment is an Model
       'run setenv "LEVTYPE" 'DLEVS
       'run setenv "DIFFTYPE" 'D
    endif
else

* INPUT Experiment is a Model
*********************************
    if( type = A )
*   CMP Experiment is an Analysis
       'run setenv "LEVTYPE" 'DLEVS
       'run setenv "DIFFTYPE" 'D

    else
*   CMP Experiment is an Model
       'run setenv "LEVTYPE" 'ALEVS
       'run setenv "DIFFTYPE" 'A
    endif
endif

* ------------------------------------------------------

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

   say 'run 1st call to getvar_2G for CMPEXP, orig EXPORT: 'OLD_EXPORT.m'  orig GC: 'OLD_GC.m
   'run getvar_2G 'OLD_EXPORT.m' 'OLD_GC.m' 'exp
        cname.numexp.m = subwrd(result,1)
        cfile.numexp.m = subwrd(result,2)
        cscal.numexp.m = subwrd(result,3)
        cdesc.numexp.m = subwrd(result,4)
         ctag.numexp.m = subwrd(result,5)
say ''
*   If OLD_EXPORT is not found, try searching for NEW_EXPORT
*   --------------------------------------------------------
    if( cname.numexp.m = 'NULL' )
      if( SUFFIX.m = '' )
          say 'run 2nd call to getvar_2G for CMPEXP,  alt EXPORT: 'NEW_EXPORT.m'  alt_GC: 'NEW_GC.m
              'run getvar_2G 'NEW_EXPORT.m' 'NEW_GC.m' 'exp
      else
          say 'run 2nd call to getvar_2G for CMPEXP,  alt EXPORT:SUFFIX 'NEW_EXPORT.m':'SUFFIX.m'  alt_GC: 'NEW_GC.m
              'run getvar_2G 'NEW_EXPORT.m':'SUFFIX.m' 'NEW_GC.m' 'exp
      endif
         cname.numexp.m = subwrd(result,1)
         cfile.numexp.m = subwrd(result,2)
         cscal.numexp.m = subwrd(result,3)
         cdesc.numexp.m = subwrd(result,4)
          ctag.numexp.m = subwrd(result,5)
say ''
          if( cname.numexp.m = 'NULL' )
              FOUND = FALSE
          endif
    endif
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

'makplot -MVAR 'qmod' -MNAME 'NAME ' -MFILE 'qfile.1' -MDESC 'qdesc.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OVAR 'cmod''numexp' -ONAME 'ctag.numexp.1' -OFILE 'cfile.numexp.1' -ODESC 'cdesc.numexp.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm' -QNAME 'qname.1

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
   exp.numexp = exp
  type.numexp = type

endwhile
 numexp = numexp - 1

'!/bin/mv HISTORY.Tmp HISTORY.T'

* ---------------------------------------------------------------------------
* Now that we have computed plots for each experiment,
* we can compute the Closeness plots to MERRA-2 and any CMPEXP ending with :V
* ---------------------------------------------------------------------------

       k  = 1
while( k <= numexp )
say 'Looping through experiments, k = 'k' CTAG = 'ctag.k.1' TYPE = 'type.k

if( ( ctag.k.1 = "MERRA-2" | type.k = V ) & cname.k.1 != 'NULL' )
     TAG   = k
     say 'Performing Closeness plots to: 'ctag.TAG.1' k = 'k

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
'run getenv "DIFFTYPE"'
             DIFFTYPE = result
'run setenv "LEVTYPE" 'DIFFTYPE'LEVS'

* Horizontal Closeness Plot (Experiment_vs_Comparison to Verification)
* --------------------------------------------------------------------
       mathparm  =  MATH
while( mathparm != 'DONE' )

       n  = 1
while( n <= numexp )

     say 'n = 'n' Testing 'qtag.1' and 'ctag.n.1' for closeness with 'ctag.TAG.1

     if( ctag.n.1 != "merra" & ctag.n.1 != "MERRA-2" & ctag.n.1 != ctag.TAG.1 & type.n != V & cname.n.1 != 'NULL' )
     say 'Closeness plot between  exp: 'qtag.1
     say '                       cexp: 'ctag.n.1
     say '                        obs: 'ctag.TAG.1
     say '              Total  numexp: 'numexp
     say ''

             flag = ""
     while ( flag = "" )
     
                  'define zobs'TAG''season' = regrid2( cmod'TAG''season',0.25,0.25,bs_p1,'lonmin','latmin' )'
                  'define zobs'n''season'   = regrid2( cmod'n''season'  ,0.25,0.25,bs_p1,'lonmin','latmin' )'
                  'define zmod'season'      = regrid2( qmod'season'     ,0.25,0.25,bs_p1,'lonmin','latmin' )'

     'closeness -CVAR 'zobs''n' -MVAR 'zmod' -OVAR 'zobs''TAG' -CNAME 'ctag.n.1' -MNAME 'NAME' -ONAME 'ctag.TAG.1' -CDESC 'cdesc.n.1' -MDESC 'qdesc.1' -ODESC 'cdesc.TAG.1' -MFILE 'qfile.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OFILE 'cfile.TAG.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm'  -QNAME 'qname.1

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
     
     'myprint -name 'OUTPUT'/'NAME''MTH'_'ctag.n.1'_closeness_'PFX''ctag.TAG.1'.'season

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
endwhile ;* END While n<=numexp Loop 

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

endif ;* END CTAG Test

k = k + 1
endwhile ;* END While_k Loop

if( cmpexp_only = TRUE ) ; return ; endif

***********************************************************************************
*                      Loop over Observational Verification Datasets
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
'makplot -MVAR 'qmod' -MNAME 'NAME ' -MFILE 'qfile.1' -MDESC 'qdesc.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OVAR 'qobs' -ONAME 'otag.1' -OFILE 'ofile.1' -ODESC 'odesc.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm' -QNAME 'qname.1
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
if( ctag.n.1 != "NULL" & ctag.n.1 != "merra" & ctag.n.1 != "MERRA-2" & type.n != V )
say 'Closeness plot between  exp: 'qtag.1
say '                       cexp: 'ctag.n.1
say '                        obs: 'otag.1
say ''
        flag = ""
while ( flag = "" )
'closeness -CVAR 'cmod''n' -MVAR 'qmod' -OVAR 'qobs' -CNAME 'ctag.n.1' -MNAME 'NAME' -ONAME 'otag.1' -CDESC 'cdesc.n.1' -MDESC 'qdesc.1' -ODESC 'odesc.1' -MFILE 'qfile.1' -MBEGDATE 'begdate' -MENDDATE 'enddate' -OFILE 'ofile.1' -OBEGDATE 'begdateo' -OENDDATE 'enddateo' -EXPID 'EXPID' -PREFIX 'PREFIX' -SEASON 'season' -OUTPUT 'OUTPUT' -CLIMATE 'climate' -GC 'GC.1' -MATH 'mathparm' -QNAME 'qname.1

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
