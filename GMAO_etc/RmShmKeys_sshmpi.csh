#!/bin/csh -f

setenv MYNAME RmShmKeys_sshmpi

                                    setenv FAILED 0
if( (! $?FVROOT ) & (! $?GEOSBIN) ) setenv FAILED 1

if ( $FAILED ) then
   env
   echo " ${MYNAME}: not all required env vars defined"
   exit 1
endif

if( $?FVROOT ) then
   set pathname = $FVROOT/bin
endif
if( $?GEOSBIN ) then
   set pathname = $GEOSBIN
endif

setenv SITE `$pathname/g5_modules site`

if( $?PBS_NODEFILE ) then
   sleep 10

   set nodes = `cat $PBS_NODEFILE | uniq`

   if ( $SITE == NCCS ) then

      echo "`hostname`: Found site ${SITE}: using sshmpi"
      foreach node ($nodes)
         echo sshmpi $node $pathname/rmshmkeyhere.sh
              sshmpi $node $pathname/rmshmkeyhere.sh &
      end

   else if ( $SITE == NAS ) then

      echo "`hostname`: Found site ${SITE}: using pbsdsh"
      set node_indices = `awk 'seen[$0] == 0 { print NR - 1 } { seen[$0]++ }' $PBS_NODEFILE`
      foreach node_index ($node_indices)
         pbsdsh -n $node_index -- $pathname/rmshmkeyhere.sh &
      end
      wait

   else

      echo "SITE: $SITE not supported for SHMEM! Contact GEOS Support"
      exit 2

   endif

   wait
endif
