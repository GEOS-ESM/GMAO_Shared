#!/bin/bash

user=$(id -un)
host=$(hostname)
segments=$(ipcs -m | awk -v user="$user" '$1 ~ /^0x/ && $3 == user { print $2 }')

if [[ -z $segments ]]
then
   echo "$host: No shared memory segments owned by $user to remove."
else
   echo "$host: Shared memory segments before cleanup:"
   ipcs -m
   echo "$host: Removing segments..."
   for seg in $segments
   do
      # Use the shmid: private segments have key 0x00000000, which cannot be
      # removed by key.
      if ! ipcrm -m "$seg"
      then
         echo "$host: Unable to remove shared memory segment $seg owned by $user"
      fi
   done
fi
