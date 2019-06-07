/* filename getarg.F
   A wrapper for extracting command line argument by index.
   It has the same name and argument as on most Unix machines. This
   file only contains definition on those machines where getarg
   does not exist for Fortran programs.
 */

#if defined(CRAY)
      subroutine getarg(ind, str)
      integer ind, ilen, ierr
      character str*(*)
      call pxfgetarg(ind, str, ilen, ierr)
      if (ilen.gt.len(str)) then
         print *, 'GETARG: given buffer (', len(str), ') is short than',
     &        ' argument ', ind, '(', ilen, ').'
         print *, '        Argument truncated.'
      else if (ierr.ne.0) then
         print *, 'GETARG: ', ind, ' is out of range.'
      endif
      return
c     end CRAY version of gerarg
      end
#endif
