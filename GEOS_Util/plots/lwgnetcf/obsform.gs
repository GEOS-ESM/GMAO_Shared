function formula (args)

'numargs  'args
 numargs = result

       n  = 1
while( n <= numargs )
       oname.n = subwrd(args,n)
       n  = n + 1
endwhile

* Formula:  lwgnetcf = (sfcupc-sfcdwnc)-(sfcup-sfcdwn)
* ----------------------------------------------------
say 'define qobs = 'oname.1' - 'oname.2' - 'oname.3' + 'oname.4
    'define qobs = 'oname.1' - 'oname.2' - 'oname.3' + 'oname.4

return
