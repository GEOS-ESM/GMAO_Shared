function pause (args)
'numargs  'args
 numargs = result
if( numargs > 0 )
 string  = subwrd(args,1)
        num  = 2
while ( num <= numargs )
        string = string' 'subwrd(args,num)
        num = num + 1
endwhile
say string
endif
say 'Hit n for Next, or c to Clear  ...'
pull flag
if( flag = 'c' )
'c'
endif
return
