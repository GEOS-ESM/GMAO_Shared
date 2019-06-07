        subroutine deteps ( eps )
c**********************************************************************
	real eps
	eps = 1.000
10	eps = 0.500 * eps
	if ( 1.000 .lt. (1.000 + eps) ) goto 10
	eps = 2.000 * eps
	return
	end
