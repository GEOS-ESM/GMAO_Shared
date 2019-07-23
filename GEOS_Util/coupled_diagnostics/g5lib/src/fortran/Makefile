SRCS=g2g.f90 varimax.f90

../futils.so: ${SRCS}
	f2py -c -m futils ${SRCS}
	mv futils.so ../../

