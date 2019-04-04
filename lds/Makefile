# INCLUDE = -I/usr/include 
# CMEX = /local/lib/matlab51/bin/cmex  -g ${INCLUDE} ${LIBRARY}

kalman: mxutils.c kalmansmooth.c
	cmex kalmansmooth.c mxutils.c ludcmp.c lubksb.c nrutil.c

