#
# @file Makefile
# @brief Hand-coded Makefile for building the Likelihood.py module
#
# @author J. Chiang
#
# $Header$
#

#CXX = /usr/local/bin/g++
#SWIG = /u1/local/bin/swig
SWIG = swig
pypath = /usr/local

CFLAGS = -w

INCLUDE = -I$(CMTPATH)/Likelihood/v1r1\
	-I$(CMTPATH)/latResponse/v0r2\
	-I$(CMTPATH)/astro/v0r4\
	-I$(CMTPATH)/optimizers/v0r1p1\
	-I$(CMTPATH)/xml/v4r3p0/\
	-I$(CMTPATH)/facilities/v2r7p2\
	-I$(CMTPATH)/f2c/v2\
	-I$(CMTPATH)/GlastPolicy/v5r2p1\
	-I$(GLAST_EXT)/cfitsio/v2410/include\
	-I$(GLAST_EXT)/CLHEP/1.8.0.0/include\
	-I$(GLAST_EXT)/xerces/1.7.0/include\
	-I$(GLAST_EXT)/xerces/1.7.0/include/xercesc

PYINCLUDE = -I$(pypath)/include/python2.2

LIBS = -L../rh72_gcc2953/ -lLikelihood\
	-L$(CMTPATH)/latResponse/v0r2/rh72_gcc2953 -llatResponse\
	-L$(CMTPATH)/astro/v0r4/rh72_gcc2953 -lastro\
	-L$(CMTPATH)/optimizers/v0r1p1/rh72_gcc2953 -loptimizers\
	-L$(CMTPATH)/xml/v4r3p0/rh72_gcc2953 -lxml\
	-L$(CMTPATH)/facilities/v2r7p2/rh72_gcc2953 -lfacilities\
	-L$(GLAST_EXT)/cfitsio/v2410/lib -lcfitsio\
	-L$(GLAST_EXT)/CLHEP/1.8.0.0/lib -lCLHEP\
	-L$(GLAST_EXT)/xerces/1.7.0/lib/ -lxerces-c1_7_0 -lpthread\
	-lg2c -lswigpy

module_base = Likelihood
module_name = _${module_base}.so

${module_name}: ${module_base}_wrap.o
	${CXX} -shared $< ${LIBS} -o ${module_name}

${module_base}_wrap.o: ${module_base}_wrap.cxx
	${CXX} -c ${CFLAGS} ${INCLUDE} ${PYINCLUDE} -o $@ $<

${module_base}_wrap.cxx: ${module_base}.i
	${SWIG} -c -python -c++ ${module_base}.i

anew:
	make clean
	make

tar_file:
	tar cvf ${MY_PROG}.tar ${SRCS} ${HFILES} ${OTHER_FILES} ${srcs}

clean:
	/bin/rm -f *~
	/bin/rm -f *.o
	/bin/rm -f *.trace
	/bin/rm -f *wrap*
	/bin/rm -f *.so
	/bin/rm -f ${module_base}.py*

very_clean:
	make clean
	/bin/rm -f core
	/bin/rm -f _${module_base}.so
