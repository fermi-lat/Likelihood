#
# @file Makefile
# @brief Hand-coded Makefile for building the Likelihood.py module
#
# @author J. Chiang
#
# $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/python/Makefile,v 1.4 2004/05/17 23:27:52 jchiang Exp $
#

SWIG = /usr/local/bin/swig -w401 -w402
#SWIG_DIR = /scratch/local
#PYTHONDIR = $(GLAST_EXT)/python/2.3/Linux
#PYINCLUDE = -I$(PYTHONDIR)/include
SWIG_DIR = /usr/local
PYINCLUDE = -I/usr/local/include/python2.3

CFLAGS = -w -g -fPIC

INCLUDE = -I$(LIKELIHOODROOT)\
	-I$(LATRESPONSEROOT)\
	-I$(ASTROROOT)\
	-I$(OPTIMIZERSROOT)\
	-I$(MAP_TOOLSROOT)\
	-I$(XMLROOT)\
	-I$(FACILITIESROOT)\
	-I$(F2CROOT)\
	-I$(CFITSIO_DIR)/include\
	-I$(CLHEP_DIR)/1.8.0.0/include\
	-I$(GLAST_EXT)/xerces/1.7.0/include\
	-I$(GLAST_EXT)/xerces/1.7.0/include/xercesc

LIBS = -L../$(BINDIR)/ -lLikelihoodLib\
	-L$(LATRESPONSEROOT)/$(BINDIR) -llatResponseLib\
	-L$(MAP_TOOLSROOT)/$(BINDIR) -lmap_toolsLib\
	-L$(ASTROROOT)/$(BINDIR) -lastro\
	-L$(OPTIMIZERSROOT)/$(BINDIR) -loptimizersLib\
	-L$(TIPROOT)/$(BINDIR) -ltipLib\
	-L$(XMLROOT)/$(BINDIR) -lxml\
	-L$(FACILITIESROOT)/$(BINDIR) -lfacilities\
	-L$(CFITSIO_DIR)/lib -lcfitsio\
	-L$(CLHEP_DIR)/1.8.0.0/lib -lCLHEP\
	-L$(GLAST_EXT)/xerces/1.7.0/lib/ -lxerces-c1_7_0 -lpthread\
	-L$(GLAST_EXT)/ROOT/v3.10.02/root/lib\
	-lCore -lCint -lTree -lMatrix -lPhysics\
	-lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript\
	-lg2c -L$(SWIG_DIR)/lib -lswigpy -fPIC -ldl

module_base = Likelihood
module_name = _${module_base}.so

${module_name}: ${module_base}_wrap.o
	${CXX} -shared $< ${LIBS} -o ${module_name}

${module_base}_wrap.o: ${module_base}_wrap.cxx
	${CXX} -c ${CFLAGS} ${INCLUDE} ${PYINCLUDE} -o $@ $<

${module_base}_wrap.cxx: ${module_base}.i
	${SWIG} -noruntime -python -c++ ${module_base}.i

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
	/bin/rm -f ${module_base}_wrap.cxx

very_clean:
	make clean
	/bin/rm -f core
	/bin/rm -f _${module_base}.so
