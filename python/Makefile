#
# @file Makefile
# @brief Hand-coded Makefile for building the Likelihood.py module
#
# @author J. Chiang
#
# $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/python/Makefile,v 1.6 2004/08/11 17:23:59 jchiang Exp $
#

#PYTHONDIR = $(GLAST_EXT)/python/2.3/Linux
#PYINCLUDE = -I$(PYTHONDIR)/include
#SWIG_DIR = /usr/local
#SWIG_DIR = /scratch/local
PYINCLUDE = -I/scratch/local/include/python2.3
SWIG = ${SWIG_DIR}/bin/swig -w401 -w402

CFLAGS = -g -w -fPIC

INCLUDE = -I$(LIKELIHOODROOT)\
	-I$(IRFINTERFACEROOT)\
	-I$(IRFLOADERROOT)\
	-I$(ASTROROOT)\
	-I$(OPTIMIZERSROOT)\
	-I$(EVTBINROOT)\
	-I$(TIPROOT)\
	-I$(MAP_TOOLSROOT)\
	-I$(XMLROOT)\
	-I$(ST_FACILITIESROOT)\
	-I$(FACILITIESROOT)\
	-I$(F2CROOT)\
	-I$(GLAST_EXT)/cfitsio/v2470/include\
	-I$(GLAST_EXT)/CLHEP/1.8.0.0/include\
	-I$(GLAST_EXT)/xerces/2.6.0/include

LIBS = -L../$(BINDIR)/ -lLikelihoodLib\
	-L$(IRFLOADERROOT)/$(BINDIR) -lirfLoaderLib\
	-L$(DC1RESPONSEROOT)/$(BINDIR) -ldc1ResponseLib\
	-L$(G25RESPONSEROOT)/$(BINDIR) -lg25ResponseLib\
	-L$(TESTRESPONSEROOT)/$(BINDIR) -ltestResponseLib\
	-L$(IRFUTILROOT)/$(BINDIR) -lirfUtilLib\
	-L$(IRFINTERFACEROOT)/$(BINDIR) -lirfInterfaceLib\
	-L$(MAP_TOOLSROOT)/$(BINDIR) -lmap_toolsLib\
	-L$(ASTROROOT)/$(BINDIR) -lastro\
	-L$(OPTIMIZERSROOT)/$(BINDIR) -loptimizersLib\
	-L$(EVTBINROOT)/$(BINDIR) -levtbinLib\
	-L$(TIPROOT)/$(BINDIR) -ltipLib\
	-L$(XMLROOT)/$(BINDIR) -lxml\
	-L$(ST_FACILITIESROOT)/$(BINDIR) -lst_facilitiesLib\
	-L$(FACILITIESROOT)/$(BINDIR) -lfacilities\
	-L$(GLAST_EXT)/cfitsio/v2470/lib -lcfitsio\
	-L$(GLAST_EXT)/CLHEP/1.8.0.0/lib -lCLHEP\
	-L$(GLAST_EXT)/xerces/2.6.0/lib/ -lxerces-c -lpthread\
	-L$(GLAST_EXT)/ROOT/v3.10.02/root/lib\
	-lCore -lCint -lTree -lMatrix -lPhysics\
	-lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript\
	-lg2c -L$(SWIG_DIR)/lib -lswigpy -fPIC 

module_base = Likelihood
module_name = _${module_base}.so

${module_name}: ${module_base}_wrap.o
	${CXX} -shared $< ${LIBS} -o ${module_name}

${module_base}_wrap.o: ${module_base}_wrap.cxx
	${CXX} -c ${CFLAGS} ${INCLUDE} ${PYINCLUDE} $<

${module_base}_wrap.cxx: ${module_base}.i
	${SWIG} -runtime -python -c++ ${module_base}.i

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
