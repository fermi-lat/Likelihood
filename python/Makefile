#
# @file Makefile
# @brief Hand-coded Makefile for building the Likelihood.py module
#
# @author J. Chiang
#
# $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/python/Makefile,v 1.10 2005/01/19 02:14:59 jchiang Exp $
#

#PYTHONDIR = $(GLAST_EXT)/python/2.3/Linux
#PYINCLUDE = $(PYTHONDIR)/include
#SWIG_DIR = /usr/local
#SWIG_DIR = /scratch/local
#PYINCLUDE = /scratch/local/include/python2.3
SWIG = ${SWIG_DIR}/bin/swig -w401 -w402

CFLAGS = -g -w -fPIC

INCLUDE = -I$(LIKELIHOODROOT)\
	-I$(IRFINTERFACEROOT)\
	-I$(IRFLOADERROOT)\
	-I$(ASTROROOT)\
	-I$(OPTIMIZERSROOT)\
	-I$(DATASUBSELECTORROOT)\
	-I$(EVTBINROOT)\
	-I$(TIPROOT)\
	-I$(MAP_TOOLSROOT)\
	-I$(XMLBASEROOT)\
	-I$(ST_FACILITIESROOT)\
	-I$(FACILITIESROOT)\
	-I$(ST_FACILITIESROOT)\
	-I$(F2CROOT)\
	-I$(GLAST_EXT)/cfitsio/v2470/include\
	-I$(GLAST_EXT)/CLHEP/1.8.0.0/include\
	-I$(GLAST_EXT)/xerces/2.6.0/include

LIBS = -L../$(BINDIR)/ -lLikelihoodLib\
	-L$(IRFLOADERROOT)/$(BINDIR) -lirfLoaderLib\
	-L$(DC1RESPONSEROOT)/$(BINDIR) -ldc1Response\
	-L$(G25RESPONSEROOT)/$(BINDIR) -lg25ResponseLib\
	-L$(TESTRESPONSEROOT)/$(BINDIR) -ltestResponseLib\
	-L$(DEVRESPONSEROOT)/$(BINDIR) -ldevResponseLib\
	-L$(IRFUTILROOT)/$(BINDIR) -lirfUtilLib\
	-L$(IRFINTERFACEROOT)/$(BINDIR) -lirfInterfaceLib\
	-L$(MAP_TOOLSROOT)/$(BINDIR) -lmap_toolsLib\
	-L$(ASTROROOT)/$(BINDIR) -lastro\
	-L$(DATASUBSELECTORROOT)/$(BINDIR) -ldataSubselectorLib\
	-L$(OPTIMIZERSROOT)/$(BINDIR) -loptimizers\
	-L$(EVTBINROOT)/$(BINDIR) -levtbinLib\
	-L$(TIPROOT)/$(BINDIR) -ltipLib\
	-L$(XMLBASEROOT)/$(BINDIR) -lxmlBase\
	-L$(ST_FACILITIESROOT)/$(BINDIR) -lst_facilitiesLib\
	-L$(FACILITIESROOT)/$(BINDIR) -lfacilities\
	-L$(GLAST_EXT)/cfitsio/v2470/lib -lcfitsio\
	-L$(GLAST_EXT)/CLHEP/1.8.0.0/lib -lCLHEP\
	-L$(GLAST_EXT)/xerces/2.6.0/lib/ -lxerces-c -lpthread\
	-L$(GLAST_EXT)/ROOT/v4.02.00/root/lib\
	-lCore -lCint -lTree -lMatrix -lPhysics\
	-lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript\
	-lg2c -L$(SWIG_DIR)/lib -lswigpy -fPIC 

module_base = Likelihood
module_name = _${module_base}.so

${module_name}: ${module_base}_wrap.o
	${CXX} -shared $< ${LIBS} -o ${module_name}

${module_base}_wrap.o: ${module_base}_wrap.cxx
	${CXX} -c ${CFLAGS} ${INCLUDE} -I${PYINCLUDE} $<

${module_base}_wrap.cxx: ${module_base}.i
	${SWIG} -python -c++ ${module_base}.i

${module_base}.i: ${module_base}.in
	python setup.py $<

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
