all:
	touch Doxyfile
	cat ${GLASTPOLICYROOT}/doc/Doxyfile Doxyfile > Doxyfile.tmp
	doxygen Doxyfile.tmp
	rm -f Doxyfile.tmp
