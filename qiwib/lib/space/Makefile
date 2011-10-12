SWIG=swig2 -w365,366,367,124 -debug-template
# -debug-tmsearch
# -debug-top 1,2,3,4

# Python include directory is OS dependent
PYTHONINCLUDE= /usr/include/python2.6

CXX=llvm-clang -x c++ -O3 
#CXX=c++ -O3 
#CXX=icpc

CXXLINK=g++ ${CXXFLAGS} 
CXXFLAGS=-Wall -fPIC -DDERIVATIVE_FIVE_POINT

swigtest: test.hh test.i
	${SWIG} -c++ -octave test.i
	mv test_wrap.cxx test_wrap.cc
	mkoctfile -o test.oct test_wrap.cc

cctest: test.o instances.o grid1d.hh grid1d.cc
	${CXXLINK} -o test test.o instances.o -lm

pymodule:  space-python_wrap.cxx
	${CXXLINK} -o _space.so -shared space-python_wrap.cxx instances.o -I${PYTHONINCLUDE}

octmodule: space-octave_wrap.cc
	mkoctfile -o space.oct space-octave_wrap.cc instances.o

test: pymodule octmodule
	./test.py
	./test.m

clean:
	rm -f *~ \#*\# *.o *.oct _*.so test *_wrap.cc *_wrap.cxx octave-core space.py *.pyc


# Intermediate targets
instances.o: grid1d.hh grid1d.cc instances.cc space.h

space-python_wrap.cxx: space-python.i instances.o
	${SWIG} -python -c++ space-python.i

space-octave_wrap.cc: space.i space-octave.i instances.o 
	${SWIG} -octave -c++ space-octave.i
	mv space-octave_wrap.cxx space-octave_wrap.cc

