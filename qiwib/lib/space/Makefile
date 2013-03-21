SWIG=swig204 -w365,366,367,368,124,362,503 
# -debug-template
# -debug-tmsearch
# -debug-top 1,2,3,4

# Python include directory is OS typedefdependent
PYTHONINCLUDE= /usr/include/python2.7

#CXX=clang -x c++ -O2
CXX=c++
#CXX=icpc

CXXLINK=g++ ${CXXFLAGS} 
CXXFLAGS=-Wall -fPIC -DDERIVATIVE_FIVE_POINT -O2 -march=native -std=c++0x -Wno-sign-compare

swigtest: test.h test.i
	${SWIG} -c++ -python test.i
	mv test_wrap.cxx test_wrap.cc
	${CXXLINK} -o _test.so -shared test_wrap.cc instances.o -I${PYTHONINCLUDE}

cctest: test.o instances.o grid1d.hh grid1d.cc
	${CXXLINK} -o test test.o instances.o -lm

cctest2: test2.o instances.o grid1d.hh grid1d.cc
	${CXXLINK} -o test2 test2.o instances.o -lm

pymodule:  space-python_wrap.cxx
	${CXXLINK} -o _space.so -shared space-python_wrap.cxx instances.o -I${PYTHONINCLUDE}

octmodule: space-octave_wrap.cc
	mkoctfile -o space.oct space-octave_wrap.cc instances.o

test: pymodule octmodule
	./test-python.py
	./test-octave.m

clean:
	rm -f *~ \#*\# *.o *.oct _*.so test *_wrap.cc *_wrap.cxx octave-core space.py *.pyc


# Intermediate targets
instances.o: grid1d.hh grid1d.cc instances.cc space.h

space-python_wrap.cxx: space-python.i instances.o
	${SWIG} -python -c++ space-python.i

space-octave_wrap.cc: space.i space-octave.i instances.o 
	${SWIG} -octave -c++ space-octave.i
	mv space-octave_wrap.cxx space-octave_wrap.cc

