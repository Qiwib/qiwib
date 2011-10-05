SWIG=swig2
#CXX=llvm-clang -x c++ -O3 
CXXLINK=g++ ${CXXFLAGS} 
CXXFLAGS=-Wall -fPIC
CXX=c++ -O3 
#CXX=icpc

cctest: test.o instances.o grid1d/grid1d.hh grid1d/grid1d.cc
	${CXXLINK} -o test test.o instances.o -lm

pymodule:  space_wrap.cxx
	${CXXLINK} -o _space.so -shared space_wrap.cxx instances.o -I /usr/include/python2.6

octmodule: space-octave_wrap.cc
	mkoctfile -o space.oct space-octave_wrap.cc instances.o

test: pymodule octmodule
	./test.py
	./test.m

clean:
	rm -f *~ \#*\# *.o


# Intermediate targets
instances.o: grid1d/grid1d.hh grid1d/grid1d.cc instances.cc space.h

space_wrap.cxx: space.i instances.o
	${SWIG} -python -c++ space.i

space-octave_wrap.cc: space.i space-octave.i instances.o 
	${SWIG} -octave -c++ space-octave.i
	mv space-octave_wrap.cxx space-octave_wrap.cc

