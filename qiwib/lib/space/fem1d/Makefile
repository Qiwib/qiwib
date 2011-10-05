MKLLAPACK=-Wl,--start-group -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_solver  -Wl,--end-group
CXXFLAGS=-O3 -Wall -I${MKLROOT}/include
LDFLAGS=-L${MKLROOT}/lib/ia32 ${MKLLAPACK}

test: FDInterval.cpp FDTest.cpp
	$(CXX) $(CXXFLAGS) ${LDFLAGS} -o $@ $^

clean:
	rm -f *~ \#*\#

