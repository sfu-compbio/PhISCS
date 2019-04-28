CXX      ?= g++
CXXFLAGS  = -O3
LDFLAGS   = 

all: clean-exe csp_maxsat clean

csp_maxsat: csp_maxsat.o
	$(CXX) csp_maxsat.o -o $@ ${LDFLAGS}

clean:
	@rm -f *.o
	
clean-exe:
	@rm -f csp_maxsat