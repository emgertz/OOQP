OOQP=../../..
OOQPINCLUDEDIR=$(OOQP)/include
OOQPLIBDIR=$(OOQP)/lib

CXX      = c++

CXXFLAGS =-O
CPPFLAGS = -I$(OOQPINCLUDEDIR)
LDFLAGS  =-L$(OOQPLIBDIR)
BLAS     = $(OOQP)/extras/blas/libblas.a
MA27LIB  = $(OOQP)/extras/MA27/libMA27.a
FLIBS    = -lg2c
LIBS     = -looqpgondzio -looqpgensparse -looqpsparse -looqpbase \
		$(BLAS) $(MA27LIB) $(FLIBS)

qpgen-sparse-random.exe : QpGenSparseRandomDriver.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o

veryclean: clean
	rm -f *.exe

distclean: veryclean

.C.o:
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
