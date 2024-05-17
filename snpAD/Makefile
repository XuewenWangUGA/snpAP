
CPP = $(CXX) -O3 -ggdb -std=c++11 -pthread -DDEBUG
OBJS = Error.o Data.o MLHet_thread.o MLHet_thread_tstv.o MLHetRB_thread.o log_add.o
HDRS = Error.h Data.h MLHet_thread.h MLHet_thread_tstv.h MLHetRB_thread.h log_add.h

OBJSCall = Error.o VCFfac.o VCFfac_GL.o snpADCall.o FastaSlurp.o log_add.o
HDRSCall = Error.h VCFfac.h VCFfac_GL.h FastaSlurp.h log_add.h

all: snpAD snpADCall snpADll snpADci

snpADCall: ${OBJSCall} ${HDRSCall}
	$(CPP) -o snpADCall ${OBJSCall} -lm -lpopt

snpADll: ${OBJS} snpADll.o ${HDRS}
	$(CPP) -o snpADll snpADll.o ${OBJS} -lnlopt -lm -lpthread -Wl,--no-as-needed -lpopt

snpADci: ${OBJS} snpADci.o ${HDRS}
	$(CPP) -o snpADci snpADci.o ${OBJS} -lnlopt -lm -lpthread -Wl,--no-as-needed -lpopt

snpAD: ${OBJS} snpAD.o ${HDRS}
	$(CPP) -o snpAD snpAD.o ${OBJS} -lnlopt -lm -lpthread -Wl,--no-as-needed -lpopt

%.o: %.cc %.h
	$(CPP) -c -o $@ $<

%.o: %.cc 
	$(CPP) -c -o $@ $<

clean: 
	rm -f ${OBJS} snpAD ${OBJSCall} snpADCall snpADll snpAD.o snpADCall.o snpADll.o snpADci.o

