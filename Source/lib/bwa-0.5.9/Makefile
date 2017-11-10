CC=			gcc-6
CXX=		g++-4.8
CFLAGS=		-g -Wall -O2 -fgnu89-inline
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_FILE_OFFSET_BITS=64
OBJS=		utils.o bwt.o bwtio.o bwtaln.o bwtgap.o is.o \
			bntseq.o bwtmisc.o bwtindex.o stdaln.o simple_dp.o \
			bwaseqio.o bwase.o bwape.o kstring.o cs2nt.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o bamlite.o
PROG=		bin/bwa libbwa.a
INCLUDES=	
LIBS=		-lm -lz -lpthread -Lbwt_gen -lbwtgen
SUBDIRS=	. bwt_gen

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

bin/bwa:lib-recur $(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

bwt.o:bwt.h
bwtio.o:bwt.h
bwtaln.o:bwt.h bwtaln.h kseq.h
bwt1away.o:bwt.h bwtaln.h
bwt2fmv.o:bwt.h
bntseq.o:bntseq.h
bwtgap.o:bwtgap.h bwtaln.h bwt.h bwt.c

bwtsw2_core.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h bwt.c
bwtsw2_aux.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h bwt.c
bwtsw2_main.o:bwtsw2.h

libbwa.a:$(OBJS)
	ar -csru $@ $(OBJS)

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal-recur
