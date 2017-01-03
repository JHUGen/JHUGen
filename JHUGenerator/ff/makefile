FC            = ifort
FFLAGS        = 

LFLAGS	      = $(FFLAGS)

DEST	      = /usr/local/lib

HDRS	      = aa.h \
		ff.h

LIBRARY	      = libff.a

MAKEFILE      = Makefile

AAOBJS	      = aaxbx.o \
		aaxcx.o \
		aaxdx.o \
		aaxex.o \
		aaxinv.o \
		aacbc.o \
		aaccc.o \
		aacinv.o \
		spence.o \
		npoin.o

OBJS	      = ff2dl2.o \
		ffabcd.o \
		ffca0.o \
		ffcb0.o \
		ffcb1.o \
		ffcb2.o \
		ffcb2p.o \
		ffcdb0.o \
		ffcc0.o \
		ffcc0p.o \
		ffcc1.o \
		ffcel2.o \
		ffcel3.o \
		ffcel4.o \
		ffcel5.o \
		ffceta.o \
		ffcli2.o \
		ffcrr.o \
		ffcxr.o \
		ffcxs3.o \
		ffcxs4.o \
		ffcxyz.o \
		ffdcc0.o \
		ffdcxs.o \
		ffdel2.o \
		ffdel3.o \
		ffdel4.o \
		ffdel5.o \
		ffdel6.o \
		ffdl2i.o \
		ffdl5p.o \
		ffdxc0.o \
		ffinit_mine.o \
		ffrcvr.o \
		fftran.o \
		ffxb0.o \
		ffxb1.o \
		ffxb2p.o \
		ffxc0.o \
		ffxc0i.o \
		ffxc0p.o \
		ffxc1.o \
		ffxd0.o \
		ffxd0h.o \
		ffxd0i.o \
		ffxd0p.o \
		ffxd1.o \
		ffxdb0.o \
		ffxdbd.o \
		ffxdi.o \
		ffxdpv.o \
		ffxe0.o \
		ffxe1.o \
		ffxf0.o \
		ffxf0h.o \
		ffxli2.o \
		ffxxyz.o

# older software can be made to run with this dummy routine
#		ffini.o

PRINT	      = pr

AASRCS	      = aaxbx.f \
		aaxcx.f \
		aaxdx.f \
		aaxex.f \
		aaxinv.f \
		aacbc.f \
		aaccc.f \
		aacinv.f \
		spence.f \
		npoin.f

SRCS	      = ff2dl2.f \
		ffabcd.f \
		ffca0.f \
		ffcb0.f \
		ffcb1.f \
		ffcb2.f \
		ffcb2p.f \
		ffcdb0.f \
		ffcc0.f \
		ffcc0p.f \
		ffcc1.f \
		ffcdb0.f \
		ffcel2.f \
		ffcel3.f \
		ffcel4.f \
		ffcel5.f \
		ffceta.f \
		ffcli2.f \
		ffcrr.f \
		ffcxr.f \
		ffcxs3.f \
		ffcxs4.f \
		ffcxyz.f \
		ffdcc0.f \
		ffdcxs.f \
		ffdel2.f \
		ffdel3.f \
		ffdel4.f \
		ffdel5.f \
		ffdel6.f \
		ffdl2i.f \
		ffdl5p.f \
		ffdxc0.f \
		ffinit_mine.f \
		ffrcvr.f \
		fftran.f \
		ffxb0.f \
		ffxb1.f \
		ffxb2p.f \
		ffxc0.f \
		ffxc0i.f \
		ffxc0p.f \
		ffxc1.f \
		ffxd0.f \
		ffxd0h.f \
		ffxd0i.f \
		ffxd0p.f \
		ffxd1.f \
		ffxdb0.f \
		ffxdbd.f \
		ffxdi.f \
		ffxdpv.f \
		ffxe0.f \
		ffxe1.f \
		ffxf0.f \
		ffxf0h.f \
		ffxli2.f \
		ffxxyz.f

all:		$(LIBRARY)

test:		npointes
		npointes > npointes.out
		diff -w npointes.out npointes.sun

npointes:	npointes.o $(LIBRARY) makefile
		$(FC) $(LFLAGS) -o npointes npointes.o $(LIBRARY)

$(LIBRARY):	$(OBJS) $(AAOBJS)
		ar cru $(LIBRARY) $(OBJS) $(AAOBJS)
		ranlib $(LIBRARY)

clean:;		rm -f $(OBJS) $(AAOBJS) $(LIBRARY) npointes.o npointes

depend:;	mkmf -f $(MAKEFILE) LIBRARY=$(LIBRARY) DEST=$(DEST)

extract:;	ar xo $(DEST)/$(LIBRARY)
		@rm -f __.SYMDEF

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		install $(LIBRARY) $(DEST)
		ranlib $(DEST)/$(LIBRARY)

library:	$(LIBRARY)

print:;		$(PRINT) $(HDRS) $(SRCS)

update:		$(DEST)/$(LIBRARY)

$(DEST)/$(LIBRARY): $(SRCS) $(HDRS) $(EXTHDRS)
		@-ar xo $(DEST)/$(LIBRARY)
		@make -f $(MAKEFILE) DEST=$(DEST) install clean
###
aaxbx.o: ff.h aa.h makefile
aaxcx.o: ff.h aa.h makefile
aaxdx.o: ff.h aa.h makefile
aaxex.o: ff.h aa.h makefile
aaxinv.o: ff.h aa.h makefile
ff2dl2.o: ff.h makefile
ffabcd.o: ff.h makefile
ffca0.o: ff.h makefile
ffcb0.o: ff.h makefile
ffcb1.o: ff.h makefile
ffcc0.o: ff.h makefile
ffcc0p.o: ff.h makefile
ffcc1.o: ff.h makefile
ffcdb0.o: ff.h makefile
ffcel2.o: ff.h makefile
ffcel3.o: ff.h makefile
ffcel4.o: ff.h makefile
ffcel5.o: ff.h makefile
ffceta.o: ff.h makefile
ffcli2.o: ff.h makefile
ffcrr.o: ff.h makefile
ffcxr.o: ff.h makefile
ffcxs3.o: ff.h makefile
ffcxs4.o: ff.h makefile
ffcxyz.o: ff.h makefile
ffdcc0.o: ff.h makefile
ffdcxs.o: ff.h makefile
ffdel2.o: ff.h makefile
ffdel3.o: ff.h makefile
ffdel4.o: ff.h makefile
ffdel5.o: ff.h makefile
ffdel6.o: ff.h makefile
ffdl2i.o: ff.h makefile
ffdl5p.o: ff.h makefile
ffdxc0.o: ff.h makefile
ffinit_mine.o: ff.h aa.h makefile
ffrcvr.o: ff.h makefile
fftran.o: ff.h makefile
ffxb0.o: ff.h makefile
ffxb1.o: ff.h makefile
ffxb2p.o: ff.h makefile
ffxc0.o: ff.h makefile
ffxc0i.o: ff.h makefile
ffxc0p.o: ff.h makefile
ffxc1.o: ff.h makefile
ffxd0.o: ff.h makefile
ffxd0h.o: ff.h makefile
ffxd0i.o: ff.h makefile
ffxd0p.o: ff.h makefile
ffxd1.o: ff.h aa.h makefile
ffxdb0.o: ff.h makefile
ffxdbd.o: ff.h makefile
ffxdi.o: ff.h makefile
ffxdpv.o: ff.h makefile
ffxe0.o: ff.h makefile
ffxe1.o: ff.h makefile
ffxf0.o: ff.h makefile
ffxf0h.o: ff.h makefile
ffxli2.o: ff.h makefile
ffxxyz.o: ff.h makefile
npoin.o: ff.h aa.h makefile
npointes.o: ff.h makefile
