FC            = gfortran
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

all:		checkcompiler $(LIBRARY)

checkcompiler:
		echo $(FC) > whichcompiler.txt.new
		cmp whichcompiler.txt.new whichcompiler.txt >& /dev/null && rm whichcompiler.txt.new || mv whichcompiler.txt.new whichcompiler.txt

test:		npointes
		npointes > npointes.out
		diff -w npointes.out npointes.sun

npointes:	npointes.o $(LIBRARY) makefile whichcompiler.txt
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
aaxbx.o: ff.h aa.h makefile whichcompiler.txt
aaxcx.o: ff.h aa.h makefile whichcompiler.txt
aaxdx.o: ff.h aa.h makefile whichcompiler.txt
aaxex.o: ff.h aa.h makefile whichcompiler.txt
aaxinv.o: ff.h aa.h makefile whichcompiler.txt
ff2dl2.o: ff.h makefile whichcompiler.txt
ffabcd.o: ff.h makefile whichcompiler.txt
ffca0.o: ff.h makefile whichcompiler.txt
ffcb0.o: ff.h makefile whichcompiler.txt
ffcb1.o: ff.h makefile whichcompiler.txt
ffcc0.o: ff.h makefile whichcompiler.txt
ffcc0p.o: ff.h makefile whichcompiler.txt
ffcc1.o: ff.h makefile whichcompiler.txt
ffcdb0.o: ff.h makefile whichcompiler.txt
ffcel2.o: ff.h makefile whichcompiler.txt
ffcel3.o: ff.h makefile whichcompiler.txt
ffcel4.o: ff.h makefile whichcompiler.txt
ffcel5.o: ff.h makefile whichcompiler.txt
ffceta.o: ff.h makefile whichcompiler.txt
ffcli2.o: ff.h makefile whichcompiler.txt
ffcrr.o: ff.h makefile whichcompiler.txt
ffcxr.o: ff.h makefile whichcompiler.txt
ffcxs3.o: ff.h makefile whichcompiler.txt
ffcxs4.o: ff.h makefile whichcompiler.txt
ffcxyz.o: ff.h makefile whichcompiler.txt
ffdcc0.o: ff.h makefile whichcompiler.txt
ffdcxs.o: ff.h makefile whichcompiler.txt
ffdel2.o: ff.h makefile whichcompiler.txt
ffdel3.o: ff.h makefile whichcompiler.txt
ffdel4.o: ff.h makefile whichcompiler.txt
ffdel5.o: ff.h makefile whichcompiler.txt
ffdel6.o: ff.h makefile whichcompiler.txt
ffdl2i.o: ff.h makefile whichcompiler.txt
ffdl5p.o: ff.h makefile whichcompiler.txt
ffdxc0.o: ff.h makefile whichcompiler.txt
ffinit_mine.o: ff.h aa.h makefile whichcompiler.txt
ffrcvr.o: ff.h makefile whichcompiler.txt
fftran.o: ff.h makefile whichcompiler.txt
ffxb0.o: ff.h makefile whichcompiler.txt
ffxb1.o: ff.h makefile whichcompiler.txt
ffxb2p.o: ff.h makefile whichcompiler.txt
ffxc0.o: ff.h makefile whichcompiler.txt
ffxc0i.o: ff.h makefile whichcompiler.txt
ffxc0p.o: ff.h makefile whichcompiler.txt
ffxc1.o: ff.h makefile whichcompiler.txt
ffxd0.o: ff.h makefile whichcompiler.txt
ffxd0h.o: ff.h makefile whichcompiler.txt
ffxd0i.o: ff.h makefile whichcompiler.txt
ffxd0p.o: ff.h makefile whichcompiler.txt
ffxd1.o: ff.h aa.h makefile whichcompiler.txt
ffxdb0.o: ff.h makefile whichcompiler.txt
ffxdbd.o: ff.h makefile whichcompiler.txt
ffxdi.o: ff.h makefile whichcompiler.txt
ffxdpv.o: ff.h makefile whichcompiler.txt
ffxe0.o: ff.h makefile whichcompiler.txt
ffxe1.o: ff.h makefile whichcompiler.txt
ffxf0.o: ff.h makefile whichcompiler.txt
ffxf0h.o: ff.h makefile whichcompiler.txt
ffxli2.o: ff.h makefile whichcompiler.txt
ffxxyz.o: ff.h makefile whichcompiler.txt
npoin.o: ff.h aa.h makefile whichcompiler.txt
npointes.o: ff.h makefile whichcompiler.txt
