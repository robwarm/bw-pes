FC = ifort 
FFLAGS =  -O2 -fp-model strict -mkl=sequential
BLD = xiar -rcv
#MKLPATH = /opt/intel/mkl/lib/intel64/
#LIBDIRS = -L. -L/usr/lib64 -L$(MKLPATH)
#LIBS = 
# BLD = /usr/bin/ar -rcv

objlibinv = inv_dp.o inv_any.o inv_sym2.o inv_sym2t1.o inv_sym.o inv_cg2.o inv_cg.o \
inv_core.o inv_mgx.o inv_mg1.o inv_mg11.o inv_mg2.o inv_mg21.o inv_mg.o \
inv_cxx.o inv_cx1.o inv_cx11.o inv_cx2.o  inv_cx21.o  \
inv_cxy.o inv_cxz.o inv_cx.o \
inv.o

objlibpes = pes0.o \
  pes_x1.o pes_y1.o pes_x2.o pes_x1y1.o  \
  pes_x2y1.o  pes_x2y1consec.o\
  pes1_xyz.o pes1c_xyz.o pes1c.o \
  pes2.o pes.o

objlibpx = px.o

objects = $(objlibinv) $(objlibpes) 
#$(objlibpx)
libraries = libpes.a 
#libpx.a

#-----------------------------------------------------------------------
# Rules.

all : $(objects) $(libraries)
obj : $(objects)
lib : $(libraries)

.PHONY : all obj lib clean realclean
.SUFFIXES : $(nil)
FORCE :

# Pattern rule to make a %.o file from a %.f90 file.
%.o : %.f90
	$(FC) -c $(FFLAGS) $<

# Here are the rules to make the libraries.
libpes.a : $(objlibinv) $(objlibpes)
	\rm -f $@
	$(BLD) $(@F) $(objlibinv) $(objlibpes)
#libpx.a : $(objlibpx)
#	\rm -f $@
#	$(BLD) $(@F) $(objlibpx)

# targets 'clean' and 'realclean' clean up the directory.
clean : FORCE
	rm -f *.o *~ core junk
realclean : clean
	rm -f *.a *.mod
# DO NOT DELETE
