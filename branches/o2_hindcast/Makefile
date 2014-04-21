OUTNAME=offtrac

#CC=gcc


#CC=gcc
#CFLAGS = -O3 -g -fopenmp -lm -lpthread --fast-math -march="athlon64" -pipe -static

# gaggle with gcc
CC = gcc
CFLAGS = -mfpmath=sse -Ofast -flto -march=native -funroll-loops -fopenmp  -g
LDFLAGS = -I/cm/shared/apps/netcdf/gcc/64/4.3.0/include -L/cm/shared/apps/netcdf/gcc/64/4.3.0/lib
# waddle with gcc
#CC=gcc
#CFLAGS = -O3 -g -fopenmp -lm -lpthread --fast-math -march="athlon64" -pipe -static
#LDFLAGS=-L/usr/local/netcdf4-gfort/lib -I/usr/local/netcdf4-gfort/include -L/usr/local/hdf5/lib -I/usr/local/hdf5/include -I/usr/local/zlib/include -L/usr/local/zlib/lib  -lhdf5 
#LDFLAGS = -L/ltraid3/ashao/waddle/local/lib -I/ltraid3/ashao/waddle/local/include

#CFLAGS = -g 
#debug for gcc
#CFLAGS = -fbounds-check -g

#CFLAGS = -ffast-math -O2 -march="athlon64" -pipe
#CFLAGS = -ffast-math -O2 -pipe

#ocean
#CFLAGS = -O1 -pipe  -pg -g
#CFLAGS = -ffast-math -O3 -march="athlon64" -pipe -g

#yucatan with icc:
#CFLAGS = -O1 -march="i686"  -pipe
#yucatan gcc debugging
#CFLAGS = -O1 -march="i686"  -pipe  -pg -Wuninitialized
#CFLAGS = -g -march="i686"  -pipe -fbounds-check -Wall
#yucatan icc debugging with idb
#CFLAGS = -march="i686"  -pipe -g -lm

#CDFFLAGS = -lm -I/usr/local/include -L/usr/local/lib -E -source_listing
#CDFFLAGS = -g -lm -lnetcdf -I/usr/local/include -L/usr/local/lib 

#CDFFLAGS = -lm -lnetcdf -I/usr/local/include -L/usr/local/lib 
#CDFFLAGS = -lm -lnetcdf -I/usr/include/netcdf-3 -L/usr/lib64 

#LDFLAGS = -lm -ftrap -lnetcdf -I/usr/local/include -L/usr/local/lib
#LDFLAGS = -lm -ftrap -lnetcdf -I/usr/include -L/usr/lib
#LDFLAGS = -lm -lnetcdf -I/usr/local/include -L/usr/local/lib
SRCDIR = src

#pendragon with icc
#CC=icc
#CFLAGS= -ip -ipo -inline-level=2 -xhost -O3 -g -openmp -lpthread 
#CFLAGS= -ip -ipo -inline-level=2 -xhost -O3 -g -lpthread 

CDFFLAGS= -lnetcdf


OFFSRC = $(SRCDIR)/offtrac.c $(SRCDIR)/read.c \
	$(SRCDIR)/initialize.c $(SRCDIR)/iocdf.c $(SRCDIR)/par_IO.c \
	$(SRCDIR)/tracadv.openmp.c $(SRCDIR)/sponge.c $(SRCDIR)/step.c \
        $(SRCDIR)/alloc_trac.c \
        $(SRCDIR)/alloc.c $(SRCDIR)/alloc_arrays.c \
	$(SRCDIR)/masks.c $(SRCDIR)/set_metrics.c\
	$(SRCDIR)/util.c $(SRCDIR)/biotic.c\
	$(SRCDIR)/oxygen.c $(SRCDIR)/phosphate.c \
	$(SRCDIR)/tracer_utilities.c 

offtrac: $(OFFSRC) $(SRCDIR)/init.h
	$(CC)  $(OFFSRC) -o $(OUTNAME) $(CDFFLAGS) $(CFLAGS) $(LDFLAGS)
	rm -f off offsf6

propre:
	rm -f *.o *.u

cleaner:
	rm -f offtrac *.o *.u 

clean:
	rm -f offtrac $(SRCDIR)/*.o *.u off

