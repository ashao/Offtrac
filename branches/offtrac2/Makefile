OUTNAME=offtrac

#CC=gcc


CC=gcc
CFLAGS = -O3 -g -fopenmp -lm -lpthread --fast-math -march="athlon64" -pipe -static
LDFLAGS=-L/usr/local/netcdf-gfort/lib -I/usr/local/netcdf-gfort/include 

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


OFFSRC = $(SRCDIR)/offtrac.c $(SRCDIR)/initialize_offtrac.c \
	$(SRCDIR)/read_netcdf.c $(SRCDIR)/util.c \
	$(SRCDIR)/iocdf.c 
 
offtrac: $(OFFSRC) $(SRCDIR)/user_options.h
	$(CC)  $(OFFSRC) -o $(OUTNAME) $(CDFFLAGS) $(CFLAGS) $(LDFLAGS)
	rm -f off offsf6

propre:
	rm -f *.o *.u

cleaner:
	rm -f offtrac *.o *.u 

clean:
	rm -f offtrac $(SRCDIR)/*.o *.u off

