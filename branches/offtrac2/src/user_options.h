// Grid options
#define MINIMUM_DEPTH	1e-6
#define NXTOT   360				//	NXTOT and NYTOT are the numbers
#define NYTOT   210				//	thickness grid points in the zonal
						//	and meridional directions of the
						//	physical domain.

#define NZ  49					//	The number of layers.

#define X0 1                 			//	X0 is the offset of the physical
						//	domain from the memory domain.
#define X1 (X0+1)     			        //    X1 is the lowest index of the
						// physical domain (as opposed to the
								// memory halo) on each processor.
#define Y0 1
#define Y1 (Y0+1)
#define NXMEM (NXTOT+2*(X0+1))
#define NYMEM (NYTOT+2*(Y0+1))

// Memory options
#define MAXLEN 100				// Maximum string length in characters

