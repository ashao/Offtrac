//HF new utility routines
#include "init.h"

void set_darray2d_zero(double **arr, int NX, int NY) {
	int x, y;
	for (x = 0; x < NX; x++)
		for (y = 0; y < NY; y++)
			arr[x][y] = 0.0;
}

void set_darray3d_zero(double ***arr, int nz, int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr[z][x][y] = 0.0;
}

void set_fix_darray3d_zero(double arr[][NXMEM][NYMEM], int nz) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NXMEM; x++)
			for (y = 0; y < NYMEM; y++)
				arr[z][x][y] = 0.0;
}

void set_fix_darray2d_zero(double arr[NXMEM][NYMEM]) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			arr[x][y] = 0.0;
}

void add_fix_darray2d(double arr1[NXMEM][NYMEM], double arr2[NXMEM][NYMEM]) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			arr1[x][y] += arr2[x][y];
}

void mult_fix_darray2d(double arr[NXMEM][NYMEM], double factor) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			arr[x][y] *= factor;
}

void add_darray3d(double ***arr1, double ***arr2, int nz, int NX, int NY) {
	int x, y, z;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] += arr2[z][x][y];
}

void mult_darray3d(double ***arr, int nz, int NX, int NY, double factor) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr[z][x][y] *= factor;
}
//BX-a
void mult_fix_darray2d_mv(double arr[NXMEM][NYMEM], double factor,
		double D[NXMEM][NYMEM], double mv) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			if (D[x][y] > MINIMUM_DEPTH) {
				arr[x][y] *= factor;
			} else {
				arr[x][y] = mv;
			}
}

void mult_darray3d_mv(double ***arr, int nz, int NX, int NY, double factor,
		double D[NXMEM][NYMEM], double mv) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				if (D[x][y] > MINIMUM_DEPTH) {
					arr[z][x][y] *= factor;
				} else {
					arr[z][x][y] = mv;
				}
}
//BX-e
/*
 void copy_fix_darray2d(double arr1[NXMEM][NYMEM], double arr2[NXMEM][NYMEM])
 {
 int x, y;
 for (x=0; x<NXMEM; x++)
 for (y=0;y<NYMEM; y++)
 arr1[x][y] = arr2[x][y];
 }
 */
void copy_darray3d(double ***arr1, double ***arr2, int nz, int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] = arr2[z][x][y];
}
//BX-a
void copy_fix_darray3d(double ***arr1, double arr2[NZ][NXMEM][NYMEM], int nz,
		int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] = arr2[z][x][y];
}
//BX-e
// begin ashao

void copy_darray3d_fix(double ***arr1, double arr2[NZ][NXMEM][NYMEM], int nz,
		int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr2[z][x][y] = arr1[z][x][y];
}
void copy_2fix_darray3d(double (*arr1)[NXMEM][NYMEM], double (*arr2)[NXMEM][NYMEM], int nz,
		int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] = arr2[z][x][y];
}
double linear_interp(double x0, double y0, double x1, double y1, double xstar) {
	double ystar;

	// If the two interpolation points are identical, set the output y to input y
	if ( (fabs(y1-y0) < .000000001) || (fabs(x1-x0) < .0000000001) ){
		 ystar=y0;
	}
	// Perform linear interpolation
	else {
		ystar = y0 + (xstar - x0) * (y1 -  y0) / (x1 - x0);
	}

	return ystar;
}

void copy_to_main_tracer_array(double ***array, int idx,
		double ****trarray)
{
	int i,j,k;
	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				trarray[idx][k][i][j] = array[k][i][j];
}

void copy_from_main_tracer_array(double ***array, int idx,
		double ****trarray)
{
	int i,j,k;
	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				array[k][i][j] = trarray[idx][k][i][j];
}

void reentrance_3d( double array[NZ][NXMEM][NYMEM] ) {
	int i,ii,j, k;
	for (k = 0; k < NZ; k++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			array[k][0][j] = array[k][nx - 1][j];
			array[k][1][j] = array[k][nx][j];
			array[k][nx + 1][j] = array[k][2][j];
			array[k][nx + 2][j] = array[k][3][j];
		}
	}
	for (i = 2; i <= nx; i++) {
		ii = (NXMEM-1) - i;
		for (k = 0; k < NZ; k++) {
			array[k][ii][ny + 1] = array[k][i][ny];
			array[k][ii][ny + 2] = array[k][i][ny - 1];
		}
	}
}

void reentrance_2d( double array[NXMEM][NYMEM]) {
	int i,ii,j;

	for (j = 0; j <= NYMEM - 1; j++) {
		array[0][j] = array[nx - 1][j];
		array[1][j] = array[nx][j];
		array[nx + 1][j] = array[2][j];
		array[nx + 2][j] = array[3][j];
	}

	for (i = 2; i <= nx; i++) {
		ii = (NXMEM-1) - i;
		array[ii][ny + 1] = array[i][ny];
		array[ii][ny + 2] = array[i][ny - 1];
	}
}

void update_average_array( double avgarray[NZ][NXMEM][NYMEM], double current[NZ][NXMEM][NYMEM], double wt)
{

	int i,j,k;
	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				avgarray[k][i][j] += current[k][i][j]*wt;


}

// end ashao
