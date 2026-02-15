/*
	File: marching_cubes_33_c.h
	Programmed by: David Vega - dvega@uc.edu.ve
	February 2026
	This library is a modified version of the library described in:
	Vega, D., Abache, J., Coll, D., A Fast and Memory Saving Marching Cubes 33
	implementation with the correct interior test, Journal of Computer Graphics
	Techniques (JCGT), vol. 8, no. 3, 1–18, 2019
*/

#ifndef marching_cubes_33_c_h
#define marching_cubes_33_c_h

/********************************** USAGE ************************************/
/*
//1. Header
	#define marching_cubes_33_c_implementation
	#include "marching_cubes_33_c.h"

//2. Read a grid file.
	_GRD* G = read_dat_file(stagbeetle832x832x494.dat);
	// The dataset stagbeetle832x832x494.dat can be downloaded at
	// https://www.cg.tuwien.ac.at/research/publications/2005/dataset-stagbeetle/
	// There are other grid reading functions that are described later.

//3. Create a MC33 structure.
	MC33 *M = create_MC33(G);

//4. calculate an isosurface with isovalue v (a float).
	surface *S = calculate_isosurface(M, v);

//5. When finished, free the memory occupied by S, M, and G.
	free_surface_memory(S);
	free_MC33(M);
	free_memory_grd(G);
*/

/******************************* CUSTOMIZING *********************************/
//You can change this defines:

//#define GRD_INTEGER // Uncomment this define for dataset with integer type
//#define GRD_TYPE_SIZE 4 // 1, 2, 4 or 8 (8 for double, if not defined GRD_INTEGER)
//#define GRD_ORTHOGONAL // If defined, the library only works with orthogonal grids.
//#define MC_NORMAL_NEG // the front and back surfaces are exchanged.
//#define DEFAULT_SURFACE_COLOR 0xFF80FF40// RGBA 0xAABBGGRR: red 64, green 255, blue 128
/*****************************************************************************/

#define MC33C_VERSION_MAJOR 5
#define MC33C_VERSION_MINOR 4

#if defined(integer_GRD)
#if GRD_TYPE_SIZE == 4
/*
GRD_data_type is the variable type of the grid data, by default it is float.
*/
typedef unsigned int GRD_data_type;
#elif GRD_TYPE_SIZE == 2
typedef unsigned short int GRD_data_type;
#elif GRD_TYPE_SIZE == 1
typedef unsigned char GRD_data_type;
#else
#error "Incorrect size of the data type. GRD_TYPE_SIZE permitted values: 1, 2 or 4."
typedef float GRD_data_type;
#endif
#elif GRD_TYPE_SIZE == 8
typedef double GRD_data_type;
#else
typedef float GRD_data_type;
#undef GRD_TYPE_SIZE
#define GRD_TYPE_SIZE 4
#endif

/*
The structure _GRD contains a function F[][][] evaluated at a grid of regularly
spaced points. N[] is the number of intervals in each dimension, so the number
of points is N[] + 1, and the grid contain (N[2] + 1)*(N[1] + 1)*(N[0] + 1)
points. L[] is the grid size in each dimension. r0[] are the coordinates of
the first grid point. d[] is the distance between adjacent points in each
dimension (can be different for each dimension), nonortho has a value of 1 when
the grid is inclined else 0. _A and A_ are the matrices that transform from
inclined to orthogonal coordinates and vice versa, respectively. If the grid is
periodic (is infinitely repeated along each dimension) the flag periodic must
be 1 (for x) + 2 (for y) + 4 (for z), else 0.

In this program, if GRD_ORTHOGONAL is defined, then nonortho, _A and A_ can be
removed from this structure, and only work with orthogonal grids. periodic
and L[] are not used here and also can be removed.
*/

typedef struct {
	GRD_data_type ***F;
	unsigned int N[3];
	double r0[3], d[3];
	float L[3]; 
#ifndef GRD_ORTHOGONAL
	float Ang[3];//not necessary
	int nonortho; //necessary if GRD_ORTHOGONAL is not defined
	double _A[3][3], A_[3][3]; //necessary if GRD_ORTHOGONAL is not defined
#endif
	int periodic;//not necessary
	int internal_data;
	char title[160];//not necessary
} _GRD;

/*The struct surface contains the data of an isosurface.
The number of points is nV, and the number of triangles is nT. The coordinates
of the points are stored in the array V[][3].
N[][3] is the array that contains the normal vectors calculated at the points
of the surface. The array color[] contains the color of each point.
T[][3] is the array that contains the sets of three point (vertex) indices
that form each triangle.*/
typedef struct {
/* A is a mask (n & A is equivalent to n % dim1).
n >> E is equivalent to n / dim1. */
	unsigned int (*T)[3];
	float (*V)[3];
	float (*N)[3];
	int *color;
	unsigned int nV, nT;
	float iso;
	union {
		void *p;
		long long ul;
		int i[2];
		short si[4];
		char c[8];
		float f[2];
		double df;
	} user; // user data
} surface;

typedef struct {
// copy of some variables of the surface structure:
	unsigned int (*T)[3];
	float (*V)[3];
	float (*N)[3];
	int *color;
	unsigned int nV, nT;
	float iso;

/*memoryfault takes a non-zero value if the system memory is insufficient
to store the surface.*/
	int memoryfault;

	unsigned int capt, capv;

// copy of some variables of the _GRD structure:
	const GRD_data_type ***F;
	float O[3], D[3], ca, cb;
	unsigned int nx, ny, nz;
	unsigned int (*store)(void *, float *);
#ifndef GRD_ORTHOGONAL
	double _A[3][3], A_[3][3];
#endif

// temporary structures that store the indices of triangle vertices:
	unsigned int **Dx, **Dy, **Ux, **Uy, **Lz;
} MC33;

extern int DefaultColorMC;
/* Example:
DefaultColorMC = 0x00FF0080
										 B G R
Red 128, green 0, blue 255.
*/

#ifndef GRD_ORTHOGONAL
void _multTSA_bf(const double (*A)[3], float *b, float *c, int t);
void _multA_bf(const double (*A)[3], float* b, float* c, int t);
extern void (*mult_Abf)(const double (*)[3], float *, float *, int);
#endif /* GRD_ORTHOGONAL */

/******************************************************************
Saves all the surface *S data (in binary format) to a "filename" file. The
return value is 0 if the call succeeds, else -1.
*/
int write_bin_s(surface *S, const char *filename);
/******************************************************************
Reads (from a "filename" file) the surface data stored in binary format.
The return value is a surface pointer if the call succeeds, else NULL.
*/
surface* read_bin_s(const char *filename);
/******************************************************************
Saves all the surface *S data (in plain text format) to a "filename" file.
The return value is 0 if the call succeeds, else -1.
*/
int write_txt_s(surface *S, const char *filename);
/******************************************************************
Saves the surface *S data (without the color) to Wavefront .obj file.
The return value is 0 if the call succeeds, else -1.
*/
int write_obj_s(surface *S, const char *filename);
/******************************************************************
Saves the surface *S data to Polygon File Format .ply file.
https://paulbourke.net/dataformats/ply/
The return value is 0 if the call succeeds, else -1.*/
int write_ply_s(surface *S, const char *filename, const char* author, const char* object);
/******************************************************************
Creates an MC33 structure from a pointer to a _GRD structure. If successful,
returns a pointer to the MC33 structure. On error, returns null pointer.
*/
MC33 *create_MC33(_GRD *G);

/******************************************************************
Calculates the isosurface (iso is the isovalue) using a MC33 structure
pointed by M. The return value is a pointer to surface. The pointer will
be NULL if there is not enough memory. The isovalue is stored in the
member iso of surface struct.*/
surface* calculate_isosurface(MC33 *M, float iso);

/******************************************************************
Return the size of surface without calculate it. The function can
calculate the number of vertices and triangles. If only the size is
required:
unsigned long long size = size_of_isosurface(M, v, 0, 0);
*/
unsigned long long size_of_isosurface(MC33 *M, float iso, unsigned int *nV, unsigned int *nT);

/******************************************************************
Releases the allocated memory occupied by MC33 structure pointed by M.
*/
void free_MC33(MC33 *M);

/******************************************************************
Releases the allocated memory pointed to by S.
*/
void free_surface_memory(surface *S);

/******************************************************************
Free the memory occupied by the _GRD structure pointed to by Z.
*/
void free_memory_grd(_GRD *Z);

/******************************************************************
Allocate memory for the grid data. Before calling this function, memory must
be allocated for the _GRD structure and values must be assigned to the number
of points in each dimension (array N of the _GRD structure). The function
returns -1 if the allocation fails. Example:
	_GRD *Z = (_GRD *)malloc(sizeof(_GRD));
	Z->N[0] = 200; Z->N[1] = 200; Z->N[2] = 100;
	if (alloc_F(Z))
		return; // memory error
// fill Z->F[][][] here.
*/
int alloc_F(_GRD* Z);

/******************************************************************
read_grd read a filename file (the file must be a output *.grd file from the
DMol program), it returns a pointer to struct _GRD that contains all the grid
data.
*/
_GRD* read_grd(const char *filename);
/* Internal binary format
*/
_GRD* read_grd_binary(const char* filename);

/******************************************************************
Reads a set of files that contain a slab of res*res scan data points, the data
points are read as unsigned short int (if order is different from 0, the bytes
of the unsigned short are exchanged). The filename must end with a number, and
the fuction read all files with end number greater or equal to filename.
(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)
*/
_GRD* read_scanfiles(const char *filename, unsigned int res, int order);

/******************************************************************
Reads a file that contains only the data points as integers of 8, 16 or 32 bits.
byte is the number of bytes of each data point (1 to 4). If the data is big endian,
byte must be negative. The vector N[3] contains the number of points in each
dimension. The size of file must be byte*N[0]*N[1]*N[2].
*/
_GRD* read_raw_file(const char *filename, unsigned int *N, int byte, int isfloat);

/******************************************************************
Reads a dat file:
http://www.cg.tuwien.ac.at/research/vis/datasets/
*/
_GRD* read_dat_file(const char *filename);

/******************************************************************
	set_data_pointer creates a _GRD struct from an external data array. data
	must be stored with the nested inner loop running from i = 0 to Nx - 1
	and the outer loop from k = 0 to Nz - 1. The data will not be erased by
	free_memory_grd function. The function returns a pointer to the created
	struct.
*/
_GRD* grid_from_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data);

/******************************************************************
	Build a grid by using a scalar function
	double fn(double x, double y, double z)
*/
_GRD* generate_grid_from_fn(
	double x_initial, double y_initial, double z_initial,
	double x_final, double y_final, double z_final,
	double x_step, double y_step, double z_step,
	double (*fn)(double x, double y, double z));

#ifdef marching_cubes_33_c_implementation

#ifndef GRD_ORTHOGONAL
/* c = Ab, A is a 3x3 upper triangular matrix. If t != 0, A is transposed. */
void _multTSA_bf(const double (*A)[3], float *b, float *c, int t) {
	if(t) {
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
		c[1] = A[0][1]*b[0] + A[1][1]*b[1];
		c[0] = A[0][0]*b[0];
	} else {
		c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		c[1] = A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][2]*b[2];
	}
}
/* Performs the multiplication of the matrix A and the vector b: c = Ab. If t != 0, A is transposed. */
void _multA_bf(const double (*A)[3], float* b, float* c, int t) {
	double u,v;
	if(t) {
		u = A[0][0]*b[0] + A[1][0]*b[1] + A[2][0]*b[2];
		v = A[0][1]*b[0] + A[1][1]*b[1] + A[2][1]*b[2];
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
	} else {
		u = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		v = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
	}
	c[0] = u;
	c[1] = v;
}

void (*mult_Abf)(const double (*)[3], float *, float *, int) = _multA_bf;

void setIdentMat3x3d(double (*A)[3]) {
	for (double *d = A[0] + 8; --d != A[0];)
		d[0] = 0.0;
	for (int i = 0; i != 3; i++)
		A[i][i] = 1.0;
}
#endif // GRD_ORTHOGONAL

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

/******************************************************************/
void free_memory_grd(_GRD *Z) {
	unsigned int k, j;
	if (Z) {
		if (Z->F) {
			if (Z->internal_data)
				for (k = 0; k <= Z->N[2]; k++) {
					if (Z->F[k]) {
						for (j = 0; j <= Z->N[1]; j++)
							free(Z->F[k][j]);
					} else
						break;
					free(Z->F[k]);
				}
			else
				for (k = 0; k <= Z->N[2]; k++)
					free(Z->F[k]);
			free(Z->F);
		}
		free(Z);
	}
}

/******************************************************************/
int alloc_F(_GRD* Z) {
	unsigned int j, k;
	Z->F = (GRD_data_type***)malloc((Z->N[2] + 1)*sizeof(void*));
	if (!Z->F)
		return -1;
	for (k = 0; k <= Z->N[2]; k++) {
		Z->F[k] = (GRD_data_type**)malloc((Z->N[1] + 1)*sizeof(void*));
		if (!Z->F[k])
			return -1;
		for (j = 0; j <= Z->N[1]; j++) {
			Z->F[k][j] = (GRD_data_type*)malloc((Z->N[0] + 1)*sizeof(GRD_data_type));
			if (!Z->F[k][j]) {
				while (j)
					free(Z->F[k][--j]);
				free(Z->F[k]);
				Z->F[k] = 0;
				return -1;
			}
		}
	}
	Z->internal_data = 1;
	return 0;
}

/******************************************************************/
_GRD* read_grd(const char *filename) {
	_GRD *Z;
	FILE *in;
	char line[128];
	unsigned int i, j, k;
	double ca, cb, sg, cg, aux1, aux2;
	int grd_xi[3], grd_ordenij;
	float Ang[3];

	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;

#ifndef GRD_ORTHOGONAL
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	in = fopen(filename,"r");
	if (!in) return 0;
	fgets(Z->title,159,in);
	fgets(line,60,in);
	fgets(line,60,in);
	sscanf(line,"%f %f %f %f %f %f",&(Z->L[0]),&(Z->L[1]),&(Z->L[2]),&Ang[0],&Ang[1],&Ang[2]);
	fgets(line,60,in);
	sscanf(line,"%d %d %d",&(Z->N[0]),&(Z->N[1]),&(Z->N[2]));
	fgets(line,60,in);
	sscanf(line,"%d %d %*d %d %*d %d %*d",&grd_ordenij,&grd_xi[0],&grd_xi[1],&grd_xi[2]);
	if (Z->N[0] < 2 || Z->N[1] < 2 || Z->N[2] < 2) return 0;
	if (grd_ordenij != 1 && grd_ordenij != 3) return 0;
	for (i = 0; i != 3; i++) {
		Z->d[i] = Z->L[i]/Z->N[i];
		Z->r0[i] = grd_xi[i]*Z->d[i];
	}

	Z->periodic = (grd_xi[0] == 0)|((grd_xi[1] == 0)<<1)|((grd_xi[2] == 0)<<2);

#ifndef GRD_ORTHOGONAL
	memcpy(Z->Ang, Ang, sizeof Ang);
	if (Ang[0] != 90 || Ang[1] != 90 || Ang[2] != 90) {
		Z->nonortho = 1;
		ca = cos(Ang[0]*(M_PI/180.0));
		cb = cos(Ang[1]*(M_PI/180.0));
		aux1 = Ang[2]*(M_PI/180.0);
		sg = sin(aux1);
		cg = cos(aux1);
		aux1 = ca - cb*cg;
		aux2 = sqrt(sg*sg + 2*ca*cb*cg - ca*ca - cb*cb);
		Z->_A[0][0] = Z->A_[0][0] = 1.0;
		Z->_A[0][1] = cg;
		Z->_A[0][2] = cb;
		Z->_A[1][1] = sg;
		Z->A_[1][1] = cb = 1.0/sg;
		Z->A_[0][1] = -cg*cb;
		Z->_A[1][2] = aux1*cb;
		Z->_A[2][2] = aux2*cb;
		aux2 = 1.0/aux2;
		Z->A_[0][2] = (cg*aux1 - ca*sg*sg)*cb*aux2;
		Z->A_[1][2] = -aux1*cb*aux2;
		Z->A_[2][2] = sg*aux2;
	} else {
		Z->nonortho = 0;
		setIdentMat3x3d(Z->A_);
		setIdentMat3x3d(Z->_A);
	}
#endif

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

	for (k = 0; k <= Z->N[2]; k++)
		if (grd_ordenij == 1) {
			for (j = 0; j <= Z->N[1]; j++)
				for (i = 0; i <= Z->N[0]; i++)
					fscanf(in,"%f",&Z->F[k][j][i]);
		} else {
			for (i = 0; i <= Z->N[0]; i++)
				for (j = 0; j <= Z->N[1]; j++)
					fscanf(in,"%f",&Z->F[k][j][i]);
		}
	fclose(in);
	return Z;
}

/******************************************************************/
_GRD* read_grd_binary(const char* filename) {
	FILE* in;
	unsigned int i, j, k;
	_GRD* Z;

	in = fopen(filename,"rb");
	if (!in)
		return 0;
	fread(&i,sizeof(int),1,in);
	if (i != 0x4452475f) // _GRD
		return 0;

	fread(&i,sizeof(int),1,in);
	if (i > 159)
		return 0;
	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;
	Z->internal_data = 1;
	fread(Z->title,i*sizeof(char),1,in);
//	Z->title[i] = '\0';
	fread(Z->N,sizeof Z->N,1,in);
	fread(Z->L,sizeof Z->L,1,in);
	fread(Z->r0,sizeof Z->r0,1,in);
	fread(Z->d,sizeof Z->d,1,in);
#ifndef GRD_ORTHOGONAL
	fread(&Z->nonortho,sizeof(int),1,in);
	if (Z->nonortho) {
		fread(Z->Ang,3*sizeof(float),1,in);
		fread(Z->_A,sizeof Z->_A,1,in);
		fread(Z->A_,sizeof Z->A_,1,in);
		mult_Abf = _multA_bf;
	} else {
		setIdentMat3x3d(Z->A_);
		setIdentMat3x3d(Z->_A);
	}
#else
	fread(&i, sizeof(int), 1, in);
	if (i)
		fseek(in, 3*sizeof(float) + 18*sizeof(double), SEEK_CUR);
#endif
	if (Z->r0[0] == 0 && Z->r0[1] == 0 && Z->r0[2] == 0)
		Z->periodic = 1;

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

	for (k = 0; k <= Z->N[2]; k++)
		for (j = 0; j <= Z->N[1]; j++)
			fread(Z->F[k][j],(Z->N[0] + 1)*sizeof(GRD_data_type),1,in);
	fclose(in);
	return Z;
}

/******************************************************************/
_GRD* read_scanfiles(const char *filename, unsigned int res, int order) {
	_GRD *Z;
	FILE *in;
	char *nm, *nm2;
	unsigned int i, j, l;
	int k = -1;
	unsigned short int n;
	GRD_data_type ***pt, **p;

	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;

#ifndef GRD_ORTHOGONAL
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	memset(Z->r0,0,sizeof Z->r0);
	for (i = 0; i != 2; i++) {
		Z->d[i] = 1;
		Z->L[i] = Z->N[i] = res - 1;
	}
	Z->d[2] = 1;
	Z->F = 0;
	l = strlen(filename) - 1;
	nm = (char*)malloc((l + 5)*sizeof(char));
	if (!nm)
		return 0;
	while (filename[l] >= '0' && filename[l] <= '9')
		l--;
	strncpy(nm,filename,l + 1);
	nm2 = nm + l + 1;
	l = atoi(filename + l + 1);
	while (1) {
		sprintf(nm2,"%-d",l++);
		in = fopen(nm,"rb");
		if (!in)
			break;
		if (!((++k)&63)) {
			pt = (GRD_data_type***)realloc(Z->F,(k + 64)*sizeof(void*));
			if(!pt)
				break;
			Z->F = pt;
		}
		Z->F[k] = (GRD_data_type**)malloc(res*sizeof(void*));
		if (!Z->F[k])
			break;
		for (j = 0; j != res; j++)
			Z->F[k][j] = (GRD_data_type*)malloc(res*sizeof(GRD_data_type));
		if (!Z->F[k][Z->N[1]])
			break;

		if (order)
			for (j = 0; j != res; j++)
				for (i = 0; i != res; i++) {
					fread(&n,sizeof(short int),1,in);
					Z->F[k][j][i] = (unsigned short int)((n>>8)|(n<<8));
				}
		else
			for (j = 0; j != res; j++)
				for (i = 0; i != res; i++) {
					fread(&n,sizeof(short int),1,in);
					Z->F[k][j][i] = n;
				}
		fclose(in);
	}
	free(nm);
	Z->L[2] = Z->N[2] = k;
	j = k>>1;
	for (i = 0; i != j; i++) {
		p = Z->F[i];
		Z->F[i] = Z->F[k - i];
		Z->F[k - i] = p;
	}
#ifndef GRD_ORTHOGONAL
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

/******************************************************************/
_GRD* read_raw_file(const char *filename, unsigned int *N, int byte, int isfloat) {
	unsigned int i, j, k;
	_GRD *Z;
	FILE *in;
	unsigned int ui = 0;
	if (isfloat) {
		if (byte != 4 && byte != 8)
			return 0;
	} else if (abs(byte) > 4 || abs(byte) == 3 || !byte)
		return 0;
	if (byte == -1) byte = 1;
	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z)
		return 0;

#ifndef GRD_ORTHOGONAL
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	memset(Z->r0,0,sizeof Z->r0);
	for (i = 0; i != 3; i++)
		Z->d[i] = 1;
	in = fopen(filename,"rb");
	if (!in)
		return 0;
	Z->L[0] = Z->N[0] = N[0] - 1;
	Z->L[1] = Z->N[1] = N[1] - 1;
	Z->L[2] = Z->N[2] = N[2] - 1;

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

#if defined(GRD_INTEGER)
	if (!isfloat && GRD_TYPE_SIZE == byte)
#else
	if (isfloat && GRD_TYPE_SIZE == byte)
#endif
	{
		byte *= N[0];
		for (k = 0; k != N[2]; k++)
			for (j = 0; j != N[1]; j++)
				fread(Z->F[k][j],byte,1,in);
	} else if (isfloat) {
		if (byte == 8) {
#if defined(GRD_INTEGER) || GRD_TYPE_SIZE == 4
			double df;
			for (k = 0; k != N[2]; k++)
				for (j = 0; j != N[1]; j++)
					for (i = 0; i != N[0]; i++) {
						fread(&df,byte,1,in);
						Z->F[k][j][i] = (GRD_data_type)df;
					}
#endif
		} else {
#if defined(GRD_INTEGER) || GRD_TYPE_SIZE == 8
			float f;
			for (k = 0; k != N[2]; k++)
				for (j = 0; j != N[1]; j++)
					for (i = 0; i != N[0]; i++) {
						fread(&f,byte,1,in);
						Z->F[k][j][i] = (GRD_data_type)f;
					}
#endif
		}
	} else if (byte < 0) {
		byte = -byte;
		for (k = 0; k != N[2]; k++)
			for (j = 0; j != N[1]; j++)
				for (i = 0; i != N[0]; i++) {
					fread(&ui,byte,1,in);
					if (byte == 2)
						Z->F[k][j][i] = (ui>>8)|((ui<<8)&0xff00);
					else //if (byte == 4)
						Z->F[k][j][i] = (ui>>24)|((ui>>8)&0xff00)|((ui<<8)&0xff0000)|(ui<<24);
				}
	} else {
		for (k = 0; k != N[2]; k++)
			for (j = 0; j != N[1]; j++)
				for (i = 0; i != N[0]; i++) {
					fread(&ui,byte,1,in);
					Z->F[k][j][i] = ui;
				}
	}
	fclose(in);
#ifndef GRD_ORTHOGONAL
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

/******************************************************************/
_GRD* read_dat_file(const char *filename) {
	_GRD *Z;
	FILE *in;
	unsigned int i, j;
	unsigned short int n, nx, ny, nz;

	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z)
		return 0;

#ifndef GRD_ORTHOGONAL
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	memset(Z->r0,0,sizeof Z->r0);
	for (i = 0; i != 3; i++)
		Z->d[i] = 1;
	in = fopen(filename,"rb");
	if (!in)
		return 0;
	fread(&nx,sizeof(short int),1,in);
	fread(&ny,sizeof(short int),1,in);
	fread(&nz,sizeof(short int),1,in);
	Z->L[0] = Z->N[0] = nx - 1;
	Z->L[1] = Z->N[1] = ny - 1;
	Z->L[2] = Z->N[2] = nz - 1;

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

	while (nz--)
		for (j = 0; j != ny; j++)
			for (i = 0; i != nx; i++) {
				fread(&n,sizeof(short int),1,in);
				Z->F[nz][j][i] = n;
			}
	fclose(in);
#ifndef GRD_ORTHOGONAL
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

/******************************************************************/
_GRD* grid_from_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data) {
	if (!data || Nx == 0 || Ny == 0 || Nz == 0)
		return 0;
	unsigned int i, j, k;
	_GRD *Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;
	Z->internal_data = 0;
	Z->F = (GRD_data_type***)malloc(Nz*sizeof(void*));
	if (!Z->F) {
		free(Z);
		return 0;
	}
	Z->N[0] = Nx - 1;
	Z->N[1] = Ny - 1;
	Z->N[2] = Nz - 1;
	for (k = 0; k < Nz; k++) {
		Z->F[k] =(GRD_data_type**)malloc(Ny*sizeof(void*));
		if (!Z->F[k]) {
			while (k)
				free(Z->F[--k]);
			free(Z->F);
			free(Z);
			return 0;
		}
		for (j = 0; j < Ny; j++)
			Z->F[k][j] = data + j*Nx;
		data += Ny*Nx;
	}
	for (i = 0; i != 3; i++) {
		Z->L[i] = Z->N[i];
		Z->d[i] = 1.0;
		Z->r0[i] = 0.0;
#ifndef GRD_ORTHOGONAL
		Z->Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_ORTHOGONAL
	Z->nonortho = 0;
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

/******************************************************************/
_GRD* generate_grid_from_fn(double xi, double yi, double zi, double xf, double yf, double zf,
	double dx, double dy, double dz, double (*fn)(double x, double y, double z)) {
	if (dx <= 0 || dy <= 0 || dz <= 0 || xi == xf || yi == yf || zi == zf)
		return 0;

	if (xi > xf) {
		double t = xi; xi = xf; xf = t;}
	if (xf - xi < dx)
		dx = xf - xi;
	if (yi > yf) {
		double t = yi; yi = yf; yf = t;}
	if (yf - yi < dy)
		dy = yf - yi;
	if (zi > zf) {
		double t = zi; zi = zf; zf = t;}
	if (zf - zi < dz)
		dz = zf - zi;
	_GRD *Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;
	Z->N[0] = (int)((xf - xi)/dx + 0.5);
	Z->N[1] = (int)((yf - yi)/dy + 0.5);
	Z->N[2] = (int)((zf - zi)/dz + 0.5);
	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}
	Z->d[0] = dx; Z->d[1] = dy; Z->d[2] = dz;
	Z->r0[0] = xi; Z->r0[1] = yi; Z->r0[2] = zi;
	if (fn) {
		unsigned int i, j, k;
		double x, y, z = zi;
		for (k = 0; k <= Z->N[2]; k++) {
			y = yi;
			for (j = 0; j <= Z->N[1]; j++) {
				x = xi;
				for (i = 0; i <= Z->N[0]; i++) {
					Z->F[k][j][i] = (GRD_data_type)fn(x,y,z);
					x += dx;
				}
				y += dy;
			}
			z += dz;
		}
	}
	for (int i = 0; i != 3; i++) {
		Z->L[i] = Z->N[i]*Z->d[i];
#ifndef GRD_ORTHOGONAL
		Z->Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_ORTHOGONAL
	Z->nonortho = 0;
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

#if defined (__SSE__) || ((defined (_M_IX86) || defined (_M_X64)) && !defined(_CHPE_ONLY_))
// https://stackoverflow.com/questions/59644197/inverse-square-root-intrinsics
// faster than 1.0f/std::sqrt, but with little accuracy.
#include <immintrin.h>
inline float invSqrt(float f) {
	__m128 temp = _mm_set_ss(f);
	temp = _mm_rsqrt_ss(temp);
	return _mm_cvtss_f32(temp);
}
#else
inline float invSqrt(float f) {
	return 1.0/sqrt(f);
}
#endif

#ifndef DEFAULT_SURFACE_COLOR
#define DEFAULT_SURFACE_COLOR 0xff5c5c5c;
#endif
int DefaultColorMC = DEFAULT_SURFACE_COLOR;

/****************** Surface managing functions ****************/

void free_surface_memory(surface *S) {
	if (S) {
		free(S->T);
		free(S->V);
		free(S->N);
		free(S->color);
		free(S);
	}
}

int write_bin_s(surface *S, const char *filename) {
	int i;
	FILE *out = fopen(filename,"wb");
	if (!out)
		return -1;
	fputs(".sup",out);

	fwrite(S->user.f + 3,sizeof(float),1,out);
	fwrite(&S->nV,sizeof(int),1,out);
	fwrite(&S->nT,sizeof(int),1,out);

	fwrite(S->T,3*S->nT*sizeof(int),1,out);
	fwrite(S->V,3*S->nV*sizeof(float),1,out);
	fwrite(S->N,3*S->nV*sizeof(float),1,out);
	i = fwrite(S->color,S->nV*sizeof(int),1,out);
	fclose(out);
	return -(i != 1);
}

surface* read_bin_s(const char *filename) {
	surface *S;
	int i;

	FILE *in = fopen(filename,"rb");
	if (!in)
		return 0;
	fread(&i,sizeof(int),1,in);
	S = (surface*)malloc(sizeof(surface));
	if (i != 0x7075732e || !S) {
		fclose(in);
		free(S);
		return 0;
	}
	fread(S->user.f + 3,sizeof(float),1,in);
	fread(&S->nV,sizeof(int),1,in);
	fread(&S->nT,sizeof(int),1,in);

	S->T = (unsigned int (*)[3])malloc(S->nT*sizeof(void*));
	S->V = (float (*)[3])malloc(S->nV*sizeof(void*));
	S->N = (float (*)[3])malloc(S->nV*sizeof(void*));
	S->color = (int *)malloc(S->nV*sizeof(void*));
	if (!(S->T && S->V && S->N && S->color)) {
		free_surface_memory(S);
		fclose(in);
		return 0;
	}
	i = !fread(S->T,3*S->nT*sizeof(int),1,in);
	i += !fread(S->V,3*S->nV*sizeof(float),1,in);
	i += !fread(S->N,3*S->nV*sizeof(float),1,in);
	i += !fread(S->color,S->nV*sizeof(int),1,in);
	fclose(in);
	if (i) {
		free_surface_memory(S);
		return 0;
	}
	return S;
}

int write_txt_s(surface *S, const char *filename) {
	FILE *out;
	unsigned int i, *t;
	float *r;

	out = fopen(filename,"w");
	if (!out)
		return -1;

	fprintf(out,"isovalue: %10.5E\n\nVERTICES:\n",S->iso);
	fprintf(out,"%d\n\n",S->nV);
	for (i = 0; i != S->nV; i++) {
		r = S->V[i];
		fprintf(out,"%9.6f %9.6f %9.6f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"\n\nTRIANGLES:\n");
	fprintf(out,"%d\n\n",S->nT);
	for (i = 0; i != S->nT; i++) {
		t = S->T[i];
		fprintf(out,"%8d %8d %8d\n",t[0],t[1],t[2]);
	}

	fprintf(out,"\n\nNORMALS:\n");
	for (i = 0; i != S->nV; i++) {
		r = S->N[i];
		fprintf(out,"%9.6f %9.6f %9.6f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"\n\nCOLORS:\n");
	for (i = 0; i != S->nV; i++)
		fprintf(out,"%d\n",S->color[i]);
	i = fprintf(out,"\nEND\n");
	fclose(out);
	return -(i < 5);
}

int write_obj_s(surface *S, const char *filename) {
	FILE *out;
	unsigned int i, *t;
	float *r;
	char s0[12], s1[12], s2[12];

	out = fopen(filename,"w");
	if (!out)
		return -1;

	fprintf(out,"# isovalue: %10.5E\n# VERTICES %d:\n", S->iso, S->nV);
	for (i = 0; i != S->nV; i++) {
		r = S->V[i];
		fprintf(out,"v %f %f %f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"# NORMALS:\n");
	for (i = 0; i != S->nV; i++) {
		r = S->N[i];
		fprintf(out,"vn %f %f %f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"# TRIANGLES %d:\n",S->nT);
	for (i = 0; i != S->nT; i++) {
	    t = S->T[i];
		sprintf(s0,"%d", t[0] + 1);
		sprintf(s1,"%d", t[1] + 1);
		sprintf(s2,"%d", t[2] + 1);
		fprintf(out,"f %s//%s %s//%s %s//%s\n", s0, s0, s1, s1, s2, s2);
	}

	i = fprintf(out,"# END");
	fclose(out);
	return -(i < 5);
}

int write_ply_s(surface *S, const char *filename, const char* author, const char* object) {
	FILE *out;
	unsigned int i, *t;
	float *r;
	unsigned char* c;
	char empty = 0;

	out = fopen(filename,"w");
	if (!out)
		return -1;

	if (!author)
		author = &empty;
	if (!object)
		object = &empty;

	fprintf(out,"ply\nformat ascii 1.0\ncomment author: %s\ncomment object: %s\n",author,object);
	fprintf(out,"element vertex %d\nproperty float x\nproperty float y\nproperty float z",S->nV);
	fprintf(out,"\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face");
	fprintf(out," %d\nproperty list uchar int vertex_index\nend_header",S->nT);

	for (i = 0; i < S->nV; i++) {
		r = S->V[i];
		c = (unsigned char*)(S->color + i);
		fprintf(out,"\n%f %f %f %d %d %d",r[0],r[1],r[2],c[0],c[1],c[2]);
	}

	for (i = 0; i < S->nT; i++) {
		t = S->T[i];
		fprintf(out,"\n3 %d %d %d",t[0],t[1],t[2]);
	}

	i = fprintf(out,"\n");
	fclose(out);
	return -(!i);
}

const unsigned short int MC33_all_tables[2310] =
{
/* INDEX
The 11 least significant bits either contain the position of the triangle pattern
in this array (cases 1, 2, 5, 8, 9, 11, and 14) or are used to calculate it. The
12th bit is used to determine whether the order of the triangle vertices (and the
normal) should be reversed. The 4 most significant bits are related to the MC33 case.

     Case              Patern position
  __________     _______________________________
  F  E  D  C  B  A  9  8  7  6  5  4  3  2  1  0
              ^
           reverse
*/
 0x0000, 0x0885, 0x0886, 0x0895, 0x0883, 0x1816, 0x089D, 0x0943,//007
 0x0884, 0x0897, 0x1814, 0x0916, 0x0891, 0x094C, 0x091F, 0x048F,//015
 0x0882, 0x089B, 0x1808, 0x0934, 0x2803, 0x3817, 0x3814, 0x0525,//023
 0x180E, 0x0928, 0x4802, 0x049D, 0x3815, 0x0541, 0x6004, 0x0110,//031
 0x0881, 0x180A, 0x0899, 0x0922, 0x1806, 0x4800, 0x0913, 0x0499,//039
 0x2802, 0x380E, 0x3811, 0x053D, 0x380D, 0x6002, 0x0521, 0x010D,//047
 0x088B, 0x0946, 0x0937, 0x0493, 0x3813, 0x6003, 0x0545, 0x012E,//055
 0x380F, 0x0531, 0x600A, 0x011C, 0x5001, 0x3001, 0x3009, 0x0087,//063
 0x0880, 0x2801, 0x1804, 0x3807, 0x088F, 0x380B, 0x092B, 0x0539,//071
 0x1802, 0x3806, 0x4806, 0x6001, 0x0925, 0x051D, 0x0495, 0x010A,//079
 0x1812, 0x380A, 0x4805, 0x6009, 0x3816, 0x5002, 0x6008, 0x3010,//087
 0x4801, 0x6000, 0x7000, 0x4003, 0x6006, 0x3004, 0x4007, 0x1010,//095
 0x0889, 0x3808, 0x093A, 0x052D, 0x0949, 0x6005, 0x0491, 0x0131,//103
 0x380C, 0x5000, 0x600B, 0x3012, 0x0549, 0x3002, 0x0119, 0x008D,//111
 0x0907, 0x0535, 0x04A1, 0x013D, 0x0529, 0x3005, 0x0140, 0x0093,//119
 0x6007, 0x3000, 0x4004, 0x1000, 0x3003, 0x2000, 0x100C, 0x007f,//127

/*
  Vertices:            Edges:               Faces:
    3 ___________2        _____2______         ____________
   /|           /|      /|           /|      /|           /|
  / |          / |     B |          A |     / |    2     / |
7/___________6/  |    /_____6_____ /  |    /___________ /  |
|   |        |   |   |   3        |   1   |   |     4  |   |
|   |        |   |   |   |        |   |   | 3 |        | 1 |     z
|   0________|___1   |   |_____0__|___|   |   |_5______|___|     |
|  /         |  /    7  /         5  /    |  /         |  /      |____y
| /          | /     | 8          | 9     | /      0   | /      /
4/___________5/      |/_____4_____|/      |/___________|/      x

*/
/*
Vertices order in triangles:
       0     1-----2
      / \     \ x /   o: Front surface
     / o \     \ /    x: Back surface
    1-----2     0
*/
/* TRIANGLE PATERNS
Each short integer corresponds to one triangle. The most significant hexadecimal
digit is set to 0 for the last triangle in the pattern. The other 3 digits are
the cube edges where the triangle vertices are located.
*/
/*position*index#vertices*/
// Case 1 (128)
/* 0*127#0*/0x0380,
/* 1*064#1*/0x0109,
/* 2*032#2*/0x021A,
/* 3*016#3*/0x0B32,
/* 4*004#5*/0x0945,
/* 5*008#4*/0x0748,
/* 6*001#7*/0x07B6,
/* 7*002#6*/0x06A5,

// Case 2 (136)
/* 0*063#01*/0x1189,0x0138,
/* 1*096#12*/0x129A,0x0092,
/* 2*048#23*/0x1B3A,0x031A,
/* 3*111#03*/0x12B0,0x0B80,
/* 4*068#15*/0x1045,0x0105,
/* 5*012#45*/0x1975,0x0987,
/* 6*119#04*/0x1340,0x0374,
/* 7*003#67*/0x1BA5,0x07B5,
/* 8*009#47*/0x1486,0x0B68,
/* 9*034#26*/0x1615,0x0216,
/*10*017#37*/0x1726,0x0732,
/*11*006#56*/0x146A,0x094A,

// Case 3.1 (160)
/* 0*123#05*/0x1945,0x0038,
/* 1*072#14*/0x1109,0x0748,
/* 2*066#16*/0x16A5,0x0109,
/* 3*036#25*/0x1945,0x021A,
/* 4*018#36*/0x16A5,0x032B,
/* 5*033#27*/0x17B6,0x021A,
/* 6*126#07*/0x17B6,0x0038,
/* 7*024#34*/0x1B32,0x0748,
/* 8*095#02*/0x121A,0x0038,
/* 9*080#13*/0x1B32,0x0109,
/*10*010#46*/0x16A5,0x0874,
/*11*005#57*/0x1945,0x07B6,

// Case 3.2 (184)
/* 0*123#05*/0x1905,0x1035,0x1453,0x0843,
/* 1*072#14*/0x1974,0x1917,0x1871,0x0081,
/* 2*066#16*/0x1605,0x1950,0x116A,0x0106,
/* 3*036#25*/0x1A45,0x1942,0x124A,0x0192,
/* 4*018#36*/0x13A5,0x1B56,0x132A,0x0B35,
/* 5*033#27*/0x11A6,0x17B2,0x1217,0x0671,
/* 6*126#07*/0x1786,0x1806,0x1B60,0x03B0,
/* 7*024#34*/0x1834,0x1324,0x1742,0x0B72,
/* 8*095#02*/0x123A,0x138A,0x11A8,0x0018,
/* 9*080#13*/0x1129,0x12B9,0x109B,0x030B,
/*10*010#46*/0x14A5,0x1A86,0x148A,0x0768,
/*11*005#57*/0x1B65,0x1794,0x17B9,0x059B,

// Case 4.1.1 (232)
/* 0*125#06*/0x16A5,0x0380,
/* 1*065#17*/0x17B6,0x0109,
/* 2*040#24*/0x121A,0x0748,
/* 3*020#35*/0x1945,0x0B32,

//The numbers in parentheses are the diagonal (interior test)
// Case 4.1.2 (240)
/* 0*(06)*125#06*/0x10A5,0x1805,0x1863,0x1685,0x16A3,0x03A0,
/* 1*(17)*065#17*/0x1796,0x11B6,0x1169,0x1970,0x17B0,0x00B1,
/* 2*(24)*040#24*/0x174A,0x17A2,0x1872,0x1A41,0x1481,0x0182,
/* 3*(35)*020#35*/0x12B5,0x1B34,0x1943,0x15B4,0x1592,0x0293,

// Case 5 (264)
/* 0*112#123*/0x1B9A,0x1930,0x09B3,
/* 1*079#023*/0x180A,0x101A,0x0B8A,
/* 2*047#013*/0x189B,0x1B12,0x0B91,
/* 3*031#012*/0x128A,0x1382,0x09A8,
/* 4*038#256*/0x1246,0x1419,0x0421,
/* 5*011#467*/0x18A5,0x1485,0x0BA8,
/* 6*110#037*/0x1026,0x1067,0x0078,
/* 7*059#015*/0x1845,0x1538,0x0135,
/* 8*014#456*/0x176A,0x1A87,0x098A,
/* 9*035#267*/0x1715,0x11B2,0x017B,
/*10*076#145*/0x1175,0x1708,0x0710,
/*11*025#347*/0x1426,0x1832,0x0248,
/*12*070#156*/0x106A,0x1460,0x010A,
/*13*055#014*/0x1149,0x1174,0x0371,
/*14*103#034*/0x142B,0x174B,0x0024,
/*15*019#367*/0x12A5,0x1532,0x0735,
/*16*050#326*/0x1635,0x136B,0x0153,
/*17*098#126*/0x1695,0x1609,0x0206,
/*18*115#045*/0x1935,0x1390,0x0753,
/*19*118#047*/0x13B6,0x1360,0x0406,
/*20*007#576*/0x19BA,0x17B4,0x0B94,
/*21*049#237*/0x17A6,0x171A,0x0317,
/*22*100#125*/0x1A25,0x1245,0x0042,
/*23*013#457*/0x1965,0x19B6,0x08B9,

// Case 6.1.1 (336)
/* 0*121#056*/0x146A,0x1A94,0x0380,
/* 1*061#016*/0x16A5,0x1189,0x0138,
/* 2*109#036*/0x16A5,0x180B,0x002B,
/* 3*124#067*/0x1BA5,0x157B,0x0038,
/* 4*093#026*/0x1615,0x1621,0x0803,
/* 5*117#046*/0x16A5,0x1034,0x0374,
/* 6*073#147*/0x18B6,0x1648,0x0109,
/* 7*067#167*/0x1BA5,0x17B5,0x0109,
/* 8*097#127*/0x129A,0x17B6,0x0092,
/* 9*062#017*/0x17B6,0x1189,0x0138,
/*10*081#137*/0x1726,0x1732,0x0109,
/*11*069#157*/0x1045,0x17B6,0x0105,
/*12*104#124*/0x129A,0x1209,0x0748,
/*13*044#245*/0x1975,0x121A,0x0879,
/*14*041#247*/0x1486,0x121A,0x08B6,
/*15*056#234*/0x131A,0x1B3A,0x0874,
/*16*087#024*/0x121A,0x1374,0x0403,
/*17*042#246*/0x1615,0x1216,0x0874,
/*18*107#035*/0x1945,0x12B0,0x0B80,
/*19*052#235*/0x1945,0x1B3A,0x031A,
/*20*022#356*/0x146A,0x194A,0x032B,
/*21*028#345*/0x1975,0x1987,0x02B3,
/*22*084#135*/0x1045,0x1510,0x0B32,
/*23*021#357*/0x1945,0x1726,0x0732,

// Case 6.1.2 (408)
/* 0*(06)*121#056*/0x136A,0x190A,0x1094,0x1804,0x1684,0x1386,0x03A0,
/* 1*(06)*061#016*/0x1895,0x11A5,0x136A,0x1591,0x13A1,0x1863,0x0856,
/* 2*(06)*109#036*/0x10A5,0x1AB6,0x102A,0x1A2B,0x186B,0x1568,0x0058,
/* 3*(06)*124#067*/0x10A5,0x1785,0x187B,0x138B,0x1A3B,0x103A,0x0058,
/* 4*(06)*093#026*/0x1685,0x1236,0x1321,0x1031,0x1501,0x1805,0x0863,
/* 5*(06)*117#046*/0x10A5,0x1456,0x1376,0x1674,0x1054,0x13A0,0x0A36,
/* 6*(17)*073#147*/0x1496,0x1948,0x1098,0x1B08,0x110B,0x161B,0x0169,
/* 7*(17)*067#167*/0x1795,0x11A5,0x11BA,0x1915,0x1097,0x1B07,0x00B1,
/* 8*(17)*097#127*/0x19A6,0x16A2,0x1B62,0x10B2,0x17B0,0x1970,0x0796,
/* 9*(17)*062#017*/0x1796,0x113B,0x1B38,0x17B8,0x1978,0x1169,0x01B6,
/*10*(17)*081#137*/0x1796,0x1730,0x1032,0x1102,0x1612,0x1916,0x0970,
/*11*(17)*069#157*/0x1165,0x1745,0x1756,0x1704,0x1B61,0x10B1,0x0B07,
/*12*(24)*104#124*/0x149A,0x1208,0x1809,0x1489,0x174A,0x127A,0x0728,
/*13*(24)*044#245*/0x19A5,0x175A,0x11A9,0x1819,0x1218,0x1728,0x027A,
/*14*(24)*041#247*/0x14A6,0x18B2,0x12B6,0x1A26,0x11A4,0x1814,0x0182,
/*15*(24)*056#234*/0x141A,0x1B7A,0x17B3,0x1873,0x1183,0x1481,0x04A7,
/*16*(24)*087#024*/0x141A,0x1014,0x1103,0x1213,0x1723,0x1A27,0x04A7,
/*17*(24)*042#246*/0x1645,0x1415,0x1746,0x1276,0x1872,0x1182,0x0814,
/*18*(35)*107#035*/0x1925,0x1B84,0x1480,0x1940,0x1290,0x1B52,0x05B4,
/*19*(35)*052#235*/0x19A5,0x1319,0x191A,0x1B5A,0x145B,0x134B,0x0439,
/*20*(35)*022#356*/0x1B6A,0x1B46,0x12BA,0x192A,0x1329,0x1439,0x034B,
/*21*(35)*028#345*/0x12B5,0x1398,0x1387,0x1B37,0x15B7,0x1925,0x0293,
/*22*(35)*084#135*/0x1B45,0x1125,0x1210,0x1320,0x1430,0x1B34,0x0B52,
/*23*(35)*021#357*/0x1265,0x1574,0x1567,0x1347,0x1943,0x1293,0x0925,

// Case 6.2 (576)
/* 0*121#056*/0x136A,0x190A,0x13A0,0x1684,0x0386,
/* 1*061#016*/0x1685,0x136A,0x1895,0x113A,0x0863,
/* 2*109#036*/0x10A5,0x1856,0x102A,0x186B,0x0058,
/* 3*124#067*/0x10A5,0x1785,0x1058,0x1A3B,0x003A,
/* 4*093#026*/0x1685,0x1236,0x1863,0x1501,0x0805,
/* 5*117#046*/0x10A5,0x1A36,0x1054,0x1376,0x03A0,
/* 6*073#147*/0x1496,0x1169,0x1B08,0x110B,0x061B,
/* 7*067#167*/0x1795,0x11BA,0x10B1,0x1097,0x0B07,
/* 8*097#127*/0x19A6,0x1796,0x10B2,0x17B0,0x0970,
/* 9*062#017*/0x1796,0x113B,0x161B,0x1978,0x0169,
/*10*081#137*/0x1796,0x1126,0x1730,0x1970,0x0916,
/*11*069#157*/0x1165,0x1704,0x1B07,0x1B61,0x00B1,
/*12*104#124*/0x149A,0x1208,0x1728,0x174A,0x027A,
/*13*044#245*/0x1A75,0x127A,0x1819,0x1218,0x0728,
/*14*041#247*/0x14A6,0x18B2,0x1182,0x11A4,0x0814,
/*15*056#234*/0x174A,0x1B7A,0x1183,0x1481,0x0A41,
/*16*087#024*/0x141A,0x1014,0x1723,0x1A27,0x04A7,
/*17*042#246*/0x1415,0x1276,0x1814,0x1872,0x0182,
/*18*107#035*/0x1925,0x1B84,0x15B4,0x1290,0x0B52,
/*19*052#235*/0x1AB5,0x1319,0x1439,0x145B,0x034B,
/*20*022#356*/0x1B46,0x192A,0x134B,0x1329,0x0439,
/*21*028#345*/0x12B5,0x1398,0x1293,0x15B7,0x0925,
/*22*084#135*/0x12B5,0x1125,0x1430,0x1B34,0x05B4,
/*23*021#357*/0x1265,0x1925,0x1734,0x1943,0x0293,

// Case 7.1 (696)
/* 0*037#257*/0x1945,0x121A,0x07B6,
/* 1*088#134*/0x12B3,0x1874,0x0109,
/* 2*026#346*/0x16A5,0x1B32,0x0874,
/* 3*091#025*/0x1945,0x121A,0x0038,
/* 4*122#057*/0x1945,0x17B6,0x0038,
/* 5*082#136*/0x16A5,0x1B32,0x0109,
/* 6*074#146*/0x16A5,0x1109,0x0748,
/* 7*094#027*/0x121A,0x17B6,0x0380,

//The characters inside of the square bracket are face test results
// Case 7.2 (720)
/* 0*037#257*[.--..+]*/0x1B65,0x19B5,0x121A,0x1794,0x0B97,
/* 1*088#134*[-..-+.]*/0x1B92,0x1874,0x1B30,0x19B0,0x0912,
/* 2*026#346*[..--.+]*/0x14A5,0x1A86,0x1B32,0x1876,0x08A4,
/* 3*091#025*[--..+.]*/0x1945,0x138A,0x1801,0x1A81,0x0A23,
/* 4*122#057*[-..-.+]*/0x1B65,0x1380,0x1794,0x1B97,0x09B5,
/* 5*082#136*[.--.+.]*/0x16A5,0x1129,0x1B92,0x1B30,0x09B0,
/* 6*074#146*[--...+]*/0x14A5,0x1876,0x1109,0x18A4,0x0A86,
/* 7*094#027*[..--+.]*/0x17B6,0x123A,0x18A3,0x1801,0x0A81,
//(760)
/* 0*037#257*[.+-..-]*/0x1A45,0x17B6,0x124A,0x1219,0x0429,
/* 1*088#134*[-..+-.]*/0x1109,0x1483,0x1243,0x12B7,0x0427,
/* 2*026#346*[..-+.-]*/0x16A5,0x1274,0x12B7,0x1483,0x0243,
/* 3*091#025*[-+..-.]*/0x1A45,0x1038,0x1219,0x1429,0x024A,
/* 4*122#057*[-..+.-]*/0x1945,0x1806,0x1678,0x103B,0x060B,
/* 5*082#136*[.+-.-.]*/0x1605,0x1095,0x116A,0x132B,0x0061,
/* 6*074#146*[-+...-]*/0x1605,0x1095,0x116A,0x1748,0x0061,
/* 7*094#027*[..-+-.]*/0x1786,0x121A,0x103B,0x160B,0x0068,
//(800)
/* 0*037#257*[.-+..-]*/0x1945,0x11A6,0x1716,0x17B2,0x0172,
/* 1*088#134*[+..--.]*/0x132B,0x1749,0x1108,0x1718,0x0179,
/* 2*026#346*[..+-.-]*/0x13A5,0x1B56,0x1748,0x135B,0x032A,
/* 3*091#025*[+-..-.]*/0x1905,0x121A,0x1350,0x1384,0x0534,
/* 4*122#057*[+..-.-]*/0x1905,0x1534,0x17B6,0x1384,0x0350,
/* 5*082#136*[.-+.-.]*/0x13A5,0x1B56,0x1109,0x132A,0x035B,
/* 6*074#146*[+-...-]*/0x16A5,0x1749,0x1179,0x1108,0x0718,
/* 7*094#027*[..+--.]*/0x11A6,0x1380,0x17B2,0x1172,0x0716,

// Case 7.3 (840)
/* 0*037#257*[.+-..+]*/0x1C65,0x1C5A,0x17C4,0x1C7B,0x1C94,0x1C19,0x1C21,0x1CA2,0x0CB6,
/* 1*088#134*[-..++.]*/0x1C74,0x1CB7,0x1C2B,0x11C9,0x1C12,0x1C09,0x1C30,0x1C83,0x0C48,
/* 2*026#346*[..-+.+]*/0x1CA5,0x1C6A,0x1C32,0x1C83,0x1C48,0x1C54,0x1C76,0x1CB7,0x0C2B,
/* 3*091#025*[-+..+.]*/0x1AC5,0x1C38,0x1C23,0x1CA2,0x1C45,0x1C94,0x1C19,0x1C01,0x0C80,
/* 4*122#057*[-..+.+]*/0x1C65,0x19C5,0x1CB6,0x1C3B,0x1C03,0x1C80,0x17C4,0x1C78,0x0C94,
/* 5*082#136*[.+-.+.]*/0x16C5,0x1C6A,0x1C95,0x1C09,0x1C30,0x1CB3,0x1C2B,0x1C12,0x0CA1,
/* 6*074#146*[-+...+]*/0x1C95,0x1C6A,0x1C10,0x1CA1,0x1C48,0x1C76,0x1C87,0x1C54,0x0C09,
/* 7*094#027*[..-++.]*/0x1C1A,0x17C6,0x1C01,0x1C80,0x1C78,0x1CB6,0x1C3B,0x1C23,0x0CA2,
//(912)
/* 0*037#257*[.-+..+]*/0x1C65,0x1CA6,0x1C59,0x11C2,0x1CB2,0x17C4,0x1C7B,0x1C94,0x0C1A,
/* 1*088#134*[+..-+.]*/0x1C2B,0x11C9,0x1C12,0x1C49,0x1C74,0x1C87,0x1C08,0x1C30,0x0CB3,
/* 2*026#346*[..+-.+]*/0x1CA5,0x1C54,0x1C76,0x1C48,0x1C2A,0x1C32,0x1CB3,0x1C6B,0x0C87,
/* 3*091#025*[+-..+.]*/0x1C45,0x12CA,0x1C84,0x1C38,0x1C23,0x1C59,0x1C1A,0x1C01,0x0C90,
/* 4*122#057*[+..-.+]*/0x1C65,0x1C03,0x1C90,0x1C59,0x1CB6,0x17C4,0x1C7B,0x1C84,0x0C38,
/* 5*082#136*[.-+.+.]*/0x1CA5,0x1C56,0x1C09,0x1C30,0x1CB3,0x1C6B,0x1C2A,0x1C12,0x0C91,
/* 6*074#146*[+-...+]*/0x1CA5,0x1C6A,0x1C54,0x1C76,0x1C87,0x1C08,0x11C9,0x1C10,0x0C49,
/* 7*094#027*[..+-+.]*/0x1CA6,0x17C6,0x1C1A,0x1C01,0x1C80,0x1C38,0x1C23,0x1CB2,0x0C7B,
//(984)
/* 0*037#257*[.++..-]*/0x1AC5,0x1CA6,0x1C94,0x1C19,0x1C21,0x1CB2,0x1C7B,0x1C67,0x0C45,
/* 1*088#134*[+..+-.]*/0x1C2B,0x11C9,0x1C49,0x1C74,0x1CB7,0x1C32,0x1C83,0x1C08,0x0C10,
/* 2*026#346*[..++.-]*/0x1CA5,0x1C56,0x1C2A,0x1C32,0x1C83,0x1C48,0x1C74,0x1CB7,0x0C6B,
/* 3*091#025*[++..-.]*/0x1AC5,0x12CA,0x1C45,0x1C84,0x1C38,0x1C03,0x1C90,0x1C19,0x0C21,
/* 4*122#057*[+..+.-]*/0x1C45,0x1CB6,0x1C3B,0x1C03,0x1C90,0x1C59,0x1C84,0x1C78,0x0C67,
/* 5*082#136*[.++.-.]*/0x16C5,0x1C2A,0x13CB,0x1C6B,0x1C95,0x1C09,0x1C10,0x1CA1,0x0C32,
/* 6*074#146*[++...-]*/0x16C5,0x1C6A,0x1C95,0x1C74,0x17C8,0x1C08,0x1C10,0x1CA1,0x0C49,
/* 7*094#027*[..++-.]*/0x1CA6,0x1C80,0x1C78,0x1C67,0x1C1A,0x1C21,0x1CB2,0x1C3B,0x0C03,

// Case 7.4.1 (1056)
/* 0*037#257*/0x1A65,0x1419,0x11B2,0x17B4,0x04B1,
/* 1*088#134*/0x174B,0x1149,0x12B1,0x11B4,0x0830,
/* 2*026#346*/0x12A5,0x1485,0x1B76,0x1832,0x0582,
/* 3*091#025*/0x1A25,0x1238,0x1458,0x1528,0x0019,
/* 4*122#057*/0x1965,0x1390,0x1B63,0x1693,0x0784,
/* 5*082#136*/0x1695,0x112A,0x1093,0x1B36,0x0396,
/* 6*074#146*/0x1495,0x176A,0x110A,0x1870,0x07A0,
/* 7*094#027*/0x17A6,0x1780,0x11A0,0x1A70,0x03B2,

// Case 7.4.2 (1096)
/* 0*(06)*037#257*/0x1465,0x1459,0x121A,0x1AB2,0x1BA6,0x17B6,0x1476,0x15A1,0x0951,
/* 1*(06)*088#134*/0x1084,0x1748,0x18B7,0x1B83,0x12B3,0x1109,0x1130,0x1123,0x0904,
/* 2*(17)*026#346*/0x16A5,0x132B,0x1B83,0x1874,0x18B7,0x1576,0x1547,0x16B2,0x0A62,
/* 3*(17)*091#025*/0x1945,0x115A,0x1038,0x1023,0x1201,0x1A21,0x1519,0x1908,0x0498,
/* 4*(24)*122#057*/0x1465,0x1380,0x1890,0x1984,0x1594,0x1647,0x1B67,0x1783,0x0B73,
/* 5*(24)*082#136*/0x16A5,0x1A95,0x19A1,0x1091,0x1312,0x1301,0x1B32,0x12A6,0x0B26,
/* 6*(35)*074#146*/0x16A5,0x1109,0x19A1,0x1A95,0x1754,0x1765,0x1874,0x1490,0x0840,
/* 7*(35)*094#027*/0x1BA6,0x1380,0x1378,0x173B,0x167B,0x1AB2,0x11A2,0x1230,0x0120,

// Case 8 (1168)
/* 0*015#0123*/0x1B9A,0x09B8,
/* 1*102#0347*/0x1426,0x0024,
/* 2*051#0145*/0x1375,0x0135,

// Case 9 (1174)
/* 0*078#0237*/0x17A6,0x1180,0x1781,0x0A71,
/* 1*039#0134*/0x1B42,0x1129,0x1492,0x04B7,
/* 2*027#0125*/0x1A35,0x13A2,0x1853,0x0584,
/* 3*114#0457*/0x1965,0x13B0,0x10B6,0x0069,

// Case 10.1.1 (1190)
/* 0*105#0356*[-.-...]*/0x146A,0x194A,0x1028,0x082B,
/* 1*060#0167*[.-.-..]*/0x17A5,0x1189,0x1381,0x07BA,
/* 2*085#0246*[....--]*/0x1625,0x1340,0x1743,0x0521,
//(1202)
/* 0*105#0356*[+.+...]*/0x1846,0x190A,0x186B,0x002A,
/* 1*060#0167*[.+.+..]*/0x1795,0x11BA,0x1789,0x03B1,
/* 2*085#0246*[....++]*/0x1015,0x1236,0x1540,0x0763,

// Case 10.1.2 (1214)
/* 0*(06)(35)*105#0356*[-.-...]*/0x126A,0x1029,0x12A9,0x1B62,0x16B8,0x1468,0x1489,0x0809,
/* 1*(06)(17)*060#0167*[.-.-..]*/0x19A5,0x1789,0x1957,0x11A9,0x1A13,0x1BA3,0x1B37,0x0387,
/* 2*(06)(24)*085#0246*[....--]*/0x1405,0x1756,0x1015,0x1021,0x1320,0x1237,0x1627,0x0745,
//(1238)
/* 0*(17)(24)*105#0356*[+.+...]*/0x146A,0x1A94,0x1904,0x1408,0x1B80,0x1B02,0x1AB2,0x0A6B,
/* 1*(35)(24)*060#0167*[.+.+..]*/0x1BA5,0x1189,0x1138,0x13B8,0x18B7,0x157B,0x115A,0x0195,
/* 2*(17)(35)*085#0246*[....++]*/0x1465,0x1156,0x1621,0x1231,0x1130,0x1403,0x1437,0x0647,

// Case 10.2 (1262)
/* 0*105#0356*[+.-...]*/0x19CA,0x1AC6,0x190C,0x102C,0x12BC,0x18CB,0x184C,0x046C,
/* 1*060#0167*[.+.-..]*/0x1C95,0x1CBA,0x1C7B,0x1C57,0x19C8,0x1C38,0x1C13,0x01CA,
/* 2*085#0246*[....-+]*/0x1C15,0x12C6,0x1C21,0x14C5,0x1C76,0x17C3,0x1C03,0x0C40,
//(1286)
/* 0*105#0356*[-.+...]*/0x1C46,0x1C2A,0x1C94,0x1CA9,0x12C0,0x1C80,0x1CB8,0x0BC6,
/* 1*060#0167*[.-.+..]*/0x1CA5,0x17C5,0x1CBA,0x1C3B,0x11C9,0x13C1,0x1C89,0x08C7,
/* 2*085#0246*[....+-]*/0x16C5,0x12C6,0x123C,0x174C,0x137C,0x10C4,0x101C,0x015C,

// Case 11 (1310)
/* 0*077#0236*/0x16B5,0x1B80,0x15B0,0x0150,
/* 1*046#0137*/0x1786,0x1189,0x1126,0x0681,
/* 2*023#0124*/0x129A,0x1974,0x1792,0x0372,
/* 3*116#0467*/0x14A5,0x13BA,0x13A4,0x0340,
/* 4*099#0345*/0x1975,0x1902,0x1927,0x072B,
/* 5*057#0156*/0x116A,0x1846,0x1861,0x0813,

// Case 14 (1334)
/* 0*113#0456*/0x176A,0x1A90,0x17A0,0x0370,
/* 1*071#0234*/0x1B1A,0x1140,0x11B4,0x074B,
/* 2*043#0135*/0x1125,0x1285,0x12B8,0x0458,
/* 3*029#0126*/0x1695,0x1236,0x1396,0x0389,
/* 4*054#0147*/0x1496,0x1139,0x1369,0x063B,
/* 5*108#0367*/0x12A5,0x1785,0x1825,0x0802,

// Case 12.1.1 (1358)
/* 0*089#0256*[-...-.]*/0x1246,0x1192,0x1429,0x0038,
/* 1*075#0235*[--....]*/0x1945,0x181A,0x1180,0x0B8A,
/* 2*045#0136*[.--...]*/0x16A5,0x1129,0x1B92,0x089B,
/* 3*053#0146*[.-...-]*/0x16A5,0x1749,0x1179,0x0371,
/* 4*030#0127*[..--..]*/0x123A,0x17B6,0x18A3,0x09A8,
/* 5*101#0346*[..-..-]*/0x16A5,0x1274,0x172B,0x0024,
/* 6*092#0267*[...--.]*/0x1715,0x17B2,0x1172,0x0380,
/* 7*120#0567*[-..-..]*/0x19BA,0x1794,0x1B97,0x0803,
/* 8*086#0247*[..-.-.]*/0x1406,0x121A,0x103B,0x060B,
/* 9*083#0245*[.-..-.]*/0x1905,0x121A,0x1350,0x0753,
/*10*058#0157*[...-.-]*/0x1345,0x17B6,0x1384,0x0135,
/*11*106#0357*[-....-]*/0x1945,0x1786,0x1068,0x0260,
//(1406)
/* 0*089#0256*[+...+.]*/0x1246,0x1384,0x1342,0x0190,
/* 1*075#0235*[++....]*/0x1A45,0x14A8,0x18AB,0x0019,
/* 2*045#0136*[.++...]*/0x16B5,0x15B9,0x112A,0x09B8,
/* 3*053#0146*[.+...+]*/0x1495,0x116A,0x1617,0x0713,
/* 4*030#0127*[..++..]*/0x1786,0x168A,0x1A89,0x023B,
/* 5*101#0346*[..+..+]*/0x14A5,0x1B76,0x1A42,0x0240,
/* 6*092#0267*[...++.]*/0x1715,0x1180,0x1817,0x0B23,
/* 7*120#0567*[+..+..]*/0x19BA,0x13B0,0x10B9,0x0784,
/* 8*086#0247*[..+.+.]*/0x11A6,0x1160,0x1064,0x03B2,
/* 9*083#0245*[.+..+.]*/0x1A35,0x123A,0x1537,0x0019,
/*10*058#0157*[...+.+]*/0x1B65,0x1B53,0x1351,0x0784,
/*11*106#0357*[+....+]*/0x1065,0x1905,0x1602,0x0784,

// Case 12.1.2 (1454)
/* 0*(06)*089#0256*[-...-.]*/0x1846,0x1948,0x1980,0x1901,0x1863,0x1362,0x1132,0x0103,
/* 1*(35)*075#0235*[--....]*/0x1AB5,0x1159,0x11A5,0x1190,0x15B4,0x14B8,0x1048,0x0094,
/* 2*(06)*045#0136*[.--...]*/0x1685,0x126A,0x1589,0x11A5,0x12B6,0x16B8,0x12A1,0x0159,
/* 3*(06)*053#0146*[.-...-]*/0x19A5,0x1A36,0x191A,0x1A13,0x1954,0x1637,0x1467,0x0456,
/* 4*(17)*030#0127*[..--..]*/0x126A,0x1738,0x137B,0x1789,0x13B2,0x1796,0x169A,0x02B6,
/* 5*(06)*101#0346*[..-..-]*/0x10A5,0x1B6A,0x1745,0x1756,0x1540,0x176B,0x1A02,0x0BA2,
/* 6*(06)*092#0267*[...--.]*/0x1015,0x1102,0x1203,0x123B,0x1058,0x1857,0x1B87,0x0B38,
/* 7*(06)*120#0567*[-..-..]*/0x13BA,0x1784,0x137B,0x1738,0x13A0,0x10A9,0x1409,0x0480,
/* 8*(24)*086#0247*[..-.-.]*/0x1B6A,0x1BA2,0x1A64,0x1B23,0x1A41,0x1140,0x1310,0x0321,
/* 9*(24)*083#0245*[.-..-.]*/0x19A5,0x1320,0x1021,0x1237,0x1019,0x127A,0x1A75,0x091A,
/*10*(17)*058#0157*[...-.-]*/0x1165,0x1456,0x1467,0x1478,0x161B,0x1B13,0x18B3,0x087B,
/*11*(35)*106#0357*[-....-]*/0x1265,0x1809,0x1894,0x1902,0x1847,0x1925,0x1756,0x0745,
//(1550)
/* 0*(17)*089#0256*[+...+.]*/0x1946,0x1312,0x1013,0x1621,0x1803,0x1961,0x1498,0x0908,
/* 1*(17)*075#0235*[++....]*/0x1945,0x115A,0x1084,0x1904,0x1B80,0x11B0,0x1AB1,0x0195,
/* 2*(24)*045#0136*[.++...]*/0x16A5,0x1195,0x1A15,0x1891,0x1281,0x1B82,0x1B26,0x02A6,
/* 3*(35)*053#0146*[.+...+]*/0x16A5,0x1764,0x1546,0x1374,0x1934,0x1139,0x119A,0x095A,
/* 4*(35)*030#0127*[..++..]*/0x12A6,0x1B26,0x19A2,0x17B6,0x1392,0x1893,0x1837,0x03B7,
/* 5*(17)*101#0346*[..+..+]*/0x16A5,0x1B2A,0x16BA,0x102B,0x1754,0x170B,0x1407,0x0765,
/* 6*(35)*092#0267*[...++.]*/0x1B25,0x178B,0x13B8,0x157B,0x1038,0x1152,0x1120,0x0230,
/* 7*(24)*120#0567*[+..+..]*/0x194A,0x1490,0x1840,0x1380,0x17A4,0x1BA7,0x1B73,0x0783,
/* 8*(35)*086#0247*[..+.+.]*/0x1BA6,0x1130,0x1231,0x1403,0x1A21,0x1B43,0x164B,0x0B2A,
/* 9*(17)*083#0245*[.+..+.]*/0x1A95,0x119A,0x1759,0x121A,0x1079,0x1370,0x1302,0x0012,
/*10*(24)*058#0157*[...+.+]*/0x1465,0x13B8,0x178B,0x1138,0x167B,0x1418,0x1514,0x0476,
/*11*(24)*106#0357*[+....+]*/0x1945,0x1765,0x1754,0x1267,0x1827,0x1028,0x1089,0x0849,

// Case 12.2 (1646)
/* 0*089#0256*[-...+.]*/0x12C6,0x1C19,0x11C0,0x13C2,0x180C,0x18C3,0x194C,0x046C,
/* 1*075#0235*[-+....]*/0x1AC5,0x119C,0x14C9,0x145C,0x1ABC,0x10C8,0x1B8C,0x01C0,
/* 2*045#0136*[.-+...]*/0x1CA5,0x156C,0x12AC,0x16BC,0x1B8C,0x11C9,0x189C,0x02C1,
/* 3*053#0146*[.-...+]*/0x1CA5,0x1AC6,0x14C5,0x16C7,0x137C,0x11C9,0x113C,0x09C4,
/* 4*030#0127*[..-+..]*/0x12CA,0x1CB6,0x13BC,0x178C,0x167C,0x189C,0x19AC,0x03C2,
/* 5*101#0346*[..-..+]*/0x1CA5,0x154C,0x176C,0x1AC6,0x140C,0x1BC2,0x102C,0x07CB,
/* 6*092#0267*[...-+.]*/0x1C15,0x123C,0x101C,0x18C3,0x180C,0x1BC7,0x157C,0x02CB,
/* 7*120#0567*[+..-..]*/0x1CBA,0x1C38,0x17C4,0x1BC7,0x1C84,0x1C90,0x1CA9,0x03C0,
/* 8*086#0247*[..+.-.]*/0x16CA,0x11C2,0x1C03,0x12CB,0x1C3B,0x1C40,0x1C64,0x0AC1,
/* 9*083#0245*[.+..-.]*/0x1AC5,0x1C19,0x11C2,0x13C0,0x10C9,0x1C37,0x1C75,0x02CA,
/*10*058#0157*[...+.-]*/0x1C45,0x17C6,0x178C,0x1BC3,0x16CB,0x113C,0x151C,0x04C8,
/*11*106#0357*[+....-]*/0x1C45,0x1C26,0x184C,0x190C,0x159C,0x102C,0x17C6,0x08C7,
//(1742)
/* 0*089#0256*[+...-.]*/0x146C,0x190C,0x184C,0x13C0,0x138C,0x11C2,0x162C,0x09C1,
/* 1*075#0235*[+-....]*/0x1C45,0x10C9,0x14C8,0x159C,0x1B8C,0x11AC,0x1ABC,0x01C0,
/* 2*045#0136*[.+-...]*/0x16C5,0x16AC,0x15C9,0x11CA,0x189C,0x12BC,0x1B8C,0x02C1,
/* 3*053#0146*[.+...-]*/0x16C5,0x16AC,0x195C,0x1A1C,0x113C,0x1C74,0x137C,0x09C4,
/* 4*030#0127*[..+-..]*/0x16CA,0x12CB,0x17BC,0x17C6,0x19AC,0x138C,0x189C,0x03C2,
/* 5*101#0346*[..+..-]*/0x16C5,0x15CA,0x1BC6,0x1AC2,0x102C,0x174C,0x140C,0x07CB,
/* 6*092#0267*[...+-.]*/0x17C5,0x1C3B,0x18C7,0x103C,0x10C8,0x121C,0x115C,0x02CB,
/* 7*120#0567*[-..+..]*/0x19CA,0x1C80,0x1C94,0x18C7,0x1C47,0x1BC3,0x1CBA,0x03C0,
/* 8*086#0247*[..-.+.]*/0x12CA,0x1CB6,0x1C23,0x1BC3,0x1C64,0x1C01,0x1C40,0x0AC1,
/* 9*083#0245*[.-..+.]*/0x19C5,0x1C1A,0x11C0,0x1C90,0x1C75,0x13C2,0x1C37,0x02CA,
/*10*058#0157*[...-.+]*/0x1C65,0x17C4,0x1BC7,0x1B6C,0x151C,0x18C3,0x113C,0x04C8,
/*11*106#0357*[-....+]*/0x1C65,0x17C4,0x194C,0x19C5,0x126C,0x180C,0x102C,0x08C7,

// Case 13.1 (1838)
/* 0*090#0257*[------]*/0x1945,0x121A,0x17B6,0x0380,
/* 1*090#0257*[++++++]*/0x1A65,0x1190,0x1B23,0x0784,

// Case 13.2 (1846)
/* 0*090#0257*[-----+]*/0x1B65,0x1B59,0x121A,0x1380,0x1794,0x0B97,
/* 1*090#0257*[----+-]*/0x1945,0x17B6,0x181A,0x1801,0x1A38,0x0A23,
/* 2*090#0257*[---+--]*/0x1945,0x121A,0x10B6,0x103B,0x1680,0x0678,
/* 3*090#0257*[--+---]*/0x1945,0x11A6,0x1038,0x1172,0x17B2,0x0167,
/* 4*090#0257*[-+----]*/0x1A45,0x17B6,0x1380,0x1429,0x1219,0x04A2,
/* 5*090#0257*[-+++++]*/0x1A65,0x1794,0x1197,0x13B2,0x1817,0x0801,
/* 6*090#0257*[+-----]*/0x1905,0x121A,0x1345,0x17B6,0x1350,0x0384,
/* 7*090#0257*[+-++++]*/0x1065,0x1590,0x11A6,0x1784,0x1B23,0x0016,
/* 8*090#0257*[++-+++]*/0x1B65,0x135A,0x1784,0x1019,0x1B53,0x0A23,
/* 9*090#0257*[+++-++]*/0x1A65,0x1190,0x1342,0x1384,0x1472,0x07B2,
/*10*090#0257*[++++-+]*/0x1A65,0x1784,0x10B9,0x103B,0x1B29,0x0219,
/*11*090#0257*[+++++-]*/0x1A45,0x18A6,0x1190,0x13B2,0x14A8,0x0678,

// Case 13.3 (1918)
/* 0*090#0257*[---+-+]*/0x1C65,0x1C59,0x121A,0x14C9,0x10C8,0x1C78,0x1C47,0x1C3B,0x16CB,0x0C03,
/* 1*090#0257*[--++--]*/0x1945,0x1CA6,0x13C0,0x11C2,0x1CB2,0x1C3B,0x1C80,0x1C78,0x17C6,0x0C1A,
/* 2*090#0257*[--+--+]*/0x1C65,0x1CA6,0x19C5,0x11AC,0x1C21,0x1CB2,0x17C4,0x1BC7,0x1C94,0x0803,
/* 3*090#0257*[---++-]*/0x1945,0x12CA,0x1CB6,0x1C3B,0x1C23,0x1C1A,0x1C01,0x1C78,0x10C8,0x0C67,
/* 4*090#0257*[--++++]*/0x1C65,0x1A6C,0x17C4,0x159C,0x194C,0x178C,0x180C,0x11AC,0x11C0,0x03B2,
/* 5*090#0257*[--+-+-]*/0x1945,0x1CA6,0x17BC,0x18C3,0x1C23,0x1CB2,0x1C67,0x1C01,0x1AC1,0x0C80,
/* 6*090#0257*[-+---+]*/0x1C65,0x1C5A,0x1CB6,0x12CA,0x17C4,0x1C7B,0x1C19,0x14C9,0x1C21,0x0038,
/* 7*090#0257*[-++---]*/0x1AC5,0x1CA6,0x1C45,0x17C6,0x1C94,0x1C19,0x1CB2,0x11C2,0x1C7B,0x0380,
/* 8*090#0257*[-+++-+]*/0x1A65,0x1C3B,0x17C4,0x18C7,0x180C,0x103C,0x1B2C,0x1C19,0x121C,0x094C,
/* 9*090#0257*[-+--+-]*/0x1AC5,0x1C45,0x17B6,0x10C8,0x14C9,0x1C19,0x1C01,0x1C38,0x1C23,0x02CA,
/*10*090#0257*[-++-++]*/0x1A65,0x119C,0x11C0,0x13C2,0x138C,0x180C,0x194C,0x17BC,0x17C4,0x0B2C,
/*11*090#0257*[-++++-]*/0x1AC5,0x1A6C,0x1C19,0x194C,0x145C,0x167C,0x180C,0x18C7,0x101C,0x0B23,
/*12*090#0257*[+-++-+]*/0x1C65,0x16CA,0x12CB,0x159C,0x121C,0x11AC,0x103C,0x10C9,0x13BC,0x0784,
/*13*090#0257*[+--+--]*/0x1C45,0x17C6,0x1C59,0x121A,0x1C84,0x1C78,0x1CB6,0x1C3B,0x1C90,0x03C0,
/*14*090#0257*[+----+]*/0x1C65,0x19C5,0x121A,0x1C38,0x17C4,0x1BC7,0x1C84,0x1C03,0x1C90,0x0CB6,
/*15*090#0257*[+-+++-]*/0x1C45,0x16CA,0x10C9,0x14C8,0x159C,0x101C,0x11AC,0x167C,0x178C,0x023B,
/*16*090#0257*[+--+++]*/0x1C65,0x12CA,0x13C2,0x159C,0x11C0,0x11AC,0x13BC,0x1B6C,0x190C,0x0784,
/*17*090#0257*[+---+-]*/0x19C5,0x12CA,0x1C45,0x17B6,0x1AC1,0x1C01,0x1C90,0x1C84,0x1C23,0x08C3,
/*28*090#0257*[+++--+]*/0x1A65,0x14C8,0x10C9,0x103C,0x138C,0x147C,0x17BC,0x121C,0x12CB,0x019C,
/*19*090#0257*[++----]*/0x1AC5,0x15C4,0x17B6,0x1C19,0x11C2,0x13C0,0x1C90,0x1CA2,0x1C84,0x0C38,
/*20*090#0257*[++-+-+]*/0x1AC5,0x1B6C,0x1C19,0x1A2C,0x121C,0x190C,0x103C,0x1BC3,0x165C,0x0784,
/*21*090#0257*[+++-+-]*/0x1AC5,0x16CA,0x12CB,0x145C,0x167C,0x17BC,0x123C,0x138C,0x14C8,0x0019,
/*22*090#0257*[++--++]*/0x1C65,0x15AC,0x17C4,0x17BC,0x1B6C,0x1A2C,0x138C,0x13C2,0x184C,0x0190,
/*23*090#0257*[++-++-]*/0x1AC5,0x1B6C,0x145C,0x1C78,0x1BC3,0x167C,0x184C,0x1A2C,0x123C,0x0019,

// Case 13.4 (2158)
/* 0*090#0257*[++---+]*/0x1C65,0x1C5A,0x1C90,0x1C03,0x1C38,0x1C84,0x1C47,0x1C7B,0x1CB6,0x1CA2,0x1C21,0x0C19,
/* 1*090#0257*[-++-+-]*/0x1AC5,0x1A6C,0x138C,0x123C,0x1B2C,0x145C,0x17BC,0x167C,0x194C,0x119C,0x101C,0x080C,
/* 2*090#0257*[--++-+]*/0x1C65,0x1CA6,0x1CB2,0x1C59,0x1C21,0x1C1A,0x1C94,0x1C47,0x1C78,0x1C80,0x1C03,0x0C3B,
/* 3*090#0257*[+--++-]*/0x1C45,0x1B6C,0x159C,0x11AC,0x101C,0x190C,0x184C,0x178C,0x167C,0x13BC,0x123C,0x0A2C,

// Case 13.5.2 (2206)
/* 0*(06)*090#0257*[-++--+]*/0x1A65,0x1784,0x1B87,0x18B3,0x1804,0x1094,0x1190,0x1310,0x1213,0x0B23,
/* 1*(17)*090#0257*[++--+-]*/0x1945,0x115A,0x17B6,0x1380,0x1302,0x1012,0x1908,0x1498,0x1951,0x0A21,
/* 2*(24)*090#0257*[+--+-+]*/0x1A65,0x19A5,0x1A91,0x1A26,0x12B6,0x13B2,0x1132,0x1031,0x1901,0x0784,
/* 3*(35)*090#0257*[--+++-]*/0x1945,0x1BA6,0x1380,0x1837,0x13B7,0x1230,0x1120,0x121A,0x12AB,0x067B,
//(2246)
/* 0*(06)*090#0257*[-++--+]*/0x1465,0x1BA6,0x1519,0x121A,0x12AB,0x15A1,0x1594,0x1476,0x17B6,0x0803,
/* 1*(17)*090#0257*[++--+-]*/0x1A65,0x13B2,0x18B3,0x1574,0x1B87,0x1B62,0x16A2,0x1756,0x1847,0x0190,
/* 2*(24)*090#0257*[+--+-+]*/0x1465,0x1459,0x121A,0x1380,0x1089,0x1849,0x1783,0x1B73,0x17B6,0x0764,
/* 3*(35)*090#0257*[--+++-]*/0x1A65,0x1190,0x1A91,0x19A5,0x1940,0x1480,0x1784,0x1574,0x1675,0x03B2,

// Case 13.5.1 (2286)
/* 0*090#0257*[-++--+]*/0x1A65,0x1380,0x1219,0x1942,0x14B2,0x07B4,
/* 1*090#0257*[++--+-]*/0x1A35,0x1584,0x17B6,0x1190,0x123A,0x0385,
/* 2*090#0257*[+--+-+]*/0x1965,0x121A,0x1B03,0x1B60,0x1690,0x0784,
/* 3*090#0257*[--+++-]*/0x1945,0x17A6,0x13B2,0x1018,0x1178,0x01A7
};

/***************** Marching cubes 33 functions ****************/

/******************************************************************
Vertices:           Faces:
    3 __________2        ___________
   /|          /|      /|          /|
  / |         / |     / |   2     / |
7/__________6/  |    /  |     4  /  |
|   |       |   |   |¯¯¯¯¯¯¯¯¯¯¯| 1 |     z
|   0_______|___1   | 3 |_______|___|     |
|  /        |  /    |  /  5     |  /      |____y
| /         | /     | /     0   | /      /
4/__________5/      |/__________|/      x


This function return a vector with all six test face results (face[6]). Each
result value is 1 if the positive face vertices are joined, -1 if the negative
vertices are joined, and 0 (unchanged) if the test must no be applied. The
return value of this function is the the sum of all six results.*/
int MC33_faceTests(int *face, int ind, const float *v) {
	if (ind&0x80)//vertex 0
	{
		face[0] = ((ind&0xCC) == 0x84? (v[0]*v[5] < v[1]*v[4]? -1: 1): 0);//0x84 = 10000100, vertices 0 and 5
		face[3] = ((ind&0x99) == 0x81? (v[0]*v[7] < v[3]*v[4]? -1: 1): 0);//0x81 = 10000001, vertices 0 and 7
		face[4] = ((ind&0xF0) == 0xA0? (v[0]*v[2] < v[1]*v[3]? -1: 1): 0);//0xA0 = 10100000, vertices 0 and 2
	}
	else
	{
		face[0] = ((ind&0xCC) == 0x48? (v[0]*v[5] < v[1]*v[4]? 1: -1): 0);//0x48 = 01001000, vertices 1 and 4
		face[3] = ((ind&0x99) == 0x18? (v[0]*v[7] < v[3]*v[4]? 1: -1): 0);//0x18 = 00011000, vertices 3 and 4
		face[4] = ((ind&0xF0) == 0x50? (v[0]*v[2] < v[1]*v[3]? 1: -1): 0);//0x50 = 01010000, vertices 1 and 3
	}
	if (ind&0x02)//vertex 6
	{
		face[1] = ((ind&0x66) == 0x42? (v[1]*v[6] < v[2]*v[5]? -1: 1): 0);//0x42 = 01000010, vertices 1 and 6
		face[2] = ((ind&0x33) == 0x12? (v[3]*v[6] < v[2]*v[7]? -1: 1): 0);//0x12 = 00010010, vertices 3 and 6
		face[5] = ((ind&0x0F) == 0x0A? (v[4]*v[6] < v[5]*v[7]? -1: 1): 0);//0x0A = 00001010, vertices 4 and 6
	}
	else
	{
		face[1] = ((ind&0x66) == 0x24? (v[1]*v[6] < v[2]*v[5]? 1: -1): 0);//0x24 = 00100100, vertices 2 and 5
		face[2] = ((ind&0x33) == 0x21? (v[3]*v[6] < v[2]*v[7]? 1: -1): 0);//0x21 = 00100001, vertices 2 and 7
		face[5] = ((ind&0x0F) == 0x05? (v[4]*v[6] < v[5]*v[7]? 1: -1): 0);//0x05 = 00000101, vertices 5 and 7
	}
	return face[0] + face[1] + face[2] + face[3] + face[4] + face[5];
}

/* Faster function for the face test, the test is applied to only one face
(int face). This function is only used for the cases 3 and 6 of MC33*/
int MC33_faceTest1(int face, const float *v) {
	switch (face) {
	case 0:
		return (v[0]*v[5] < v[1]*v[4]? 0x48: 0x84);
	case 1:
		return (v[1]*v[6] < v[2]*v[5]? 0x24: 0x42);
	case 2:
		return (v[3]*v[6] < v[2]*v[7]? 0x21: 0x12);
	case 3:
		return (v[0]*v[7] < v[3]*v[4]? 0x18: 0x81);
	case 4:
		return (v[0]*v[2] < v[1]*v[3]? 0x50: 0xA0);
	default:
		return (v[4]*v[6] < v[5]*v[7]? 0x05: 0x0A);
	}
}

// an ugly signbit for float type with
// warning: dereferencing type-punned pointer will break strict-aliasing rules
// Silence dereferencing type-punned pointer warning in GCC
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

inline unsigned int signbf(float x) {
	return ((*(unsigned int*)(void*)(&x))&0x80000000);
}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/******************************************************************
Interior test function. If the test is positive, the function returns a value
different from 0. The integer i must be 0 to test if the vertices 0 and 6 are
joined. 1 for vertices 1 and 7, 2 for vertices 2 and 4, and 3 for 3 and 5.
For case 13, the integer flag13 must be 1, and the function returns 2 if one
of the vertices 0, 1, 2 or 3 is joined to the center point of the cube (case
13.5.2), returns 1 if one of the vertices 4, 5, 6 or 7 is joined to the
center point of the cube (case 13.5.2 too), and it returns 0 if the vertices
are no joined (case 13.5.1)*/
int MC33_interiorTest(int i, int flag13, const float *v) {
	//Signs of cube vertices were changed to use signbit function in calc_isosurface
	//A0 = -v[0], B0 = -v[1], C0 = -v[2], D0 = -v[3]
	//A1 = -v[4], B1 = -v[5], C1 = -v[6], D1 = -v[7]
	//But the function still works
	float At = v[4] - v[0], Bt = v[5] - v[1], Ct = v[6] - v[2], Dt = v[7] - v[3];
	float t = At*Ct - Bt*Dt; // the "a" value.
	if (signbf(t)) {
		if (i&0x01)
			return 0;
	} else {
		if (!(i&0x01) || t == 0)
			return 0;
	}
	t = 0.5f*(v[3]*Bt - v[2]*At + v[1]*Dt - v[0]*Ct)/t; // t = -b/2a
	if (t > 0 && t < 1) {
		At = v[0] + At*t;
		Bt = v[1] + Bt*t;
		Ct = v[2] + Ct*t;
		Dt = v[3] + Dt*t;
		Ct *= At;
		Dt *= Bt;
		if (i&0x01) {
			if (Ct < Dt && signbf(Dt) == 0)
				return (signbf(Bt) == signbf(v[i])) + flag13;
		} else {
			if (Ct > Dt && signbf(Ct) == 0)
				return (signbf(At) == signbf(v[i])) + flag13;
		}
	}
	return 0;
}

unsigned int MC33_fail_mem_VN(MC33 *M) {
	M->nV = 0;
	M->memoryfault = 1;
	return 0;
}

void MC33_fail_mem_T(MC33 *M) {
	M->nT = 0;
	M->capt >>= 1;
	M->memoryfault = 1;
	return;
}

/******************************************************************
Assign memory for the vertex r[3], normal n[3]. The return value is the new
vertex label.
*/
unsigned int MC33_spn0(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	for (int i = 0; i != 3; i++)
		p[i] = *(r++);
	// now r points to normal coordinates
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
unsigned int MC33_spnA(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	for (int i = 0; i != 3; i++)
		p[i] = *(r++)*M->D[i] + M->O[i];
	// now r points to normal coordinates
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
unsigned int MC33_spnB(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	for (int i = 0; i != 3; i++)
		p[i] = *(r++)*M->D[i] + M->O[i];
	// now r points to normal coordinates
	r[0] *= M->ca; // normal[0]
	r[1] *= M->cb; // normal[1]
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
#ifndef GRD_ORTHOGONAL
unsigned int MC33_spnC(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	mult_Abf(M->_A,r,r,0);
	for (int i = 0; i != 3; i++)
		p[i] = *(r++) + M->O[i];
	// now r points to normal coordinates
	mult_Abf(M->A_,r,r,1);
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
#endif

/******************************************************************
Auxiliary function that calculates the normal if a vertex
of the cube lies on the isosurface.
*/
unsigned int MC33_surfint(MC33 *M, unsigned int x, unsigned int y, unsigned int z, float *r) {
	r[0] = x; r[1] = y; r[2] = z;
	if (x == 0)
		r[3] = M->F[z][y][0] - M->F[z][y][1];
	else if (x == M->nx)
		r[3] = M->F[z][y][x - 1] - M->F[z][y][x];
	else
		r[3] = 0.5f*(M->F[z][y][x - 1] - M->F[z][y][x + 1]);
	if (y == 0)
		r[4] = M->F[z][0][x] - M->F[z][1][x];
	else if (y == M->ny)
		r[4] = M->F[z][y - 1][x] - M->F[z][y][x];
	else
		r[4] = 0.5f*(M->F[z][y - 1][x] - M->F[z][y + 1][x]);
	if (z == 0)
		r[5] = M->F[0][y][x] - M->F[1][y][x];
	else if (z == M->nz)
		r[5] = M->F[z - 1][y][x] - M->F[z][y][x];
	else
		r[5] = 0.5f*(M->F[z - 1][y][x] - M->F[z + 1][y][x]);
	return M->store(M, r);
}

/******************************************************************
This function find the MC33 case (using the index i, and the face and interior
tests). The correct triangle pattern is obtained from the arrays contained in
the MC33_LookUpTable.h file. The necessary vertices (intersection points) are
also calculated here (using trilinear interpolation).
       _____2_____
     /|          /|
   11 |<-3     10 |
   /____6_____ /  1     z
  |   |       |   |     |
  |   |_____0_|___|     |____y
  7  /        5  /     /
  | 8         | 9     x
  |/____4_____|/

The temporary matrices: M->Lz, M->Dx, M->Dy, M->Ux and M->Uy are filled
and used here.*/
#define FF 0xFFFFFFFF
void MC33_findCase(MC33 *M, unsigned int x, unsigned int y, unsigned int z, unsigned int i, const float *v) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	unsigned int ti[3];//for vertex indices of a triangle
	union { // memory saving
		int f[6];//for the face tests
		float r[6];//for intercept and normal coordinates
	} u;
	const unsigned short int *pcase = MC33_all_tables;
	unsigned int c, m, k, n;
	float t;
	if (i&0x80) {
		c = pcase[i^0xFF];
		m = (c&0x800) == 0;
		n = !m;
	} else {
		c = pcase[i];
		n = (c&0x800) == 0;
		m = !n;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // case 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&MC33_faceTest1(k>>2, v)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (MC33_interiorTest(k, 0, v)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&MC33_faceTest1(k%6, v))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (MC33_interiorTest(k/6, 0, v)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (u.f[4] + u.f[5] < 0? (u.f[0] + u.f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (u.f[4] + u.f[5] < 0? 983: (u.f[0] + u.f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (MC33_interiorTest(k>>1, 0, v)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2:
				if (k == 2? MC33_interiorTest(0, 0, v): MC33_interiorTest(0, 0, v)||MC33_interiorTest(k? 1: 3, 0, v))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (u.f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k == 2? MC33_interiorTest(1, 0, v): MC33_interiorTest(2, 0, v)||MC33_interiorTest(k? 3: 1, 0, v))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2: //12.1
				pcase += (MC33_interiorTest((0xDA010C>>(k<<1))&3, 0, v)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (u.f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (MC33_interiorTest((0xA7B7E5>>(k<<1))&3, 0, v)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(MC33_faceTests(u.f, 165, v))) {
			case 0:
				k = ((u.f[1] < 0)<<1)|(u.f[5] < 0);
				if (u.f[0]*u.f[1] == u.f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = MC33_interiorTest(k, 1, v); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((u.f[0] < 0? u.f[2] > 0: 12 + (u.f[2] < 0)) + (u.f[1] < 0? u.f[3] < 0: 6 + (u.f[3] > 0)));
				if (u.f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*u.f[0] + 4*u.f[1] + 3*u.f[2] + 2*u.f[3] + u.f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*u.f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) { // the vertices r[3] and normals (r + 3)[3] are calculated here
				case 0:
					if (z || x)
						ti[--k] = p[0] = M->Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = M->Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = M->Dx[y][0];
							else if (y? signbf(M->iso - M->F[0][y - 1][0]): 0)
								p[0] = M->Dy[y - 1][0];
							else
								p[0] = MC33_surfint(M,0,y,0,u.r);
						} else if (v[1] == 0) {
							if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = (p[1] != FF? p[1]: MC33_surfint(M,0,y + 1,0,u.r));
						} else {
							t = v[0]/(v[0] - v[1]);
							u.r[0] = u.r[2] = 0;
							u.r[1] = y + t;
							u.r[3] = (v[4] - v[0])*(1 - t) + (v[5] - v[1])*t;
							u.r[4] = v[1] - v[0];
							u.r[5] = (v[3] - v[0])*(1 - t) + (v[2] - v[1])*t;
							p[0] = M->store(M, u.r);
						}
						M->Dy[y][0] = ti[--k] = p[0];
					}
					break;
				case 1:
					if (x)
						ti[--k] = p[1] = M->Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = M->Dy[y][0];
							//else if (z && signbf(v[5]))
							//	p[1] = M->Dx[y + 1][0];
							else if (z && y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][0]): 0)
								p[1] = M->Dy[y + 1][0];
							else if (z? signbf(M->iso - M->F[z - 1][y + 1][0]): 0) {
								ti[--k] = p[1] = M->Lz[y + 1][0]; // value of previous slice
								break;
							} else
								p[1] = MC33_surfint(M, 0, y + 1, z, u.r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = (p[2] != FF? p[2]: MC33_surfint(M, 0, y + 1, z + 1, u.r));
						} else {
							t = v[1]/(v[1] - v[2]);
							u.r[0] = 0; u.r[1] = y + 1;
							u.r[2] = z + t;
							u.r[3] = (v[5] - v[1])*(1 - t) + (v[6] - v[2])*t;
							u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[z][y][0] - M->F[z][y + 2][0])*(1 - t)
										+ (M->F[z + 1][y][0] - M->F[z + 1][y + 2][0])*t):
										(v[1] - v[0])*(1 - t) + (v[2] - v[3])*t);
							u.r[5] = v[2] - v[1];
							p[1] = M->store(M, u.r);
						}
						M->Lz[y + 1][0] = ti[--k] = p[1];
					}
					break;
				case 2:
					if (x)
						ti[--k] = p[2] = M->Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = M->Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = M->Ux[y][0];
							else if (y? signbf(M->iso - M->F[z + 1][y - 1][0]): 0)
								p[2] = M->Uy[y - 1][0];
							else
								p[2] = MC33_surfint(M,0,y,z + 1,u.r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = (p[1] != FF? p[1]: MC33_surfint(M,0,y + 1,z + 1,u.r));
						} else {
							t = v[3]/(v[3] - v[2]);
							u.r[0] = 0; u.r[2] = z + 1;
							u.r[1] = y + t;
							u.r[3] = (v[7] - v[3])*(1 - t) + (v[6] - v[2])*t;
							u.r[4] = v[2] - v[3];
							u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][y][0] - M->F[z + 2][y][0])*(1 - t)
										+ (M->F[z][y + 1][0] - M->F[z + 2][y + 1][0])*t):
										(v[3] - v[0])*(1 - t) + (v[2] - v[1])*t);
							p[2] = M->store(M, u.r);
						}
						M->Uy[y][0] = ti[--k] = p[2];
					}
					break;
				case 3:
					if (y || x)
						ti[--k] = p[3] = M->Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = M->Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = M->Dx[0][0];
							else if (z? signbf(M->iso - M->F[z - 1][0][0]): 0) {
								ti[--k] = p[3] = M->Lz[0][0]; // value of previous slice
								break;
							} else
								p[3] = MC33_surfint(M,0,0,z,u.r);
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else
								p[3] = (p[11] != FF? p[11]: MC33_surfint(M,0,0,z + 1,u.r));
						} else {
							t = v[0]/(v[0] - v[3]);
							u.r[0] = u.r[1] = 0;
							u.r[2] = z + t;
							u.r[3] = (v[4] - v[0])*(1 - t) + (v[7] - v[3])*t;
							u.r[4] = (v[1] - v[0])*(1 - t) + (v[2] - v[3])*t;
							u.r[5] = v[3] - v[0];
							p[3] = M->store(M,u.r);
						}
						M->Lz[0][0] = ti[--k] = p[3];
					}
					break;
				case 4:
					if (z)
						ti[--k] = p[4] = M->Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[0]))
								p[4] = M->Dx[y][x];
							else if (y && signbf(v[7]))
								p[4] = M->Lz[y][x + 1];
							else if (y? signbf(M->iso - M->F[0][y - 1][x + 1]): 0)
								p[4] = M->Dy[y - 1][x + 1];
							else if (y && x + 1 < M->nx? signbf(M->iso - M->F[0][y][x + 2]): 0)
								p[4] = M->Dx[y][x + 1];
							else
								p[4] = MC33_surfint(M,x + 1,y,0,u.r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[4] = p[5];
							else
								p[4] = (p[9] != FF? p[9]: MC33_surfint(M,x + 1,y + 1,0,u.r));
						} else {
							t = v[4]/(v[4] - v[5]);
							u.r[0] = x + 1; u.r[2] = 0;
							u.r[1] = y + t;
							u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[0][y][x] - M->F[0][y][x + 2])*(1 - t)
										+ (M->F[0][y + 1][x] - M->F[0][y + 1][x + 2])*t):
										(v[4] - v[0])*(1 - t) + (v[5] - v[1])*t);
							u.r[4] = v[5] - v[4];
							u.r[5] = (v[7] - v[4])*(1 - t) + (v[6] - v[5])*t;
							p[4] = M->store(M,u.r);
						}
						M->Dy[y][x + 1] = ti[--k] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (z) {
							if (signbf(v[4]))
								p[5] = p[4] = M->Dy[y][x + 1];
							else if (signbf(v[1]))
								p[5] = p[9] = M->Dx[y + 1][x];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z][y + 1][x + 2]): 0)
								p[5] = M->Dx[y + 1][x + 1];
							else if (y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][x + 1]): 0)
								p[5] = M->Dy[y + 1][x + 1];
							else if (signbf(M->iso - M->F[z - 1][y + 1][x + 1])) {
								ti[--k] = p[5] = M->Lz[y + 1][x + 1]; // value of previous slice
								break;
							} else
								p[5] = MC33_surfint(M,x + 1,y + 1,z,u.r);
						} else
							p[5] = MC33_surfint(M,x + 1,y + 1,0,u.r);
					} else if (v[6] == 0)
						p[5] = MC33_surfint(M,x + 1,y + 1,z + 1,u.r);
					else {
						t = v[5]/(v[5] - v[6]);
						u.r[0] = x + 1; u.r[1] = y + 1;
						u.r[2] = z + t;
						u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[z][y + 1][x] - M->F[z][y + 1][x + 2])*(1 - t)
									+ (M->F[z + 1][y + 1][x] - M->F[z + 1][y + 1][x + 2])*t):
									(v[5] - v[1])*(1 - t) + (v[6] - v[2])*t);
						u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[z][y][x + 1] - M->F[z][y + 2][x + 1])*(1 - t)
									+ (M->F[z + 1][y][x + 1] - M->F[z + 1][y + 2][x + 1])*t):
									(v[5] - v[4])*(1 - t) + (v[6] - v[7])*t);
						u.r[5] = v[6] - v[5];
						p[5] = M->store(M,u.r);
					}
					M->Lz[y + 1][x + 1] = ti[--k] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = M->Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = M->Lz[y][x + 1];
							else if (signbf(M->iso - M->F[z + 1][y - 1][x + 1]))
								p[6] = M->Uy[y - 1][x + 1];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z + 1][y][x + 2]): 0)
								p[6] = M->Ux[y][x + 1];
							else
								p[6] = MC33_surfint(M,x + 1,y,z + 1,u.r);
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = MC33_surfint(M,x + 1,0,z + 1,u.r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[6] = p[5];
						else
							p[6] = (p[10] == FF? MC33_surfint(M,x + 1,y + 1,z + 1,u.r): p[10]);
					} else {
						t = v[7]/(v[7] - v[6]);
						u.r[0] = x + 1;
						u.r[1] = y + t; u.r[2] = z + 1;
						u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[z + 1][y][x] - M->F[z + 1][y][x + 2])*(1 - t)
									+ (M->F[z + 1][y + 1][x] - M->F[z + 1][y + 1][x + 2])*t):
									(v[7] - v[3])*(1 - t) + (v[6] - v[2])*t);
						u.r[4] = v[6] - v[7];
						u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][y][x + 1] - M->F[z + 2][y][x + 1])*(1 - t)
										+ (M->F[z][y + 1][x + 1] - M->F[z + 2][y + 1][x + 1])*t):
						(v[7] - v[4])*(1 - t) + (v[6] - v[5])*t);
						p[6] = M->store(M,u.r);
					}
					M->Uy[y][x + 1] = ti[--k] = p[6];
					break;
				case 7:
					if (y)
						ti[--k] = p[7] = M->Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = M->Dx[0][x];
							else if (z && signbf(v[5]))
								p[7] = M->Dy[0][x + 1];
							else if (z && x + 1 < M->nx? signbf(M->iso - M->F[z][0][x + 2]): 0)
								p[7] = M->Dx[0][x + 1];
							else if (z? signbf(M->iso - M->F[z - 1][0][x + 1]): 0) {
								ti[--k] = p[7] = M->Lz[0][x + 1]; // value of previous slice
								break;
							} else
								p[7] = MC33_surfint(M,x + 1,0,z,u.r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else
								p[7] = (p[11] != FF? p[11]: MC33_surfint(M,x + 1,0,z + 1,u.r));
						} else {
							t = v[4]/(v[4] - v[7]);
							u.r[0] = x + 1; u.r[1] = 0;
							u.r[2] = z + t;
							u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[z][0][x] - M->F[z][0][x + 2])*(1 - t)
										+ (M->F[z + 1][0][x] - M->F[z + 1][0][x + 2])*t):
										(v[4] - v[0])*(1 - t) + (v[7] - v[3])*t);
							u.r[4] = (v[5] - v[4])*(1 - t) + (v[6] - v[7])*t;
							u.r[5] = v[7] - v[4];
							p[7] = M->store(M,u.r);
						}
						M->Lz[0][x + 1] = ti[--k] = p[7];
					}
					break;
				case 8:
					if (z || y)
						ti[--k] = p[8] = M->Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = M->Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = M->Dy[0][x];
							else if (x? signbf(M->iso - M->F[0][0][x - 1]): 0)
								p[8] = M->Dx[0][x - 1];
							else
								p[8] = MC33_surfint(M,x,0,0,u.r);
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else
								p[8] = (p[7] != FF? p[7]: MC33_surfint(M,x + 1,0,0,u.r));
						} else {
							t = v[0]/(v[0] - v[4]);
							u.r[1] = u.r[2] = 0;
							u.r[0] = x + t;
							u.r[3] = v[4] - v[0];
							u.r[4] = (v[1] - v[0])*(1 - t) + (v[5] - v[4])*t;
							u.r[5] = (v[3] - v[0])*(1 - t) + (v[7] - v[4])*t;
							p[8] = M->store(M,u.r);
						}
						M->Dx[0][x] = ti[--k] = p[8];
					}
					break;
				case 9:
					if (z)
						ti[--k] = p[9] = M->Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							//else if (p[1] != FF)
							//	p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = M->Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = M->Lz[y + 1][x];
							else if (x? signbf(M->iso - M->F[0][y + 1][x - 1]): 0)
								p[9] = M->Dx[y + 1][x - 1];
							else
								p[9] = MC33_surfint(M,x,y + 1,0,u.r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else
								p[9] = (p[4] != FF? p[4]: MC33_surfint(M,x + 1,y + 1,0,u.r));
						} else {
							t = v[1]/(v[1] - v[5]);
							u.r[1] = y + 1; u.r[2] = 0;
							u.r[0] = x + t;
							u.r[3] = v[5] - v[1];
							u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[0][y][x] - M->F[0][y + 2][x])*(1 - t)
										+ (M->F[0][y][x + 1] - M->F[0][y + 2][x + 1])*t):
										(v[1] - v[0])*(1 - t) + (v[5] - v[4])*t);
							u.r[5] = (v[2] - v[1])*(1 - t) + (v[6] - v[5])*t;
							p[9] = M->store(M,u.r);
						}
						M->Dx[y + 1][x] = ti[--k] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = M->Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = M->Uy[y][x];
							else if (signbf(M->iso - M->F[z + 1][y + 1][x - 1]))
								p[10] = M->Ux[y + 1][x - 1];
							else
								p[10] = MC33_surfint(M,x,y + 1,z + 1,u.r);
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = MC33_surfint(M,0,y + 1,z + 1,u.r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[10] = p[5];
						else
							p[10] = (p[6] != FF? p[6]: MC33_surfint(M,x + 1,y + 1,z + 1,u.r));
					} else {
						t = v[2]/(v[2] - v[6]);
						u.r[0] = x + t;
						u.r[1] = y + 1; u.r[2] = z + 1;
						u.r[3] = v[6] - v[2];
						u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[z + 1][y][x] - M->F[z + 1][y + 2][x])*(1 - t)
									+ (M->F[z + 1][y][x + 1] - M->F[z + 1][y + 2][x + 1])*t):
									(v[2] - v[3])*(1 - t) + (v[6] - v[7])*t);
						u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][y + 1][x] - M->F[z + 2][y + 1][x])*(1 - t)
									+ (M->F[z][y + 1][x + 1] - M->F[z + 2][y + 1][x + 1])*t):
									(v[2] - v[1])*(1 - t) + (v[6] - v[5])*t);
						p[10] = M->store(M,u.r);
					}
					M->Ux[y + 1][x] = ti[--k] = p[10];
					break;
				case 11:
					if (y)
						ti[--k] = p[11] = M->Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = M->Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = M->Uy[0][x];
							else if (x? signbf(M->iso - M->F[z + 1][0][x - 1]): 0)
								p[11] = M->Ux[0][x - 1];
							else
								p[11] = MC33_surfint(M,x,0,z + 1,u.r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else
								p[11] = (p[7] != FF? p[7]: MC33_surfint(M,x + 1,0,z + 1,u.r));
						} else {
							t = v[3]/(v[3] - v[7]);
							u.r[1] = 0; u.r[2] = z + 1;
							u.r[0] = x + t;
							u.r[3] = v[7] - v[3];
							u.r[4] = (v[2] - v[3])*(1 - t) + (v[6] - v[7])*t;
							u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][0][x] - M->F[z + 2][0][x])*(1 - t)
										+ (M->F[z][0][x + 1] - M->F[z + 2][0][x + 1])*t):
										(v[3] - v[0])*(1 - t) + (v[7] - v[4])*t);
							p[11] = M->store(M,u.r);
						}
						M->Ux[0][x] = ti[--k] = p[11];
					}
				break;
				default:
					u.r[0] = x + 0.5f; u.r[1] = y + 0.5f; u.r[2] = z + 0.5f;
					u.r[3] = v[4] + v[5] + v[6] + v[7] - v[0] - v[1] - v[2] - v[3];
					u.r[4] = v[1] + v[2] + v[5] + v[6] - v[0] - v[3] - v[4] - v[7];
					u.r[5] = v[2] + v[3] + v[6] + v[7] - v[0] - v[1] - v[4] - v[5];
					ti[--k] = p[12] = M->store(M,u.r);
				}
			} else
				ti[--k] = p[c];//now ti contains the vertex indices of the triangle
		}
		if (ti[0] != ti[1] && ti[0] != ti[2] && ti[1] != ti[2]) { //to avoid zero area triangles
			if (M->nT == M->capt) {
				unsigned int (*pt)[3] = M->T;
				M->capt <<= 1;
				M->T = (unsigned int (*)[3])realloc(pt,M->capt*3*sizeof(int));
				if (!M->T) {
					M->T = pt;
					MC33_fail_mem_T(M);
				}
			}
			unsigned int *vp = M->T[M->nT++];
#ifndef MC_NORMAL_NEG
			*vp = ti[n]; *(++vp) = ti[m]; *(++vp) = ti[2];
#else
			*vp = ti[m]; *(++vp) = ti[n]; *(++vp) = ti[2];
#endif
		}
	}
}

void MC33_Case_count(MC33 *M, unsigned int x, unsigned int y, unsigned int z, unsigned int i, const float *v) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	union {
		int f[6];
		unsigned int ti[3];
	} u;
	const unsigned short int *pcase = MC33_all_tables;
	unsigned int c, m, k;
	if (i&0x80) {
		c = pcase[i^0xFF];
		m = (c&0x800) == 0;
	} else {
		c = pcase[i];
		m = (c&0x800) != 0;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // cases 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&MC33_faceTest1(k>>2, v)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (MC33_interiorTest(k, 0, v)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&MC33_faceTest1(k%6, v))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (MC33_interiorTest(k/6, 0, v)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (u.f[4] + u.f[5] < 0? (u.f[0] + u.f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (u.f[4] + u.f[5] < 0? 983: (u.f[0] + u.f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (MC33_interiorTest(k>>1, 0, v)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2:
				if (k == 2? MC33_interiorTest(0, 0, v): MC33_interiorTest(0, 0, v)||MC33_interiorTest(k? 1: 3, 0, v))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (u.f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k == 2? MC33_interiorTest(1, 0, v): MC33_interiorTest(2, 0, v)||MC33_interiorTest(k? 3: 1, 0, v))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2: //12.1
				pcase += (MC33_interiorTest((0xDA010C>>(k<<1))&3, 0, v)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (u.f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (MC33_interiorTest((0xA7B7E5>>(k<<1))&3, 0, v)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(MC33_faceTests(u.f, 165, v))) {
			case 0:
				k = ((u.f[1] < 0)<<1)|(u.f[5] < 0);
				if (u.f[0]*u.f[1] == u.f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = MC33_interiorTest(k, 1, v); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 - 40*c + (c? 10: 6)*k;
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((u.f[0] < 0? u.f[2] > 0: 12 + (u.f[2] < 0)) + (u.f[1] < 0? u.f[3] < 0: 6 + (u.f[3] > 0)));
				if (u.f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*u.f[0] + 4*u.f[1] + 3*u.f[2] + 2*u.f[3] + u.f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*u.f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) {
				case 0:
					if (z || x)
						p[0] = M->Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = M->Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = M->Dx[y][0];
							else if (y? signbf(M->iso - M->F[0][y - 1][0]): 0)
								p[0] = M->Dy[y - 1][0];
							else
								p[0] = M->nV++;
						} else if (v[1] == 0) {
							if (p[1] != FF)
								p[0] = p[1];
							else if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = M->nV++;
						} else
							p[0] = M->nV++;

						M->Dy[y][0] = p[0];
					}
					break;
				case 1:
					if (x)
						p[1] = M->Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = M->Dy[y][0];
							else if (z && signbf(v[5]))
								p[1] = M->Dx[y + 1][0];
							else if (z && y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][0]): 0)
								p[1] = M->Dy[y + 1][0];
							else if (z? signbf(M->iso - M->F[z - 1][y + 1][0]): 0) {
								p[1] = M->Lz[y + 1][0];
								break;
							} else
								p[1] = M->nV++;
						} else if (v[2] == 0) {
							if (p[2] != FF)
								p[1] = p[2];
							else if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = M->nV++;
						} else
							p[1] = M->nV++;
						M->Lz[y + 1][0] = p[1];
					}
					break;
				case 2:
					if (x)
						p[2] = M->Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = M->Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = M->Ux[y][0];
							else if (y? signbf(M->iso - M->F[z + 1][y - 1][0]): 0)
								p[2] = M->Uy[y - 1][0];
							else
								p[2] = M->nV++;
						} else if (v[2] == 0) {
							if (p[1] != FF)
								p[2] = p[1];
							else if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = M->nV++;
						} else
							p[2] = M->nV++;
						M->Uy[y][0] = p[2];
					}
					break;
				case 3:
					if (y || x)
						p[3] = M->Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = M->Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = M->Dx[0][0];
							else if (z? signbf(M->iso - M->F[z - 1][0][0]): 0) {
								p[3] = M->Lz[0][0];
								break;
							} else
								p[3] = M->nV++;
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else if (p[11] != FF)
								p[3] = p[11];
							else
								p[3] = M->nV++;
						} else
							p[3] = M->nV++;
						M->Lz[0][0] = p[3];
					}
					break;
				case 4:
					if (z)
						p[4] = M->Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[0]))
								p[4] = M->Dx[y][x];
							else if (y && signbf(v[7]))
								p[4] = M->Lz[y][x + 1];
							else if (y? signbf(M->iso - M->F[0][y - 1][x + 1]): 0)
								p[4] = M->Dy[y - 1][x + 1];
							else if (y && x + 1 < M->nx? signbf(M->iso - M->F[0][y][x + 2]): 0)
								p[4] = M->Dx[y][x + 1];
							else
								p[4] = M->nV++;
						} else if (v[5] == 0) {
							if (p[9] != FF)
								p[4] = p[9];
							else if (p[5] != FF)
								p[4] = p[5];
							else
								p[4] = M->nV++;
						} else
							p[4] = M->nV++;
						M->Dy[y][x + 1] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (z) {
							if (signbf(v[4]))
								p[5] = p[4] = M->Dy[y][x + 1];
							else if (signbf(v[1]))
								p[5] = p[9] = M->Dx[y + 1][x];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z][y + 1][x + 2]): 0)
								p[5] = M->Dx[y + 1][x + 1];
							else if (y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][x + 1]): 0)
								p[5] = M->Dy[y + 1][x + 1];
							else if (signbf(M->iso - M->F[z - 1][y + 1][x + 1])) {
								p[5] = M->Lz[y + 1][x + 1]; // value of previous slice
								break;
							} else
								p[5] = M->nV++;
						} else
							p[5] = M->nV++;
					} else
						p[5] = M->nV++;
					M->Lz[y + 1][x + 1] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = M->Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = M->Lz[y][x + 1];
							else if (signbf(M->iso - M->F[z + 1][y - 1][x + 1]))
								p[6] = M->Uy[y - 1][x + 1];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z + 1][y][x + 2]): 0)
								p[6] = M->Ux[y][x + 1];
							else
								p[6] = M->nV++;
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = M->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[6] = (p[10] == FF? M->nV++: p[10]);
						else
							p[6] = p[5];
					} else
						p[6] = M->nV++;
					M->Uy[y][x + 1] = p[6];
					break;
				case 7:
					if (y)
						p[7] = M->Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = M->Dx[0][x];
							else if (z && signbf(v[5]))
								p[7] = M->Dy[0][x + 1];
							else if (z && x + 1 < M->nx? signbf(M->iso - M->F[z][0][x + 2]): 0)
								p[7] = M->Dx[0][x + 1];
							else if (z? signbf(M->iso - M->F[z - 1][0][x + 1]): 0) {
								p[7] = M->Lz[0][x + 1];
								break;
							} else
								p[7] = M->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else if (p[11] != FF)
								p[7] = p[11];
							else
								p[7] = M->nV++;
						} else
							p[7] = M->nV++;
						M->Lz[0][x + 1] = p[7];
					}
					break;
				case 8:
					if (z || y)
						p[8] = M->Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = M->Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = M->Dy[0][x];
							else if (x? signbf(M->iso - M->F[0][0][x - 1]): 0)
								p[8] = M->Dx[0][x - 1];
							else
								p[8] = M->nV++;
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else if (p[7] != FF)
								p[8] = p[7];
							else
								p[8] = M->nV++;
						} else
							p[8] = M->nV++;
						M->Dx[0][x] = p[8];
					}
					break;
				case 9:
					if (z)
						p[9] = M->Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							else if (p[1] != FF)
								p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = M->Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = M->Lz[y + 1][x];
							else if (x? signbf(M->iso - M->F[0][y + 1][x - 1]): 0)
								p[9] = M->Dx[y + 1][x - 1];
							else
								p[9] = M->nV++;
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else if (p[4] != FF)
								p[9] = p[4];
							else
								p[9] = M->nV++;
						} else
							p[9] = M->nV++;
						M->Dx[y + 1][x] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = M->Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = M->Uy[y][x];
							else if (signbf(M->iso - M->F[z + 1][y + 1][x - 1]))
								p[10] = M->Ux[y + 1][x - 1];
							else
								p[10] = M->nV++;
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = M->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[10] = (p[6] == FF? M->nV++: p[6]);
						else
							p[10] = p[5];
					} else
						p[10] = M->nV++;
					M->Ux[y + 1][x] = p[10];
					break;
				case 11:
					if (y)
						p[11] = M->Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = M->Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = M->Uy[0][x];
							else if (x? signbf(M->iso - M->F[z + 1][0][x - 1]): 0)
								p[11] = M->Ux[0][x - 1];
							else
								p[11] = M->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else if (p[7] != FF)
								p[11] = p[7];
							else
								p[11] = M->nV++;
						} else
							p[11] = M->nV++;
						M->Ux[0][x] = p[11];
					}
				break;
				default:
					p[12] = M->nV++;
				}
			}
			u.ti[--k] = p[c];
		}
		if (u.ti[0] != u.ti[1] && u.ti[0] != u.ti[2] && u.ti[1] != u.ti[2])
			M->nT++;
	}
}
#undef FF

void MC33_freeTemp_O_N(MC33 *M) {
	free(M->Dx); free(M->Ux); free(M->Dy); free(M->Uy);
	free(M->Lz);
	free(M);
}

void free_MC33(MC33 *M) {
	unsigned int y;
	if (M) {
		for (y = 0; y != M->ny; y++) {
			free(M->Dx[y]); free(M->Ux[y]); free(M->Dy[y]); free(M->Uy[y]);
			free(M->Lz[y]);
		}
		free(M->Dx[M->ny]); free(M->Ux[M->ny]);
		free(M->Lz[M->ny]);
		MC33_freeTemp_O_N(M);
	}
}

MC33 *create_MC33(_GRD* G) {
	unsigned int x, y;
	MC33 *M;
	if (!G)
		return 0;
	M = (MC33*)malloc(sizeof(MC33));
	if (!M)
		return 0;
	M->nx = G->N[0];
	M->ny = G->N[1];
	M->nz = G->N[2];

#ifndef GRD_ORTHOGONAL
	if (G->nonortho) {
		M->store = MC33_spnC;
		for (int j = 0; j != 3; j++)
			for (int i = 0; i != 3; i++) {
				M->_A[j][i] = G->_A[j][i]*G->d[i]; // true transformation matrices
				M->A_[j][i] = G->A_[j][i]/G->d[j];
			}
	} else
#endif
	if (G->d[0] != G->d[1] || G->d[1] != G->d[2]) {
		M->ca = G->d[2]/G->d[0];
		M->cb = G->d[2]/G->d[1];
		M->store = MC33_spnB;
	} else
		M->store = (G->d[0] == 1 && G->r0[0] == 0 && G->r0[1] == 0 && G->r0[2] == 0? MC33_spn0: MC33_spnA);

	for (int j = 0; j != 3; j++) {
		M->O[j] = G->r0[j];
		M->D[j] = G->d[j];
	}
	M->Lz = (unsigned int**)malloc((M->ny + 1)*sizeof(int*));//edges 1, 3, 5 (only write) and 7
	M->Dy = (unsigned int**)malloc(M->ny*sizeof(int*));//edges 0 (only read) and 4
	M->Uy = (unsigned int**)malloc(M->ny*sizeof(int*));//edges 2 and 6 (only write)
	M->Dx = (unsigned int**)malloc((M->ny + 1)*sizeof(int*));//edges 8 and 9
	M->Ux = (unsigned int**)malloc((M->ny + 1)*sizeof(int*));//edges 10 (only write) and 11
	if (!M->Ux) {
		MC33_freeTemp_O_N(M);
		return 0;
	}
	M->F = (const GRD_data_type***)G->F;
	x = M->nx*sizeof(int);
	for (y = 0; y != M->ny; y++) {
		M->Dx[y] = (unsigned int*)malloc(x);
		M->Ux[y] = (unsigned int*)malloc(x);
		M->Lz[y] = (unsigned int*)malloc(x + sizeof(int));
		M->Dy[y] = (unsigned int*)malloc(x + sizeof(int));
		M->Uy[y] = (unsigned int*)malloc(x + sizeof(int));
	}
	if (M->Uy[y - 1]) {
		M->Dx[y] = (unsigned int*)malloc(x);
		M->Ux[y] = (unsigned int*)malloc(x);
		M->Lz[y] = (unsigned int*)malloc(x + sizeof(int));
		if (M->Lz[y])
			return M;
	} else
		M->Dx[y] = M->Ux[y] = M->Lz[y] = 0;
	free_MC33(M);
	return 0;
}

surface* calculate_isosurface(MC33 *M, float iso) {
	unsigned int x, y, z, Nx = M->nx;
	float Vt[12];
	float *v1 = Vt, *v2 = Vt + 4;
	const GRD_data_type ***F = M->F, **F0, **F1, *V00, *V01, *V11, *V10;
	surface *S = (surface*)malloc(sizeof(surface));
	if (!S)
		return 0;
	M->nT = M->nV = 0;
	M->memoryfault = 0;
	M->capt = M->capv = 4096;
	M->T = (unsigned int(*)[3])malloc(3*4096*sizeof(int));
	M->N = (float(*)[3])malloc(3*4096*sizeof(float));
	M->V = (float(*)[3])malloc(3*4096*sizeof(float));
	M->iso = iso;
	if (M->V)
		for (z = 0; z != M->nz; z++) {
			F0 = *F;
			F1 = *(++F);
			for (y = 0; y != M->ny; y++) {
				V00 = *F0;
				V01 = *(++F0);
				V10 = *F1;
				V11 = *(++F1);
				v2[0] = iso - *V00;//the difference was inverted to use signbit function
				v2[1] = iso - *V01;
				v2[2] = iso - *V11;
				v2[3] = iso - *V10;
				//the eight least significant bits of i correspond to vertex indices. (x...x01234567)
				//If the bit is 1 then the vertex value is greater than zero.
				unsigned int i = signbf(v2[3]) != 0;
				if (signbf(v2[2])) i |= 2;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[0])) i |= 8;
				for (x = 0; x != Nx; x++) {
					{float *P = v1; v1 = v2; v2 = P;}//v1 and v2 are exchanged
					v2[0] = iso - *(++V00);
					v2[1] = iso - *(++V01);
					v2[2] = iso - *(++V11);
					v2[3] = iso - *(++V10);
					i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
					if (signbf(v2[2])) i |= 2;
					if (signbf(v2[1])) i |= 4;
					if (signbf(v2[0])) i |= 8;
					if (i && i^0xFF) {
						if (v1 > v2) {float *t = v2; float *s = t + 8; *s = *t; *(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t);}
						MC33_findCase(M,x,y,z,i,v1);
					}
				}
			}
			{unsigned int** P = M->Dx; M->Dx = M->Ux; M->Ux = P;}//M->Dx and M->Ux are exchanged
			{unsigned int** P = M->Dy; M->Dy = M->Uy; M->Uy = P;}//M->Dy and M->Uy are exchanged
		}
	else
		M->memoryfault = 1;
	if (M->nV) {
		M->color = (int*)malloc(M->nV*sizeof(int));
		memcpy(S, M, offsetof(MC33, iso));
		if (M->color) {
			*M->color = DefaultColorMC;
			while (--M->nV)
				*(++M->color) = DefaultColorMC;
		} else
			M->memoryfault = 1;
	} else {
		free(M->V); free(M->N); free(M->T);
		memset(S, 0, sizeof(surface));
	}
	if (M->memoryfault) {
		free_surface_memory(S);
		return 0;
	}
	S->iso = iso;
	return S;
}

// modified from calculate_isosurface function
unsigned long long size_of_isosurface(MC33 *M, float iso, unsigned int *nV, unsigned int *nT) {
	unsigned int x, y, z, Nx = M->nx;
	float Vt[12];
	float *v1 = Vt, *v2 = Vt + 4;
	const GRD_data_type ***F = M->F, **F0, **F1, *V00, *V01, *V11, *V10;
	M->nT = M->nV = 0;
	M->iso = iso;
	for (z = 0; z != M->nz; z++) {
		F0 = *F;
		F1 = *(++F);
		for (y = 0; y != M->ny; y++) {
			V00 = *F0;
			V01 = *(++F0);
			V10 = *F1;
			V11 = *(++F1);
			v2[0] = iso - *V00;
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			unsigned int i = signbf(v2[3]) != 0;
			if (signbf(v2[2])) i |= 2;
			if (signbf(v2[1])) i |= 4;
			if (signbf(v2[0])) i |= 8;
			for (x = 0; x != Nx; x++) {
				{float *P = v1; v1 = v2; v2 = P;}
				v2[0] = iso - *(++V00);
				v2[1] = iso - *(++V01);
				v2[2] = iso - *(++V11);
				v2[3] = iso - *(++V10);
				i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
				if (signbf(v2[2])) i |= 2;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[0])) i |= 8;
				if (i && i^0xFF) {
					if (v1 > v2) {float *t = v2; float *s = t + 8; *s = *t; *(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t);}
					MC33_Case_count(M,x,y,z,i,v1);
				}
			}
		}
		{unsigned int** P = M->Dx; M->Dx = M->Ux; M->Ux = P;}
		{unsigned int** P = M->Dy; M->Dy = M->Uy; M->Uy = P;}
	}
	if (nV)
		*nV = M->nV;
	if (nT)
		*nT = M->nT;
	// number of vertices * (size of vertex and normal + size of color ) + number of triangle * size of triangle + size of struct surface
	return M->nV * (6 * sizeof(float) + sizeof(int)) + M->nT * (3 * sizeof(int)) + sizeof(surface);
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif /*marching_cubes_33_c_implementation*/
#endif /*marching_cubes_33_c_h*/

