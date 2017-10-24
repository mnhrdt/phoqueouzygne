/************************************************************************
 * utilities to read and write GRD files headers and data               *
 * files can be read and written on a row-by-row basis                  *
 * modified from readgrd                        		                *
 ************************************************************************/
/************************************************************************
 * Creator: Rob J. Mellors, San Deigo State University                  *
 * Date   : December 18, 2007                                           *
 *                                                                      *
 * Modification history:                                                *
 *   some slight format changes - Dec 18, 2007 - RJM 		            *
 ************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"gmtsar.h"
# include"lib_functions.h"

/*----------------------------------------------------------------------*/
int create_full_GMT_binary_hdr(int nx, int ny, double x_inc, double y_inc, 
				double xmin, double ymin, double zmin, double zmax, 
				char *command, char *remark, FILE *outfile)
{
struct GMT_binary hdr;

	hdr.nx = nx;
	hdr.ny = ny;
	hdr.x_inc = x_inc;
	hdr.y_inc = y_inc;
	strcpy(hdr.command, command);
	strcpy(hdr.remark, remark);

	hdr.node_offset = 1;
	hdr.x_min = xmin;
	hdr.x_max = hdr.nx * hdr.x_inc + hdr.x_min;
	hdr.y_min = ymin;
	hdr.y_max = hdr.ny * hdr.y_inc + hdr.y_min;
	hdr.z_min = zmin;
	hdr.z_max = zmax;
	hdr.z_scale_factor = 1.0;
	hdr.z_add_offset = 0.0;

	strcpy(hdr.x_units, "");
	strcpy(hdr.y_units, "");
	strcpy(hdr.z_units, "");

	strcpy(hdr.title, "");

	if (verbose) fprintf(stderr," creating GMT binary header %d %d %lf %lf %lf %lf %lf %lf %s %s \n", 
		hdr.nx, hdr.ny, hdr.x_inc, hdr.y_inc, hdr.x_min, hdr.y_min, hdr.z_min, hdr.z_max, hdr.command, hdr.remark);

	write_GMT_binary_hdr(hdr, outfile);

	return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------*/
int create_GMT_binary_hdr(int nx, int ny, double x_inc, double y_inc, char *command, char *remark, FILE *outfile)
{
struct GMT_binary hdr;

	hdr.nx = nx;
	hdr.ny = ny;
	hdr.x_inc = x_inc;
	hdr.y_inc = y_inc;
	strcpy(hdr.command, command);
	strcpy(hdr.remark, remark);

	hdr.node_offset = 1;
	hdr.x_min = 0.0;
	hdr.x_max = hdr.nx * hdr.x_inc + hdr.x_min;
	hdr.y_min = 0.0;
	hdr.y_max = hdr.ny * hdr.y_inc + hdr.y_min;
	hdr.z_min = 0.0;
	hdr.z_max = 100.0;	/* default */
	hdr.z_scale_factor = 1.0;
	hdr.z_add_offset = 0.0;

	strcpy(hdr.x_units, "");
	strcpy(hdr.y_units, "");
	strcpy(hdr.z_units, "");

	strcpy(hdr.title, "");

	if (verbose) fprintf(stderr," creating GMT binary header %d %d %lf %lf %s %s \n", hdr.nx, hdr.ny, hdr.x_inc, hdr.y_inc, hdr.command, hdr.remark);

	write_GMT_binary_hdr(hdr, outfile);

	return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------*/
/* write out each variable separately so no problems with 64 bit 	*/
/* (size of struct may vary on 32 versus 64 due to alignment)		*/

int write_GMT_binary_hdr(struct GMT_binary hdr, FILE *outfile)
{
	fwrite(&hdr.nx, sizeof(int), 1, outfile);
	fwrite(&hdr.ny, sizeof(int), 1, outfile);
	fwrite(&hdr.node_offset, sizeof(int), 1, outfile);

	fwrite(&hdr.x_min, sizeof(double), 1, outfile);
	fwrite(&hdr.x_max, sizeof(double), 1, outfile);
	fwrite(&hdr.y_min, sizeof(double), 1, outfile);
	fwrite(&hdr.y_max, sizeof(double), 1, outfile);
	fwrite(&hdr.z_min, sizeof(double), 1, outfile);
	fwrite(&hdr.z_max, sizeof(double), 1, outfile);
	fwrite(&hdr.x_inc, sizeof(double), 1, outfile);
	fwrite(&hdr.y_inc, sizeof(double), 1, outfile);
	fwrite(&hdr.z_scale_factor, sizeof(double), 1, outfile);
	fwrite(&hdr.z_add_offset, sizeof(double), 1, outfile);

	fwrite(&hdr.x_units, sizeof(char), 80, outfile);
	fwrite(&hdr.y_units, sizeof(char), 80, outfile);
	fwrite(&hdr.z_units, sizeof(char), 80, outfile);
	fwrite(&hdr.title, sizeof(char), 80, outfile);
	fwrite(&hdr.command, sizeof(char), 320, outfile);
	fwrite(&hdr.remark, sizeof(char), 160, outfile);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------*/
int read_GMT_binary_header(FILE *in, struct GMT_binary *hdr)
{
	rewind(in);

	fread(&(hdr->nx), sizeof(int), 1, in);
	fread(&(hdr->ny), sizeof(int), 1, in);
	fread(&(hdr->node_offset), sizeof(int), 1, in);

	fread(&(hdr->x_min), sizeof(double), 1, in);
	fread(&(hdr->x_max), sizeof(double), 1, in);
	fread(&(hdr->y_min), sizeof(double), 1, in);
	fread(&(hdr->y_max), sizeof(double), 1, in);
	fread(&(hdr->z_min), sizeof(double), 1, in);
	fread(&(hdr->z_max), sizeof(double), 1, in);
	fread(&(hdr->x_inc), sizeof(double), 1, in);
	fread(&(hdr->y_inc), sizeof(double), 1, in);
	fread(&(hdr->z_scale_factor), sizeof(double), 1, in);
	fread(&(hdr->z_add_offset), sizeof(double), 1, in);

	fread(&(hdr->x_units), sizeof(char), 80, in);
	fread(&(hdr->y_units), sizeof(char), 80, in);
	fread(&(hdr->z_units), sizeof(char), 80, in);
	fread(&(hdr->title), sizeof(char), 80, in);
	fread(&(hdr->command), sizeof(char), 320, in);
	fread(&(hdr->remark), sizeof(char), 160, in);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------*/
int write_GMT_binary_float(struct GMT_binary hdr, float *data, FILE *outfile)
{
	fwrite(data, sizeof(float), hdr.nx * hdr.ny, outfile);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------*/
int read_GMT_binary_dimensions(FILE *in, int *xdim, int *ydim)
{
	rewind(in);

	fread(xdim, sizeof(int), 1, in);
	fread(ydim, sizeof(int), 1, in);

	if (verbose) fprintf(stderr, "read_GMT_binary: xdim %d ydim %d\n", *xdim, *ydim);
	rewind(in);

	/* skip to end of header (use 892 since sizeof (struct ..) would give 896 under 64-bit */
	fseek(in, 892U, SEEK_SET);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------*/
