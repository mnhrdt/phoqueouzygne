/************************************************************************
 * utilities to read and write GRD files headers and data               *
 * modified from readgrd                        		                *
 ************************************************************************/
/************************************************************************
 * Creator: Rob J. Mellors, San Deigo State University                  *
 * Date   : December 18, 2007                                           *
 *                                                                      *
 * Modification history:                                                *
 *   some slight format changes - Dec 18, 2007 - RJM 		            *
 ************************************************************************/

# include "gmt.h"
# include "gmtsar.h"
# include"lib_functions.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

/*----------------------------------------------------------------------------*/
int readgrdsize_(int *nx, int *ny, char *file)
/*
 nx;                 number of x points 	
 ny;                 number of y points 	
 filein;           filename of input file 	
 */ 

{
  int argc2 = 1;
  char *argv2[2] = {"readgrdsize",0};
  struct GRD_HEADER grd;

  if (verbose) fprintf(stderr," readgrdsize: reading GRD file %s\n", file);

/* Initialize with default values */
	
  GMT_begin (argc2,argv2);
  GMT_grd_init(&grd, argc2, argv2, FALSE);
  if (GMT_read_grd_info (file, &grd)) die("Error opening file ", file); 
  *nx = grd.nx;
  *ny = grd.ny;

  return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------------*/
int readgrdinfo_(int *nx, int *ny, double *x0, double *y0, double *x_inc, double *y_inc, char *file)
/*
 nx;                 number of x points 	
 ny;                 number of y points 	
 filein;           filename of input file 	
 */ 

{
  int argc2 = 1;
  char *argv2[2] = {"readgrdinfo",0};
  struct GRD_HEADER grd;

  if (verbose) fprintf(stderr," readgrdsize: reading GRD file %s\n", file);

/* Initialize with default values */
	  
  GMT_begin (argc2,argv2);
  GMT_grd_init(&grd, argc2, argv2, FALSE);
  if (GMT_read_grd_info (file, &grd)) die("Error opening file ", file); 
  *nx = grd.nx;
  *ny = grd.ny;
  *x0 = grd.x_min;
  *y0 = grd.y_min;
  *x_inc = grd.x_inc;
  *y_inc = grd.y_inc;

  return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------------*/
int initialize_grdfile_(int nx, int ny, double rlt0, double rln0, double dlt, double dln, char *title, char *file)
{
  int argc2 = 1;
  char *argv2[2] = {"initialize_grdfile_",0};
  struct GRD_HEADER grd;

  if (verbose) fprintf(stderr," initialize_grdfile: initializing GRD file %s\n",file);
	fprintf(stderr," %d %d %lf %lf %lf %lf %s %s\n",nx, ny, rlt0, rln0, dlt, dln, title, file);

   	GMT_begin (argc2,argv2);
   	GMT_grd_init(&grd, argc2, argv2, FALSE);

	if (verbose) fprintf(stderr," initialize_grdfile: initializing GRD file %s\n", file);
	make_grd_header(&grd, nx, ny, rlt0, rln0, dlt, dln, title);

   	if (verbose) {
		fprintf(stderr,"initialize_grdfile\n");
		fprintf(stderr," file %s title %s offset %d\n", file, grd.title, grd.node_offset);
		fprintf(stderr," nx %d ny %d \n",grd.nx, grd.ny);
		fprintf(stderr," xmin %6.2lf xmax %6.2lf xinc %6.2lf \n",grd.x_min, grd.x_max,grd.x_inc);
		fprintf(stderr," ymin %6.2lf ymax %6.2lf yinc %6.2lf \n",grd.y_min, grd.y_max,grd.y_inc);
		}

/*  write the file header but no data */
	  
	GMT_write_grd_info(file, &grd);

	if (verbose) fprintf(stderr,"done with initialize_grdfile\n");
 	return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------------*/
int make_grd_header(struct GRD_HEADER *grd, int nx, int ny, double rlt0, double rln0, double dlt, double dln, char *title)
/*   
     nx              number of x points 	
     ny              number of y points 	
     rlt0            starting latitude 		
     rln0            starting longitude		
     dlt             latitude spacing 	
     dln             longitude spacing 		
     title           title 			
*/

{
  double rdum;

	rdum = -999999;

/* Initialize with default values */
	  
   	GMT_grd_init(grd, 0, NULL, FALSE);

/* Calculate header parameters */
	  
	grd->nx = nx;
	grd->ny = ny;

   	grd->x_max = rln0 + nx * dln;
   	grd->x_min = rln0;

   	if (grd->x_max < grd->x_min) {
     		grd->x_min = grd->x_max;
     		grd->x_max = rln0;
     		}

   	grd->x_inc = fabs((double) dln);

   	grd->y_max = rlt0 + ny * dlt;
   	grd->y_min = rlt0;

   	if (grd->y_max < grd->y_min) {
     		grd->y_min = grd->y_max;
     		grd->y_max = rlt0;
     		}

   	grd->y_inc = fabs((double) dlt);

/* update the header using values passed */
	  
   	grd->node_offset = TRUE;
   	strncpy(grd->title, title, GRD_TITLE_LEN); 

   	if (verbose) {
		fprintf(stderr,"make_grd_header\n");
		fprintf(stderr," title %s offset %d\n",grd->title, grd->node_offset);
		fprintf(stderr," nx %d ny %d \n",grd->nx, grd->ny);
		fprintf(stderr," xmin %lf xmax %lf xinc %lf \n",grd->x_min, grd->x_max,grd->x_inc);
		fprintf(stderr," ymin %lf ymax %lf yinc %lf \n",grd->y_min, grd->y_max,grd->y_inc);
		}

 	return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------------*/
int read_grd_row(char *file, int nx, int ny, int row_no, int n_rows, float *data)
/* netcdf (gmt) file is assumed to exist with a header already                */
/* read n_rows of data to file beginning at row row_no                        */
/*										                                      */

  {
  int ncid;
  int status;
  int z_id;
  size_t start[2];
  size_t count[2];

	if (verbose) fprintf(stderr,"read_grd_row %d nx %d ny %d\n",row_no, nx, ny);

/* netcdf dimensions */
	  
	start[0] = (size_t) row_no;		/* start y */
	start[1] = (size_t) 0;			/* start x */
	count[0] = (size_t) n_rows;		/* end y */
	count[1] = (size_t) nx;			/* end x */

	status = nc_open(file, NC_NOWRITE, &ncid);
	if (status !=0 ) die("error opening ",file);

	status = nc_get_vara_float(ncid, z_id, start, count, data);
	if (status !=0 ) die("error reading ",file);

	status = nc_close(ncid);
	if (status !=0 ) die("error closing ",file);

 	return(EXIT_SUCCESS);
}
/*----------------------------------------------------------------------------*/
int write_grd_row(char *file, int nx, int ny, int row_no, int n_rows, float *data)
/* netcdf (gmt) file is assumed to exist with a header already                */
/* write n_rows of data to file beginning at row row_no                       */
/* data is in *data                                                           */
/* no checking is done to verify that *data is large enough                   */
/*                                                                            */

  {
  int ncid;
  int status;
  int z_id;
  size_t start[2];
  size_t count[2];

	if (debug) fprintf(stderr,"add_grd_row %d nx %d ny %d\n",row_no, nx, ny);

/* netcdf dimensions */
	  
	start[0] = (size_t) row_no;		/* start y */
	start[1] = (size_t) 0;			/* start x */

	count[0] = (size_t) n_rows;		/* # y to change */
	count[1] = (size_t) nx;			/* end x */

	status = nc_open(file, NC_WRITE, &ncid);
	if (status !=0 ) die("error opening for write  ",file);

	status = nc_inq_varid (ncid, "z", &z_id);
	if (status !=0 ) die("error reading ",file);

	status = nc_put_vara_float(ncid, z_id, start, count, data);
	if (status !=0 ) die("error writing ",file);
		
	status = nc_close(ncid);
	if (status !=0 ) die("error closing ",file);

 	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------*/
