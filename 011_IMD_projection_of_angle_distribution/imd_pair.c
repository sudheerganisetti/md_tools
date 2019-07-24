
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_pair -- calculate pair distibution functions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef PAIR
#define PAIR
#endif

#define MAIN

#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-a<nnn>] [-e<nnn>] [-v] [-p paramter-file]\n",progname); 
  exit(1); 
}


/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  int tablesize;
  int i,j,k;
  printf ("BECAREFUL: I made some changes in the original code. I beleive the pair \n \
distribution function implemented in IMD is not a correct one so I followed the \n \
approach from the following sources\n \
1) http://homepage.univie.ac.at/franz.vesely/simsp/dx/node22.html\n \
2) Jochen Zausch@2009, Dissertation:Dynamic, Rheology and Critical Properties of \
Colloidal Fluid Mixtures \n");

  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();

  tablesize = slots*ntypes*ntypes*sizeof(real);
  histogram = (real *) malloc(tablesize);
  histogram2 = (real *) malloc(tablesize);
  if (NULL==histogram) error("Cannot allocate memory for histograms.");
  if (NULL==histogram2) error("Cannot allocate memory for histograms2.");
  
  hist_dim.x = slots;
  hist_dim.y = ntypes;
  hist_dim.z = ntypes;

  for (i=0; i<slots; ++i)
    for (j=0; j<ntypes; ++j)
      for (k=0; k<ntypes; ++k)
      {
	*PTR_3D_V(histogram,i,j,k,hist_dim) = 0.0;
        *PTR_3D_V(histogram2,i,j,k,hist_dim) = 0.0;
      }

  r2_cut = SQR(r_max);

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Calculate the distances */
  do_work(do_cell_pair);

  /* Output results */
  write_data();

  return 0;

}


/******************************************************************************
*
*  write_data writes histogram to *.pair file
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int i,j,k;
  real r;
  real f;

  if (-1==restart)
    sprintf(fname,"%s.pair",infilename);
  else
    sprintf(fname,"%s.%05d.pair",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open histograms file.");

  for (i=1; i<slots; ++i) {
    r = ((float) i / slots * 180);
    fprintf(out,"%f ", r);
    f = 2;
    for (j=0; j<ntypes; ++j)
      for (k=j; k<ntypes; ++k)
      {
	fprintf(out,"%f ",*PTR_3D_V(histogram,i,j,k,hist_dim)/f);
        //fprintf(out,"%f ",*PTR_3D_V(histogram2,i,j,k,hist_dim)/f);
      }
//#ifdef TWOD
//    f = natoms * 2 * 3.14159265 * r;
//#else
//    f = 4* 3.14159265 * SQR(r);
//#endif
//f=0.0;    
//f = box_x.x*box_y.y*box_z.z/(4*4*3.14159265 * SQR(r)*natoms_ptype*natoms_qtype);
//    for (j=0; j<ntypes; ++j)
//      for (k=j; k<ntypes; ++k)
//	{
//	fprintf(out,"%f ",*PTR_3D_V(histogram,i,j,k,hist_dim)*f);
//	}
    fprintf(out,"\n");
  }
  fclose(out);
}


/******************************************************************************
*
*  do_cell_pair calulates the distances for atoms in two cells
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i,j,k;
  int temp;
  vektor d;
  real radius;
  int p_typ,q_typ;
  real projection;
  real phi;
  int ang;
  /* For each atom in first cell */
  for (i = 0;i < p->n; ++i) 
    /* For each atom in neighbouring cell */
    /* Nasty little trick: If p==q, use only rest of atoms */
    for (j = ((p==q) ? i+1 : 0); j < q->n; ++j) {
      
      /* Calculate distance */
      d.x =q->ort[j].x - p->ort[i].x + pbc.x;
      d.y =q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
      d.z =q->ort[j].z - p->ort[i].z + pbc.z;
#endif
      radius = sqrt( (double)(SPROD(d,d)) );
      
      projection = d.z;
      
      if(radius < r_max)
      {
        if(projection<0) {projection=-1*projection;}
        phi = acos( projection / radius )*180/3.141592654;  // result is in radians
        ang = (int) ( slots * phi / 180 );           // results is in integers slots
//      projection = d.z;
//      if (projection<0) {projection =-1*projection;}
//      k     = (int) ( slots * (projection - r_min) / (r_max - r_min));
        p_typ = p->sorte[i];
        q_typ = q->sorte[j];

        if (q_typ > p_typ) {temp = p_typ; p_typ = q_typ; q_typ = temp;}
      
        if ((ang>0) && (ang<slots))
        {
          *PTR_3D_V(histogram, ang , q_typ, p_typ, hist_dim) += 1;
        }
        
      }
}
}
