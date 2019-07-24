/* *******************************************************************************************************************
 * Sudheer Ganisetti, Sa 26. Mär 16:14:30 CET 2016
 * This is my own code written to compute/see the orientation of bonds in the space 
 * This in not generalised for all conditions, no guarentee for other cases, use it on your own risk
 * Assumptions: 
 * 1) SiO2 system of 72000 atoms; Si_type=0, O_type=1
 * 2) box will be devided into cells only based on box_xx, box_yy, box_zz (so carefull if it is an orthogonal)
 * 3) Periodic Boundary Conditions in all directions
 * 4) 
 * Theory:
 * suppose "r = ix+jy+kz"  is a bond between two atoms
 * calculated the orientation of this bond using two angles alpha,beta
 * alpha = angle between x-axis and component of r onto xy plane
 * alpha = acos(x/sqrt(x**2+y**2))
 * beta = angle between x-axis and component of r onto xz plane
 * beta = acos(x/sqrt(x**2+z**2))
 * the box is devided into cells of lenth '~ 2*rcut'
 * atoms in each cell interacts with neighbouring cells (to avoid the double counting, reduced the neighbour cells from 26 to 13)
 * ******************************************************************************************************************* */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define num_atoms 72000
#define nnl_slots 72000
#define MAXCELLS_X 30
#define MAXCELLS_Y 30
#define MAXCELLS_Z 30
#define MAXATOMS_INACELL 1000

int main(int argc, char **argv)
{

  if(argc != 4)
  {
    printf("usage:./this_binary  rcut slots chkpt_file.chkpt \n");
    exit(0);
  }
  double rcut = atof(argv[1]);
  int slots = atoi(argv[2]);
  char *chkptfile1 = argv[3];
  
  char Output1[255];
  sprintf(Output1,"%s.ProjAngSud",chkptfile1);
  
  FILE *FPchkptfile1;
  FILE *FPOutput1;
  
  static int dump_id,dump_type;
  static double dump_mass,dump_posx,dump_posy,dump_posz,dump_velx,dump_vely,dump_velz;
  
  static double head1_XX,head1_XY,head1_XZ,head1_YX,head1_YY,head1_YZ,head1_ZX,head1_ZY,head1_ZZ;
  static int id1[num_atoms],type1[num_atoms];
  static double mass1[num_atoms],posx1[num_atoms],posy1[num_atoms],posz1[num_atoms],velx1[num_atoms],vely1[num_atoms],velz1[num_atoms];
  static int cell1[num_atoms];
  static int cell1_atoms_count[MAXCELLS_X][MAXCELLS_Y][MAXCELLS_Z];			// this gives the maximum number of atoms in that cell
  static int cell1_atoms_list[MAXCELLS_X][MAXCELLS_Y][MAXCELLS_Z][MAXATOMS_INACELL];    // list of atom global ids in a cell
  static int cell1_max_atoms_count[100];
  static int proj_angle0[1000][1000],proj_angle1[1000][1000],proj_angle2[1000][1000];
  
  char line1[200];
  
  static int i,j,k;

  /* Initialization */
  for(i=0;i<MAXCELLS_X;i++)
  {
    for(j=0;j<MAXCELLS_Y;j++)
    {
      for(k=0;k<MAXCELLS_Z;k++)
      { 
	cell1_atoms_count[i][j][k]=0;
      }
    }
  }
  
  for (i=0;i<slots;i++)
    for(j=0;j<2*slots;j++)
    {
      proj_angle0[i][j]=0;
      proj_angle1[i][j]=0;
      proj_angle2[i][j]=0;
    }

  /* reading atoms from chkptfile1*/ 
  FPchkptfile1=fopen(chkptfile1,"r");
  while(fgets(line1, sizeof(line1),FPchkptfile1) != NULL){
        if(line1[0] != '#')
	{
	   sscanf(line1,"%d %d %lf %lf %lf %lf %lf %lf %lf",&dump_id, &dump_type, &dump_mass, &dump_posx, &dump_posy, &dump_posz, &dump_velx, &dump_vely, &dump_velz);
           id1[dump_id]=dump_id;
           type1[dump_id]=dump_type;
           mass1[dump_id]=dump_mass;
           posx1[dump_id]=dump_posx;
           posy1[dump_id]=dump_posy;
           posz1[dump_id]=dump_posz;
           velx1[dump_id]=dump_velx;
           vely1[dump_id]=dump_vely;
           velz1[dump_id]=dump_velz;
           
           
           //printf("%d %d %lf %lf %lf %lf\n",dump_id,type1[dump_id],mass1[dump_id],posx1[dump_id],posy1[dump_id],posz1[dump_id]);fflush(stdout);
        }
        else if(line1[1] == 'X')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head1_XX, &head1_XY, &head1_XZ);
        }
        else if(line1[1] == 'Y')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head1_YX, &head1_YY, &head1_YZ);
        }
        else if(line1[1] == 'Z')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head1_ZX, &head1_ZY, &head1_ZZ);
        }
  }
  fclose(FPchkptfile1); 
  /* ********************************************* */

  /* making sample box into small volumes (cells) based on rcut */
  static int no_cellsX,no_cellsY,no_cellsZ;
  static double cellX_width,cellY_width,cellZ_width;
  no_cellsX=head1_XX/(2*rcut); // total no of cells along X 
  no_cellsY=head1_YY/(2*rcut); // total no of cells along Y
  no_cellsZ=head1_ZZ/(2*rcut); // total no of cells along Z
  
  cellX_width=head1_XX/no_cellsX;  // each cell width along X
  cellY_width=head1_YY/no_cellsY;  // each cell width along Y
  cellZ_width=head1_ZZ/no_cellsZ;  // each cell width along Z
  
  /* assign the atoms into respective cell_number */
  static int calc_cellX,calc_cellY,calc_cellZ,cell_number;
  for(i=0;i<num_atoms;i++)
  {
     calc_cellX=posx1[i]/cellX_width;  // cell number along X the atom belongs to (i.e 0 to no_cellsX-1)
     calc_cellY=posy1[i]/cellY_width;
     calc_cellZ=posz1[i]/cellZ_width;
     
     cell1_atoms_list[calc_cellX][calc_cellY][calc_cellZ][cell1_atoms_count[calc_cellX][calc_cellY][calc_cellZ]]=i;
     cell1_atoms_count[calc_cellX][calc_cellY][calc_cellZ]=cell1_atoms_count[calc_cellX][calc_cellY][calc_cellZ]+1;
  }
  
  /* main loop for interactions of each atom in one cell with other atoms in the same cell and also with the neighbouring cells */
  static int atom1,atom2,atom1_id,atom2_id;
  static double pos_x1,pos_y1,pos_z1,pos_x2,pos_y2,pos_z2;
  static int new_i,new_j,new_k;
  static int i1,j1;
  static int radius,radius2;
  static double dx,dy,dz;
  static double temp1,temp2,temp3,temp4;
  static double alpha,beta;
  static int alpha_slot,beta_slot;
  
  for(k=0;k<no_cellsZ;k++)  // cell belongs to first atom among interacting two atoms 
  {
    for(j=0;j<no_cellsY;j++)
    {
      for(i=0;i<no_cellsX;i++)
      {   
	for(atom1=0;atom1<cell1_atoms_count[i][j][k];atom1++)  // all atoms in the cell
	{
	  atom1_id=cell1_atoms_list[i][j][k][atom1];
	  pos_x1=posx1[atom1_id];
	  pos_y1=posy1[atom1_id];
	  pos_z1=posz1[atom1_id];
	  
	  // below loop is for new_k =k
	  new_k=k;
	  for(i1=0;i1<=1;i1++)
	  {
	    for(j1=-1;j1<=1;j1++)
	    {
	      if((i1+j1)>=0)
	      {
		new_i=i+i1;			// dealing PBC
		new_j=j+j1;			// dealing PBC
		if(new_i<0)			// dealing PBC
		  new_i=no_cellsX-1;		// dealing PBC
		if(new_j<0)			// dealing PBC
		  new_j=no_cellsY-1;		// dealing PBC
		
		for(atom2=0;atom2<cell1_atoms_count[new_i][new_j][new_k];atom2++)
		{
		  atom2_id=cell1_atoms_list[new_i][new_j][new_k][atom2];
		  if(atom1_id != atom2_id)
		  {
		    pos_x2=posx1[atom2_id];
		    pos_y2=posy1[atom2_id];
		    pos_z2=posz1[atom2_id];
		    
		    dx=pos_x1-pos_x2;
		    dy=pos_y1-pos_y2;
		    dz=pos_z1-pos_z2;
		    radius2=(dx*dx+dy*dy+dz*dz);
		    radius=sqrt(radius2);
		    if(radius < rcut)
		    {
		      if(dy<0)
		      {
			dx=-1*dx;
			dy=-1*dy;
			dz=-1*dz;
		      }
		      temp1=dx*dx+dy*dy;
		      temp2=sqrt(temp1);
		      alpha = acos(dx/temp2)*180/3.141592654;
		      alpha_slot=(int)(slots * alpha / 180);
		      
		      temp3=dx*dx+dz*dz;
		      temp4=sqrt(temp3);
		      beta = acos(dx/temp4)*180/3.141592654;
		      if(dz<0)
			beta=360-beta;
		      beta_slot=(int)(slots * beta / 180);
		      if((type1[atom1_id]==0)&&(type1[atom2_id]==0))
			proj_angle0[alpha_slot][beta_slot]=proj_angle0[alpha_slot][beta_slot]+1;
		      if(((type1[atom1_id]==0)&&(type1[atom2_id]==1)) || ((type1[atom1_id]==1)&&(type1[atom2_id]==0)))
			proj_angle1[alpha_slot][beta_slot]=proj_angle1[alpha_slot][beta_slot]+1;
		      if((type1[atom1_id]==1)&&(type1[atom2_id]==1))
			proj_angle2[alpha_slot][beta_slot]=proj_angle2[alpha_slot][beta_slot]+1;
		      
		    }
		      
		  }		  		  
		}
	      }   
	    }
	  }
	  
	  // below loop is for new_k=k+1
	  new_k=k+1;
	  for(i1=-1;i1<=1;i1++)
	  {
	    for(j1=-1;j1<=1;j1++)
	    {
		new_i=i+i1;			// dealing PBC
		new_j=j+j1;			// dealing PBC
		if(new_i<0)			// dealing PBC
		  new_i=no_cellsX-1;		// dealing PBC
		if(new_j<0)			// dealing PBC
		  new_j=no_cellsY-1;		// dealing PBC
		
		for(atom2=0;atom2<cell1_atoms_count[new_i][new_j][new_k];atom2++)
		{
		  atom2_id=cell1_atoms_list[new_i][new_j][new_k][atom2];
		  if(atom1_id != atom2_id)
		  {
		    pos_x2=posx1[atom2_id];
		    pos_y2=posy1[atom2_id];
		    pos_z2=posz1[atom2_id];
		    
		    dx=pos_x1-pos_x2;
		    dy=pos_y1-pos_y2;
		    dz=pos_z1-pos_z2;
		    radius2=(dx*dx+dy*dy+dz*dz);
		    radius=sqrt(radius2);
		    if(radius < rcut)
		    {
		      if(dy<0)
		      {
			dx=-1*dx;
			dy=-1*dy;
			dz=-1*dz;
		      }
		      temp1=dx*dx+dy*dy;
		      temp2=sqrt(temp1);
		      alpha = acos(dx/temp2)*180/3.141592654;
		      alpha_slot=(int)(slots * alpha / 180);
		      
		      temp3=dx*dx+dz*dz;
		      temp4=sqrt(temp3);
		      beta = acos(dx/temp4)*180/3.141592654;
		      if(dz<0)
			beta=360-beta;
		      beta_slot=(int)(slots * beta / 180);
		      if((type1[atom1_id]==0)&&(type1[atom2_id]==0))
			proj_angle0[alpha_slot][beta_slot]=proj_angle0[alpha_slot][beta_slot]+1;
		      if(((type1[atom1_id]==0)&&(type1[atom2_id]==1)) || ((type1[atom1_id]==1)&&(type1[atom2_id]==0)))
			proj_angle1[alpha_slot][beta_slot]=proj_angle1[alpha_slot][beta_slot]+1;
		      if((type1[atom1_id]==1)&&(type1[atom2_id]==1))
			proj_angle2[alpha_slot][beta_slot]=proj_angle2[alpha_slot][beta_slot]+1;
		      
		    }
		  }	  
		} 
	    }
	  }
	}
      }
    }
  }

  /* writing out the data */
  int i_slot,j_slot;
  FPOutput1=fopen(Output1,"w");
  fprintf(FPOutput1,"#alpha=x-axis to (ix,jy,0); beta=x-axis to (ix,0,kz) \n");
  fprintf(FPOutput1,"# %s %s %s %s \n",argv[0],argv[1],argv[2],argv[3]);
  for(i=0;i<slots;i++)
  {
      i_slot=i*180/slots;
      for(j=0;j<2*slots;j++)
      {
	j_slot=j*180/slots;
	fprintf(FPOutput1,"%d %d %d %d %d\n",i_slot,j_slot,proj_angle0[i][j],proj_angle1[i][j],proj_angle2[i][j]);
      }
  }
  
  
  fclose(FPOutput1);
 
  printf ("Sudheer, it is succeeful \n");
  return 0;
}
