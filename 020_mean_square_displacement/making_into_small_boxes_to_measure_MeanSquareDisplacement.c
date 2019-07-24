/* *******************************************************************************************************************
 * Sudheer Ganisetti, Mo 21. Mär 21:58:04 CET 2016
 * A small piece of code to calculate the Mean Square Displacement (MSD)
 * How I am calculating the MSD
 * 1) Take reference sample and make it into voxels (of given size i.e VoxelSizeX,VoxelSizeY, VoxelSizeZ)
 * 2) Take the Current sample
 * 3) Deal with periodic boundary conditions
 * 4) calculate SD of each atom
 * 5) MSD = Average all SD's in a reference voxel
 * 6) Output SD and MSD to CurChkpt
 *********************************************************************************************************************/
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

int main(int argc, char **argv)
{
  
  if(argc != 6)
  {
    printf("usage:./this_binary VoxelSizeX VoxelSizeY VoxelSizeZ ReferenceChkpt CurrentChkpt \n");
    exit(0);
  }

  double VoxelSizeX = atof(argv[1]);
  double VoxelSizeY = atof(argv[2]);
  double VoxelSizeZ = atof(argv[3]);
  char *RefChkpt1 = argv[4];
  char *CurChkpt1 = argv[5];
  printf("%lf %lf %lf \n",VoxelSizeX,VoxelSizeY,VoxelSizeZ);

  char OutChkpt1[255];
  sprintf(OutChkpt1,"%s.MSD",CurChkpt1);

  FILE *FPRefChkpt1;
  FILE *FPCurChkpt1;
  FILE *FPOutChkpt1;
 
  char line1[200];
  char line2[200];


  int 		 i;
  static int     dump_id=0,dump_type=0;
  static double  dump_mass=0.0,dump_posx=0.0,dump_posy=0.0,dump_posz=0.0,dump_velx=0.0,dump_vely=0.0,dump_velz=0.0;
  
  static double  head1_XX, head1_XY, head1_XZ, head1_YX, head1_YY, head1_YZ, head1_ZX, head1_ZY, head1_ZZ;
  static double  head2_XX, head2_XY, head2_XZ, head2_YX, head2_YY, head2_YZ, head2_ZX, head2_ZY, head2_ZZ;
  
  static int     id1[num_atoms];
  static int     type1[num_atoms];
  static double  mass1[num_atoms];
  static double  posx1[num_atoms],posy1[num_atoms],posz1[num_atoms];
  static double  velx1[num_atoms],vely1[num_atoms],velz1[num_atoms];
  static double  MSD[num_atoms],SD[num_atoms];
  static int     voxel[num_atoms],AtomsInVoxel[num_atoms];
  static double  MSD_voxel[num_atoms];

  static int     id2[num_atoms];
  static int     type2[num_atoms];
  static double  mass2[num_atoms];
  static double  posx2[num_atoms],posy2[num_atoms],posz2[num_atoms];
  static double  velx2[num_atoms],vely2[num_atoms],velz2[num_atoms];

  /* ********************************************* */
  /* Initialization                                */
  for(i=0;i<num_atoms;i++)
  {
    AtomsInVoxel[i]=0;
    MSD_voxel[i]=0;
  }
  /* ----------------------------------------------*/

  /* ********************************************* */
  /*  Reading Reference Chkpt                      */
  FPRefChkpt1=fopen(RefChkpt1,"r");
  while(fgets(line1, sizeof(line1),FPRefChkpt1) != NULL){
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
  fclose(FPRefChkpt1);
  /* ----------------------------------------------*/

  /* ********************************************* */
  /*  Reading Current Chkpt                        */
  FPCurChkpt1=fopen(CurChkpt1,"r");
  while(fgets(line1, sizeof(line1),FPCurChkpt1) != NULL){
          if(line1[0] != '#')
        {
                sscanf(line1,"%d %d %lf %lf %lf %lf %lf %lf %lf",&dump_id, &dump_type, &dump_mass, &dump_posx, &dump_posy, &dump_posz, &dump_velx, &dump_vely, &dump_velz);
                id2[dump_id]=dump_id;
                type2[dump_id]=dump_type;
                mass2[dump_id]=dump_mass;
                posx2[dump_id]=dump_posx;
                posy2[dump_id]=dump_posy;
                posz2[dump_id]=dump_posz;
                velx2[dump_id]=dump_velx;
                vely2[dump_id]=dump_vely;
                velz2[dump_id]=dump_velz;
                //printf("%d %d %lf %lf %lf %lf\n",dump_id,type1[dump_id],mass1[dump_id],posx1[dump_id],posy1[dump_id],posz1[dump_id]);fflush(stdout);
        }
        else if(line1[1] == 'X')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head2_XX, &head2_XY, &head2_XZ);
        }
        else if(line1[1] == 'Y')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head2_YX, &head2_YY, &head2_YZ);
        }
        else if(line1[1] == 'Z')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head2_ZX, &head2_ZY, &head2_ZZ);
        }
  }
  fclose(FPCurChkpt1);
  /* ----------------------------------------------*/

  /* ********************************************* */
  /*  Calculating MSD                              */
  int XSlices,YSlices,ZSlices;
  int TotalVoxels;
  double XSlices_temp,YSlices_temp,ZSlices_temp;
  double radius1,radius2;
  double new_posx2,new_posy2,new_posz2;

  XSlices_temp=head1_XX/VoxelSizeX;
  YSlices_temp=head1_YY/VoxelSizeY;
  ZSlices_temp=head1_ZZ/VoxelSizeZ;

  XSlices=XSlices_temp;
  YSlices=YSlices_temp;
  ZSlices=ZSlices_temp;

  if(XSlices < XSlices_temp)
  XSlices=XSlices+1;
 
  if(YSlices < YSlices_temp)
  YSlices=YSlices+1;

  if(ZSlices < ZSlices_temp)
  ZSlices=ZSlices+1;

  TotalVoxels=XSlices*YSlices*ZSlices;
  /* calculate square displacesment */
  for(i=0;i<num_atoms;i++)
  { 
    /* Dealing with periodic boundary conditions */
    if((posx1[i]<10)&&(posx2[i]>(head2_XX-10)))
      new_posx2=posx2[i]-head2_XX;
    else if((posx1[i]>(head1_XX-10))&&(posx2[i]<10))
      new_posx2=head2_XX-posx2[i];
    else
      new_posx2=posx2[i];
      
    if((posy1[i]<10)&&(posy2[i]>(head2_YY-10)))
      new_posy2=posy2[i]-head2_YY;
    else if((posy1[i]>(head1_YY-10))&&(posy2[i]<10))
      new_posy2=head2_YY-posy2[i];
    else
      new_posy2=posy2[i];
    
    if((posz1[i]<10)&&(posz2[i]>(head2_ZZ-10)))
      new_posz2=posz2[i]-head2_ZZ;
    else if((posz1[i]>(head1_ZZ-10))&&(posz2[i]<10))
      new_posz2=head2_ZZ-posz2[i];
    else
      new_posz2=posz2[i];
  
    /* ------------------------------------------ */
    
    radius1=sqrt((posx1[i]*posx1[i])+(posy1[i]*posy1[i])+(posz1[i]*posz1[i]));
    radius2=sqrt((new_posx2*new_posx2)+(new_posy2*new_posy2)+(new_posz2*new_posz2));
    SD[i]=(radius2-radius1)*(radius2-radius1);
  }

  /* now calculate mean of these square displacements over Voxel */
  int posx1_voxel;
  int posy1_voxel;
  int posz1_voxel;
  int VoxelNumber;

  for(i=0;i<num_atoms;i++)
  {
    posx1_voxel=posx1[i]/VoxelSizeX;
    posy1_voxel=posy1[i]/VoxelSizeY;
    posz1_voxel=posz1[i]/VoxelSizeZ;
    VoxelNumber=posz1_voxel*XSlices*YSlices+posy1_voxel*XSlices+posx1_voxel;
    voxel[i]=VoxelNumber;
    MSD_voxel[VoxelNumber]=MSD_voxel[VoxelNumber]+SD[i];
    AtomsInVoxel[VoxelNumber]=AtomsInVoxel[VoxelNumber]+1;
  }

  for(i=0;i<num_atoms;i++)
  {
   MSD[i]=MSD_voxel[voxel[i]]/AtomsInVoxel[voxel[i]];
  }
  /* ----------------------------------------------*/

  /* ********************************************* */
  /*  Writing out MSD to Current Chkpt             */
  FPOutChkpt1= fopen(OutChkpt1,"w");
  fprintf(FPOutChkpt1,"#F A 1 1 1 3 3 2 \n#C number type mass x y z vx vy vz SD MSD\n");
  fprintf(FPOutChkpt1,"#X   %lf %lf %lf \n",head2_XX,head2_XY,head2_XZ);
  fprintf(FPOutChkpt1,"#Y   %lf %lf %lf \n",head2_YX,head2_YY,head2_YZ);
  fprintf(FPOutChkpt1,"#Z   %lf %lf %lf \n",head2_ZX,head2_ZY,head2_ZZ);
  fprintf(FPOutChkpt1,"#E MeanSquareDisp Over Voxel; VoxelSize: %s x %s x %s; Ref: %s \n",argv[1],argv[2],argv[3],RefChkpt1);

  for(i=0;i<num_atoms;i++)
  {
     fprintf(FPOutChkpt1,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",id2[i],type2[i],mass2[i],posx2[i],posy2[i],posz2[i],velx2[i],vely2[i],velz2[i],SD[i],MSD[i]);
  }
  
  /* ----------------------------------------------*/


  printf ("Sudheer, it is succeeful \n");
  return 0;
}








  
