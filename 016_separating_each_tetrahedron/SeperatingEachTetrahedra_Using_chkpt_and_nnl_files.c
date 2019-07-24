/* *******************************************************************************************************************
 * Sudheer Ganisetti, Thu Jan 28 09:38:40 CET 2016 
 * A small piece of code to get selected atoms from a chkpt file which is useful for analysing Silica stucture
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
  
  if(argc != 3)
  {
    printf("usage:./this_binary  chkpt_file.chkpt nnlfile \n");
  }
  
  char *chkptfile1 = argv[1];
  char *inputdatafile1 = argv[2];
    
  FILE *chkptfile1fp;
  FILE *inputdatafile1fp;
  
  char outputchkptfile1[255];
  FILE *outputchkptfile1fp;
  
  
  char line1[200];
  char line2[200];

  int i,j;
  int atom_count;
  static int inputdatafile1_atom;
  static int selected_atoms[num_atoms];
  int chkpt_int;
  
  static double  head_X11, head_X12, head_X13, head_Y11, head_Y12, head_Y13, head_Z11, head_Z12, head_Z13;
  static int     dump_id=0,dump_type=0;
  static double  dump_mass=0.0,dump_posx=0.0,dump_posy=0.0,dump_posz=0.0,dump_velx=0.0,dump_vely=0.0,dump_velz=0.0;
  static int     dump_nnl0=0,dump_nnl1=0,dump_nnl2=0,dump_nnl3=0,dump_nnl4=0,dump_nnl5=0,dump_nnl6=0,dump_nnl7=0,dump_nnl8=0,dump_nnl9=0;
  
  
  static int     id1[num_atoms];
  static int     type1[num_atoms];
  static double  mass1[num_atoms];
  static double  posx1[num_atoms],posy1[num_atoms],posz1[num_atoms];
  static double  velx1[num_atoms],vely1[num_atoms],velz1[num_atoms];
  static int     nnl[num_atoms][11];
  static int     selected_si_atoms[num_atoms];

  for(i=0;i<num_atoms;i++)
  {
    selected_si_atoms[i]=0;
  }
    
  /* reading atoms from chkptfile1*/ 
  chkptfile1fp=fopen(chkptfile1,"r");
  while(fgets(line2, sizeof(line2),chkptfile1fp) != NULL){
        if(line2[0] != '#')
        {
                sscanf(line2,"%d %d %lf %lf %lf %lf %lf %lf %lf",&dump_id, &dump_type, &dump_mass, &dump_posx, &dump_posy, &dump_posz, &dump_velx, &dump_vely, &dump_velz);
                id1[dump_id]=dump_id;
                type1[dump_id]=dump_type;
                mass1[dump_id]=dump_mass;
                posx1[dump_id]=dump_posx;
                posy1[dump_id]=dump_posy;
                posz1[dump_id]=dump_posz;
                velx1[dump_id]=dump_velx;
                vely1[dump_id]=dump_vely;
                velz1[dump_id]=dump_velz;
                
                if(dump_type == 0)
                {
                 selected_si_atoms[dump_id]=1;
                }
                //printf("%d %d %lf %lf %lf %lf\n",dump_id,type1[dump_id],mass1[dump_id],posx1[dump_id],posy1[dump_id],posz1[dump_id]);fflush(stdout);
        }
        else if(line2[1] == 'X')
        {
          sscanf(line2+2,"%lf %lf %lf ",&head_X11, &head_X12, &head_X13);
        }
        else if(line2[1] == 'Y')
        {
          sscanf(line2+2,"%lf %lf %lf ",&head_Y11, &head_Y12, &head_Y13);
        }
        else if(line2[1] == 'Z')
        {
          sscanf(line2+2,"%lf %lf %lf ",&head_Z11, &head_Z12, &head_Z13);
        }
  }
  fclose(chkptfile1fp);
  
  /* ********************************************* */

  
  /* ********************************************* */
  /* reading atoms from inputdatafile1 */
  atom_count=0;
  inputdatafile1fp=fopen(inputdatafile1,"r");
  while(fgets(line1, sizeof(line1),inputdatafile1fp) != NULL){
        if(line1[0] != '#')
        {
          sscanf(line1,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld",&dump_id,&dump_nnl0,&dump_nnl1,&dump_nnl2,&dump_nnl3,&dump_nnl4,&dump_nnl5,&dump_nnl6,&dump_nnl7,&dump_nnl8,&dump_nnl9);
          nnl[dump_id][0]=dump_id;
          nnl[dump_id][1]=dump_nnl0;
          nnl[dump_id][2]=dump_nnl1;
          nnl[dump_id][3]=dump_nnl2;
          nnl[dump_id][4]=dump_nnl3;
          nnl[dump_id][5]=dump_nnl4;
          nnl[dump_id][6]=dump_nnl5;
          nnl[dump_id][7]=dump_nnl6;
          nnl[dump_id][8]=dump_nnl7;
          nnl[dump_id][9]=dump_nnl8;
          nnl[dump_id][10]=dump_nnl9;
         
          //printf("%d\n",selected_atoms[atom_count]);
          atom_count=atom_count+1;
        }
  }
  fclose(inputdatafile1fp);       
  
  /* ********************************************* */
  
  
  for(i=0;i<atom_count;i++)
  {
    if(selected_si_atoms[i]==1)
    {
      sprintf(outputchkptfile1,"%s.Si.%05d",chkptfile1,i);

  
      /* Now writing all neighbouring atoms of the seleccted atoms into outputfile1 */
      outputchkptfile1fp = fopen(outputchkptfile1,"w");
      fprintf(outputchkptfile1fp,"#F A 1 1 1 3 0 0 \n#C number type mass x y z\n");
      fprintf(outputchkptfile1fp,"#X   %lf %lf %lf\n",head_X11,head_X12,head_X13);
      fprintf(outputchkptfile1fp,"#Y   %lf %lf %lf\n",head_Y11,head_Y12,head_Y13);
      fprintf(outputchkptfile1fp,"#Z   %lf %lf %lf\n",head_Z11,head_Z12,head_Z13);
      fprintf(outputchkptfile1fp,"#E Separating Each Tetrahedra;script:SeperatingEachTetrahedra_Using_chkpt_and_nnl_files.c \n");
      for(j=0;j<11;j++)
      {
       if(nnl[i][j] != -1)
       {
          fprintf(outputchkptfile1fp,"%d %d %lf %lf %lf %lf %lf %lf %lf\n", id1[nnl[i][j]],type1[nnl[i][j]],mass1[nnl[i][j]],posx1[nnl[i][j]],posy1[nnl[i][j]],posz1[nnl[i][j]],velx1[nnl[i][j]],vely1[nnl[i][j]],velz1[nnl[i][j]]);
       }
      }
      fclose(outputchkptfile1fp);
    }
  }
  /* ********************************************* */
  
  
  printf ("Sudheer, it is succeeful \n");
  return 0;
}
