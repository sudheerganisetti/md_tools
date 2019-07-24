/* *******************************************************************************************************************
 * Sudheer Ganisetti, Thu Jan 28 09:38:40 CET 2016 
 * A small piece of code to get selected atoms from a chkpt file which is useful for analysing Silica stucture
 *********************************************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define num_atoms 72000
#define nnl_slots 72000

int main(int argc, char **argv)
{
  
  if(argc != 3)
  {
    printf("usage:./this_binary  selectedAtoms_file.data  chkpt_file.chkpt \n");
  }
  char *inputdatafile1 = argv[1];
  char *chkptfile1 = argv[2];
  
  FILE *inputdatafile1fp;
  FILE *chkptfile1fp;

  char outputchkptfile1[255];
  sprintf(outputchkptfile1,"%s.SelectedAtoms",chkptfile1);
  FILE *outputchkptfile1fp;
    
  char line1[48];
  char line2[200];

  int i,j;
  int atom_count;
  static int inputdatafile1_atom;
  static int selected_atoms[num_atoms];
  int chkpt_int;
  
  static double  head_X11, head_X12, head_X13, head_Y11, head_Y12, head_Y13, head_Z11, head_Z12, head_Z13;
  static int     dump_id=0,dump_type=0,dump_nn=0,dump_nb=0,dump_sb=0,dump_bb=0,dump_bs=0;
  static double  dump_mass=0.0,dump_posx=0.0,dump_posy=0.0,dump_posz=0.0;
  
  static int     id1[num_atoms];
  static int     type1[num_atoms];
  static double  mass1[num_atoms];
  static double  posx1[num_atoms],posy1[num_atoms],posz1[num_atoms];
  static short int nn1[num_atoms],nb1[num_atoms],sb1[num_atoms],bb1[num_atoms],bs1[num_atoms];


    
  /* ********************************************* */
  
  
  /* reading atoms from inputdatafile1 */
  atom_count=0;
  inputdatafile1fp=fopen(inputdatafile1,"r");
  while(fgets(line1, sizeof(line1),inputdatafile1fp) != NULL){
        if(line1[0] != '#')
        {
          sscanf(line1,"%ld",&inputdatafile1_atom);
          selected_atoms[atom_count]=inputdatafile1_atom;
          //printf("%d\n",selected_atoms[atom_count]);
          atom_count=atom_count+1;
        }
  }
  fclose(inputdatafile1fp);       
  
  /* ********************************************* */
  
  
  /* reading atoms from chkptfile1*/ 
  chkptfile1fp=fopen(chkptfile1,"r");
  while(fgets(line2, sizeof(line2),chkptfile1fp) != NULL){
        if(line2[0] != '#')
        {
                sscanf(line2,"%d %d %lf %lf %lf %lf %d %d %d %d %d",&dump_id, &dump_type, &dump_mass, &dump_posx, &dump_posy, &dump_posz, &dump_nn, &dump_nb, &dump_sb, &dump_bb, &dump_bs);
                id1[dump_id]=dump_id;
                type1[dump_id]=dump_type;
                mass1[dump_id]=dump_mass;
                posx1[dump_id]=dump_posx;
                posy1[dump_id]=dump_posy;
                posz1[dump_id]=dump_posz;
                nn1[dump_id]=dump_nn;
                nb1[dump_id]=dump_nb;
                sb1[dump_id]=dump_sb;
                bb1[dump_id]=dump_bb;
                bs1[dump_id]=dump_bs;
                
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
  
  
  /* Now writing all neighbouring atoms of the seleccted atoms into outputfile1 */
  outputchkptfile1fp = fopen(outputchkptfile1,"w");
  fprintf(outputchkptfile1fp,"#F A 1 1 1 3 5 0 \n#C number type mass x y z nn new_bonds switched_bonds broken_bonds bond_status \n");
  fprintf(outputchkptfile1fp,"#X   %lf %lf %lf\n",head_X11,head_X12,head_X13);
  fprintf(outputchkptfile1fp,"#Y   %lf %lf %lf\n",head_Y11,head_Y12,head_Y13);
  fprintf(outputchkptfile1fp,"#Z   %lf %lf %lf\n",head_Z11,head_Z12,head_Z13);
  fprintf(outputchkptfile1fp,"#E Selected Atoms;script:ReadingAtomNumbersFromAFile_And_CollectingThemFrom_chkpt_File.c \n");
  
  for(i=0;i<atom_count;i++)
  {
   
     fprintf(outputchkptfile1fp,"%d %d %lf %lf %lf %lf %d %d %d %d %d \n",id1[selected_atoms[i]],type1[selected_atoms[i]],mass1[selected_atoms[i]],posx1[selected_atoms[i]],posy1[selected_atoms[i]],posz1[selected_atoms[i]],nn1[selected_atoms[i]],nb1[selected_atoms[i]],sb1[selected_atoms[i]],bb1[selected_atoms[i]],bs1[selected_atoms[i]]);
   
  }
  fclose(outputchkptfile1fp);
  
  /* ********************************************* */
  
  
  printf ("Sudheer, it is succeeful \n");
  return 0;
}
