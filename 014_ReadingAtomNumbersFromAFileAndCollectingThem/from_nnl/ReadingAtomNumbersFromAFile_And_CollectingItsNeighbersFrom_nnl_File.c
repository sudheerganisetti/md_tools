/* *******************************************************************************************************************
 * Sudheer Ganisetti, Thu Jan 28 09:38:40 CET 2016 
 * A small piece of code to get all neighbouring atoms of selected atoms which is useful for analysing Silica stucture
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
    printf("usage:./this_binary  selectedAtoms_file.data  nearest_neighbours_list_file.nnl \n");
  }
  char *inputdatafile1 = argv[1];
  char *nnlfile1 = argv[2];
  
  FILE *inputdatafile1fp;
  FILE *nnlfile1fp;

  char outputdatafile1[255];
  sprintf(outputdatafile1,"%s.nnlOfSelectedAtoms",inputdatafile1);
  FILE *outputdatafile1fp;
    
  char line1[48];
  char line2[200];

  int i,j;
  int atom_count;
  static int inputdatafile1_atom;
  static int selected_atoms[num_atoms];
  static int final_selection_of_atoms[num_atoms];
  int inputdatafile1_total_atoms;
  int nnl_int;
  
  static long     nnl[num_atoms][11];
  static long     dump_nnl,dump_nnl1,dump_nnl2,dump_nnl3,dump_nnl4,dump_nnl5,dump_nnl6,dump_nnl7,dump_nnl8,dump_nnl9,dump_nnl10;

  /* Initialising all "nnl[][]" values with -1 */
  for(i=0;i<num_atoms;i++)
  {
   for(j=0;j<11;j++)
   {
     nnl[i][j]=-1;
   }
   final_selection_of_atoms[i]=0;
  }
    
  /* ********************************************* */
  
  
  /* reading atoms from inputdatafile1 */
  atom_count=0;
  inputdatafile1fp=fopen(inputdatafile1,"r");
  while(fgets(line1, sizeof(line1),inputdatafile1fp) != NULL){
        if(line1[0] != '#')
        {
          sscanf(line1,"%ld",&inputdatafile1_atom);
          selected_atoms[atom_count]=inputdatafile1_atom;
          atom_count=atom_count+1;
        }
  }
  fclose(inputdatafile1fp);       
  inputdatafile1_total_atoms=atom_count;
  
  /* ********************************************* */
  
  
  /* reading atoms from nnlfile1*/ 
  nnl_int=0;
  nnlfile1fp=fopen(nnlfile1,"r");
  while(fgets(line2, sizeof(line2),nnlfile1fp) != NULL){
        if(line2[0] != '#')
        {
                sscanf(line2,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld",&dump_nnl, &dump_nnl1, &dump_nnl2, &dump_nnl3, &dump_nnl4, &dump_nnl5, &dump_nnl6, &dump_nnl7, &dump_nnl8, &dump_nnl9, &dump_nnl10);
                
                if(dump_nnl != -1)
                {
                nnl[dump_nnl][0]=dump_nnl;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl1 != -1)
                {
                nnl[dump_nnl][1]=dump_nnl1;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl2 != -1)
                {
                nnl[dump_nnl][2]=dump_nnl2;
                //nnl_int=nnl_int+1;
                }
                
                if(dump_nnl3 != -1)
                {
                nnl[dump_nnl][3]=dump_nnl3;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl4 != -1)
                {
                nnl[dump_nnl][4]=dump_nnl4;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl5 != -1)
                {
                nnl[dump_nnl][5]=dump_nnl5;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl6 != -1)
                {
                nnl[dump_nnl][6]=dump_nnl6;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl7 != -1)
                {
                nnl[dump_nnl][7]=dump_nnl7;
                //nnl_int=nnl_int+1;
                }
                
                if(dump_nnl8 != -1)
                {
                nnl[dump_nnl][8]=dump_nnl8;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl9 != -1)
                {
                nnl[dump_nnl][9]=dump_nnl9;
                //nnl_int=nnl_int+1;
                }

                if(dump_nnl10 != -1)
                {
                nnl[dump_nnl][10]=dump_nnl10;
                //nnl_int=nnl_int+1;
                }

        }
  }
  fclose(nnlfile1fp);
  
  /* ********************************************* */
  
  
  /* Selecting all nnl atoms (including the atoms from file 1) */
  for(i=0;i<atom_count;i++)
  {
    for(j=0;j<11;j++)
    {
      if(nnl[selected_atoms[i]][j] != -1)
      {
        final_selection_of_atoms[nnl[selected_atoms[i]][j]]=1;
      }
    }
  }
  
  /* ********************************************* */
  
  
  /* Now writing all neighbouring atoms of the seleccted atoms into outputfile1 */
  outputdatafile1fp = fopen(outputdatafile1,"w");
  for(i=0;i<num_atoms;i++)
  {
   if(final_selection_of_atoms[i] == 1)
   {
     fprintf(outputdatafile1fp,"%d \n",i);
   }
  }
  //for(i=0;i<atom_count;i++)
  //{
  //  for(j=0;j<11;j++)
  //  {
  //    if(nnl[selected_atoms[i]][j] != -1)
  //    {
  //      fprintf(outputdatafile1fp,"%d \n",nnl[selected_atoms[i]][j]);
  //    }
  //  }
  //}
  fclose(outputdatafile1fp);
  
  /* ********************************************* */
  
  
  printf ("Sudheer, it is succeeful \n");
  return 0;
}
