/* *******************************************************************************************************************
 * Sudheer Ganisetti, Wed Feb  3 17:40:22 CET 2016 
 * A small piece of code to compare file1 and file2 then print the common atoms from both the files
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
    printf("usage:./this_binary  file1.data  file2.chkpt \n");
  }
  char *file1 = argv[1];
  char *file2 = argv[2];

  FILE *file1fp;
  FILE *file2fp;

  char line1[50];
  char line2[50];

  int i,j;
  static int atom1[num_atoms];
  static int atom2[num_atoms];
  int dummy_atom;


  /* Initialisation */
  for(i=0;i<=num_atoms;i++)
  {
    atom1[i]=-1;
    atom2[i]=-1;
  }

  /* reading atoms in file1 */
  file1fp=fopen(file1,"r");
  while(fgets(line1,sizeof(line1),file1fp) != NULL){
	sscanf(line1,"%ld",&dummy_atom);
        atom1[dummy_atom]=1;
  }

  /* reading atoms in file2 */
  file2fp=fopen(file2,"r");
  while(fgets(line2,sizeof(line2),file2fp) != NULL){
        sscanf(line2,"%d",&dummy_atom);
        atom2[dummy_atom]=1;
  }

  /* Now check for the common atoms in the two files and print them */
  for(i=0;i<num_atoms;i++)
  {
     if((atom1[i]==1)&&(atom2[i]==1))
     {
	printf("%ld \n",i);
     }
  }


  return 0;
}

