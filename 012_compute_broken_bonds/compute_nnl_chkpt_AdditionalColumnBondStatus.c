/* To compute the broken bonds information, new bonds information, switched bonds information
   Two chkpt files and corresponding nnl files must be provided input to this code
   It it will compare the two nnl information and based on that, it will write, how many bonds
   were broken, how many new bonds were formed and how many bonds were switched into the second 
   chkpt file */
/* An additional column (bond_status) is printed with the following values to be used for color coding 
 * 0 - no changes in bond
 * 1 - existing bonds were broken i.e lost bonds
 * 2 - additional bonds were formed i.e bond number was increased
 * 3 - n number of bonds were broken and n number of new bonds were formed 
14Sep2015 
By Erik's suggestion I changed the following numbering for the bond status
 * 0 - no changes in bond
 * 1 - new bonds were formed
 * 2 - no of bonds broken = no of bonds formed
 * 3 - bonds were broken
05Oct2015
By mistake switched bonds column is printing with zeros, however broken bonds, new bonds and bond status are correct
Now i enable the switched bonds column aswell
02Dec2015
Now I am increasing the neighbors list checking count from 6 neighbours to 8 neighbours
i.e If the number of nearest neighbours are more in the new configuration ==> New Bond (It doesnt matter, whether bonds are switched or not)
    If the number of nearest neighbours are less in the new configuration ==> Broken Bond (It doesnt matter, whether bonds are switched or not)
    If the number of nearest neighbours are same but neighbouring atoms are not same ==> Switched Bond
    
    
    broken_bonds    = number of broken bonds    (it does not depend on the bond_status)
    new_bonds       = number of new bonds       (it does not depend on the bond_status)
    switched_bonds  = number of switched bonds  (it does not depend on the bond_status)
    bond_status     = I didnt change
     * 0 - no changes in bond
     * 1 - new bonds were formed
     * 2 - no of bonds broken = no of bonds formed
     * 3 - bonds were broken
      
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define num_atoms 72000
#define MAX_NEIGHB 6
int main(int argc, char **argv)
{
char *inputdatafile1 = argv[1];
char *inputdatafile2 = argv[2];
char *nnlfile1 = argv[3];
char *nnlfile2 = argv[4];

FILE *fp1;
FILE *fp2;
FILE *fp3;
FILE *fp4;
FILE *fpout1;
//FILE *fpout2;
char outputdatafile1[255];
//char outputdatafile2[255];
sprintf(outputdatafile1,"%s.nnlchkpt",inputdatafile2);
//sprintf(outputdatafile2,"%s.nnlchkpt",inputdatafile2);

char line1[200];
char line2[200];
char line3[200];
char line4[200];

// header variables
double  head_X11, head_X12, head_X13, head_Y11, head_Y12, head_Y13, head_Z11, head_Z12, head_Z13;
double  head_X21, head_X22, head_X23, head_Y21, head_Y22, head_Y23, head_Z21, head_Z22, head_Z23;

int     dump_id=0,dump_type=0;
double  dump_mass=0.0,dump_posx=0.0,dump_posy=0.0,dump_posz=0.0;
static int     id1[num_atoms];
static int     id2[num_atoms];
static int     type1[num_atoms];
static int     type2[num_atoms];
static double  mass1[num_atoms];
static double  mass2[num_atoms];
static double  posx1[num_atoms],posy1[num_atoms],posz1[num_atoms];
static double  posx2[num_atoms],posy2[num_atoms],posz2[num_atoms];
static int     tot_nna[num_atoms],tot_nnb[num_atoms];
static int     broken_bonds[num_atoms],new_bonds[num_atoms],bonds_switched[num_atoms],bond_status[num_atoms];

int     i,type,mass,ni;
int     na_neighbours;
//int     printing1[num_atoms];
//int     printing2[num_atoms];
int	broken_bonds_count,new_bonds_count,count_tot_nn;

int     dump_na1=-1,dump_na2=1,dump_na3=-1,dump_na4=-1,dump_na5=-1,dump_na6=-1,dump_na7=-1,dump_na8=-1;
int     dump_nb1=-1,dump_nb2=-1,dump_nb3=-1,dump_nb4=-1,dump_nb5=-1,dump_nb6=-1,dump_nb7=-1,dump_nab=-1;
static long na[num_atoms],na1[num_atoms],na2[num_atoms],na3[num_atoms],na4[num_atoms],na5[num_atoms],na6[num_atoms],na7[num_atoms],na8[num_atoms];
static long nb[num_atoms],nb1[num_atoms],nb2[num_atoms],nb3[num_atoms],nb4[num_atoms],nb5[num_atoms],nb6[num_atoms],nb7[num_atoms],nb8[num_atoms];


for(i=0;i<num_atoms;i++)
{
id1[i]=0;
type1[i]=0;
mass1[i]=0;
posx1[i]=0;
posy1[i]=0;
posz1[i]=0;
id2[i]=0;
type2[i]=0;
mass2[i]=0;
posx2[i]=0;
posy2[i]=0;
posz2[i]=0;
na[i]=na1[i]=na2[i]=na3[i]=na4[i]=na5[i]=na6[i]=na7[i]=na8[i]=-1;
nb[i]=nb1[i]=nb2[i]=nb3[i]=nb4[i]=nb5[i]=nb6[i]=nb7[i]=nb8[i]=-1;
broken_bonds[i]=0;
new_bonds[i]=0;
bonds_switched[i]=0;
bond_status[i]=0;
//printing1[i]=0;
//printing2[i]=0;
}



fp1 = fopen(inputdatafile1,"r");
while(fgets(line1, sizeof(line1), fp1) != NULL){
	if(line1[0] != '#')
	{
		sscanf(line1,"%d %d %lf %lf %lf %lf",&dump_id, &dump_type, &dump_mass, &dump_posx, &dump_posy, &dump_posz);
		id1[dump_id]=dump_id;
                type1[dump_id]=dump_type;
                mass1[dump_id]=dump_mass;
                posx1[dump_id]=dump_posx;
                posy1[dump_id]=dump_posy;
                posz1[dump_id]=dump_posz;
                //printf("%d %d %lf %lf %lf %lf\n",dump_id,type1[dump_id],mass1[dump_id],posx1[dump_id],posy1[dump_id],posz1[dump_id]);fflush(stdout);
	}
	else if(line1[1] == 'X')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head_X11, &head_X12, &head_X13);
        }
        else if(line1[1] == 'Y')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head_Y11, &head_Y12, &head_Y13);
        }
        else if(line1[1] == 'Z')
        {
          sscanf(line1+2,"%lf %lf %lf ",&head_Z11, &head_Z12, &head_Z13);
        }
}
fclose(fp1);


fp2 = fopen(inputdatafile2,"r");
while(fgets(line2, sizeof(line2), fp2) != NULL){
        if(line2[0] != '#')
        {
                sscanf(line2,"%d %d %lf %lf %lf %lf",&dump_id,&dump_type,&dump_mass,&dump_posx,&dump_posy,&dump_posz);
                id2[dump_id]=dump_id;
                type2[dump_id]=dump_type;
                mass2[dump_id]=dump_mass;
                posx2[dump_id]=dump_posx;
                posy2[dump_id]=dump_posy;
                posz2[dump_id]=dump_posz;
                //printf("%d %d %lf %lf %lf %lf\n",dump_id,type2[dump_id],mass2[dump_id],posx2[dump_id],posy2[dump_id],posz2[dump_id]);fflush(stdout);
//                no_of_atoms2++;
        }
        else if(line2[1] == 'X')
        {
          sscanf(line2+2,"%lf %lf %lf ",&head_X21, &head_X22, &head_X23);
        }
        else if(line2[1] == 'Y')
        {
          sscanf(line2+2,"%lf %lf %lf ",&head_Y21, &head_Y22, &head_Y23);
        }
        else if(line2[1] == 'Z')
        {
          sscanf(line2+2,"%lf %lf %lf ",&head_Z21, &head_Z22, &head_Z23);
        }
        
}
fclose(fp2);


fp3 = fopen(nnlfile1,"r");
while(fgets(line3, sizeof(line3), fp3) != NULL){
        if(line3[0] != '#')
        {
                sscanf(line3,"%d %d %d %d %d %d %d %d %d",&dump_id, &dump_na1, &dump_na2, &dump_na3, &dump_na4, &dump_na5, &dump_na6, &dump_na7, &dump_na8);
                na[dump_id]=dump_id;
                na1[dump_id]=dump_na1;
                na2[dump_id]=dump_na2;
                na3[dump_id]=dump_na3;
                na4[dump_id]=dump_na4;
                na5[dump_id]=dump_na5;
                na6[dump_id]=dump_na6;
                na7[dump_id]=dump_na7;
                na8[dump_id]=dump_na8;
                //printf("%d %d %d %d\n",na1[dump_id],na2[dump_id],na3[dump_id],na4[dump_id]);
                //printf("na= %d %d %d %d %d %d %d \n",na[dump_id],na1[dump_id],na2[dump_id],na3[dump_id],na4[dump_id],na5[dump_id],na6[dump_id]);
        }
}
fclose(fp3);


fp4 = fopen(nnlfile2,"r");
while(fgets(line4, sizeof(line4), fp4) != NULL){
        if(line4[0] != '#')
        {
                sscanf(line4,"%d %d %d %d %d %d %d %d %d",&dump_id, &dump_na1, &dump_na2, &dump_na3, &dump_na4, &dump_na5, &dump_na6, &dump_na7, &dump_na8);
                nb[dump_id]=dump_id;
                nb1[dump_id]=dump_na1;
                nb2[dump_id]=dump_na2;
                nb3[dump_id]=dump_na3;
                nb4[dump_id]=dump_na4;
                nb5[dump_id]=dump_na5;
                nb6[dump_id]=dump_na6;
                nb7[dump_id]=dump_na7;
                nb8[dump_id]=dump_na8;
                
        }
}
fclose(fp4);

/* Counting total number of nearest neighbours for each atom in each file */
for(i=0;i<num_atoms;i++)
{
    count_tot_nn =0;
    
    if(na1[i] != -1)
    count_tot_nn++;
    if(na2[i] != -1)
    count_tot_nn++;
    if(na3[i] != -1)
    count_tot_nn++;
    if(na4[i] != -1)
    count_tot_nn++;
    if(na5[i] != -1)
    count_tot_nn++;
    if(na6[i] != -1)
    count_tot_nn++;
    if(na7[i] != -1)
    count_tot_nn++;
    if(na8[i] != -1)
    count_tot_nn++;
    
    tot_nna[i]=count_tot_nn;
}

for(i=0;i<num_atoms;i++)
{
    count_tot_nn =0;
    
    if(nb1[i] != -1)
    count_tot_nn++;
    if(nb2[i] != -1)
    count_tot_nn++;
    if(nb3[i] != -1)
    count_tot_nn++;
    if(nb4[i] != -1)
    count_tot_nn++;
    if(nb5[i] != -1)
    count_tot_nn++;
    if(nb6[i] != -1)
    count_tot_nn++;
    if(nb7[i] != -1)
    count_tot_nn++;
    if(nb8[i] != -1)
    count_tot_nn++;
    
    tot_nnb[i]=count_tot_nn;
}



for(i=0;i< num_atoms;i++)
{
  na_neighbours =0;
  // checking "na atom" nearest neighbours in file 1 are equal to "nb atom" nearest neighours in file 2, 
  
  broken_bonds_count=0;
  new_bonds_count=0;
  
  if(na1[i] != -1)
  if( (na1[i]==nb1[i]) || (na1[i]==nb2[i]) || (na1[i]==nb3[i]) || (na1[i]==nb4[i]) || (na1[i]==nb5[i]) || (na1[i]==nb6[i]) || (na1[i]==nb7[i]) ||(na1[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  if(na2[i] != -1)
  if( (na2[i]==nb1[i]) || (na2[i]==nb2[i]) || (na2[i]==nb3[i]) || (na2[i]==nb4[i]) || (na2[i]==nb5[i]) || (na2[i]==nb6[i]) || (na2[i]==nb7[i]) ||(na2[i]==nb8[i]) )
  {
  }
  else
  {
     broken_bonds_count++;
  }
  
  if(na3[i] != -1)
  if( (na3[i]==nb1[i]) || (na3[i]==nb2[i]) || (na3[i]==nb3[i]) || (na3[i]==nb4[i]) || (na3[i]==nb5[i]) || (na3[i]==nb6[i]) || (na3[i]==nb7[i]) ||(na3[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  if(na4[i] != -1)
  if( (na4[i]==nb1[i]) || (na4[i]==nb2[i]) || (na4[i]==nb3[i]) || (na4[i]==nb4[i]) || (na4[i]==nb5[i]) || (na4[i]==nb6[i]) || (na4[i]==nb7[i]) ||(na4[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  if(na5[i] != -1)
  if( (na5[i]==nb1[i]) || (na5[i]==nb2[i]) || (na5[i]==nb3[i]) || (na5[i]==nb4[i]) || (na5[i]==nb5[i]) || (na5[i]==nb6[i]) || (na5[i]==nb7[i]) ||(na5[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  if(na6[i] != -1)
  if( (na6[i]==nb1[i]) || (na6[i]==nb2[i]) || (na6[i]==nb3[i]) || (na6[i]==nb4[i]) || (na6[i]==nb5[i]) || (na6[i]==nb6[i]) || (na6[i]==nb7[i]) ||(na6[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  if(na7[i] != -1)
  if( (na7[i]==nb1[i]) || (na7[i]==nb2[i]) || (na7[i]==nb3[i]) || (na7[i]==nb4[i]) || (na7[i]==nb5[i]) || (na7[i]==nb6[i]) || (na7[i]==nb7[i]) ||(na7[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  if(na8[i] != -1)
  if( (na8[i]==nb1[i]) || (na8[i]==nb2[i]) || (na8[i]==nb3[i]) || (na8[i]==nb4[i]) || (na8[i]==nb5[i]) || (na8[i]==nb6[i]) || (na8[i]==nb7[i]) ||(na8[i]==nb8[i]) )
  {
  }
  else
  {
    broken_bonds_count++;
  }
  
  
  
  // Searching for whether the same nearest neighbors exist or changed, if changed print only the new nearest neighbor  
  
  if(nb1[i] != -1)
  if( (nb1[i]==na1[i]) || (nb1[i]==na2[i]) || (nb1[i]==na3[i]) || (nb1[i]==na4[i]) || (nb1[i]==na5[i]) ||(nb1[i]==na6[i]) || (nb1[i]==na7[i]) ||(nb1[i]==na8[i]) )
  {
  }
  else
  {
    new_bonds_count++;
  }
    
  if(nb2[i] != -1)
  if( (nb2[i]==na1[i]) || (nb2[i]==na2[i]) || (nb2[i]==na3[i]) || (nb2[i]==na4[i]) || (nb2[i]==na5[i]) ||(nb2[i]==na6[i]) || (nb2[i]==na7[i]) ||(nb2[i]==na8[i]) )
  {
  }
  else
  {
      new_bonds_count++;
  }     
  
  if(nb3[i] != -1)
  if( (nb3[i]==na1[i]) || (nb3[i]==na2[i]) || (nb3[i]==na3[i]) || (nb3[i]==na4[i]) || (nb3[i]==na5[i]) ||(nb3[i]==na6[i]) || (nb3[i]==na7[i]) ||(nb3[i]==na8[i]) )
  {
  }
  else
  {
      new_bonds_count++;
  }
  
  if(nb4[i] != -1)
  if( (nb4[i]==na1[i]) || (nb4[i]==na2[i]) || (nb4[i]==na3[i]) || (nb4[i]==na4[i]) || (nb4[i]==na5[i]) ||(nb4[i]==na6[i]) || (nb4[i]==na7[i]) ||(nb4[i]==na8[i]) )
  {
  }
  else
  {
      new_bonds_count++;
  }
  
  if(nb5[i] != -1)
  if( (nb5[i]==na1[i]) || (nb5[i]==na2[i]) || (nb5[i]==na3[i]) || (nb5[i]==na4[i]) || (nb5[i]==na5[i]) ||(nb5[i]==na6[i]) || (nb5[i]==na7[i]) ||(nb5[i]==na8[i]) )
  {
  }
  else
  {
      new_bonds_count++;
  }  
    
  if(nb6[i] != -1)
  if( (nb6[i]==na1[i]) || (nb6[i]==na2[i]) || (nb6[i]==na3[i]) || (nb6[i]==na4[i]) || (nb6[i]==na5[i]) ||(nb6[i]==na6[i]) || (nb6[i]==na7[i]) ||(nb6[i]==na8[i]) )
  {
  }
  else
  {
     new_bonds_count++;
  }
    
  if(nb7[i] != -1)
  if( (nb7[i]==na1[i]) || (nb7[i]==na2[i]) || (nb7[i]==na3[i]) || (nb7[i]==na4[i]) || (nb7[i]==na5[i]) ||(nb7[i]==na6[i]) || (nb7[i]==na7[i]) ||(nb7[i]==na8[i]) )
  {
  }
  else
  {
     new_bonds_count++;
  }
    
  if(nb8[i] != -1)
  if( (nb8[i]==na1[i]) || (nb8[i]==na2[i]) || (nb8[i]==na3[i]) || (nb8[i]==na4[i]) || (nb8[i]==na5[i]) ||(nb8[i]==na6[i]) || (nb8[i]==na7[i]) ||(nb8[i]==na8[i]))
  {
  }
  else
  {
     new_bonds_count++;
  }
  
  
  
  broken_bonds[i]=broken_bonds_count;
  new_bonds[i]=new_bonds_count;
  if((broken_bonds[i]>0)&&(new_bonds[i]>0))
  {
    if(broken_bonds[i]<=new_bonds[i])
    {
      bonds_switched[i]=broken_bonds[i];
    }
    else
    {
      bonds_switched[i]=new_bonds[i];
    }
    
  }


  if((broken_bonds[i]-new_bonds[i])>0)
  {
    bond_status[i]=3;
  }
  else if((broken_bonds[i]-new_bonds[i])<0)
  {
    bond_status[i]=1;
  }
  else if((broken_bonds[i] == new_bonds[i]) && broken_bonds[i] != 0)
  {
    bond_status[i]=2;
  }
  
  /*
  for(ni=0;ni<MAX_NEIGHB;ni++)
  {
   if((broken_bonds_count == ni) && (new_bonds_count >= ni))
   {
     bonds_switched[i] = ni;
     ni=MAX_NEIGHB;
   }
  }

  if(broken_bonds[i] != 0)
    bond_status[i]=1;
  if(new_bonds[i] != 0)
    bond_status[i]=2;
  if(bonds_switched[i] != 0)
    bond_status[i]=3;
  */
  
  }
  

  
  
  fpout1 = fopen(outputdatafile1,"w");
  fprintf(fpout1,"#F A 1 1 1 3 5 0 \n#C number type mass x y z tot_nn broken_bonds new_bonds switched_bonds bond_status\n");
  fprintf(fpout1,"#X   %lf %lf %lf\n",head_X21,head_X22,head_X23);
  fprintf(fpout1,"#Y   %lf %lf %lf\n",head_Y21,head_Y22,head_Y23);
  fprintf(fpout1,"#Z   %lf %lf %lf\n",head_Z21,head_Z22,head_Z23);
  fprintf(fpout1,"#E \n");
  for(i=0;i<num_atoms;i++)
  {

      fprintf(fpout1,"%d  %d  %lf  %lf  %lf  %lf  %d  %d  %d  %d %d \n",id2[nb[i]],type2[nb[i]],mass2[nb[i]],posx2[nb[i]],posy2[nb[i]],posz2[nb[i]],tot_nnb[i],broken_bonds[i],new_bonds[i],bonds_switched[i],bond_status[i]);

  }
  fclose(fpout1); 
  
//printf ("Sudheer, it is succeeful \n");

return 0;
}



