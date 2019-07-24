#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define num_atoms 72000
#define nnl_slots 72000
int main(int argc, char **argv)
{

int i;
static int nnl_int=0;
int total_nnl_int;
char *inputdatafile1 = argv[1];
char *nnlfile1 = argv[2];


FILE *inputdatafile1fp;
FILE *nnlfile1fp;


char outputdatafile1[255];
sprintf(outputdatafile1,"%s.modchkpt",inputdatafile1);
FILE *outputdatafile1fp;


char line1[200];
char line2[200];

static int     id1[num_atoms];
static int     type1[num_atoms];
static double  mass1[num_atoms];
static double  posx1[num_atoms],posy1[num_atoms],posz1[num_atoms];

//int nnl_slots=4*num_atoms;
static long     nnl[nnl_slots];
static long     dump_nnl,dump_nnl1,dump_nnl2,dump_nnl3,dump_nnl4,dump_nnl5,dump_nnl6,dump_nnl7,dump_nnl8,dump_nnl9,dump_nnl10;

static double  head_X11, head_X12, head_X13, head_Y11, head_Y12, head_Y13, head_Z11, head_Z12, head_Z13;
static int     dump_id=0,dump_type=0;
static double  dump_mass=0.0,dump_posx=0.0,dump_posy=0.0,dump_posz=0.0;

for(nnl_int=0;nnl_int<nnl_slots;nnl_int++)
{
nnl[nnl_int]=-1;
}


inputdatafile1fp = fopen(inputdatafile1,"r");
while(fgets(line1, sizeof(line1), inputdatafile1fp) != NULL){
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
fclose(inputdatafile1fp);

nnl_int=0;
nnlfile1fp=fopen(nnlfile1,"r");
while(fgets(line2, sizeof(line2),nnlfile1fp) != NULL){
        if(line2[0] != '#')
        {
                sscanf(line2,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld",&dump_nnl, &dump_nnl1, &dump_nnl2, &dump_nnl3, &dump_nnl4, &dump_nnl5, &dump_nnl6, &dump_nnl7, &dump_nnl8, &dump_nnl9, &dump_nnl10);
		if(dump_nnl != -1)
		{
		nnl[nnl_int]=dump_nnl;
		nnl_int=nnl_int+1;
		}

                if(dump_nnl1 != -1)
                {
                nnl[nnl_int]=dump_nnl1;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl2 != -1)
                {
                nnl[nnl_int]=dump_nnl2;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl3 != -1)
                {
                nnl[nnl_int]=dump_nnl3;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl4 != -1)
                {
                nnl[nnl_int]=dump_nnl4;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl5 != -1)
                {
                nnl[nnl_int]=dump_nnl5;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl6 != -1)
                {
                nnl[nnl_int]=dump_nnl6;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl7 != -1)
                {
                nnl[nnl_int]=dump_nnl7;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl8 != -1)
                {
                nnl[nnl_int]=dump_nnl8;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl9 != -1)
                {
                nnl[nnl_int]=dump_nnl9;
                nnl_int=nnl_int+1;
                }

                if(dump_nnl10 != -1)
                {
                nnl[nnl_int]=dump_nnl10;
                nnl_int=nnl_int+1;
                }

		/*nnl[nnl_int+1]=dump_nnl1;
                nnl[nnl_int+2]=dump_nnl2;
                nnl[nnl_int+3]=dump_nnl3;
                nnl[nnl_int+4]=dump_nnl4;
                nnl[nnl_int+5]=dump_nnl5;
                nnl[nnl_int+6]=dump_nnl6;
                nnl[nnl_int+7]=dump_nnl7;
                nnl[nnl_int+8]=dump_nnl8;
                nnl[nnl_int+9]=dump_nnl9;
                nnl[nnl_int+10]=dump_nnl10;
                
                nnl_int=nnl_int+11; */


	}
}
fclose(nnlfile1fp);
total_nnl_int=nnl_int;

for(nnl_int=0;nnl_int<total_nnl_int;nnl_int++)
{
if(type1[nnl[nnl_int]]==0)
type1[nnl[nnl_int]]=2;
else if(type1[nnl[nnl_int]]==1)
type1[nnl[nnl_int]]=3;

}


  outputdatafile1fp = fopen(outputdatafile1,"w");
  fprintf(outputdatafile1fp,"#F A 1 1 1 3 0 0 \n#C number type mass x y z\n");
  fprintf(outputdatafile1fp,"#X   %lf %lf %lf\n",head_X11,head_X12,head_X13);
  fprintf(outputdatafile1fp,"#Y   %lf %lf %lf\n",head_Y11,head_Y12,head_Y13);
  fprintf(outputdatafile1fp,"#Z   %lf %lf %lf\n",head_Z11,head_Z12,head_Z13);
  fprintf(outputdatafile1fp,"#E The atom type (Originally Si=0, O=1) of selected atoms are modified (to Si=2, O=3) \n");
  for(i=0;i<num_atoms;i++)
  {

      fprintf(outputdatafile1fp,"%d  %d  %lf  %lf  %lf  %lf \n",id1[i],type1[i],mass1[i],posx1[i],posy1[i],posz1[i]);

  }
  fclose(outputdatafile1fp);


printf ("Sudheer, it is succeeful \n");


return 0;
}

