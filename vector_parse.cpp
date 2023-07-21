#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define axis_a 13.7052176656990365
#define axis_b 13.705218
#define axis_c 9.8739316857179276
#define o_number 48
#define si_number 24
#define alpha  90.0
#define beta   90.0
#define gamma  120.0
#define atom_number 72
#define arch   57.29578
#define pi     3.1415926535

int main(void)
{
  double thita,fi,aa,bb,r[500],r_sum,d,d_save,d_save1;
	double x[500],y[500],z[500],a[500],b[500],c[500];
	double x0[500],y0[500],z0[500],a0[500],b0[500],c0[500];
	double xc[500],yc[500],zc[500],ac[500],bc[500],cc[500];
	double ac_d[500],bc_d[500],cc_d[500],a0_d[500],b0_d[500],c0_d[500];
	float distance_a0_try,distance_b0_try,distance_c0_try;
	float distance_a0[150][150],distance_b0[150][150],distance_c0[150][150],overall_distance_c[150][150],vx[150][150],vy[150][150],vz[150][150];
	float overall_distance[150][150],projection[200][10],module_xyz;
	float va[150][150],vb[150][150],vc[150][150],dot_multiply;
	float overall_distance_try,d1,d2,d3,module_vector[180];
	float v_a[180],v_b[180],v_c[180],module[180],module_multiply;
	double v_x[180],v_y[180],v_z[180],c_bending[180],c_rocking[180];
	double bending_x[180],bending_y[180],bending_z[180];
	double rocking_x[180],rocking_y[180],rocking_z[180];
	double bond_length[180],cos_thita[180],r_save_o[20],bond_angle_o[180];
	float neighbor_distance[150][150],neighbor_distance_o[150][150];
	int i,ii,ii1,ii2,ii3,iii,n,tag_a,tag_b,tag_c,k,encut,box_judge;
	int i1,i2,i3,i4,i5,i6,i7,i8,index,step[20],s,s1,s2,s3,s4,s5,s6,s7,s8;
	int neighbor_number,neighbor_index_o[180][5],neighbor_index[180][5];
	int o_index,si_index,index_o,index_si,r4_n,o_step[20];
	int index_save1,index_save,r_save[20],si1,si2;
	char l1,l2,l3,jump[1024];
	FILE *fp,*gp,*hp,*jp;
	fp=fopen("POSCAR.phon","rb");
	printf("Reading POSCAR...\n");
	for(i=1;i<=8;i++)
	{
	  fgets(jump,512,fp);
	}
	for(i=1;i<=atom_number;i++)
    {
	  fscanf(fp,"%lf%lf%lf",&a0_d[i],&b0_d[i],&c0_d[i]);
	  fgets(jump,512,fp);
	  a0[i]=a0_d[i]*axis_a;
	  b0[i]=b0_d[i]*axis_b;
	  c0[i]=c0_d[i]*axis_c;
	  x0[i]=a0[i]+b0[i]*cos(gamma/arch)+c0[i]*cos(beta/arch);
	  y0[i]=b0[i]*sin(gamma/arch);
	  z0[i]=c0[i]*sin(beta/arch);
	  module[i]=sqrt(a0[i]*a0[i]+b0[i]*b0[i]+c0[i]*c0[i]);
	}
	fclose(fp);
	printf("starting calculating distances among atoms...\n");
	for(i=1;i<=atom_number;i++)
    {
	  for(ii=1;ii<=atom_number;ii++)
	  {
        if(i==ii)
	      continue;
	    if(i!=ii)
	    {
	      distance_a0[i][ii]=fabs(a0[i]-a0[ii]);
	      distance_b0[i][ii]=fabs(b0[i]-b0[ii]);
	      distance_c0[i][ii]=fabs(c0[i]-c0[ii]);
	      va[i][ii]=a0[i]-a0[ii];
	      vb[i][ii]=b0[i]-b0[ii];
	      vc[i][ii]=c0[i]-c0[ii];
	      overall_distance[i][ii]=sqrt(distance_a0[i][ii]*distance_a0[i][ii]+distance_b0[i][ii]*distance_b0[i][ii]+distance_c0[i][ii]*distance_c0[i][ii]);
	    for(d1=-1.0;d1<=1.0;d1++)
	    {
		distance_a0_try=fabs(d1*axis_a+a0[i]-a0[ii]);
		for(d2=-1.0;d2<=1.0;d2++)
		{
		  distance_b0_try=fabs(d2*axis_b+b0[i]-b0[ii]);
		  for(d3=-1.0;d3<=1.0;d3++)
		  {
		    distance_c0_try=fabs(d3*axis_c+c0[i]-c0[ii]);
		    overall_distance_try=sqrt(distance_a0_try*distance_a0_try+distance_b0_try*distance_b0_try+distance_c0_try*distance_c0_try);
		    if(overall_distance_try<overall_distance[i][ii])
		    {
		      overall_distance[i][ii]=overall_distance_try;
		      va[i][ii]=d1*axis_a+a0[i]-a0[ii];
		      vb[i][ii]=d2*axis_b+b0[i]-b0[ii];
		      vc[i][ii]=d3*axis_c+c0[i]-c0[ii];
		    }
		  //printf("%d %d %lf %lf %lf\n",i,ii,distance_a0[i][ii],distance_b0[i][ii],distance_c0[i][ii]);
		  }
		}
	    }
	    }
	    //overall_distance_c[i][ii]=sqrt(distance_x0[i][ii]*distance_x0[i][ii]+distance_y0[i][ii]*distance_y0[i][ii]+distance_z0[i][ii]*distance_z0[i][ii]);
	    printf("%d %d  %f\n",i,ii,overall_distance[i][ii]);
	    vx[i][ii]=va[i][ii]+vb[i][ii]*cos(gamma/arch)+vc[i][ii]*cos(beta/arch);
	    vy[i][ii]=vb[i][ii]*sin(gamma/arch);
	    vz[i][ii]=vc[i][ii]*sin(beta/arch);
	    overall_distance_c[i][ii]=sqrt(vx[i][ii]*vx[i][ii]+vy[i][ii]*vy[i][ii]+vz[i][ii]*vz[i][ii]);
	  }
	}
	printf("Starting classifying atom neighbors...\n");
	for(i=1;i<=atom_number;i++)
	{
	  neighbor_number=0;
	  for(ii=1;ii<=atom_number;ii++)
	  {
	    //printf("%d  %d\n",i,ii);
            if(i==ii)
	      continue;
	    if(i!=ii)
	    {
	      if(i<=o_number&&ii>o_number)
		{
		  if(neighbor_number<2)
		    {
		      neighbor_number=neighbor_number+1;
		      neighbor_distance[i][neighbor_number]=overall_distance[i][ii];
		      neighbor_index[i][neighbor_number]=ii;
		      //printf("Initial %d %d %d \n",i,neighbor_number,neighbor_index[i][neighbor_number]);
		    }
		  else
		    {
		      d_save=overall_distance[i][ii];
		      index_save=ii;
		      d_save1=overall_distance[i][ii];
		      index_save1=ii;
		      for(iii=1;iii<=2;iii++)
			{
			  if(neighbor_distance[i][iii]>d_save)
			    {
			      d_save1=neighbor_distance[i][iii];
			      index_save1=neighbor_index[i][iii];
			      neighbor_distance[i][iii]=d_save;
			      neighbor_index[i][iii]=index_save;
			      d_save=d_save1;
			      index_save=index_save1;
			      //printf(" %d %d %d \n",i,iii,neighbor_index[i][iii]);
			      //printf("%d  ",neighbor_index_o[i][iii]);
			    }
			}
		    }
		}
	      if(i>o_number&&ii<=o_number)
		{
		  if(neighbor_number<4)
		    {
		      neighbor_number=neighbor_number+1;
		      neighbor_distance_o[i][neighbor_number]=overall_distance[i][ii];
		      neighbor_index_o[i][neighbor_number]=ii;
		    }
		  else
		    {
		      d_save=overall_distance[i][ii];
		      index_save=ii;
		      d_save1=overall_distance[i][ii];
		      index_save1=ii;
		      for(iii=1;iii<=4;iii++)
			{
			  if(neighbor_distance_o[i][iii]>d_save)
			    {
			      d_save1=neighbor_distance_o[i][iii];
			      index_save1=neighbor_index_o[i][iii];
			      neighbor_distance_o[i][iii]=d_save;
			      neighbor_index_o[i][iii]=index_save;
			      d_save=d_save1;
			      index_save=index_save1;
			    }
			}
		    }
		}
	    }
	  }
	}
	for(i=o_number+1;i<=atom_number;i++)
	{
	  neighbor_number=0;
	  for(ii=o_number+1;ii<=atom_number;ii++)
	  {
	    if(i==ii)
	      continue;
	    else
	    {
	    for(i2=1;i2<=4;i2++)
	    {
		for(i3=1;i3<=4;i3++)
		{
		  if(neighbor_index_o[i][i2]==neighbor_index_o[ii][i3])
		  {
		    //printf("%d %d  ",i,ii);
		    if(neighbor_number<4)
		    {
		      neighbor_number=neighbor_number+1;
		      neighbor_distance[i][neighbor_number]=overall_distance[i][ii];
		      neighbor_index[i][neighbor_number]=ii;
		    }
		    else
		    {
		      d_save=overall_distance[i][ii];
		      index_save=ii;
		      d_save1=overall_distance[i][ii];
		      index_save1=ii;
		    for(iii=1;iii<=4;iii++)
		    {
			if(neighbor_distance[i][iii]>=d_save)
			{
			  d_save1=neighbor_distance[i][iii];
			  index_save1=neighbor_index[i][iii];
			  neighbor_distance[i][iii]=d_save;
			  neighbor_index[i][iii]=index_save;
			  d_save=d_save1;
			  index_save=index_save1;
			  //printf(" %d %d %d \n",i,iii,neighbor_index[i][iii]);
			  //printf("%d  ",neighbor_index_o[i][iii]);
			}
		    }
		    }
		  }
		}
	    }
	    }
	  }
	}
	printf("\nReading vibrational vectors...\n");
	hp=fopen("VECTR","rb");
	for(i=1;i<=atom_number;i++)
	{
	  fscanf(hp,"   %d %f  %f  %f\n",&index,&v_a[i],&v_b[i],&v_c[i]);
	  v_x[i]=v_a[i]+v_b[i]*cos(gamma/arch)+v_c[i]*cos(beta/arch);
	  v_y[i]=v_b[i]*sin(gamma/arch);
	  v_z[i]=v_c[i]*sin(beta/arch);
	  module_vector[i]=sqrt(v_x[i]*v_x[i]+v_y[i]*v_y[i]+v_z[i]*v_z[i]);
	  fgets(jump,512,hp);
	  fgets(jump,512,hp);
	  printf("   %d %f  %f  %f\n",index,v_a[i],v_b[i],v_c[i]);
	}
	fclose(hp);
	for(i=1;i<=o_number;i++)
	{
	  printf("\n%d ",i);
	  for(ii=1;ii<=2;ii++)
	  {
	    printf(" %d %d",ii,neighbor_index[i][ii]);
	  }
	}
	for(i=o_number+1;i<=atom_number;i++)
	{
	  printf("\n%d ",i);
	  for(ii=1;ii<=4;ii++)
	  {
	    printf(" %d_%d_%d",ii,neighbor_index[i][ii],neighbor_index_o[i][ii]);
	  }
	}
	for(i=1;i<=o_number;i++)
	{
	    si1=neighbor_index[i][1];
	    si2=neighbor_index[i][2];
	    //bond_length[i]=sqrt((x0[si1]-x0[si2])*(x0[si1]-x0[si2])+(y0[si1]-y0[si2])*(y0[si1]-y0[si2])+(z0[si1]-z0[si2])*(z0[si1]-z0[si2]));
	    bending_x[i]=vx[i][si1]/overall_distance[i][si1]+vx[i][si2]/overall_distance[i][si2];
	    bending_y[i]=vy[i][si1]/overall_distance[i][si1]+vy[i][si2]/overall_distance[i][si2];
	    bending_z[i]=vz[i][si1]/overall_distance[i][si1]+vz[i][si2]/overall_distance[i][si2];
	    rocking_x[i]=vy[i][si1]*vz[i][si2]-vz[i][si1]*vy[i][si2];
	    rocking_y[i]=vz[i][si1]*vx[i][si2]-vx[i][si1]*vz[i][si2];
	    rocking_z[i]=vx[i][si1]*vy[i][si2]-vy[i][si1]*vx[i][si2];
	    dot_multiply=v_x[i]*bending_x[i]+v_y[i]*bending_y[i]+v_z[i]*bending_z[i];
	    module_multiply=module_vector[i]*sqrt(bending_x[i]*bending_x[i]+bending_y[i]*bending_y[i]+bending_z[i]*bending_z[i]);
	    c_bending[i]=dot_multiply/module_multiply;
	    dot_multiply=v_x[i]*rocking_x[i]+v_y[i]*rocking_y[i]+v_z[i]*rocking_z[i];
	    module_multiply=module_vector[i]*sqrt(rocking_x[i]*rocking_x[i]+rocking_y[i]*rocking_y[i]+rocking_z[i]*rocking_z[i]);
	    c_rocking[i]=dot_multiply/module_multiply;
	    dot_multiply=vx[i][si1]*vx[i][si2]+vy[i][si1]*vy[i][si2]+vz[i][si1]*vz[i][si2];
	    module_multiply=overall_distance_c[i][si1]*overall_distance_c[i][si2];
	    cos_thita[i]=dot_multiply/module_multiply;
	    if(cos_thita[i]<-1.00&&cos_thita[i]>-1.01)
	      cos_thita[i]=-0.99999999999999;
	    printf("%d %f %f\n",i,module_multiply,dot_multiply);
	    bond_angle_o[i]=acos(cos_thita[i])*180/pi;
	}
	gp=fopen("EDI_vector_512.dat","wb");
	fprintf(gp,"index  bond_angle  c_bending  c_rocking module vector\n\n");
	for(i=1;i<=o_number;i++)
	{
	  si1=neighbor_index[i][1];
	  si2=neighbor_index[i][2];
	  fprintf(gp,"\n%d  %f  %f  %f  %f",i,bond_angle_o[i],c_bending[i],c_rocking[i],module_vector[i]);
	}
	for(i=o_number+1;i<=atom_number;i++)
	{
	  fprintf(gp,"\n%d            %f",i,module_vector[i]);
	}
	fclose(gp);
	return 0;
}
