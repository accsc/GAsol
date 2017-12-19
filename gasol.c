/********************************************************************

GAsol rev6

Alvaro Cortes and Lucia Fusani
GlaxoSmithKline 2017

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors 
may be used to endorse or promote products derived from this software without 
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY 
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Compilation instructions:

gcc gasol.c -o gasol -O3 -lm -fopenmp

To select the number of threads please use:

export OMP_NUM_THREADS=12 (bash)
setenv OMP_NUM_THREADS 12 (csh/tcsh)

Otherwise, OpenMP will determine the number of CPUs autmatically and adjust the number
of threads accordingly.

Rev6 - Added water numbers check 
       Fixed bug in water-water distance threshold check

Rev5 - Added ratio parameter to allow very low density waters (ratio=population/radius)
       Renamed program to "GAsol" or Genetic Algorithm for SOLvent placement

Rev4 - Fixed bugs related to ligand parsing
       Added molarity parameter

Rev3 - First version based on population pre-calculation for speed
       Added option to read a PDB file with a ligand to set centre of the sphere

Please send any feedback to alvaro.x.cortes@gsk.com


********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>


void show_usage()
{

        fprintf(stderr,"Usage: gasol <params>\n\n");
        fprintf(stderr,"Valid params:\n");
        fprintf(stderr,"\t-d or --dx <DX file>. Set the grid filename.\n");
        fprintf(stderr,"\t-x <value>. Set the center of the spatial filter (optional).\n");
        fprintf(stderr,"\t-y <value>. Set the center of the spatial filter (optional).\n");
        fprintf(stderr,"\t-z <value>. Set the center of the spatial filter (optional).\n");
        fprintf(stderr,"\t-r or --radius <value>. Set the radius of the spatial filter (optional).\n");
        fprintf(stderr,"\t-t or --threshold <value>. Set the minimum density threshold (default 5.0).\n");
        fprintf(stderr,"\t-p or --population <value>. Set the population for the genetic algorithm (default number of genes x 10).\n");
        fprintf(stderr,"\t-i or --iterations <value>. Set the number of generations for the genetic algorithm (default 5000).\n");
        fprintf(stderr,"\t-l or --ligand <PDB file>. Use a ligand to set the center of spatial filter (optional).\n");
        fprintf(stderr,"\t-m or --molarity <value>. Define concentration in 3D-RISM calculation (default 55.5 M).\n");
        fprintf(stderr,"\t-s or --seed <value>. Seed for the random number generator (Default -1 for current time).\n");
        fprintf(stderr,"\t-c or --ratio <value>. Minimum ratio of density/radius to consider a water site (Default: 0.15).\n");

}


/* Generates a random number between 0 and 1*/
double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}

/* Random integer in an interval */
unsigned int rand_interval(unsigned int min, unsigned int max)
{
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    do
    {
        r = rand();
    } while (r >= limit);

    return min + (r / buckets);
}

int get_COM_ligand(char *name, float *x, float *y, float *z)
{
	FILE *in = NULL;
	char *buffer_line = NULL;
	float com[3], cx = 0, cy =0 , cz = 0;
	char tmp_number[20];
	int n_atoms = 0;

	com[0] = com[1] = com[2] - 0;
	if( (in = fopen(name,"rb")) == NULL)
	{
		fprintf(stderr,"Cannot open ligand file %s\n",name);
		fflush(stderr);
		return -1;
	}

        if ((buffer_line = (char *) calloc(sizeof(char),1024)) == NULL)
        {
                fprintf(stderr,"Error allocating memory\n");
                fflush(stderr);
                fclose(in);
                return -2;
        }


        while( (fgets(buffer_line,1023,in)) != NULL)
        {

		if( (buffer_line[0] == 'A' && buffer_line[1] == 'T' && buffer_line[2] == 'O' && buffer_line[3] == 'M') || ( buffer_line[0] == 'H'  && buffer_line[1] == 'E' && buffer_line[2] == 'T' && buffer_line[3] == 'A' && buffer_line[4] == 'T' && buffer_line[5] == 'M' ))
		{
			cx = cy = cz = 0.0f;
			strncpy(tmp_number,&buffer_line[30],8);
			tmp_number[8] = 0;
			cx = atof(tmp_number);
        	        strncpy(tmp_number,&buffer_line[38],8);
                        tmp_number[8] = 0;
                	cy = atof(tmp_number);
                	strncpy(tmp_number,&buffer_line[46],8);
                        tmp_number[8] = 0;
                	cz = atof(tmp_number);
			com[0] += cx;
			com[1] += cy;
			com[2] += cz;
			n_atoms++;
		}


	}

	if( n_atoms == 0)
	{
		fprintf(stderr,"PDB file contains no atoms\n");
		fflush(stderr);
		return -1;
	}

	com[0] /= (float) n_atoms;
	com[1] /= (float) n_atoms;
	com[2] /= (float) n_atoms;
	*x = com[0];
	*y = com[1];
	*z = com[2];

	free(buffer_line);
	fclose(in);

	return 0;
}


int main( int argc, char *argv[])
{

	FILE *dx_in = NULL;
	char *buffer_line = NULL;
	int lines = 0, nx = 0, ny = 0, nz = 0, i = 0, j = 0, k = 0, l =0, current_delta = 0, points = 0;
	float min[3], delta[3];
	char *token = NULL;
	int current_token = 0, state_read = 0, current_point = 0;
	float *data = NULL, *g = NULL, *max_g = NULL;
	char **line_tokens = NULL;
	int *x_index = NULL, *y_index = NULL, *z_index = NULL, n_points = 0;
	int **population = NULL, **offspring = NULL, *best_ever = NULL;
	double *fitness = NULL, *new_fitness = NULL, current_g = 0;
	int max_ind = 0, nind = 0, ngen = 0, max_ngen = 0, ig = 0, jg = 0, kg = 0, lg = 0;
        double best_fitness = -9999.0f, w_sum = 0, all_g = 0, full_integ = 0.0;
	float *points_radii = NULL;
	int **forbid_combination = NULL;
	float current_distance = 0.0f, current_population = 0.0f;
        double w[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 = 0;
	double d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0, d6 = 0;
	double penalty = 0, e = 0;
	float dx = 0.0f, dy = 0.0f, dz = 0.0f, dr = 0.0f;
	float concentration = 55.0;

	int t1 = -1, t2 = -1, t3 = -1, best_t = -1, bad_w = 0;
	int c1 = 0, c2 = 0, cxpoint1 = 0, cxpoint2 = 0;

	float x = 0.0f, y = 0.0f, z =0.0f;
	char c;

	char *grid_name = NULL, *ligand_name = NULL;
	float select_x = 0, select_y = 0, select_z = 0, select_radius = 0, threshold = 5.0;
	float min_ratio = 0.15;
	int seed = -1;

        max_ngen = 10000;


         while (1)
         {
           static struct option long_options[] =
             {
               {"dx",   required_argument, 0, 'd'},
               {"x",   required_argument, 0, 'x'},
               {"y",  required_argument, 0, 'y'},
               {"z",  required_argument, 0, 'z'},
               {"radius",  required_argument, 0, 'r'},
               {"threshold",    required_argument, 0, 't'},
               {"population",    required_argument, 0, 'p'},
               {"iterations",    required_argument, 0, 'i'},
               {"ligand",    required_argument, 0, 'l'},
               {"molarity",    required_argument, 0, 'm'},
               {"seed",    required_argument, 0, 's'},
               {"ratio",    required_argument, 0, 'c'},
               {0, 0, 0, 0}
             };
           int option_index = 0;

           c = getopt_long (argc, argv, "d:x:y:z:t:r:p:i:l:m:s:c:",
                            long_options, &option_index);
           if (c == -1)
             break;

           switch (c)
             {
             case 0:
               /* If this option set a flag, do nothing else now. */
               if (long_options[option_index].flag != 0)
                 break;
               printf ("option %s", long_options[option_index].name);
               if (optarg)
                 printf (" with arg %s", optarg);
               printf ("\n");
               break;
             case 'c':
                min_ratio = atof(optarg);
                break;
	     case 'd':
		grid_name = optarg;
		break;
	     case 'x':
                select_x = atof(optarg);
		break;
             case 'y':
                select_y = atof(optarg);
                break;
             case 'z':
                select_z = atof(optarg);
                break;
             case 'm':
                concentration = atof(optarg);
                break;
             case 'r':
                select_radius = atof(optarg);
                select_radius *= select_radius;
                break;
             case 't':
                threshold = atof(optarg);
                break;
             case 'p':
                max_ind = atoi(optarg);
		break;
             case 'i':
                max_ngen = atoi(optarg);
                break;
	     case 'l':
                ligand_name = optarg;
                break;
	     case 's':
                seed = atoi(optarg);
		break;
             case '?':
               /* getopt_long already printed an error message. */
               break;

             default:
               abort ();
             }
         }


        if( seed == -1)
                srand(time(NULL));
        else
                srand(seed);


        fprintf(stderr,"GAsol rev6 - A program to place water molecules using density grids\n");
        fprintf(stderr,"By Alvaro Cortes and Lucia Fusani. 2017 GlaxoSmithKline\n");
        fflush(stderr);

	
	if( grid_name == NULL)
	{
		fprintf(stderr,"No grid file provided\n");
		fflush(stderr);
                show_usage();
		exit(-1);
	}

        if( (dx_in = fopen(grid_name,"rb")) == NULL)
	{
		fprintf(stderr,"Cannot open the DX file %s\n",grid_name);
		fflush(stderr);
		show_usage();
		exit(-1);
	}

	
	if ((buffer_line = (char *) calloc(sizeof(char),1024)) == NULL)
	{
		fprintf(stderr,"Error allocating memory\n");
		fflush(stderr);
		fclose(dx_in);
		exit(-2);
	}
	

	if( (line_tokens = (char **) calloc(sizeof(char *), 20)) == NULL)
	{
                fprintf(stderr,"Error allocating memory\n");
                fflush(stderr);
                fclose(dx_in);
		free(buffer_line);
                exit(-2);
	}

	if( ligand_name != NULL)
	{
		if ( get_COM_ligand(ligand_name, &select_x, &select_y, &select_z) != 0)
		{
			fprintf(stderr,"Error reading PDB file\n");
			fflush(stderr);
			exit(-1);
		}
	}

	for( i = 0; i < 20; i ++)
	{
		line_tokens[i] = (char *) calloc(sizeof(char), 30);
	}

	/* This ad-hoc grid reader is not very general. It might fail with badly */
        /* formatted DX files. Need to check with different input types */
	fprintf(stderr,"Reading grid file %s ... ",grid_name);
	fflush(stderr);
        while( (fgets(buffer_line,1023,dx_in)) != NULL)
        {	
		current_token = 0;
		token = strtok(buffer_line, " ");
		strcpy(line_tokens[current_token],token);
   		while( token != NULL ) 
   		{
			current_token++;
/*      			printf( "%i %s\n", current_token-1,token );*/
      			token = strtok(NULL, " ");
			if ( token != NULL)
			strncpy(line_tokens[current_token],token,29);
   		}

	
		if( state_read != 0){

                        for(i = 0; i < current_token; i++)
                        {
                          data[current_point] = atof(line_tokens[i]);
                          ++current_point;
                        }

		if( current_point >= (nx*ny*nz))
			break;


		}else{

			if( current_token > 3 && strcmp(line_tokens[0],"object") == 0 && strcmp(line_tokens[1],"1") == 0 && strcmp(line_tokens[3],"gridpositions") == 0)
			{
				nx = atoi(line_tokens[5]);
				ny = atoi(line_tokens[6]);
				nz = atoi(line_tokens[7]);
				if( (data = (float *) calloc(sizeof(float),nx*ny*nz)) == NULL)
				{
					fprintf(stderr,"Error allocating memory\n");
					fflush(stderr);
					exit(-2);
				}
			}else if(current_token > 2 && strcmp(line_tokens[0],"origin") == 0){
				min[0] = atof(line_tokens[1]);
				min[1] = atof(line_tokens[2]);
				min[2] = atof(line_tokens[3]);
			}else if( strcmp(line_tokens[0],"delta") == 0){
				delta[current_delta] = atof(line_tokens[1+current_delta]);
				++current_delta;
		 	}else if(strcmp(line_tokens[0],"object") == 0 && strcmp(line_tokens[1],"3") == 0 && strcmp(line_tokens[3],"array") == 0){
				points = atoi(line_tokens[9]);
#ifdef DEBUG_OVERKILL
				fprintf(stderr,"%i %i\n",points,nx*ny*nz);
				fflush(stderr);
#endif
				state_read = 1;
			}

		}
        }

        fclose(dx_in);
	fprintf(stderr," done\n");
	fflush(stderr);

	if( current_point != nx*ny*nz)
	{
		fprintf(stderr,"Grid size and total number of points do not match\n");
		fflush(stderr);
		exit(-3);
	}
	fprintf(stderr,"Grid dimensions: %i,%i,%i - Total points: %i\n",nx,ny,nz,nx*ny*nz);
	fprintf(stderr,"Grid origin: %f,%f,%f - Grid spacing: %f,%f,%f\n",min[0],min[1],min[2],delta[0],delta[1],delta[2]);


	if( select_radius > 0)
	{
		fprintf(stderr,"Spatial filter is on. Solutions centered at %f,%f,%f with radius %f\n",select_x,select_y,select_z,sqrtf(select_radius));
		fflush(stderr);
	}

	/* Count the number of grid points with g(r) values above the threshold value */
	/* and optionally use the distance threshold if defined */
	l = 0;
        for( i = 0; i < nx; i++)
        {
                for(j = 0; j < ny; j++)
                {
                        for( k = 0; k < nz; k++)
                        {
				if( select_radius > 0)
				{
	 				dx = (min[0] + delta[0]*i) - select_x;
	 				dy = (min[1] + delta[1]*j) - select_y;
	 				dz = (min[2] + delta[2]*k) - select_z;
					dr = (dx*dx) + (dy*dy) + (dz*dz);
					if( data[l] >= threshold && dr <= select_radius)
					    ++n_points;
				}else{
					if( data[l] >= threshold)
						++n_points;
				}
				++l;
			}
		}
	}
	for( i = 0; i < nx*ny*nz; i++)
	{
		if( data[i] >= threshold)
		 ++n_points;
	}

	/* Select and store information for the grid points that meet the requirements */
	x_index = (int *) calloc(sizeof(int),n_points);
	y_index = (int *) calloc(sizeof(int),n_points);
	z_index = (int *) calloc(sizeof(int),n_points);
	g = (float *) calloc(sizeof(float),n_points);
	max_g = (float *) calloc(sizeof(float),n_points);
	
	current_point = 0;
	l = 0;
	all_g = 0;
	/* First pass - Select the points above the treshold and optionally use the distance
           threshold too */
	fprintf(stderr,"Selecting grid points and calculating populations ...");
	fflush(stderr);
        full_integ = 0.0;
	for( i = 0; i < nx; i++)
	{
		for(j = 0; j < ny; j++)
		{
			for( k = 0; k < nz; k++)
			{

                                if( select_radius > 0)
                                {
                                        dx = (min[0] + delta[0]*i) - select_x;
                                        dy = (min[1] + delta[1]*j) - select_y;
                                        dz = (min[2] + delta[2]*k) - select_z;
                                        dr = (dx*dx) + (dy*dy) + (dz*dz);
					if( dr <= select_radius)
						full_integ += data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
					if( data[l] >= threshold && dr <= select_radius) 
					{
						x_index[current_point] = i;
						y_index[current_point] = j;
						z_index[current_point] = k;
						g[current_point] = data[l] *delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
						all_g = all_g + data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
						++current_point;
					}
				}else{

					full_integ += data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;

                                        if( data[l] >= threshold)
                                        {
                                                x_index[current_point] = i;
                                                y_index[current_point] = j;
                                                z_index[current_point] = k;
                                                g[current_point] = data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
                                                all_g = all_g + data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
                                                ++current_point;
                                        }
				}
				++l;
			}
		}
	}
    	points_radii = (float *) calloc(sizeof(float),current_point);
	forbid_combination = (int **) calloc(sizeof(int *), current_point);
	for( i = 0; i < current_point; i++)
	{
		forbid_combination[i] = (int *) calloc(sizeof(int),current_point);
	}
	all_g = 0;
	/* Calculate sphere with unit population for all the candidate sites */
        #pragma omp parallel for \
                    private(ig,jg,kg,l,lg,current_distance,current_population,i,j,k,dx,dy,dz,dr) \
                    shared(x_index,y_index,z_index,data,delta,nx,ny,nz,max_g,points_radii ) \
		    reduction(+:all_g) \
                    schedule(static)
        for( l = 0; l < current_point; l++)
        {
                ig = x_index[l];
                jg = y_index[l];
                kg = z_index[l];
		lg = 0;
                current_distance = 0.9;
		current_population = data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
                while( current_distance < 2.6 && current_population < 1.0)
                {
			current_distance += 0.1;
			current_population = data[l]*delta[0]*delta[1]*delta[2] * concentration * 6.0221415E-4;
			lg = 0;
                        for( i = 0; i < nx; i++)
                        {
                                for(j = 0; j < ny; j++)
                                {
                                        for( k = 0; k < nz; k++)
                                        {

                                        dx = (min[0] + delta[0]*i) - (min[0] + delta[0]*ig);
                                        dy = (min[1] + delta[1]*j) - (min[1] + delta[1]*jg);
                                        dz = (min[2] + delta[2]*k) - (min[2] + delta[2]*kg);
                                        dr = (dx*dx) + (dy*dy) + (dz*dz);
                                        if( dr <= (current_distance*current_distance) && !(i == ig && j == jg && k == kg))
                                        {
						current_population += data[lg]*delta[0]*delta[1]*delta[2]*concentration*6.0221415E-4;

                                        }
					++lg;
					}
                                }
                        }
                }

		max_g[l] = current_population/current_distance;
		all_g += current_population/current_distance;
		/* Bugfix - 15/08/2017. If population is not close to 1 at 2.5 A */
		/* Very low density waters filter */
                if(max_g[l] < min_ratio)
			points_radii[l] = 999.9; /* Basically eliminates this molecule from all solutions. TODO: I should remove this site instead of making it invisible */
		else
			points_radii[l] = 1.5; /*current_distance/current_population; */

#ifdef DEBUG
	fprintf(stderr,"%i point %f density %f distance\n",l,current_population,current_distance);
	fflush(stderr);
#endif

        }

        for( l = 0; l < current_point; l++)
        {
            ig = x_index[l];
            jg = y_index[l];
            kg = z_index[l];
		
			for( lg = l+1; lg < current_point; lg++)
			{
                
                dx = (min[0] + delta[0]*x_index[lg]) - (min[0] + delta[0]*ig);
                dy = (min[1] + delta[1]*y_index[lg]) - (min[1] + delta[1]*jg);
                dz = (min[2] + delta[2]*z_index[lg]) - (min[2] + delta[2]*kg);
                dr = (dx*dx) + (dy*dy) + (dz*dz);
                /*if( dr <= (points_radii[l]*points_radii[l]) || dr <= (points_radii[lg]*points_radii[lg]))*/
                if( dr <= ((points_radii[l]+points_radii[lg])*(points_radii[l]+points_radii[lg])))
                {
					forbid_combination[l][lg] = 1;
					forbid_combination[lg][l] = 1;
				}
            }
        }

		

	fprintf(stderr," done\n");
        fprintf(stderr,"Number of predicted water molecules in the site: %f\n",full_integ);
	fflush(stderr);

	fprintf(stderr,"Chromosomes of %i genes\n",current_point);
	if( max_ind == 0)
   	  max_ind = current_point * 10;

        /* Initialize random population */
        population = (int **) calloc(sizeof(int *), max_ind);
        offspring = (int **) calloc(sizeof(int *), max_ind+1);
        fitness = (double *) calloc(sizeof(double), max_ind);
        new_fitness = (double *) calloc(sizeof(double), max_ind);
        best_ever = (int *) calloc(sizeof(int), current_point);

        for( nind = 0; nind < max_ind; nind++)
        {
                fitness[nind] = -999.9;
                population[nind] = (int *) calloc(sizeof(int),current_point+1);
                offspring[nind] = (int *) calloc(sizeof(int),current_point+1);
                for( i = 0; i < current_point; i++)
                {
                    population[nind][i] = (int) rand() % 2;
                }

        }

        /* Evolve population */
        best_fitness = -9999.9;

	fprintf(stderr,"Running Genetic algorithm for %i generations with a population of %i individuals\n",max_ngen,max_ind);
	fflush(stderr);

	/* Generations loop */
	/* TODO: check for convergence!! */
        for( ngen = 0; ngen < max_ngen; ngen++)
        {
		/* Only parallelized for individuals */
        #pragma omp parallel for \
                    private(nind,i,j,k,current_g,d1,d2,p2,w_sum,e,penalty,p6,d6, bad_w, dx, dy, dz) \
                    shared(g, forbid_combination,max_ind, fitness, lines, current_point, all_g, x_index, y_index, z_index, points_radii ) \
                    schedule(dynamic) 
                for( nind = 0 ; nind < max_ind; nind++) /* Evaluate population */
                {
                   if( fitness[nind] < -1)
                   {
                        d1 = d2 = 0.0;
			d2 = 1.0;
                        p2 = 0.0001;

			current_g = 0;
			bad_w = 0;
			for( j = 0; j < current_point; j++)
			{
				if( population[nind][j] == 1)
				{
				  current_g += max_g[j]; /* Add normalized population for on bits */

				  /* Check if two water molecules are overlapping in this solution */
				  for( k = j+1; k < current_point; k++)
				  {
					if( population[nind][k] == 1 && forbid_combination[j][k] == 1)
					{
						d2 = 0.0;
					 	bad_w += 1.0;
					}

				  }

				}
			}

			d1 = current_g / all_g;
                        w_sum = 0;
                        e = 1.0;
			p2 = (double) bad_w / (double) current_point;
                        penalty = 1.0;
                        for ( j = 0; j < 2; j++)
                        {
                                w_sum += w[j];
                        }
                        penalty *= p2;
                        e *= powf(d1,w[0]);
                        e *= powf(d2,w[1]);
                        e = powf(e,1.0/ (double) w_sum);
                        penalty = p2-0.0001;
                        e -= penalty;
                        fitness[nind] = e;
/*                        fitness[nind] = (d1*d2)-(p2-0.0001);*/

		   } /* Fitness */
                } /* Individuals loop */

		/* Update best solution */
                for( i = 0; i < max_ind; ++i)
                {
                        if ( fitness[i] > best_fitness)
                        {
#ifdef DEBUG
                                fprintf(stderr,"New Fitness: %f\n",fitness[i]);
#endif
                                best_fitness = fitness[i];
                                for( j = 0; j < current_point; j++)
                                {
                                 best_ever[j] = population[i][j];
#ifdef DEBUG
                                 fprintf(stderr,"%i, ", best_ever[j]);
#endif
                                }
#ifdef DEBUG
                                fprintf(stderr,"\n");
				fflush(stderr);
#endif
                        }
                }
                fprintf(stderr,"Iteration %i. Best fitness so far: %f\r",ngen,best_fitness);
                fflush(stderr);



                /* Selection round with 3 inidivuals at the same time */
                for( i = 0; i < max_ind; ++i)
                {
                        t1 = rand_interval(0,max_ind-1);
                        t2 = rand_interval(0,max_ind-1);
                        t3 = rand_interval(0,max_ind-1);
                        best_t = t1;
                        if( fitness[t1] > fitness[t2])
                            best_t = t1;
                        else
                            best_t = t2;
                        if( fitness[best_t] < fitness[t3])
                            best_t = t3;

                        for( j = 0; j < current_point; j++)
                        {
                                offspring[i][j] = population[best_t][j];
                                new_fitness[i] = fitness[t3];
                        }

                }

                /* Crossover */
                for( c1 = 1; c1 < max_ind; c1 = c1 + 2)
                {
                        c2 = c1 - 1;
                        if( r2() < 0.5)
                        {
                                new_fitness[c1] = -2;
                                new_fitness[c2] = -2;
                                cxpoint1 = rand_interval(0, current_point);
                                cxpoint2 = rand_interval(0, current_point);
                                if (cxpoint2 >= cxpoint1)
                                    cxpoint2 += 1;
                                else{
                                    j = cxpoint1;
                                    cxpoint1 = cxpoint2;
                                    cxpoint2 = j;
                                }

                                for (i = cxpoint1; i < cxpoint2; i++)
                                {
                                    j = offspring[c2][i];
                                    offspring[c2][i] = offspring[c1][i];
                                    offspring[c1][i] = j;
                                }
                        }

                }

                /* Mutation */
                for( c1 = 0; c1 < max_ind; c1++)
                {
                        if( r2() < 0.2)
                        {
                                new_fitness[c1] = -2;
                                for( i = 0; i < current_point; i++)
                                {
                                        if( r2() < 0.05)
                                        {
                                                if( offspring[c1][i] == 0)
                                                    offspring[c1][i] = 1;
                                                else
                                                     offspring[c1][i] = 0;
                                        }
                                }
                        }
                }

		/* "Puberty" */
		/* Kill the parents and promote the children */
                for( c1 = 0; c1 < max_ind; c1++)
                {
                        for( i = 0; i < current_point; i++)
                        {
                              	population[c1][i] = offspring[c1][i];
                                fitness[c1] = new_fitness[c1];
                        }
                }
	}

	/* Print best solution in the form of a PDB */
	fprintf(stderr,"\nBest solution fitness: %f\n",best_fitness);
	printf("MODEL 1\n");
	printf("REMARK Fitness: %f\n",best_fitness);
        printf("REMARK Generations: %i\n",max_ngen);
        printf("REMARK Individuals: %i\n",max_ind);
        printf("REMARK Centre x,y,z: %f,%f,%f\n",select_x,select_y,select_z);
	printf("REMARK Grid: %s\n",grid_name);
	printf("REMARK Radius: %f\n",sqrtf(select_radius));
	printf("REMARK Concentration: %f\n",concentration);
	printf("REMARK Ramdom_seed: %i\n",seed);
	i = 0;
        for( j = 0; j < current_point; j++)
        {
		if( best_ever[j] == 1)
		{
			++i;
			x = min[0] + delta[0]*(x_index[j]);
			y = min[1] + delta[1]*(y_index[j]);
			z = min[2] + delta[2]*(z_index[j]);
		printf("ATOM      1  %s   HOH A%4i    %8.3f%8.3f%8.3f  1.00 %2.7f\n", "O",j,x,y,z,max_g[j]);
		}
        }
	printf("TER\n");
	printf("ENDMDL\n");
	fflush(stdout);
	if( ceil(full_integ) < i)
	{
		fprintf(stderr,"Warning: current solution has %i more water molecule/s than the integral of g(r)\n",i- (int) ceil(full_integ));
		fprintf(stderr,"It may mean that some waters are partially occupied sites or false positives\n");
		fprintf(stderr,"Please Increase the ratio threshold with --ratio\n");
		fflush(stderr);
	}
        
	/* Save the whales and free the mallocs! */
	free(buffer_line);
	for(i = 0; i < 20; i++)
		free(line_tokens[i]);
	free(line_tokens);

	for( i = 0; i < current_point; i++)
		free(forbid_combination[i]); 
	
	free(forbid_combination);
	free(x_index); free(y_index); free(z_index);
	free(g); free(max_g);

	free(fitness); free(best_ever); free(new_fitness);
	for( i = 0; i < nind; i++)
	{
		free(population[i]);
		free(offspring[i]);
	}

	free(population); free(offspring); free(points_radii);
	free(data);
	exit(0);
}
