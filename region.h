/* Copyright (c) 2013 Government of Canada  */
#ifndef regionDEF
#define regionDEF
#include <vector>
#include "domain.h"
#include "utility.h"
#include "geometry.h"


//controls local organization of domains for the simulation
class region
{
public:
	int Regid;
	vector<domain *> Dom;
	int dim;
	int num_domain;
	cart_vector global_size;
	int num_comms;
	double **recv_val;
	int *recv_val_size;
	int **recv_coords_array;

	region()
	{
		this->num_domain = 0;
	}
	~region()
	{}
	void Adapt_Domains()
	{
		unsigned int i;
		for (i = 0; i < Dom.size(); i++)
		{
			Dom[i]->Adapt_Domain();
		}
	}
	void Unadapt_Domains()
	{
		unsigned int i;
		for (i = 0; i < Dom.size(); i++)
		{
			Dom[i]->Unadapt_Domain();
		}
	}
	void Adapt_Parallel(int dom_id, cart_vector coord, int dx)
	{
		int i;
		domain *c_dom=NULL;
		for (i = 0; i<num_domain; i++)//search for the domain at info[m]
		{
			if (Dom[i]->DP.dom_id == dom_id)
			{
				c_dom = Dom[i];
				break;
			}
		}
		if (c_dom != NULL)
		{
			//	c_dom->Split_To_Coordinate(coord, (c_dom->Get_Uni_Size() - 1)*dx * 2 + 1);
			if (dx == 0)//if the parallel neighbour is splitting, set the mesh to be one resolution level different.
				dx = 2;
			if (dx > 0)
				c_dom->Split_To_Coordinate(coord, (c_dom->Get_Uni_Size() - 1)*dx * 2 + 1);
			else
				c_dom->Split_To_Coordinate(coord, (int)pow(2, MESH_PARAMS.uniRes) + 1);//(c_dom->Get_Uni_Size() - 1)*1 + 1);
//				c_dom->Split_To_Coordinate(coord, (c_dom->Get_Uni_Size() - 1)*1 * 2 + 1);
			
		}
	}

	void Seed_Adapt_Domains(adapt_geometry **ad_geo, int num_geom)
	{
		unsigned int i;

		for (i = 0; i < Dom.size(); i++)
		{
			Dom[i]->Seed_Adapt_Domain(ad_geo, num_geom);
		}
	}
/*	void Init_Create_Domain(int dom_id, int numComms)
	{
		char inputline[BUFSIZ];
		char filename[BUFSIZ];

	}*/
	void Create_Domain(char *domainFilename, int time, int dom_id, int numComms)
	{
		char inputline[BUFSIZ];
		char filename[BUFSIZ];
		domain *t_Dom;
		dom_params DP;
		FILE *fp;
		sprintf(filename, "./restart_%d/region_%d/%s", time, this->Regid, domainFilename);
		fp = fopen(filename, "r");
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d", &(DP.dim));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d", &(DP.numComms));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d %d %d", &(DP.maxRes), &(DP.uniRes), &(DP.buffer));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d %d", &(DP.num_fields),&(DP.num_aux_fields));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d  %d  %d", &(DP.O.Vx), &(DP.O.Vy), &(DP.O.Vz));
		fclose(fp);
		
		DP.global_size.Vx = this->global_size.Vx;
		DP.global_size.Vy = this->global_size.Vy;
		DP.global_size.Vz = this->global_size.Vz;
		DP.dom_id = dom_id;
		DP.size = (int)pow(2.0, DP.uniRes + DP.maxRes - 1) + 1;
		t_Dom = new domain(DP);
		Dom.push_back(t_Dom);
		num_domain++;
		Load_Domain_From_File(filename, t_Dom);
	}
	void Load_Domain_From_File(char *domainFile, domain *t_Dom)
	{
		/*
		char inputline[BUFSIZ];
		dom_params DP;
		subdomain *tmp_sub;
		cart_vector origin;
		int N, size;
		FILE *fp;
		fp = fopen(domainFile, "r");

		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d", &(DP.dim));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d", &(DP.numComms));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d %d %d", &(DP.maxRes), &(DP.uniRes), &(DP.buffer));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d %d", &(DP.num_fields), &(DP.num_aux_fields));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d  %d  %d", &(DP.O.Vx), &(DP.O.Vy), &(DP.O.Vz));
		while (fgets(inputline, BUFSIZ, fp))
		{
			//subdomain origin
			sscanf(inputline, "%d %d %d", &(origin.Vx), &(origin.Vy), &(origin.Vz));

			//subdomain size
			fgets(inputline, BUFSIZ, fp);
			sscanf(inputline, "%d", &size);
			fgets(inputline, BUFSIZ, fp);
			sscanf(inputline, "%d", &N);//Field N
			//field header
			fgets(inputline, BUFSIZ, fp);//Field Data line
			tmp_sub = new subdomain(NULL, origin, size, N, DP.num_fields, DP.num_aux_fields, DP.dim, DP.buffer,DP.numComms);
			tmp_sub->fld->Load_Field_Data_From_File(fp);
			
			t_Dom->Insert_Subdomain(tmp_sub);
			delete tmp_sub;
//			break;
		}
		
		fclose(fp);	//close domain file
		*/
	}
	void Get_Field_List(list<field*> *fields)
	{
		int i;
		fields->clear();
		for (i = 0; i<this->num_domain; i++)
		{
			Dom[i]->Add_Fields_To_List(fields);
		}
	}
	// ***************** COMMUNICATION FUNCTIONS ************/
	void Get_Comm_Coords(int **recv_coords, int *recv_coords_size, int comm)
	{
		int *t_recv;
		int i, j, start;
		int *tmp_coord_buffer;
		(*recv_coords_size) = 0;
		//delete this comms array for values and coordinates
		if (recv_val[comm] != NULL)
			delete recv_val[comm];
		if (recv_coords_array[comm] != NULL)
			delete recv_coords_array[comm];
		for (i = 0; i < this->num_domain; i++)
		{
			(*recv_coords_size) += this->Dom[i]->Assemble_Buffers(comm);
		}
		
		t_recv = new int[(*recv_coords_size)];

		// **** copy coords into region request list
		start = 0;
		for (i = 0; i<num_domain; i++)
		{
			tmp_coord_buffer = Dom[i]->Retrieve_External_Coordinate_Buffers(comm);

			for (j = 0; j < Dom[i]->coord_buffer_size[1][comm]; j++)
			{
				t_recv[j + start] = tmp_coord_buffer[j];
			}
			start = start + j;
		}
		(*recv_coords) = t_recv;
		recv_val_size[comm] = (*recv_coords_size) / (this->dim + 1);
		recv_val[comm] = new double[recv_val_size[comm]];
		recv_coords_array[comm] = new int[(*recv_coords_size)];
		for (i = 0; i<(*recv_coords_size); i++)
		{
			recv_coords_array[comm][i] = t_recv[i];
		}
		
	}
	void Local_Communicate(int comm)
	{
		int i;
		for (i = 0; i < this->num_domain; i++)
		{
			Dom[i]->Fill_Buffers(comm);
		}
	}
	void Fill_Buffer(double *fill, int fill_size, int *info, int info_size, int comm)
	{
		int i;
		domain *c_dom=NULL;
		int m = 2, n = 0;
		int j = 2;
		int c = 0;
		if (this->dim == 2)
		{
			while (j<info_size)
			{
				if (j == m)
				{
					for (i = 0; i<num_domain; i++)//search for the domain at info[m]
					{
						if (Dom[i]->DP.dom_id == info[m])
						{
							c_dom = Dom[i];
							break;
						}
					}
					m = m + 2 + info[j + 1] * (this->dim + 1);
					j += 2;
				}
				else
				{
					fill[n] = c_dom->Retrieve_Value(info[j], info[j + 1], info[j + 2]);
					n++;
					j += (this->dim + 1);
				}
				c++;
			}
		}
		else if (this->dim == 3)
		{
			while (j<info_size)
			{
				if (j == m)
				{
					for (i = 0; i<num_domain; i++)//search for the domain at info[m]
					{
						if (Dom[i]->DP.dom_id == info[m])
						{
							c_dom = Dom[i];
							break;
						}
					}
					m = m + 2 + info[j + 1] * (this->dim + 1);//4 is the number of coordinate information variables, i.e. x,y,z, id
					j += 2;
				}
				else
				{
					fill[n] = c_dom->Retrieve_Value(info[j], info[j + 1], info[j + 2], info[j + 3]);
					n++;
					j += (this->dim + 1);
				}
				c++;
			}
		}
	}
	void Load_Values_Into_Domains(int comm)
	{
		int i, j, n, k;
		n = 0;
		for (i = 0; i<num_domain; i++)
		{
			k = 0;
			for (j = 0; j < Dom[i]->coord_buffer_size[1][comm] / (this->dim + 1); j++)//generalize for number of dimensions, currently 2 is d  : changed 2 to 3 to account for the id : Now changed to d+1
			{
				Dom[i]->coord_buffer_values[1][comm][k] = recv_val[comm][n];
				k++;
				n++;
			}
//			Dom[i]->load_requested_values(comm);
			//Dom[i]->outputDomain_ASCII(1);//remove
		}

	}

	// **************** OUTPUT FUNCTIONS **********/
	void Output_Load(int time, FILE *fp,int myid)
	{
		int i;
		for (i = 0; i < this->Dom.size(); i++)
		{
			Dom[i]->Output_Load(time,fp,myid);
		}
	}
	void Output_VTK(int time, int lowres)
	{
		unsigned int i;
		char directory[BUFSIZ];
		sprintf(directory, "./restart_%d/region_%d", time, this->Regid);
		for (i = 0; i < this->Dom.size(); i++)
		{
			Dom[i]->Output_VTK_3D(time, directory,lowres);
		}
		fflush(stdout);
	}
	void Output_Restart(int time)
	{
		unsigned int i;
		char directory[BUFSIZ];
		#ifdef _WIN32
				sprintf(directory, "restart_%d\\region_%d", time, this->Regid);
		#elif __linux
				sprintf(directory, "./restart_%d/region_%d", time, this->Regid);
		#elif __unix
				sprintf(directory, "./restart_%d/region_%d", time, this->Regid);
		#endif
		for (i = 0; i < this->Dom.size(); i++)
		{
			this->Dom[i]->Output_Restart(time, directory);
		}
	}
};
#endif