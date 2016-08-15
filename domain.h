#ifndef domainDEF
#define domainDEF

#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include <stack>
#include "utility.h"
#include "subdomain.h"
#include "geometry.h"
#include "mesh.h"
#include "mpi.h"
//#include "field.h"
//#include "utility.h"

using namespace std;

class dom_params
{
public:
	//geometry parameters
	int dim;					//dimensions
	int maxRes;					//number of subdivides
	int uniRes;					//size of uniform mesh
	int size;					//size domain with dx=1
	int buffer;					//width of buffer in fields
	cart_vector O;				//origin of domain
	cart_vector global_size;	//global_size of simulation

	//Simulation parameters
	int num_fields;				//number of values in fields
	int num_aux_fields;

	//Communication Parameters
	int dom_id;
	int numComms;


	dom_params()
	{}
	~dom_params()
	{}
};

class domain
{
	list<subdomain *>::iterator itRes;
	list<subdomain *> *subdomainRes;
	subdomain *root;	
	//communication stuff
	vector <int> *comm_fld;	//vector sorting field communications by comm
	//communication storage : 0 inner domain,1 extra domain

public:
	int **coord_buffer[2];
	double **coord_buffer_values[2];
	dom_params DP;

	int *coord_buffer_size[2];
	// House Keeping Functions
	domain(dom_params DP)
	{
		int i;
		this->DP = DP;
		int dx;
		subdomainRes = new list<subdomain *>[DP.maxRes];
		cart_vector grid_size;
		if (MESH_PARAMS.periodic.Vx){
			grid_size.Vx = (int)pow(2.0, DP.uniRes);
/*			if (this->DP.O.Vx + this->DP.size == MESH_PARAMS.global_size.Vx)//check if domain outer boundary touches global size
				grid_size.Vx = (int)pow(2.0, DP.uniRes);
			else
				grid_size.Vx = (int)pow(2.0, DP.uniRes) + 1;
				*/
		}
		else{
			if (this->DP.O.Vx + this->DP.size == MESH_PARAMS.global_size.Vx)//check if domain outer boundary touches global size
				grid_size.Vx = (int)pow(2.0, DP.uniRes) + 1;
			else
				grid_size.Vx = (int)pow(2.0, DP.uniRes);
		}

		if (MESH_PARAMS.periodic.Vy)
		{
			grid_size.Vy = (int)pow(2.0, DP.uniRes);
		}
		else
		{
			if (this->DP.O.Vy + this->DP.size == MESH_PARAMS.global_size.Vy)//check if domain outer boundary touches global size
				grid_size.Vy= (int)pow(2.0, DP.uniRes) + 1;
			else
				grid_size.Vy = (int)pow(2.0, DP.uniRes);

		}
		if (DP.dim == 3)
		{
			if (MESH_PARAMS.periodic.Vz)
			{
				grid_size.Vz = (int)pow(2.0, DP.uniRes);
			}
			else
			{
				if (this->DP.O.Vz + this->DP.size == MESH_PARAMS.global_size.Vz)//check if domain outer boundary touches global size
					grid_size.Vz = (int)pow(2.0, DP.uniRes) + 1;
				else
					grid_size.Vz = (int)pow(2.0, DP.uniRes);
			}
		}
		else
		{
			grid_size.Vz = 0;
		}
		dx = (int)pow(2.0, DP.maxRes-1);
		root = new subdomain(NULL, DP.O, DP.size, grid_size, DP.num_fields,DP.num_aux_fields, DP.dim, DP.buffer,DP.numComms,dx);
		subdomainRes[0].push_front(root);
		root->me = subdomainRes[0].begin();
		coord_buffer[0] = new int*[DP.numComms];
		coord_buffer_values[0] = new double*[DP.numComms];
		coord_buffer_size[0] = new int[DP.numComms];
		coord_buffer[1] = new int*[DP.numComms];
		coord_buffer_values[1] = new double*[DP.numComms];
		coord_buffer_size[1] = new int[DP.numComms];
		
		for (i = 0; i < DP.numComms; i++)
		{
			coord_buffer[0][i] = NULL;
			coord_buffer_values[0][i] = NULL;
			coord_buffer_size[0][i] = 0;
			coord_buffer[1][i] = NULL;
			coord_buffer_values[1][i] = NULL;
			coord_buffer_size[1][i] = 0;
		}
		
	}
	~domain()
	{
		int i;
		if (root != NULL)
			delete root;
		if (subdomainRes != NULL)
			for (i = 0; i < DP.maxRes; i++)
				subdomainRes[i].clear();
			delete subdomainRes;
	}
	dom_params Get_Domain_Parameters()
	{
		return this->DP;
	}
	void add_fields_to_list(list<field*> *fields)
	{
		int i;
		field *Field;
		for (i = 0; i<DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				Field = (*itRes)->fld;
				fields->push_back(Field);
				itRes++;
			}
		}
	}
	void Insert_Subdomain(subdomain *sub)
	{
		subdomain *it;
		//Split to size
		Split_To_Size(sub->origin, sub->size);
		//Find Subdomain position by navigating the tree
		it = this->root->Find_Subdomain_Containing_Origin(sub->origin);
		*it = *sub;
	}
	void Add_Fields_To_List(list<field*> *fields)
	{
		int i;
		list<subdomain *>::iterator itRes;
		field *Field;
		for (i = 0; i<DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				Field = (*itRes)->fld;
				fields->push_back(Field);
				itRes++;
			}
		}
	}
	//Adaption Functions
	bool Seed_Adapt_Domain(adapt_geometry **ad_geo, int num_geom)
	{
		int c = 1;

		while (c != 0)
		{
			c = 0;
			Seed_Calculate_Split_Bool(ad_geo,num_geom);
			c = Split_Domain();
		}
		return 1;
	}
	void Seed_Calculate_Split_Bool(adapt_geometry **ad_geo, int num_geom)
	{
		int i;
		subdomain *tmp;

		for (i = 0; i < DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				tmp = (*itRes);
				tmp->fld->Seed_Set_Split_Status(ad_geo, num_geom);
				itRes++;
			}
		}
	}
	bool Adapt_Domain()
	{
		int c = 1;//adaption counter, check adaption for as long as a change was made
		//split
		while (c != 0)
		{
			c = 0;
			Calculate_Split_Bool();
			c = Split_Domain();
		}
		return 1;
	}
	bool Unadapt_Domain()
	{
		Calculate_Split_Bool();

		Unsplit_Domain();
		return 1;
	}
	void Calculate_Split_Bool()
	{
		int i;
		subdomain *tmp;
		for (i = 0; i < DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				tmp = (*itRes);
				tmp->fld->Set_Split_Status();
				itRes++;
			}
		}
	}
	int Split_Domain()
	{
		int x, y, z, inc;
		cart_vector target;
		int k;
		int c = 0;//split counter
		subdomain *tmp;
		//split for check split
		for (k = DP.maxRes - 2; k >= 0; k--)
		{
			itRes = subdomainRes[k].begin();
			while (itRes != subdomainRes[k].end())
			{
				tmp = (*itRes);
				itRes++;
				if (tmp->Check_Split())
				{
					Split_And_Push(tmp);
					c++;
				}
			}
		}
		
		//split checking for neighbour variances
		
		for (k = 0; k < DP.maxRes; k++)
		{
			itRes = subdomainRes[k].begin();
			while (itRes != subdomainRes[k].end())
			{
				tmp = (*itRes);
				itRes++;
//				inc = (tmp->size - 1) / 2 + 1;
				inc = (tmp->size - 1) / 2;//changed due to grid_size change
				if (this->DP.dim == 3)
				{
					for (z = -1; z < 3; z++)for (y = -1; y < 3; y++)for (x = -1; x < 3; x++)
					{
						target.Vx = tmp->origin.Vx + x*inc;
						target.Vy = tmp->origin.Vy + y*inc;
						target.Vz = tmp->origin.Vz + z*inc;
						
						//						if (tmp->Check_Split())
						//						{
						//							c += Split_To_Size(target, tmp->size);
						//						}
						//						else
						c += Split_To_Size(target, (tmp->size - 1) * 2 + 1);
					}
				}
				else if( this->DP.dim == 2)
				{
					for (y = -1; y < 3; y++)for (x = -1; x < 3; x++)
					{
						target.Vx = tmp->origin.Vx + x*inc;
						target.Vy = tmp->origin.Vy + y*inc;
						target.Vz = 0;
						/*
						if (tmp->Check_Split())
						{
							c += Split_To_Size(target, tmp->size);
						}
						else
						*/
							c += Split_To_Size(target, (tmp->size - 1) * 2 + 1);
					}
				}
			}
		}
		
		
		return c;
	}
	int Unsplit_Domain()
	{
		int k;
		int c = 0,cs;
		int x, y, z;
		cart_vector O;
		
		subdomain *tmp;
		subdomain *tmpN;
		for (k = DP.maxRes - 1; k > 0; k--)
		{
			itRes = subdomainRes[k].begin();
			while (itRes != subdomainRes[k].end())
			{
				tmp = (*itRes)->P;
				if (tmp == NULL)
				{
					printf("hmmmm: Problem with Unsplit\n");
					while (1);
				}
				itRes++;
				cs = 0;
				if (DP.dim == 3)
				{
					//looping over children parent, if any children want to split 
					for (z = 0; z < 2; z++)for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
					{
						if (tmp->C[x][y][z] != NULL)
						{
							if (tmp->C[x][y][z]->C[0][0][0] != NULL)
								cs++;
							else
								if (tmp->C[x][y][z]->Check_Split())
									cs++;
						}
					}
					if (cs == 0)//looping over neighbours if primary indicates split
					{
						for (z = -1; z < 2; z++)for (y = -1; y < 2; y++)for (x = -1; x < 2; x++)
						{
							O.Vx = tmp->origin.Vx + x*(tmp->size - 1);
							O.Vy = tmp->origin.Vy + y*(tmp->size - 1);
							O.Vz = tmp->origin.Vz + z*(tmp->size - 1);
							tmpN = this->root->Find_Subdomain(O , tmp->size);
							if (tmpN != NULL)
							{
								if (tmpN->Has_Grand_Child())
									cs++;
							}
						}
					}
					if (cs == 0)// Criteria met, unsplit
					{
						//push parent back on list at level below
						subdomainRes[k - 1].push_front(tmp);		
						tmp->me = subdomainRes[k - 1].begin();
						//looping over children to remove them from list
						for (z = 0; z<2; z++)for (y = 0; y<2; y++)for (x = 0; x<2; x++)
						{
							if (tmp->C[x][y][z] != NULL)
							{
								if (itRes == tmp->C[x][y][z]->me) //is the current iterator set the child
									itRes++;
								subdomainRes[k].erase(tmp->C[x][y][z]->me);	//remove my iterator from current list
							}
						}
						tmp->Unsplit();	//unsplit domain
						c++;
					}
				}
				else if (DP.dim == 2)
				{
					for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
					{
						if (tmp->C[x][y][0] != NULL)
						{
							if (tmp->C[x][y][0]->C[0][0][0] != NULL)
								cs++;
							else
								if (tmp->C[x][y][0]->Check_Split())
									cs++;
						}
					}
					if (cs == 0)//looping over neighbours if primary indicates split
					{
						for (y = -1; y < 2; y++)for (x = -1; x < 2; x++)
						{
							O.Vx = tmp->origin.Vx + x*(tmp->size - 1);
							O.Vy = tmp->origin.Vy + y*(tmp->size - 1);
							O.Vz = 0;
							tmpN = this->root->Find_Subdomain(O, tmp->size);
							if (tmpN != NULL)
							{
								if (tmpN->Has_Grand_Child())
									cs++;
							}
						}
					}
					if (cs == 0)// Criteria met, unsplit
					{
						//push parent back on list at level below
						subdomainRes[k - 1].push_front(tmp);
						tmp->me = subdomainRes[k - 1].begin();
						//looping over children to remove them from list
						for (y = 0; y<2; y++)for (x = 0; x<2; x++)
						{
							if (tmp->C[x][y][0] != NULL)
							{
								//is the current iterator set to the child
								if (itRes == tmp->C[x][y][0]->me) 
									itRes++;
								//remove my iterator from current list
								subdomainRes[k].erase(tmp->C[x][y][0]->me);
							}
						}
						tmp->Unsplit();	//unsplit domain
						c++;
					}
				}
			}
		}
		return c;
	}
	bool Split_And_Push(subdomain *toSplit)
	{
		int clevel;
		int x, y, z;
		toSplit->split();
		clevel = (int)(log((double)(DP.size - 1) / (double)(toSplit->size - 1)) / log(2.0) + 0.1);
		subdomainRes[clevel].erase(toSplit->me);
		if (this->DP.dim == 3)
		{
			for (z = 0; z < 2; z++)for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			{
				subdomainRes[clevel + 1].push_front(toSplit->C[x][y][z]);
				toSplit->C[x][y][z]->me = subdomainRes[clevel + 1].begin();
			}
		}
		else if (this->DP.dim == 2)
		{
			for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			{
				subdomainRes[clevel + 1].push_front(toSplit->C[x][y][0]);
				toSplit->C[x][y][0]->me = subdomainRes[clevel + 1].begin();
			}
		}
		
		return 1;
	}
	int Split_To_Size(cart_vector O, int size)
	{
		int splits = 0;
		subdomain *sd;
		sd = this->root->Find_Subdomain_Containing_Origin(O);
		if (sd != NULL)
		{
			if (sd->size > size)
			{
				Split_And_Push(sd);
				splits++;
				splits += Split_To_Size(O, size);
			}
			return splits;
		}
		else
		{
			return 0;
		}

		/*
		if (sd != NULL)
		{
			if (sd->size > size)
			{
				Split_And_Push(sd);

				return 1;
			}
			else
				return 0;
		}
		else
		{
			return 0;
		}
		*/
	}
	int Get_Uni_Size()
	{
		return (int)pow(2.0, DP.uniRes) + 1;
	}
	void Split_To_Coordinate(cart_vector coord, int size)
	{
		subdomain *it;

		//Find Subdomain position by navigating the tree
		it = root->Find_Subdomain_With_Coordinate(coord);

		if (it != NULL)
		{
			while (!(it->size <= size))
			{
				Split_And_Push(it);
				it = root->Find_Subdomain_With_Coordinate(coord);
			}
		}
		else
		{
			printf("Error, not in domain received : NULL\n");
			int xt, yt, zt;
			xt = this->DP.O.Vx + this->DP.size;
			yt = this->DP.O.Vy + this->DP.size;
			zt = this->DP.O.Vz + this->DP.size;
			printf("Domain coords and size : %d - %d , %d - %d , %d - %d : %d\n", this->DP.O.Vx, xt, this->DP.O.Vy, yt, this->DP.O.Vz, zt, this->DP.size);
		}
//		Split_Domain();//split domain to ensure level variance??
	}
	//Communication Functions
	void Fill_Buffers(int comm)
	{
		int x, y, z;
		int k;
		unsigned int i;
		int num_int,num_ext;
		field *fld;
		//k = (comm - 1) % DP.maxRes;
		if (this->DP.dim == 3)
		{
			//first calculate the needed values in each field and push them into coord_buffer_values
			for (k = DP.maxRes - 1; k >= 0; k--)
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					fld = (*itRes)->fld;
					for (i = 0; i < fld->coord_request[comm].size(); i++)
					{
						x = fld->coord_request[comm][i].coord.Vx;
						y = fld->coord_request[comm][i].coord.Vy;
						z = fld->coord_request[comm][i].coord.Vz;
						coord_buffer_values[0][comm][fld->coord_request[comm][i].index] = fld->Retrieve_Value(x, y,z, fld->coord_request[comm][i].id);
					}
					itRes++;
				}
			}
			//Second iterate over fields to load values from coord_buffer_values into the fields
			num_int = 0;
			num_ext = 0;
			k = (comm - 1) % DP.maxRes;
//			for (k = DP.maxRes - 1; k >= 0; k--)
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					fld = (*itRes)->fld;
					fld->Load_Coordinate_Requests(this->coord_buffer_values, &num_int, &num_ext, this->DP.O, this->DP.size, comm,this->DP.global_size);
					itRes++;
				}
			}
		}
		else if (this->DP.dim == 2)
		{
			//first calculate the needed values in each field and push them into coord_buffer_values
			for (k = DP.maxRes - 1; k >= 0; k--)
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					fld = (*itRes)->fld;
					for (i = 0; i < fld->coord_request[comm].size(); i++)
					{
						x = fld->coord_request[comm][i].coord.Vx;
						y = fld->coord_request[comm][i].coord.Vy;
						coord_buffer_values[0][comm][fld->coord_request[comm][i].index] = fld->Retrieve_Value(x, y, fld->coord_request[comm][i].id);
					}
					itRes++;
				}
			}
			//Second iterate over fields to load values from coord_buffer_values into the fields
			num_int = 0;
			num_ext = 0;
			k = (comm - 1) % DP.maxRes;
//			for (k = DP.maxRes - 1; k >= 0; k--)
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					fld = (*itRes)->fld;
					fld->Load_Coordinate_Requests(this->coord_buffer_values, &num_int, &num_ext, this->DP.O, this->DP.size, comm,this->DP.global_size);
					itRes++;
				}
			}
			
		}
	}
	int Assemble_Buffers(int comm)
	{
		int k,i;
		coordinate coord;
		subdomain *sd;
		int num_vals[2];//0 inner domain, 1 extra domain, 
		subdomain *tmp;
		if (comm == 0)
			Assemble_Adaption_Comm();
		else
		{
			//reset communication arrays
			coord_buffer_size[0][comm] = 0;
			coord_buffer_size[1][comm] = 0;
			if (coord_buffer_values[0][comm] != NULL)
				delete coord_buffer_values[0][comm];
			if (coord_buffer[0][comm] != NULL)
				delete coord_buffer[0][comm];
			if (coord_buffer_values[1][comm] != NULL)
				delete coord_buffer_values[1][comm];
			if (coord_buffer[1][comm] != NULL)
				delete coord_buffer[1][comm];
			//Count communication sizes
			num_vals[0] = 0;
			num_vals[1] = 0;
			
			int num_int=0, num_ext=0;
			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					(*itRes)->fld->coord_request[comm].clear();
					itRes++;
				}
			}
			k = (comm - 1)%DP.maxRes;
			itRes = subdomainRes[k].begin();
			while (itRes != subdomainRes[k].end())
			{
				tmp = (*itRes);
				tmp->fld->Count_Extra_And_Local_Domain(&(num_int), &(num_ext), this->DP.O, this->DP.size, comm, MESH_PARAMS.global_size);
				num_vals[0] = num_int;
				num_vals[1] = num_ext;
				itRes++;
			}
			coord_buffer_size[0][comm] = num_vals[0] * (this->DP.dim + 1);
			coord_buffer[0][comm] = new int[coord_buffer_size[0][comm]];
			coord_buffer_values[0][comm] = new double[num_vals[0]];
			coord_buffer_size[1][comm] = num_vals[1] * (this->DP.dim + 1);
			coord_buffer[1][comm] = new int[coord_buffer_size[1][comm]];
			coord_buffer_values[1][comm] = new double[num_vals[1]];
			num_int = 0;
			num_ext = 0;
			k = (comm - 1) % DP.maxRes;
			itRes = subdomainRes[k].begin();
			
			while (itRes != subdomainRes[k].end())
			{
				tmp = (*itRes);
				tmp->fld->Get_Coordinate_Requests(coord_buffer, &num_int, &num_ext, this->DP.O, this->DP.size, comm, MESH_PARAMS.global_size);
				itRes++;
			}
			//Partition local requests into fields
			if (this->DP.dim == 3)
			{
				for (i = 0; i < coord_buffer_size[0][comm]; i += this->DP.dim + 1)
				{
					coord.coord.Vx = coord_buffer[0][comm][i];
					coord.coord.Vy = coord_buffer[0][comm][i + 1];
					coord.coord.Vz = coord_buffer[0][comm][i + 2];
					coord.id = coord_buffer[0][comm][i + 3];
					coord.index = i / (this->DP.dim + 1);//index into the values array
					sd = this->root->Find_Subdomain_With_Coordinate(coord.coord);
					if (sd == NULL)
					{
						printf("Could not find coordinate: %d %d\n", coord.coord.Vx, coord.coord.Vy);
						printf("domain.h : Assemble_buffers\n");
						fflush(stdout);
						exit(0);

					}
					else{
						sd->fld->coord_request[comm].push_back(coord);
					}
				}
			}
			else if (this->DP.dim == 2)
			{
				for (i = 0; i < coord_buffer_size[0][comm]; i += this->DP.dim + 1)
				{
					coord.coord.Vx = coord_buffer[0][comm][i];
					coord.coord.Vy = coord_buffer[0][comm][i + 1];
					coord.coord.Vz = 0;
					coord.id = coord_buffer[0][comm][i + 2];
					coord.index = i / (this->DP.dim + 1);//index into the values array

					sd = this->root->Find_Subdomain_With_Coordinate(coord.coord);
					if (sd == NULL)
					{
						printf("Could not find coordinate: %d %d\n",coord.coord.Vx,coord.coord.Vy);
						printf("domain.h : Assemble_buffers\n");
						fflush(stdout);
						exit(0);
					}
					else{
						sd->fld->coord_request[comm].push_back(coord);
					}
				}
			}	
		}
		return coord_buffer_size[1][comm];
	}
	void Assemble_Adaption_Comm()
	{
		int i, j, k,q;
		subdomain *tmp;
		//Count the number of coordinates
		int num_subD = 0;
		int x, y, z;
		int c = 0;
		int N, B,dx;


		//set coord buffer size
		coord_buffer_size[1][0] = 0;

		//delete coord_buffer_values
		if (coord_buffer_values[1][0] != NULL)
			delete coord_buffer_values[1][0];
		//delete coord_buffer
		if (coord_buffer[1][0] != NULL)
			delete coord_buffer[1][0];
		if (this->DP.dim == 3)
		{
			
			int xend, yend,zend;
			if (DP.O.Vx + DP.size == MESH_PARAMS.global_size.Vx)
				xend = DP.O.Vx + DP.size;
			else
				xend = DP.O.Vx + DP.size - 1;
			if (DP.O.Vy + DP.size == MESH_PARAMS.global_size.Vy)
				yend = DP.O.Vy + DP.size;
			else
				yend = DP.O.Vy + DP.size - 1;
			if (DP.O.Vz + DP.size == MESH_PARAMS.global_size.Vz)
				zend = DP.O.Vz + DP.size;
			else
				zend = DP.O.Vz + DP.size - 1;
			
			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					tmp = (*itRes);
					
					if (tmp->edge_subdomain(this->DP.O, this->DP.size))//check if edge subdomain
					{
					
						B = tmp->fld->B;
						dx = tmp->fld->dx;
						
						for (q = -B; q < tmp->fld->grid_size.Vz + B; q++)for (j = -B; j < tmp->fld->grid_size.Vy + B; j++)for (i = -B; i < tmp->fld->grid_size.Vx + B; i++)
						{
						
							if (MESH_PARAMS.periodic.Vx)
								x = (i*dx + tmp->origin.Vx + this->DP.global_size.Vx) % this->DP.global_size.Vx;//periodic
							else
								x = (i*dx + tmp->origin.Vx);
							if (MESH_PARAMS.periodic.Vy)
								y = (j*dx + tmp->origin.Vy + this->DP.global_size.Vy) % this->DP.global_size.Vy;//periodic
							else
								y = (j*dx + tmp->origin.Vy);
							if (MESH_PARAMS.periodic.Vz)
								z = (q*dx + tmp->origin.Vz + this->DP.global_size.Vz) % this->DP.global_size.Vz;//periodic
							else
								z = (q*dx + tmp->origin.Vz);

							if (!(x >= this->DP.O.Vx && x < xend && y >= this->DP.O.Vy && y < yend && z >= this->DP.O.Vz && z < zend))
							{
								if (x >= 0 && x < MESH_PARAMS.global_size.Vx && y >= 0 && y < MESH_PARAMS.global_size.Vy && z >= 0 && z < MESH_PARAMS.global_size.Vz)
								{
									num_subD++;
									if (tmp->Check_Split())
									{
										num_subD++;
									}
								}
							}
							
						}
						
					}
					itRes++;
				}
			}
			
			coord_buffer_size[1][0] = num_subD*(this->DP.dim + 1);
			//allocate coord_buffer
			coord_buffer[1][0] = new int[coord_buffer_size[1][0]];
			//allocate coord_buffer_values, not needed for comm 0 since we send no values and piggy back on the coordinate sync
			coord_buffer_values[1][0] = new double[coord_buffer_size[1][0] / (this->DP.dim + 1)];
			//load coordinates into coord_buffer
			c = 0;
			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					tmp = (*itRes);
					if (tmp->edge_subdomain(this->DP.O, this->DP.size))//check if edge subdomain
					{
						B = tmp->fld->B;
						dx = tmp->fld->dx;

						for (q = -B; q < tmp->fld->grid_size.Vz + B; q++)for (j = -B; j < tmp->fld->grid_size.Vy + B; j++)for (i = -B; i < tmp->fld->grid_size.Vx + B; i++)
						{
							if (MESH_PARAMS.periodic.Vx)
								x = (i*dx + tmp->origin.Vx + this->DP.global_size.Vx) % this->DP.global_size.Vx;//periodic
							else
								x = (i*dx + tmp->origin.Vx);
							if (MESH_PARAMS.periodic.Vy)
								y = (j*dx + tmp->origin.Vy + this->DP.global_size.Vy) % this->DP.global_size.Vy;//periodic
							else
								y = (j*dx + tmp->origin.Vy);
							if (MESH_PARAMS.periodic.Vz)
								z = (q*dx + tmp->origin.Vz + this->DP.global_size.Vz) % this->DP.global_size.Vz;//periodic
							else
								z = (q*dx + tmp->origin.Vz);

							if (!(x >= this->DP.O.Vx && x < xend && y >= this->DP.O.Vy && y < yend && z >= this->DP.O.Vz && z < zend))
							{
								if (x >= 0 && x < MESH_PARAMS.global_size.Vx && y >= 0 && y < MESH_PARAMS.global_size.Vy && z >= 0 && z < MESH_PARAMS.global_size.Vz)
								{
									coord_buffer[1][0][c] = x;
									coord_buffer[1][0][c + 1] = y;
									coord_buffer[1][0][c + 2] = z;
									coord_buffer[1][0][c + 3] = tmp->fld->dx;
									c += 4;
									if (tmp->Check_Split())
									{
										coord_buffer[1][0][c] = x;
										coord_buffer[1][0][c + 1] = y;
										coord_buffer[1][0][c + 2] = z;
										coord_buffer[1][0][c + 3] = 0;
										c += 4;
									}
								}
							}
						}
					}
					itRes++;
				}
			}
/*			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					tmp = (*itRes);
					if (tmp->edge_subdomain(this->DP.O, this->DP.size))//check if edge subdomain
					{
						N = tmp->fld->N;
						B = tmp->fld->B;
						dx = tmp->fld->dx;
						for (q = -B; q < N + B; q++)for (j = -B; j < N + B; j++)for (i = -B; i < N + B; i++)
						{
							x = (i*dx + tmp->origin.Vx + this->DP.global_size.Vx) % this->DP.global_size.Vx;
							y = (j*dx + tmp->origin.Vy + this->DP.global_size.Vy) % this->DP.global_size.Vy;
							z = (q*dx + tmp->origin.Vz + this->DP.global_size.Vz) % this->DP.global_size.Vz;
							if (!(x >= this->DP.O.Vx && x <= this->DP.O.Vx + this->DP.size && y >= this->DP.O.Vy && y <= this->DP.O.Vy && z >= this->DP.O.Vz && z <= this->DP.O.Vz + this->DP.size))
							{
								num_subD++;
							}
						}
					}
					itRes++;
				}
			}
			coord_buffer_size[1][0] = num_subD*(this->DP.dim + 1);
			//allocate coord_buffer
			coord_buffer[1][0] = new int[coord_buffer_size[1][0]];
			//allocate coord_buffer_values, not needed for comm 0 since we send no values and piggy back on the coordinate sync
			coord_buffer_values[1][0] = new double[coord_buffer_size[1][0] / (this->DP.dim + 1)];
			//load coordinates into coord_buffer
			c = 0;
			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					tmp = (*itRes);
					if (tmp->edge_subdomain(this->DP.O, this->DP.size))//check if edge subdomain
					{
						N = tmp->fld->N;
						B = tmp->fld->B;
						dx = tmp->fld->dx;
						for (q = -B; q < N + B; q++)for (j = -B; j < N + B; j++)for (i = -B; i < N + B; i++)
						{
							x = (i*dx + tmp->origin.Vx + this->DP.global_size.Vx) % this->DP.global_size.Vx;
							y = (j*dx + tmp->origin.Vy + this->DP.global_size.Vy) % this->DP.global_size.Vy;
							z = (q*dx + tmp->origin.Vz + this->DP.global_size.Vz) % this->DP.global_size.Vz;
							if (!(x >= this->DP.O.Vx && x <= this->DP.O.Vx + this->DP.size && y >= this->DP.O.Vy && y <= this->DP.O.Vy && z >= this->DP.O.Vz && z <= this->DP.O.Vz + this->DP.size))
							{
								coord_buffer[1][0][c] = x;
								coord_buffer[1][0][c + 1] = y;
								coord_buffer[1][0][c + 2] = z;
								coord_buffer[1][0][c + 3] = tmp->fld->dx;
								c += 4;
							}
						}
					}
					itRes++;
				}
			}
			*/
		}
		else if (this->DP.dim == 2)
		{
			int xend, yend;
			if (DP.O.Vx + DP.size == MESH_PARAMS.global_size.Vx)
				xend = DP.O.Vx + DP.size;
			else
				xend = DP.O.Vx + DP.size - 1;
			if (DP.O.Vy + DP.size == MESH_PARAMS.global_size.Vy)
				yend = DP.O.Vy + DP.size;
			else
				yend = DP.O.Vy + DP.size - 1;
			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{
					tmp = (*itRes);
					if (tmp->edge_subdomain(this->DP.O, this->DP.size))//check if edge subdomain
					{
						B = tmp->fld->B;
						dx = tmp->fld->dx;
						
						for (j = -B; j < tmp->fld->grid_size.Vy + B; j++)for (i = -B; i < tmp->fld->grid_size.Vx + B; i++)
						{
							if (MESH_PARAMS.periodic.Vx)
								x = (i*dx + tmp->origin.Vx + this->DP.global_size.Vx) % this->DP.global_size.Vx;//periodic
							else
								x = (i*dx + tmp->origin.Vx);
							if (MESH_PARAMS.periodic.Vy)
								y = (j*dx + tmp->origin.Vy + this->DP.global_size.Vy) % this->DP.global_size.Vy;//periodic
							else
								y = (j*dx + tmp->origin.Vy);
							if (!(x >= this->DP.O.Vx && x < xend && y >= this->DP.O.Vy && y < yend))
							{
								if (x >= 0 && x < MESH_PARAMS.global_size.Vx && y >= 0 && y < MESH_PARAMS.global_size.Vy)
								{
									num_subD++;
									if (tmp->Check_Split())
									{
										num_subD++;
									}
								}
							}
						}
					}
					itRes++;
				}
			}
			
			coord_buffer_size[1][0] = num_subD*(this->DP.dim + 1);
			//allocate coord_buffer
			coord_buffer[1][0] = new int[coord_buffer_size[1][0]];
			//allocate coord_buffer_values, not needed for comm 0 since we send no values and piggy back on the coordinate sync
			coord_buffer_values[1][0] = new double[coord_buffer_size[1][0] / (this->DP.dim + 1)];
			//load coordinates into coord_buffer
			c = 0;
			for (k = DP.maxRes - 1; k >= 0; k--)//check for edge sub domains
			{
				itRes = subdomainRes[k].begin();
				while (itRes != subdomainRes[k].end())
				{

					tmp = (*itRes);
					if (tmp->edge_subdomain(this->DP.O, this->DP.size))//check if edge subdomain
					{
						B = tmp->fld->B;
						dx = tmp->fld->dx;

						for (j = -B; j < tmp->fld->grid_size.Vy + B; j++)for (i = -B; i < tmp->fld->grid_size.Vx + B; i++)
						{
							if (MESH_PARAMS.periodic.Vx)
								x = (i*dx + tmp->origin.Vx + this->DP.global_size.Vx) % this->DP.global_size.Vx;//periodic
							else
								x = (i*dx + tmp->origin.Vx);
							if (MESH_PARAMS.periodic.Vy)
								y = (j*dx + tmp->origin.Vy + this->DP.global_size.Vy) % this->DP.global_size.Vy;//periodic
							else
								y = (j*dx + tmp->origin.Vy);
							if (!(x >= this->DP.O.Vx && x < xend && y >= this->DP.O.Vy && y < yend))
							{
								if (x >= 0 && x < MESH_PARAMS.global_size.Vx && y >= 0 && y <MESH_PARAMS.global_size.Vy)
								{
									coord_buffer[1][0][c] = x;
									coord_buffer[1][0][c + 1] = y;
									coord_buffer[1][0][c + 2] = tmp->fld->dx;
									c += 3;
									if (tmp->Check_Split())
									{
										coord_buffer[1][0][c] = x;
										coord_buffer[1][0][c + 1] = y;
										coord_buffer[1][0][c + 2] = 0;
										c += 3;
									}
								}
							}
						}
					}
					itRes++;
				}
			}
		}
	}
	int * Retrieve_External_Coordinate_Buffers(int comm)
	{
		return coord_buffer[1][comm];
	}
	double Retrieve_Value(int xi, int yi, int id)
	{
//		double value;
		cart_vector coord;
		subdomain *sd;
		coord.Vx = xi;
		coord.Vy = yi;
//		value = root->Retrieve_Value(xi, yi, id);
		sd = root->Find_Subdomain_With_Coordinate(coord);
		return sd->fld->Retrieve_Value(xi, yi, id);

//		return value;
	}
	double Retrieve_Value(int xi, int yi, int zi, int id)
	{
		cart_vector coord;
		subdomain *sd;
		coord.Vx = xi;
		coord.Vy = yi;
		coord.Vz = zi;
		//		value = root->Retrieve_Value(xi, yi, id);
		sd = root->Find_Subdomain_With_Coordinate(coord);
		return sd->fld->Retrieve_Value(xi, yi, zi,id);

		/*
		double value;

		value = root->Retrieve_Value(xi, yi, zi, id);

		return value;
		*/
	}
	//Output Functions
	void Output_Restart(int time, char *directory)
	{
		int i;
		char filename[BUFSIZ];
		FILE *fp;
		sprintf(filename, "%s/domain_%d_%d_%d.dat", directory, this->DP.O.Vx,this->DP.O.Vy,this->DP.O.Vz);
		fp = fopen(filename, "w");
		fprintf(fp, "%d	:	Number of dimensions (2 or 3)\n",this->DP.dim);
		fprintf(fp, "%d	:	Number of Communications\n", this->DP.numComms);
		fprintf(fp, "%d %d %d	:	maxRes, uniRes, buffer size\n", this->DP.maxRes, this->DP.uniRes, this->DP.buffer);
		fprintf(fp, "%d %d	:	Number of simulation fields, auxilliary fields\n",this->DP.num_fields,this->DP.num_aux_fields);
		fprintf(fp, "%d %d %d	:	Domain Origin\n",this->DP.O.Vx,this->DP.O.Vy,this->DP.O.Vz);
		for (i = 0; i < this->DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				(*itRes)->Output_Restart(fp);
				itRes++;
			}
		}

		fclose(fp);
	}
	void Output_VTK_3D(int time, char *outdirectory, int lowres)
	{
		int V;
		char filename[BUFSIZ];
		FILE *fp;
		int i;
		int np = 0;
		int nc = 0;
		int np2 = 0;
		
		if (lowres > MESH_PARAMS.uniRes)
		{
			lowres = MESH_PARAMS.uniRes;
		}
//		if (lowres == 0)
//			sprintf(filename, "%s/domain_%d_%d_%d.vtk", outdirectory, this->DP.O.Vx, this->DP.O.Vy, this->DP.O.Vz);
//		else
		sprintf(filename, "%s/domain_%d_res_%d_%d_%d.vtk", outdirectory,lowres, this->DP.O.Vx, this->DP.O.Vy, this->DP.O.Vz);
		fp = fopen(filename, "w");
		//write tvk header
		fprintf(fp, "# vtk DataFile Version 3.0\n");
		fprintf(fp, "vtk output\n");
		fprintf(fp, "ASCII\n");
		fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

		//build output lists
	
		//points
		//cells
		//cell type
		//point_data NUM
		//SCALARS Phase double
		//LOOPUP_TABLE default
		//...data...
		for (i = 0; i<DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				np += (*itRes)->fld->vtk_Get_Num_Points(lowres);
				nc += (*itRes)->fld->vtk_Get_Num_Cells(lowres);
				itRes++;
			}
		}
		//points
		fprintf(fp, "POINTS %d FLOAT\n", np);
		for (i = 0; i<DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				(*itRes)->fld->vtk_Output_Points(fp,lowres);
				itRes++;
			}
		}
		//cells
		np2 = 0;
		fprintf(fp, "CELLS %d %d\n", nc, nc * 9);
		for (i = 0; i<DP.maxRes; i++)
		{
			itRes = subdomainRes[i].begin();
			while (itRes != subdomainRes[i].end())
			{
				np2+=(*itRes)->fld->vtk_Output_Cells(np2,fp,lowres);
				itRes++;
			}
		}
		fprintf(fp, "CELL_TYPES %d\n", nc);
		for (i = 0; i<nc; i++)
			fprintf(fp, "11\n");
		//point_data NUM
		//SCALARS Phase double
		//LOOPUP_TABLE default
		//...data...
		fprintf(fp, "POINT_DATA %d\n", np);

		for (V = 0; V < DP.num_fields; V++)
		{
			fprintf(fp, "SCALARS Field_%d double\n", V);
			fprintf(fp, "LOOKUP_TABLE default\n");
			for (i = 0; i < DP.maxRes; i++)
			{
				itRes = subdomainRes[i].begin();
				while (itRes != subdomainRes[i].end())
				{
					(*itRes)->fld->vtk_Output_Point_Data(V,fp,lowres);
					itRes++;
				}
			}
		}
		fclose(fp);
	}

};

#endif
