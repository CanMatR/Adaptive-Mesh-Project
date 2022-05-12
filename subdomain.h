#ifndef subdomainDEF
#define subdomainDEF

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <list>
#include <cmath>
#include "field.h"
#include "utility.h"
#include "mesh.h"

using namespace std;


class subdomain
{
public:
	subdomain *P;
	subdomain *C[2][2][2];
	cart_vector origin;
	int size;
	int dim;
	list<subdomain *>::iterator me;

	field *fld;
	subdomain(subdomain *Parent, cart_vector origin, int size, cart_vector grid_size, int num_fields,int num_fields_aux,int dim,int buffer,int num_comms,int dx)
	{
		int x, y, z;
		this->P = Parent;
		this->origin = origin;
		this->size = size;
		this->dim = dim;
//		dx = (size - 1) / (N-1);
		fld = new field(origin, grid_size, buffer, dx, num_fields,num_fields_aux, dim,num_comms);
		//Set all children to NULL, for 2D still set z to NULL, but don't use them in simulation
		for (z = 0; z < 2; z++)for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			C[x][y][z] = NULL;
	}
	subdomain(subdomain *Parent, cart_vector origin, int size, int dim)
	{
		int x, y, z;
		this->P = Parent;
		this->origin = origin;
		this->size = size;
		this->dim = dim;
		this->fld = NULL;
		//Set all children to NULL, for 2D still set z to NULL, but don't use them in simulation
		for (z = 0; z < 2; z++)for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			C[x][y][z] = NULL;

	}
	~subdomain()
	{
		int x, y, z;
		for (z = 0; z < 2; z++)for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			if (C[x][y][z] != NULL)
				delete C[x][y][z];
		delete fld;
	}
	subdomain & operator=(const subdomain &rhs)
	{
		this->origin = rhs.origin;
		this->size = rhs.size;
		this->dim = rhs.dim;
		(*this->fld) = *(rhs.fld);
		return *this;
	}
	//OUTPUT
	void Output_Restart(FILE *fp)
	{
		fprintf(fp,"%d %d %d	:	Subdomain Origin\n",this->origin.Vx,this->origin.Vy, this->origin.Vz);
		fprintf(fp, "%d	:	Subdomain Size\n",this->size);
		this->fld->Output_Restart(fp);
	}

	//Adaption Functions
	bool Check_Split()
	{
		if (fld->Get_Split_Status())
			return 1;
		else 
			return 0;
	}
	bool Has_Grand_Child()
	{
		int xc, yc, zc;
		int xcg, ycg, zcg;
		if (this->dim == 3)
		{
			for (zc = 0; zc < 2; zc++) for (yc = 0; yc < 2; yc++)for (xc = 0; xc < 2; xc++)//looping over children
			{
				if (C[xc][yc][zc] != NULL)
				{
					for (zcg = 0; zcg < 2; zcg++) for (ycg = 0; ycg < 2; ycg++)for (xcg = 0; xcg < 2; xcg++)//looping over children's children
					{
						if (C[xc][yc][zc]->C[xcg][ycg][zcg] != NULL)
						{
							return 1;
						}
					}
				}
			}
		}
		else if (this->dim == 2)
		{
			for (yc = 0; yc < 2; yc++)for (xc = 0; xc < 2; xc++)//looping over children
			{
				if (C[xc][yc][0] != NULL)
				{
					for (ycg = 0; ycg < 2; ycg++)for (xcg = 0; xcg < 2; xcg++)//looping over children's children
					{
						if (C[xc][yc][0]->C[xcg][ycg][0] != NULL)
						{
							return 1;
						}
					}
				}
			}
		}
		return 0;
	}
	subdomain *Find_Subdomain_Containing_Origin(cart_vector Orig)
	{
		int x, y, z;
		if (this->dim == 3)
		{
			if (Orig.Vx >= this->origin.Vx 
				&& Orig.Vx < this->origin.Vx + this->size - 1 
				&& Orig.Vy >= this->origin.Vy
				&& Orig.Vy < this->origin.Vy + this->size - 1
				&& Orig.Vz >= this->origin.Vz
				&& Orig.Vz < this->origin.Vz + this->size - 1) 
				//if the coordinate is in the subdomain search children
			{
				//if the subdomain has no children return itself
				if (C[0][0][0] == NULL)
					return this;
				for (z = 0; z<2; z++) for (y = 0; y<2; y++)for (x = 0; x<2; x++)
				{
					if (C[x][y][z]->Find_Subdomain_Containing_Origin(Orig) != NULL)
					{
						return C[x][y][z]->Find_Subdomain_Containing_Origin(Orig);
					}
				}
			}
			return NULL;
		}
		else if (this->dim == 2)
		{
			if (Orig.Vx >= this->origin.Vx
				&& Orig.Vx < this->origin.Vx + this->size - 1
				&& Orig.Vy >= this->origin.Vy
				&& Orig.Vy < this->origin.Vy + this->size - 1)
				//if the coordinate is in the subdomain search children
			{
				//if the subdomain has no children return itself
				if (C[0][0][0] == NULL)
					return this;
				for (y = 0; y<2; y++)for (x = 0; x<2; x++)
				{
					if (C[x][y][0]->Find_Subdomain_Containing_Origin(Orig) != NULL)
					{
						return C[x][y][0]->Find_Subdomain_Containing_Origin(Orig);
					}
				}
			}
			return NULL;
		}
		return NULL;
	}
	subdomain *Find_Subdomain(cart_vector Orig, int size)
	{
		int xc, yc, zc;
		subdomain *tmp;
		if (Orig.Vx >= this->origin.Vx 
			&& Orig.Vx < this->origin.Vx + this->size-1   
			&& Orig.Vy >= this->origin.Vy
			&& Orig.Vy < this->origin.Vy + this->size-1   
			&& Orig.Vz >= this->origin.Vz 
			&& Orig.Vz < this->origin.Vz + this->size-1)
			//if the coordinate lies in current subdomain
		{
			if (   Orig.Vx == this->origin.Vx 
				&& Orig.Vy == this->origin.Vy 
				&& Orig.Vz == this->origin.Vz 
				&& size == this->size)
				//does current domain coordinates and size match the subdomain i am looking for
			{
				return this;	//return me
			}
			else
			{
				//does current domain have children
				if (C[0][0][0] != NULL)
				{
					if (this->dim == 3)
					{
						for (zc = 0; zc < 2; zc++) for (yc = 0; yc < 2; yc++)for (xc = 0; xc < 2; xc++)
						{
							//set the domain i am looking for to the return of the call to this function
							tmp = C[xc][yc][zc]->Find_Subdomain(Orig, size);	
							//if its the domain i am looking for or i am successful in finding the domain
							if (tmp != NULL)
								return tmp;
						}
					}
					else if (this->dim == 2)
					{
						for (yc = 0; yc < 2; yc++)for (xc = 0; xc < 2; xc++)
						{
							//set the domain i am looking for to the return of the call to this function
							tmp = C[xc][yc][0]->Find_Subdomain(Orig,size);	
							//if its the domain i am looking for or i am successful in finding the domain
							if (tmp != NULL)	
								return tmp;
						}
					}
				}
				else
				{
					return NULL;
				}
			}
		}
		return NULL;
	}
	subdomain *Find_Subdomain_With_Coordinate(cart_vector coord)
	{
		int x, y, z;
		subdomain *tmp;
		int xend, yend, zend;
		if (this->origin.Vx + this->size == MESH_PARAMS.global_size.Vx && !MESH_PARAMS.periodic.Vx)
			xend = this->origin.Vx+this->size;
		else
			xend = this->origin.Vx+this->size-1;

		if (this->origin.Vy + this->size == MESH_PARAMS.global_size.Vy && !MESH_PARAMS.periodic.Vy)
			yend = this->origin.Vy + this->size;
		else
			yend = this->origin.Vy + this->size - 1;

		if (this->origin.Vz + this->size == MESH_PARAMS.global_size.Vz && !MESH_PARAMS.periodic.Vz)
			zend = this->origin.Vz + this->size;
		else
			zend = this->origin.Vz + this->size - 1;

		if (coord.Vx >= this->origin.Vx
			&& coord.Vx < xend
			&& coord.Vy >= this->origin.Vy
			&& coord.Vy < yend
			&& coord.Vz >= this->origin.Vz
			&& coord.Vz < zend )
			//if the coordinate lies in current domain
		{
			if (this->C[0][0][0] == NULL)//if i'm a leaf
			{
				return this;
			}
			else
			{
				if (this->dim == 3)
				{
					for (z = 0; z < 2; z++)	for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
					{
						tmp = this->C[x][y][z]->Find_Subdomain_With_Coordinate(coord);
						if (tmp != NULL)
						{
							return tmp;
						}
					}
				}
				else if (this->dim == 2)
				{
					
					for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
					{
						tmp = this->C[x][y][0]->Find_Subdomain_With_Coordinate(coord);
						if (tmp !=NULL)
						{
							return tmp;
						}
					}
				}
			}
		}
		return NULL;
	}
	bool split()
	{
		int x, y, z;
		cart_vector O_tmp;
		//check to make sure the domain calling this function does not have children
		if (C[0][0][0] != NULL)
			return 0;

		//creating new children for the domain
		if (this->dim == 3)
		{
			for (z = 0; z < 2; z++)for (y = 0; y < 2; y++) for (x = 0; x < 2; x++)
			{
				O_tmp.Vx = origin.Vx + x*((this->size-1) / 2);
				O_tmp.Vy = origin.Vy + y*((this->size-1) / 2);
				O_tmp.Vz = origin.Vz + z*((this->size-1) / 2);
				C[x][y][z] = new subdomain(this, O_tmp, (this->size-1)/2+1, this->dim);
				C[x][y][z]->fld = new field(O_tmp,this->fld);//creates the new field based off the parent's field
			}
		}
		else if (this->dim == 2)
		{
			for (y = 0; y < 2; y++) for (x = 0; x < 2; x++)
			{
				O_tmp.Vx = origin.Vx + x*((this->size - 1) / 2);
				O_tmp.Vy = origin.Vy + y*((this->size - 1) / 2);
				O_tmp.Vz = 0;

				C[x][y][0] = new subdomain(this, O_tmp, (this->size - 1) / 2 + 1, this->dim);
				C[x][y][0]->fld = new field(O_tmp, this->fld);//creates the new field based off the parent's field
			}
		}
		return 1;
	}
	bool Unsplit()
	{
		int x, y, z;
		//check to make sure the domain calling this function to unsplit actually has childre
		if (C[0][0][0] == NULL)
		{
			return 0;
		}
		if (this->dim == 3)
		{
			for (z = 0; z < 2; z++)for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			{
				this->fld->Merge(this->C[x][y][z]->fld);
				delete C[x][y][z];
				C[x][y][z] = NULL;
			}
		}
		else if (this->dim == 2)
		{
			for (y = 0; y < 2; y++)for (x = 0; x < 2; x++)
			{
				this->fld->Merge(this->C[x][y][0]->fld);
				delete C[x][y][0];
				C[x][y][0] = NULL;
			}
		}
		return 1;
	}
	bool edge_subdomain(cart_vector origin, int size)
	{
		if (origin.Vx == this->origin.Vx || origin.Vx + size == this->origin.Vx + this->size)
		{
			return 1;
		}
		if (origin.Vy == this->origin.Vy || origin.Vy + size == this->origin.Vy + this->size)
			return 1;
		if (this->dim == 3)
		{
			if (origin.Vz == this->origin.Vz || origin.Vz + size == this->origin.Vz + this->size)
				return 1;
		}
		return 0;
	}
	//Communication Functions
	bool In_Subdomain(int x, int y)
	{
		if (x >= this->origin.Vx && x <= this->origin.Vx + this->size - 1 && y >= this->origin.Vy && y <= this->origin.Vy + this->size - 1)
			return 1;
		else
			return 0;
	}
	bool In_Subdomain(int x, int y, int z)
	{
		
		if (x >= this->origin.Vx && x <= this->origin.Vx + this->size - 1 && y >= this->origin.Vy && y <= this->origin.Vy + this->size - 1 && z >= this->origin.Vz && z <= this->origin.Vz + this->size - 1)
			return 1;
		else
			return 0;
	}
	double Retrieve_Value(int xi, int yi, int id)
	{
		// recursively find the subdomain that has point x,y,z with the field id and ask the field for the value.
		subdomain *Cur;
		Cur = this;

		while (Cur->C[0][0][0] != NULL)
		{
			if (Cur->C[0][0][0]->In_Subdomain(xi, yi))
			{
				Cur = Cur->C[0][0][0];
			}
			else if(Cur->C[0][1][0]->In_Subdomain(xi, yi))
			{
				Cur = Cur->C[0][1][0];
			}
			else if(Cur->C[1][0][0]->In_Subdomain(xi, yi))
			{
				Cur = Cur->C[1][0][0];
			}
			else
			{
				Cur = Cur->C[1][1][0];
			}
			
		}
		return Cur->fld->Retrieve_Value(xi, yi, id);
	}
	double Retrieve_Value(int xi, int yi, int zi, int id)
	{
		// recursively find the subdomain that has point x,y,z with the field id and ask the field for the value.
		subdomain *Cur;
		Cur = this;

		while (Cur->C[0][0][0] != NULL)
		{
			if (Cur->C[0][0][0]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[0][0][0];
			}
			else if (Cur->C[0][0][1]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[0][0][1];
			}
			else if (Cur->C[0][1][0]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[0][1][0];
			}
			else if (Cur->C[0][1][1]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[0][1][1];
			}
			else if (Cur->C[1][0][0]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[1][0][0];
			}
			else if (Cur->C[1][1][0]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[1][1][0];
			}
			else if (Cur->C[1][0][1]->In_Subdomain(xi, yi, zi))
			{
				Cur = Cur->C[1][0][1];
			}
			else
			{
				Cur = Cur->C[1][1][1];
			}
		}
		return Cur->fld->Retrieve_Value(xi, yi, zi, id);
	}
	/**** OUTPUT FUNCTIONS ****/
};
#endif
