#ifndef fieldDEF
#define fieldDEF

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <list>
#include <stack>
#include <cmath>
#include <vector>
#include "utility.h"
#include "geometry.h"
#include "mesh.h"

using namespace std;

class communications
{
public:
	vector <int> field;
	communications()
	{}
	~communications()
	{}
	void push_comm(int field)
	{
		this->field.push_back(field);
	}
};
vector <communications> COMMS;
class coordinate
{
public:
	cart_vector coord;
	int id;
	int index;
};

class field
{
public:
	bool split_flag;
	int dx;		//mesh spacing of field
	int B;		//buffer size
	int N;		//size of field without buffer
	cart_vector grid_size;
	int num_comms;
	int nf;		//number of fields
	int nf_a;	//number of auxilliary fields
	int dim;
	double **data;
	double **data_aux;
	cart_vector origin;
	vector <coordinate> *coord_request;
	double *coord_request2;

	field(cart_vector origin, cart_vector grid_size, int buffer, int dx, int num_fields, int num_fields_aux, int dim,int num_comms)
	{
		this->origin = origin;
//		this->N = N;
		this->grid_size = grid_size;
		this->B = buffer;
		this->dx = dx;
		this->nf = num_fields;
		this->nf_a = num_fields_aux;
		this->dim = dim;
		this->num_comms = num_comms;
		data = new double*[nf];
		data_aux = new double*[nf_a];
		coord_request = NULL;
		coord_request = new vector<coordinate>[num_comms];
		coord_request2 = NULL;
		int i,ind,x,y,z;
		for (i = 0; i < nf; i++)
		{
			if (dim == 3)
			{
				data[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)*(grid_size.Vz + 2 * B)];
				for (z = 0; z < this->grid_size.Vz + 2 * this->B; z++)for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vz + 2 * this->B; x++)
				{
					ind = index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					data[i][ind] = 0.0;
				}
			}
			else if (dim == 2)
			{
				data[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)];
				for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					ind = index(x, y, this->grid_size.Vx + 2 * this->B);
					data[i][ind] = 0.0;
				}
			}
			else
			{
				printf("Simulation only presently supports 2 and 3 dimensions.\n  Field declaration of %d dimensions\n",this->dim);
			}
		}
		for (i = 0; i < nf_a; i++)
		{
			if (dim == 3)
				data_aux[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)*(grid_size.Vz + 2 * B)];
			else if (dim == 2)
				data_aux[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)];
			else
			{
				printf("Simulation only presently supports 2 and 3 dimensions.\n  Field declaration of %d dimensions\n", this->dim);
			}
		}
	}
	field(cart_vector origin, int N, int buffer, int dx, int num_fields, int num_fields_aux, int dim, int num_comms)
	{
		this->origin = origin;
				this->N = N;
		//this->grid_size = grid_size;
		this->B = buffer;
		this->dx = dx;
		this->nf = num_fields;
		this->nf_a = num_fields_aux;
		this->dim = dim;
		this->num_comms = num_comms;
		data = new double*[nf];
		data_aux = new double*[nf_a];
		coord_request = NULL;
		coord_request = new vector<coordinate>[num_comms];
		coord_request2 = NULL;
		int i, ind, x, y, z;
		for (i = 0; i < nf; i++)
		{
			if (dim == 3)
			{
				data[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)*(grid_size.Vz + 2 * B)];
				for (z = 0; z < this->grid_size.Vz + 2 * this->B; z++)for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vz + 2 * this->B; x++)
				{
					ind = index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					data[i][ind] = 0.0;
				}
			}
			else if (dim == 2)
			{
				data[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)];
				for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					ind = index(x, y, this->grid_size.Vx + 2 * this->B);
					data[i][ind] = 0.0;
				}
			}
			else
			{
				printf("Simulation only presently supports 2 and 3 dimensions.\n  Field declaration of %d dimensions\n", this->dim);
			}
		}
		for (i = 0; i < nf_a; i++)
		{
			if (dim == 3)
				data_aux[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)*(grid_size.Vz + 2 * B)];
			else if (dim == 2)
				data_aux[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)];
			else
			{
				printf("Simulation only presently supports 2 and 3 dimensions.\n  Field declaration of %d dimensions\n", this->dim);
			}
		}
	}
	field(cart_vector origin, field *fld)//creates a field from the parent fld assuming a split
	{
		int x, y, z;
		int i,ind;
		this->origin = origin;
		this->N = fld->N;
		int Dom_size = (int)pow(2, MESH_PARAMS.maxRes + MESH_PARAMS.uniRes - 1) + 1;	
		this->B = fld->B;
		this->dx = fld->dx/2;

		//set grid size here
		if (MESH_PARAMS.periodic.Vx)
			this->grid_size.Vx = fld->grid_size.Vx;
		else
		{
			if (fld->origin.Vx + (fld->grid_size.Vx-1)*fld->dx+1 == MESH_PARAMS.global_size.Vx)
			{
				if (fld->origin.Vx == origin.Vx)
					this->grid_size.Vx = fld->grid_size.Vx-1;
				else
					this->grid_size.Vx = fld->grid_size.Vx;
			}
			else
				this->grid_size.Vx = fld->grid_size.Vx;
		}
		if (MESH_PARAMS.periodic.Vy)
			this->grid_size.Vy = fld->grid_size.Vy;
		else
		{
			if (fld->origin.Vy + (fld->grid_size.Vy-1)*fld->dx+1 == MESH_PARAMS.global_size.Vy)
			{
				if (fld->origin.Vy == origin.Vy)
					this->grid_size.Vy = fld->grid_size.Vy - 1;
				else
					this->grid_size.Vy = fld->grid_size.Vy;
			}
			else
				this->grid_size.Vy = fld->grid_size.Vy;
		}
		if (MESH_PARAMS.DIMENSION == 3)
		{
			if (MESH_PARAMS.periodic.Vz)
				this->grid_size.Vz = fld->grid_size.Vz;
			else
			{
				if (fld->origin.Vz + (fld->grid_size.Vz -1)*fld->dx+1 == MESH_PARAMS.global_size.Vz)
				{
					if (fld->origin.Vz == origin.Vz)
						this->grid_size.Vz = fld->grid_size.Vz - 1;
					else
						this->grid_size.Vz = fld->grid_size.Vz;
				}
				else
					this->grid_size.Vz = fld->grid_size.Vz;
			}
		}
		//End set grid size
		this->nf = fld->nf;
		this->nf_a = fld->nf_a;
		this->dim = fld->dim;
		this->num_comms = fld->num_comms;
		data = new double*[nf];
		data_aux = new double*[nf_a];
		coord_request = NULL;
		coord_request = new vector<coordinate>[this->num_comms];
		coord_request2 = NULL;
		//Allocate memory for the fields and auxilliary fields
		for (i = 0; i < nf; i++)
		{
			if (dim == 3)
			  {
				  
			    data[i] = new (nothrow) double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)*(grid_size.Vz + 2 * B)];
			    if(data[i] == NULL)
			      {
				printf("Failed to allocate memory, exiting\n");
				exit(0);
			      }
			  }
			else if (dim == 2)
				data[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)];
			else
			{
				printf("Simulation only presently supports 2 and 3 dimensions.\n  Field declaration of %d dimensions\n", this->dim);
			}
		}
		for (i = 0; i < nf_a; i++)
		{
			if (dim == 3)
				data_aux[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)*(grid_size.Vz + 2 * B)];
			else if (dim == 2)
				data_aux[i] = new double[(grid_size.Vx + 2 * B)*(grid_size.Vy + 2 * B)];
			else
			{
				printf("Simulation only presently supports 2 and 3 dimensions.\n  Field declaration of %d dimensions\n", this->dim);
			}
		}
		//Copy relevant data from parent to child for data fields
		for (i = 0; i < nf; i++)
		{
			if (dim == 3)
			{
				for (z = 0; z < this->grid_size.Vz+this->B*2; z++)for (y = 0; y < this->grid_size.Vy+this->B*2; y++)for (x = 0; x < this->grid_size.Vx+this->B*2; x++)
				{
					ind = index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					data[i][ind] = fld->Retrieve_Value(this->origin.Vx + (x-this->B)*this->dx, this->origin.Vy + (y-this->B)*this->dx, this->origin.Vz + (z-this->B)*this->dx, i);
				}
			}
			else if (dim == 2)
			{
				for (y = 0; y < this->grid_size.Vy + this->B * 2; y++)for (x = 0; x < this->grid_size.Vx + this->B * 2; x++)
				{
					ind = index(x, y, this->grid_size.Vx + 2 * this->B);
					data[i][ind] = fld->Retrieve_Value(this->origin.Vx + (x-this->B)*this->dx, this->origin.Vy + (y-this->B)*this->dx, i);
				}
			}
		}
	}
	~field()
	{
		int i;
		if (data != NULL)
		{
			for (i = 0; i < this->nf; i++)
			{
				if (data[i] != NULL)
					delete data[i];
			}
			delete data;
		}
		if (data_aux != NULL)
		{
			for (i = 0; i < this->nf_a; i++)
			{
				if (data_aux[i] != NULL)
					delete data_aux[i];
			}
			delete data_aux;
		}
		if (coord_request != NULL)
		{
			delete[] coord_request;
		}
	}
	field & operator=(const field &rhs)
	{
		int i, x, y, z,ind;
		this->origin = rhs.origin;
		this->dx = rhs.dx;
		this->B = rhs.B;
		this->N = rhs.N;
		this->nf = rhs.nf;
		this->nf_a = rhs.nf_a;
		this->dim = rhs.dim;
		this->num_comms = rhs.num_comms;
		
		if (this->dim == 3)
		{
			for (i = 0; i < this->nf; i++)
			{
				for (z = this->B; z < this->N + this->B; z++)for (y = this->B; y < this->N + this->B; y++)for (x = this->B; x < this->N + this->B; x++)
				{
					ind = index(x, y, z, this->N + this->B * 2, this->N + this->B * 2);
					this->data[i][ind] = rhs.data[i][ind];
				}
			}
		}
		else if (this->dim == 2)
		{
			for (i = 0; i < this->nf; i++)
			{
				for (y = this->B; y < this->N + this->B; y++)for (x = this->B; x < this->N + this->B; x++)
				{
					ind = index(x, y, this->N+this->B*2);
					this->data[i][ind] = rhs.data[i][ind];
				}
			}
		}
		return *this;
	}
	void Load_Field_Data_From_File(FILE *fp)
	{
		int i, x, y, z;
		int ind;
		double tmpV;
		char inputline[BUFSIZ];
		if (this->dim == 3)
		{
			for (i = 0; i < this->nf; i++)
			{
				for (z = this->B; z < this->N +this->B; z++)for (y = this->B; y < this->N + this->B; y++)for (x = this->B; x < this->N + this->B; x++)
				{
					ind = index(x, y, z, this->N + 2 * this->B, this->N + 2 * this->B);
					fgets(inputline, BUFSIZ, fp);
					sscanf(inputline, "%lf", &tmpV);
					this->data[i][ind] = tmpV;
				}
			}
		}
		else if (this->dim == 2)
		{
			for (i = 0; i < this->nf; i++)
			{
				for (y = this->B; y < this->N + this->B; y++)for (x = this->B; x < this->N + this->B; x++)
				{
					ind = index(x, y, this->N+2*this->B);
					fgets(inputline, BUFSIZ, fp);
					sscanf(inputline, "%lf", &tmpV);
					this->data[i][ind] = tmpV;
				}
			}
		}
	}

	void Merge(field *fld)//Merges a given field into it's parent field.
	{
		int x, y, i,ind;
		int xs, ys;
		int xe, ye;
		int dx_rat;
		dx_rat = this->dx / fld->dx;
		xs = (fld->origin.Vx - this->origin.Vx) / this->dx;
		xe = ceil(((float)(fld->origin.Vx+fld->grid_size.Vx*fld->dx-this->origin.Vx)) / (float) this->dx);
		ys = (fld->origin.Vy - this->origin.Vy) / this->dx;
		ye = ceil(((float)(fld->origin.Vy + fld->grid_size.Vy*fld->dx - this->origin.Vy)) / (float)this->dx);
		if (this->dim == 3)
		{
			int z,zs, ze;
			zs = (fld->origin.Vz - this->origin.Vz) / this->dx;
			ze = ceil(((fld->origin.Vz + fld->grid_size.Vz*fld->dx - this->origin.Vz)) / (float)this->dx);
			for (i = 0; i < fld->nf; i++)
			{
				for (z = zs; z < ze; z++)for (y = ys; y < ye; y++)for (x = xs; x < xe; x++)
				{
					ind = index(x + this->B, y + this->B, z + this->B, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					data[i][ind] = fld->Retrieve_Value(this->origin.Vx + x*this->dx, this->origin.Vy + y*this->dx, this->origin.Vz + z*this->dx, i);
				}
			}
		}
		else if (this->dim == 2)
		{

			for (i = 0; i < fld->nf; i++)
			{
				for (y = ys; y < ye; y++)for (x = xs; x < xe; x++)
				{
					ind = index(x + this->B, y + this->B, this->grid_size.Vx+ 2 * this->B);
					data[i][ind] = fld->Retrieve_Value(this->origin.Vx + x*this->dx, this->origin.Vy + y*this->dx, i);
				}
			}
		}
	}
	bool Get_Split_Status()
	{
		return this->split_flag;
	}
	void Seed_Set_Split_Status(adapt_geometry **ad_geo, int num_geom)
	{
		int x, y, z,i;
		int counter;
		this->split_flag = 0;
		double Nx = (this->grid_size.Vx+1)*this->dx;
		double Ny = (this->grid_size.Vy+1)*this->dx;
		double Nz = (this->grid_size.Vz+1)*this->dx;
		if (this->dim == 3)
		{
			counter = 0;

			for (i = 0; i < num_geom; i++)
			{
				if (ad_geo[i]->check_contained_in(this->origin.Vx, this->origin.Vy, this->origin.Vz, Nx,Ny,Nz))
					counter++;
				for (z = 0; z <= this->grid_size.Vz; z++)for (y = 0; y <= this->grid_size.Vy; y++)for (x = 0; x <= this->grid_size.Vx; x++)
				{
					if (ad_geo[i]->check_inside(x*this->dx + this->origin.Vx, y*this->dx + this->origin.Vy, z*this->dx + this->origin.Vz))
						counter++;
				}
			}
			if (counter > 0)
				this->split_flag = 1;

			
		}
		else if (dim == 2)
		{
			counter = 0;
			for (i = 0; i < num_geom; i++)
			{
				if (ad_geo[i]->check_contained_in(this->origin.Vx, this->origin.Vy, Nx,Ny))
					counter++;

				for (y = 0; y < this->grid_size.Vy; y++)for (x = 0; x < this->grid_size.Vx; x++)
				{
					if (ad_geo[i]->check_inside(x*this->dx + this->origin.Vx, y*this->dx + this->origin.Vy))
						counter++;					
				}
			}
			if (counter > 0)
				this->split_flag = 1;
		}
	}
	void Set_Split_Status()
	{
		int xi, yi, zi, V;
		double counter;
		this->split_flag = 0;
		int xp, xm, yp, ym, zp, zm, ind;
		int start, endx, endy, endz;
		start = this->B;
		endx = this->grid_size.Vx + this->B - 1;
		endy = this->grid_size.Vy + this->B - 1;
		endz = this->grid_size.Vz + this->B - 1;
		double n2 = (start - endx)*(start - endy);
		double n3 = (endx - start)*(endy-start)*(endz-start);
		double prefactor[5];
		prefactor[0] = 1.0;
		prefactor[1] = 0.0;
		prefactor[2] = 0.0;
		prefactor[3] = 0.0;
		prefactor[4] = 0.0;

		if (this->dim == 3)
		{
			counter = 0.0;
			for (zi = start; zi < endz; zi++)for (yi = start; yi < endy; yi++)for (xi = start; xi < endx; xi++)
			{
				ind = index(xi, yi,zi, grid_size.Vx + 2 * this->B, grid_size.Vy + 2 * this->B);
				xp = index(xi+1, yi, zi, grid_size.Vx + 2 * this->B, grid_size.Vx + 2 * this->B);
				yp = index(xi, yi+1, zi, grid_size.Vx + 2 * this->B, grid_size.Vx + 2 * this->B);
				zp = index(xi, yi, zi+1, grid_size.Vx + 2 * this->B, grid_size.Vx + 2 * this->B);

				for (V = 0; V < this->nf; V++)
				{
					counter += prefactor[V] * (fabs(data[V][xp] - data[V][ind])
						+ fabs(data[V][yp] - data[V][ind])
						+ fabs(data[V][zp] - data[V][ind])
						) / (this->dx*n3);

				}
			}
			if (counter > 1E-4)
			{
				this->split_flag = 1;
			}

		}
		else if (this->dim == 2)
		{
			counter = 0.0;
			for (yi = start; yi < endy; yi++)for (xi = start; xi < endx; xi++)
			{
				ind = index(xi, yi, grid_size.Vx + 2 * this->B);
				xp = index(xi + 1, yi, grid_size.Vx + 2 * this->B);
				//				xm = index(xi - 1, yi, N+2*this->B);
				yp = index(xi, yi + 1, grid_size.Vx + 2 * this->B);
				//				ym = index(xi, yi - 1, N+2*this->B);
				// loop through all fields
				//				for (V = 0; V < 1; V++)
				for (V = 0; V < this->nf; V++)
				{
				  counter += prefactor[V] * (fabs(data[V][xp]-data[V][ind])+fabs(data[V][yp]-data[V][ind]))/(this->dx*n2);
				  //					counter += prefactor[V] * fabs(data[0][ind] - 1)*(fabs(data[V][xp] - data[V][ind])
				  //	+ fabs(data[V][yp] - data[V][ind])
				  //	) / (this->dx*n2);
				}
			}
			if (counter > 1E-4)
				this->split_flag = 1;
		}
	}
	/*
	inline double Retrieve_Value(int Gx, int Gy, int Gz, int V)
	{
		double Val;
		int xi, yi, zi;
		int rx, ry, rz;
		double xp, yp, zp;
		rx = Gx - this->origin.Vx;
		ry = Gy - this->origin.Vy;
		rz = Gz - this->origin.Vz;
		xi = rx / dx + this->B;
		yi = ry / dx + this->B;
		zi = rz / dx + this->B;
		int ind[8];
		//int Na = this->N+2*this->B;//array size with buffer
		int Nx = this->grid_size.Vx + 2 * this->B;
		int Ny = this->grid_size.Vy + 2 * this->B;

		ind[0] = index(xi, yi, zi, Nx, Ny);

		if (rx%dx == 0 && ry%dx == 0 && rz%dx == 0)
		{
			return data[V][ind[0]];
		}
		else
		{
			ind[1] = index(xi + 1, yi, zi,  Nx, Ny);
			ind[2] = index(xi, yi + 1, zi,  Nx, Ny);
			ind[3] = index(xi, yi, zi + 1,  Nx, Ny);
			ind[4] = index(xi, yi + 1, zi + 1,  Nx, Ny);
			ind[5] = index(xi + 1, yi, zi + 1,  Nx, Ny);
			ind[6] = index(xi + 1, yi + 1, zi,  Nx, Ny);
			ind[7] = index(xi + 1, yi + 1, zi + 1, Nx, Ny);
			xp = ((float)(rx%dx)) / (float)dx;
			yp = ((float)(ry%dx)) / (float)dx;
			zp = ((float)(rz%dx)) / (float)dx;

			//		linear interpolation function for 3D
			Val =
				(
				(
				(
				-data[V][ind[0]] + data[V][ind[1]] + data[V][ind[2]] - data[V][ind[4]] + data[V][ind[3]] - data[V][ind[5]] + data[V][ind[7]] - data[V][ind[6]]
				)*zp
				+ data[V][ind[0]] - data[V][ind[2]] - data[V][ind[1]] + data[V][ind[6]]
				)*yp
				+
				(
				-data[V][ind[1]] + data[V][ind[0]] + data[V][ind[5]] - data[V][ind[3]]
				)*zp
				+
				data[V][ind[1]] - data[V][ind[0]]
				)*xp
				+
				(
				(
				data[V][ind[0]] - data[V][ind[2]] - data[V][ind[3]] + data[V][ind[4]]
				)*zp
				- data[V][ind[0]] + data[V][ind[2]]
				)*yp
				+
				(
				data[V][ind[3]] - data[V][ind[0]]
				)*zp
				+ data[V][ind[0]];
			return Val;
		}
	}*/
	inline double Retrieve_Value(int Gx, int Gy, int Gz, int V)
	{
		double Val;
		int xi, yi,zi;
		int rx, ry,rz;
		double xp, yp,zp;
		rx = Gx - this->origin.Vx + this->B*dx;
		ry = Gy - this->origin.Vy + this->B*dx;
		rz = Gz - this->origin.Vz + this->B*dx;
		xi = rx / dx;
		yi = ry / dx;
		zi = rz / dx;

		int ind[8];
		int Nx = this->grid_size.Vx + 2 * this->B;
		int Ny = this->grid_size.Vy + 2 * this->B;
		/*
		//for model with single id field
		if(V==3)//special case for id field, remove for general modelling
		  {
		    int xx,yy,zz,nn;
		    for(zz=0;zz<2;zz++)for(yy=0;yy<2;yy++)for(xx=0;xx<2;xx++)
		    {
		      nn = index(xi+xx,yi+yy,zi+zz,Nx,Ny);
		      if(data[3][nn] != 0)
			 return data[V][nn];
		      else
			return 0.0;
		      //                        ind[0] = index(xi, yi, zi, Nx, Ny);
		      //                        return data[V][ind[0]];
		    }
		  }
		  */
		if (rx%dx == 0 && ry%dx == 0 && rz%dx ==0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			return data[V][ind[0]];
		}
		else if (rx%dx == 0 && ry%dx == 0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi, yi, zi+1, Nx, Ny);
			return 0.5*(data[V][ind[1]] + data[V][ind[0]]);
		}
		else if (rx%dx == 0 && rz%dx == 0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi, yi+1, zi, Nx, Ny);
			return 0.5*(data[V][ind[1]] + data[V][ind[0]]);
		}
		else if (ry%dx == 0 && rz%dx == 0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi+1, yi, zi, Nx, Ny);
			return 0.5*(data[V][ind[1]] + data[V][ind[0]]);
		}
		else if (rx%dx == 0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi, yi+1, zi, Nx, Ny);
			ind[2] = index(xi, yi, zi+1, Nx, Ny);
			ind[3] = index(xi, yi+1, zi+1, Nx, Ny);
			return 0.25*(data[V][ind[0]] + data[V][ind[1]] + +data[V][ind[2]] + data[V][ind[3]]);
		}
		else if (ry%dx == 0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi+1, yi, zi, Nx, Ny);
			ind[2] = index(xi, yi, zi + 1, Nx, Ny);
			ind[3] = index(xi+1, yi, zi + 1, Nx, Ny);
			return 0.25*(data[V][ind[0]] + data[V][ind[1]] + +data[V][ind[2]] + data[V][ind[3]]);
		}
		else if (rz%dx == 0)
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi, yi + 1, zi, Nx, Ny);
			ind[2] = index(xi+1, yi, zi, Nx, Ny);
			ind[3] = index(xi+1, yi + 1, zi, Nx, Ny);
			return 0.25*(data[V][ind[0]] + data[V][ind[1]] + +data[V][ind[2]] + data[V][ind[3]]);
		}
		else
		{
			ind[0] = index(xi, yi, zi, Nx, Ny);
			ind[1] = index(xi, yi + 1, zi, Nx, Ny);
			ind[2] = index(xi + 1, yi, zi, Nx, Ny);
			ind[3] = index(xi + 1, yi + 1, zi, Nx, Ny);
			ind[4] = index(xi, yi, zi+1, Nx, Ny);
			ind[5] = index(xi, yi + 1, zi+1, Nx, Ny);
			ind[6] = index(xi + 1, yi, zi+1, Nx, Ny);
			ind[7] = index(xi + 1, yi + 1, zi+1, Nx, Ny);
			return 0.125*(data[V][ind[0]] + data[V][ind[1]] + +data[V][ind[2]] + data[V][ind[3]] + data[V][ind[4]] + data[V][ind[5]] + +data[V][ind[6]] + data[V][ind[7]]);
		}
	}
	inline double Retrieve_Value(int Gx, int Gy, int V)
	{
		double Val;
		int xi, yi;
		int rx, ry;
		double xp, yp;
		rx = Gx - this->origin.Vx+this->B*dx;
		ry = Gy - this->origin.Vy+this->B*dx;
		xi = rx / dx;
		yi = ry / dx;
		int ind[4];
		int Nx = this->grid_size.Vx + 2 * this->B;
		if (rx%dx == 0 && ry%dx == 0 )
		{
			ind[0] = index(xi, yi, Nx);
			return data[V][ind[0]];
		}
		else if (rx%dx == 0)
		{
			ind[0] = index(xi, yi, Nx);
			ind[1] = index(xi, yi+1, Nx);
			return 0.5*(data[V][ind[1]] + data[V][ind[0]]);
		}
		else if (ry%dx == 0)
		{
			ind[0] = index(xi, yi, Nx);
			ind[1] = index(xi+1, yi, Nx);
			return 0.5*(data[V][ind[1]] + data[V][ind[0]]);
		}
		else
		{
			ind[0] = index(xi, yi, Nx);
			ind[1] = index(xi + 1, yi, Nx);
			ind[2] = index(xi, yi+1, Nx);
			ind[3] = index(xi + 1, yi+1, Nx);
			return 0.25*(data[V][ind[1]] + data[V][ind[0]] + data[V][ind[2]] + data[V][ind[3]]);
		}
	}
	// *************** Buffer fill functions ****************
	void Count_Extra_And_Local_Domain(int *num_int, int *num_ext, cart_vector origin, int size,int comm,cart_vector global)
	{
		unsigned int i;
		int x, y,z;
		int xt, yt, zt;
		int comm_level;
		int xend, yend, zend;
		if (this->dim == 3)
		{
			//translate comm, correcting for maxRes
			comm_level = (comm - 1) / MESH_PARAMS.maxRes + 1;
			
			if (origin.Vx + size == MESH_PARAMS.global_size.Vx)
				xend = origin.Vx + size;
			else
				xend = origin.Vx + size - 1;
			if (origin.Vy + size == MESH_PARAMS.global_size.Vy)
				yend = origin.Vy + size;
			else
				yend = origin.Vy + size - 1;
			if (origin.Vz + size == MESH_PARAMS.global_size.Vz)
				zend = origin.Vz + size;
			else
				zend = origin.Vz + size - 1;


			for (i = 0; i < COMMS[comm_level].field.size(); i++)
			{
				for (z = 0; z < this->grid_size.Vz + 2 * this->B; z++)for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					if (!(x >= this->B && x < this->grid_size.Vx + this->B && y >= this->B && y < this->grid_size.Vy + this->B && z >= this->B && z < this->grid_size.Vz + this->B))
					{
						if (MESH_PARAMS.periodic.Vx)
							xt = ((x - this->B)*this->dx + this->origin.Vx + global.Vx) % global.Vx;//periodic
						else
							xt = ((x - this->B)*this->dx + this->origin.Vx);
						if (MESH_PARAMS.periodic.Vy)
							yt = ((y - this->B)*this->dx + this->origin.Vy + global.Vy) % global.Vy;
						else
							yt = ((y - this->B)*this->dx + this->origin.Vy);
						if (MESH_PARAMS.periodic.Vz)
							zt = ((z - this->B)*this->dx + this->origin.Vz + global.Vz) % global.Vz;
						else
							zt = ((z - this->B)*this->dx + this->origin.Vz);
						if (xt >= origin.Vx && xt < xend && yt >= origin.Vy && yt < yend && zt >= origin.Vz && zt < zend)
						{
							(*num_int)++;
						}
						else if (xt >= 0 && xt < MESH_PARAMS.global_size.Vx && yt >= 0 && yt < MESH_PARAMS.global_size.Vy && zt >= 0 && zt < MESH_PARAMS.global_size.Vz)//external but contained in another domain
						{
							(*num_ext)++;
						}
						else//set by boundary condition, ignore
						{
						}
					}
				}
			}
		}
		else if (this->dim == 2)
		{
			comm_level = (comm-1)/MESH_PARAMS.maxRes+1;
			if (origin.Vx + size == MESH_PARAMS.global_size.Vx)
				xend = origin.Vx + size;
			else
				xend = origin.Vx + size - 1;
			if (origin.Vy + size == MESH_PARAMS.global_size.Vy)
				yend = origin.Vy + size;
			else
				yend = origin.Vy + size - 1;

			for (i = 0; i < COMMS[comm_level].field.size(); i++)
			{
				for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					if (!(x >= this->B && x < this->grid_size.Vx + this->B && y >= this->B && y < this->grid_size.Vy + this->B))
					{
						if (MESH_PARAMS.periodic.Vx)
							xt = ((x - this->B)*this->dx + this->origin.Vx + global.Vx) % global.Vx;//periodic
						else
							xt = ((x - this->B)*this->dx + this->origin.Vx );
						if (MESH_PARAMS.periodic.Vy)
							yt = ((y - this->B)*this->dx + this->origin.Vy + global.Vy) % global.Vy;
						else
							yt = ((y - this->B)*this->dx + this->origin.Vy);

						if (xt >= origin.Vx && xt < xend && yt >= origin.Vy && yt < yend)
						{
							(*num_int)++;
						}
						else if (xt >= 0 && xt < MESH_PARAMS.global_size.Vx && yt >= 0 && yt < MESH_PARAMS.global_size.Vy)//external but contained in another domain
						{
							(*num_ext)++;
						}
						else
						{
						}
					}
				}
			}

		}
	}
	void Get_Coordinate_Requests(int ***coord_buffer, int *num_int, int *num_ext, cart_vector origin, int size,int comm,cart_vector global)
	{
		unsigned int i;
		int x, y, z;
		int xt, yt, zt;
		int comm_level;
		int xend, yend, zend;

		if (this->dim == 3)
		{
			
			comm_level = (comm - 1) / MESH_PARAMS.maxRes + 1;
			if (origin.Vx + size == MESH_PARAMS.global_size.Vx)
				xend = origin.Vx + size;
			else
				xend = origin.Vx + size - 1;
			if (origin.Vy + size == MESH_PARAMS.global_size.Vy)
				yend = origin.Vy + size;
			else
				yend = origin.Vy + size - 1;
			if (origin.Vz + size == MESH_PARAMS.global_size.Vz)
				zend = origin.Vz + size;
			else
				zend = origin.Vz + size - 1;

			for (i = 0; i < COMMS[comm_level].field.size(); i++)
			{
				for (z = 0; z < this->grid_size.Vz + 2 * this->B; z++)for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					if (!(x >= this->B && x < this->grid_size.Vx + this->B && y >= this->B && y < this->grid_size.Vy + this->B && z >= this->B && z < this->grid_size.Vz + this->B))
					{
						if (MESH_PARAMS.periodic.Vx)
							xt = ((x - this->B)*this->dx + this->origin.Vx + global.Vx) % global.Vx;//periodic
						else
							xt = ((x - this->B)*this->dx + this->origin.Vx);
						if (MESH_PARAMS.periodic.Vy)
							yt = ((y - this->B)*this->dx + this->origin.Vy + global.Vy) % global.Vy;
						else
							yt = ((y - this->B)*this->dx + this->origin.Vy);
						if (MESH_PARAMS.periodic.Vz)
							zt = ((z - this->B)*this->dx + this->origin.Vz + global.Vz) % global.Vz;
						else
							zt = ((z - this->B)*this->dx + this->origin.Vz);
						if (xt >= origin.Vx && xt < xend && yt >= origin.Vy && yt < yend && zt >= origin.Vz && zt < zend)
						{
							coord_buffer[0][comm][(*num_int)] = xt;
							coord_buffer[0][comm][(*num_int) + 1] = yt;
							coord_buffer[0][comm][(*num_int) + 2] = zt;
							coord_buffer[0][comm][(*num_int) + 3] = COMMS[comm_level].field[i];
							(*num_int) += 4;
						}
						else if (xt >= 0 && xt < MESH_PARAMS.global_size.Vx && yt >= 0 && yt < MESH_PARAMS.global_size.Vy && zt >= 0 && zt < MESH_PARAMS.global_size.Vz)//external but contained in another domain
						{
							coord_buffer[1][comm][(*num_ext)] = xt;
							coord_buffer[1][comm][(*num_ext) + 1] = yt;
							coord_buffer[1][comm][(*num_ext) + 2] = zt;
							coord_buffer[1][comm][(*num_ext) + 3] = COMMS[comm_level].field[i];
							(*num_ext) += 4;
						}
					}
				}
			}
		}
		else if (this->dim == 2)
		{
			comm_level = (comm - 1) / MESH_PARAMS.maxRes + 1;
			if (origin.Vx + size == MESH_PARAMS.global_size.Vx)
				xend = origin.Vx + size;
			else
				xend = origin.Vx + size - 1;
			if (origin.Vy + size == MESH_PARAMS.global_size.Vy)
				yend = origin.Vy + size;
			else
				yend = origin.Vy + size - 1;

			for (i = 0; i < COMMS[comm_level].field.size(); i++)
			{
				for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					if (!(x >= this->B && x < this->grid_size.Vx + this->B && y >= this->B && y < this->grid_size.Vy + this->B))
					{
						if (MESH_PARAMS.periodic.Vx)
							xt = ((x - this->B)*this->dx + this->origin.Vx + global.Vx) % global.Vx;//periodic
						else
							xt = ((x - this->B)*this->dx + this->origin.Vx);
						if (MESH_PARAMS.periodic.Vy)
							yt = ((y - this->B)*this->dx + this->origin.Vy + global.Vy) % global.Vy;
						else
							yt = ((y - this->B)*this->dx + this->origin.Vy);
						if (xt >= origin.Vx && xt < xend && yt >= origin.Vy && yt < yend)
						{
							coord_buffer[0][comm][(*num_int)] = xt;
							coord_buffer[0][comm][(*num_int) + 1] = yt;
							coord_buffer[0][comm][(*num_int) + 2] = COMMS[comm_level].field[i];
							(*num_int) += 3;
						}
						else if (xt >= 0 && xt < MESH_PARAMS.global_size.Vx && yt >= 0 && yt < MESH_PARAMS.global_size.Vy)//external but contained in another domain
						{
							coord_buffer[1][comm][(*num_ext)] = xt;
							coord_buffer[1][comm][(*num_ext) + 1] = yt;
							coord_buffer[1][comm][(*num_ext) + 2] = COMMS[comm_level].field[i];
							(*num_ext) += 3;
						}
						
					}
				}
			}
		}
	}
	void Load_Coordinate_Requests(double ***coord_buffer_values, int *num_int, int *num_ext, cart_vector origin, int size, int comm, cart_vector global)
	{
		unsigned int i;
		int x, y, z;
		int ind;
		int xt, yt, zt;
		int comm_level;
		int xend, yend, zend;
		if (this->dim == 3)
		{
			comm_level = (comm - 1) / MESH_PARAMS.maxRes + 1;
			if (origin.Vx + size == MESH_PARAMS.global_size.Vx)
				xend = origin.Vx + size;
			else
				xend = origin.Vx + size - 1;
			if (origin.Vy + size == MESH_PARAMS.global_size.Vy)
				yend = origin.Vy + size;
			else
				yend = origin.Vy + size - 1;
			if (origin.Vz + size == MESH_PARAMS.global_size.Vz)
				zend = origin.Vz + size;
			else
				zend = origin.Vz + size - 1;

			for (i = 0; i < COMMS[comm_level].field.size(); i++)
			{
				for (z = 0; z < this->grid_size.Vz + 2 * this->B; z++)for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					if (!(x >= this->B && x < this->grid_size.Vx + this->B && y >= this->B && y < this->grid_size.Vy + this->B && z >= this->B && z < this->grid_size.Vz + this->B))
					{
						ind = index(x, y, z, this->B * 2 + this->grid_size.Vx, this->B * 2 + this->grid_size.Vy);
						if (MESH_PARAMS.periodic.Vx)
							xt = ((x - this->B)*this->dx + this->origin.Vx + global.Vx) % global.Vx;//periodic
						else
							xt = ((x - this->B)*this->dx + this->origin.Vx);
						if (MESH_PARAMS.periodic.Vy)
							yt = ((y - this->B)*this->dx + this->origin.Vy + global.Vy) % global.Vy;
						else
							yt = ((y - this->B)*this->dx + this->origin.Vy);
						if (MESH_PARAMS.periodic.Vz)
							zt = ((z - this->B)*this->dx + this->origin.Vz + global.Vz) % global.Vz;
						else
							zt = ((z - this->B)*this->dx + this->origin.Vz);

						if (xt >= origin.Vx && xt < xend && yt >= origin.Vy && yt < yend && zt >= origin.Vz && zt < zend)
						{
							this->data[COMMS[comm_level].field[i]][ind] = coord_buffer_values[0][comm][*num_int];
							(*num_int)++;
						}
						else if (xt >= 0 && xt < MESH_PARAMS.global_size.Vx && yt >= 0 && yt < MESH_PARAMS.global_size.Vy && zt >= 0 && zt < MESH_PARAMS.global_size.Vz)//external but contained in another domain
						{
							this->data[COMMS[comm_level].field[i]][ind] = coord_buffer_values[1][comm][*num_ext];
							(*num_ext)++;
						}
					}
				}
			}
		}
		else if (this->dim == 2)
		{
			comm_level = (comm - 1) / MESH_PARAMS.maxRes + 1;
			if (origin.Vx + size == MESH_PARAMS.global_size.Vx)
				xend = origin.Vx + size;
			else
				xend = origin.Vx + size - 1;
			if (origin.Vy + size == MESH_PARAMS.global_size.Vy)
				yend = origin.Vy + size;
			else
				yend = origin.Vy + size - 1;

			for (i = 0; i < COMMS[comm_level].field.size(); i++)
			{
				for (y = 0; y < this->grid_size.Vy + 2 * this->B; y++)for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
				{
					if (!(x >= this->B && x < this->grid_size.Vx + this->B && y >= this->B && y < this->grid_size.Vy + this->B))
					{
						ind = index(x, y, this->B * 2 + this->grid_size.Vx);
						if (MESH_PARAMS.periodic.Vx)
							xt = ((x - this->B)*this->dx + this->origin.Vx + global.Vx) % global.Vx;//periodic
						else
							xt = ((x - this->B)*this->dx + this->origin.Vx);
						if (MESH_PARAMS.periodic.Vy)
							yt = ((y - this->B)*this->dx + this->origin.Vy + global.Vy) % global.Vy;
						else
							yt = ((y - this->B)*this->dx + this->origin.Vy);
						if (xt >= origin.Vx && xt < xend && yt >= origin.Vy && yt < yend)
						{
							this->data[COMMS[comm_level].field[i]][ind] = coord_buffer_values[0][comm][*num_int];
							(*num_int)++;
						}
						else if (xt >= 0 && xt < MESH_PARAMS.global_size.Vx && yt >= 0 && yt < MESH_PARAMS.global_size.Vy)//external but contained in another domain
						{
							this->data[COMMS[comm_level].field[i]][ind] = coord_buffer_values[1][comm][*num_ext];
							(*num_ext)++;
						}
					}
				}
			}
		}
	}
	// *********** OUTPUT FUNCTIONS ****************
	void Output_Restart(FILE *fp)
	{
		int i, x, y, z;
		int ind;
		fprintf(fp,"%d	%d  %d:	Field N\n",this->grid_size.Vx,this->grid_size.Vy,this->grid_size.Vz);
		fprintf(fp,"#Field Data\n");
		for (i = 0; i < this->nf; i++)
		{
			if (this->dim == 3)
			{
				for (z = this->B; z < this->grid_size.Vz + this->B; z++)for (y = this->B; y < this->grid_size.Vy + this->B; y++)for (x = this->B; x < this->grid_size.Vx + this->B; x++)
				{
					ind = index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					fprintf(fp, "%f\n", this->data[i][ind]);
				}
			}
			else if (this->dim == 2)
			{
				for (y = this->B; y < this->grid_size.Vy + this->B; y++)for (x = this->B; x < this->grid_size.Vx + this->B; x++)
				{
					ind = index(x, y, this->grid_size.Vx + 2 * this->B);
					fprintf(fp, "%f\n", this->data[i][ind]);
				}
			}
		}
	}
	int vtk_Get_Num_Points(int lowres)
	{
		if (lowres == 0)
		{
			if (this->dim == 3)
			{
				return (this->grid_size.Vx + this->B * 2)* (this->grid_size.Vy + this->B * 2)* (this->grid_size.Vz + this->B * 2);
			}
			else if (this->dim == 2)
			{
				return (this->grid_size.Vx + this->B * 2)*(this->grid_size.Vy + this->B * 2) * 2;
			}
		}
		else
		{
			int x_size,y_size,z_size;
			int inc, cx, cy, cz;
			if (this->dim == 3)
			{
				cx = 1+ceil(grid_size.Vx / (pow(2, lowres)));
				cy = 1+ceil(grid_size.Vy / (pow(2, lowres)));
				cz = 1+ceil(grid_size.Vz / (pow(2, lowres)));
				return cx*cy*cz;
//				return (1 + ceil((double)grid_size.Vx / (pow(2, lowres))))*(2 + ceil((double)grid_size.Vy / (pow(2, lowres))))*(2 + ceil((double)grid_size.Vz / (pow(2, lowres))));
//				return 4 * 4 * 4;
			}
			else if (this->dim == 2)
			{
				cx = 1+ceil(grid_size.Vx / (pow(2, lowres)));
				cy = 1+ceil(grid_size.Vy / (pow(2, lowres)));
				return 2 * cx*cy;
//				return 2*(1 + ceil((double)grid_size.Vx / (pow(2, lowres))))*(2 + ceil((double)grid_size.Vy / (pow(2, lowres))));
//				return 4 * 4 * 2;
			}

		}
	}
	int vtk_Get_Num_Cells(int lowres)
	{
		if (lowres == 0)
		{
			if (this->dim == 3)
			{
				return ((this->grid_size.Vx + this->B * 2) - 1)*((this->grid_size.Vy + this->B * 2) - 1)*((this->grid_size.Vz + this->B * 2) - 1);
			}
			else if (this->dim == 2)
			{
				return ((this->grid_size.Vx + this->B * 2) - 1)*((this->grid_size.Vy + this->B * 2) - 1);
			}
		}
		else
		{
			int cx, cy, cz;
			if (this->dim == 3)
			{
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				cz = ceil(grid_size.Vz / (pow(2, lowres)));
				return cx*cy*cz;

//				return 3*3*3;
//				return (0 + ceil((double)grid_size.Vx / (pow(2, lowres))))*(2-1 + ceil((double)grid_size.Vy / (pow(2, lowres))))*(2 + ceil((double)grid_size.Vz / (pow(2, lowres))));
			}
			else if (this->dim == 2)
			{
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				return cx*cy;

//				return (0 + ceil((double)grid_size.Vx / (pow(2, lowres))))*(2-1 + ceil((double)grid_size.Vy / (pow(2, lowres))));
//				return 3*3;
			}
		}
	}
	void vtk_Output_Points(FILE *fp,int lowres)
	{
		int x, y, z;
		if (lowres == 0)
		{
			if (this->dim == 3)
			{
				for (z = 0; z < this->grid_size.Vz + 2 * B; z++)for (y = 0; y < this->grid_size.Vy + 2 * B; y++)
				{
					for (x = 0; x < this->grid_size.Vx + 2 * B; x++)
					{
						fprintf(fp, "%6.1d %6.1d %6.1d ", this->origin.Vx + (x - B)*dx, this->origin.Vy + (y - B)*dx, this->origin.Vz + (z - B)*dx);
					}
					fprintf(fp, "\n");
				}
			}
			else if (this->dim == 2)
			{
				for (z = 0; z < 2; z++)for (y = 0; y < this->grid_size.Vy + 2 * B; y++)
				{
					for (x = 0; x < this->grid_size.Vx + 2 * B; x++)
					{
						fprintf(fp, "%6.1d %6.1d %6.1d ", this->origin.Vx + (x - B)*dx, this->origin.Vy + (y - B)*dx, this->origin.Vz + z + this->dx);
					}
					fprintf(fp, "\n");
				}
			}
		}
		else
		{
			if (this->dim == 3)
			{
				int indx[2000], indy[2000], indz[2000];
				int inc,cx, cy, cz;
				inc = pow(2,lowres);
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				cz = ceil(grid_size.Vz / (pow(2, lowres)));

				for (int i = 0; i < cx; i++)
					indx[i] = i*inc;
				indx[cx] = grid_size.Vx +  B - 1;

				for (int i = 0; i < cy; i++)
					indy[i] = i*inc;
				indy[cy] = grid_size.Vy + B - 1;

				for (int i = 0; i < cz; i++)
					indz[i] = i*inc;
				indz[cz] = grid_size.Vz + B - 1;

				/*
				indx[0] = -B;
				indx[1] = 0;
				indx[2] = grid_size.Vx-1;
				indx[3] = grid_size.Vx-1+B;
				indy[0] = -B;
				indy[1] = 0;
				indy[2] = grid_size.Vy - 1;
				indy[3] = grid_size.Vy - 1 + B;
				indz[0] = -B;
				indz[1] = 0;
				indz[2] = grid_size.Vz - 1;
				indz[3] = grid_size.Vz - 1 + B;
				*/

				for (z = 0; z < cz+1; z++)for (y = 0; y < cy+1; y++)
				{
					for (x = 0; x < cx+1; x++)
					{
						fprintf(fp, "%6.1d %6.1d %6.1d ", this->origin.Vx + (indx[x])*dx, this->origin.Vy + (indy[y])*dx, this->origin.Vz + (indz[z])*dx);
					}
					fprintf(fp, "\n");
				}				
			}
			else if (this->dim == 2)
			{
				int indx[2000], indy[2000];
				int inc, cx, cy;
				inc = pow(2, lowres);
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
								
				for (int i = 0; i < cx; i++)
					indx[i] = i*inc;
				indx[cx] = grid_size.Vx +  B - 1;

				for (int i = 0; i < cy; i++)
					indy[i] = i*inc;
				indy[cy] = grid_size.Vy +  B - 1;

				/*
				int indx[4],indy[4];
				indx[0] = -B;
				indx[1] = 0;
				indx[2] = grid_size.Vx - 1;
				indx[3] = grid_size.Vx - 1 + B;
				indy[0] = -B;
				indy[1] = 0;
				indy[2] = grid_size.Vy - 1;
				indy[3] = grid_size.Vy - 1 + B;
				*/

				for (z = 0; z < 2; z++)for (y = 0; y < cy+1; y++)
				{
					for (x = 0; x < cx+1; x++)
					{
						fprintf(fp, "%6.1d %6.1d %6.1d ", this->origin.Vx + (indx[x])*dx, this->origin.Vy + (indy[y])*dx, this->origin.Vz + z + this->dx);
					}
					fprintf(fp, "\n");
				}

			}
		}
	}
	void vtk_Output_Point_Data(int V, FILE *fp,int lowres)
	{
		int x, y, z,ind;
		if (lowres == 0)
		{
			if (this->dim == 3)
			{
				for (z = 0; z < this->grid_size.Vz + 2 * B; z++)for (y = 0; y < this->grid_size.Vy + 2 * B; y++)
				{
					for (x = 0; x < this->grid_size.Vx + 2 * B; x++)
					{
						ind = index(x, y, z, this->grid_size.Vx + this->B * 2, this->grid_size.Vy + this->B * 2);
						fprintf(fp, "%1.4f ", data[V][ind]);
					}
					fprintf(fp, "\n");
				}
			}
			else if (this->dim == 2)
			{
				for (z = 0; z < 2; z++)for (y = 0; y < this->grid_size.Vy + 2 * B; y++)
				{
					for (x = 0; x < this->grid_size.Vx + 2 * B; x++)
					{
						ind = index(x, y, this->grid_size.Vx + this->B * 2);
						fprintf(fp, "%1.4f ", data[V][ind]);
					}
					fprintf(fp, "\n");
				}
			}
		}
		else
		{
			if (this->dim == 3)
			{
				int indx[2000], indy[2000], indz[2000];
				int inc, cx, cy, cz;
				inc = pow(2, lowres);
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				cz = ceil(grid_size.Vz / (pow(2, lowres)));

				for (int i = 0; i < cx; i++)
					indx[i] = B+i*inc;
				indx[cx] = grid_size.Vx + 2 * B - 1;

				for (int i = 0; i < cy; i++)
					indy[i] = B+i*inc;
				indy[cy] = grid_size.Vy + 2 * B - 1;

				for (int i = 0; i < cz; i++)
					indz[i] = B+i*inc;
				indz[cz] = grid_size.Vz + 2 * B - 1;
				/*
				int indx[4], indy[4], indz[4];
				indx[0] = 0;
				indx[1] = B;
				indx[2] = grid_size.Vx - 1+B;
				indx[3] = grid_size.Vx - 1 + 2*B;
				indy[0] = 0;
				indy[1] = B;
				indy[2] = grid_size.Vy - 1+B;
				indy[3] = grid_size.Vy - 1 + 2*B;
				indz[0] = 0;
				indz[1] = B;
				indz[2] = grid_size.Vz - 1+B;
				indz[3] = grid_size.Vz - 1 + 2*B;
				*/
				for (z = 0; z < cz+1; z++)for (y = 0; y < cy+1; y++)
				{
					for (x = 0; x < cx+1; x++)
					{
						ind = index(indx[x], indy[y], indz[z], this->grid_size.Vx + this->B * 2, this->grid_size.Vy + this->B * 2);
						fprintf(fp, "%1.4f ", data[V][ind]);
					}
					fprintf(fp, "\n");
				}
			}
			else if (this->dim == 2)
			{
				int indx[2000], indy[2000];
				int inc, cx, cy;
				inc = pow(2, lowres);
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				
				for (int i = 0; i < cx; i++)
					indx[i] = B + i*inc;
				indx[cx] = grid_size.Vx + 2 * B - 1;

				for (int i = 0; i < cy; i++)
					indy[i] = B + i*inc;
				indy[cy] = grid_size.Vy + 2 * B - 1;

				/*
				int indx[4], indy[4];
				indx[0] = 0;
				indx[1] = B;
				indx[2] = grid_size.Vx - 1+B;
				indx[3] = grid_size.Vx - 1 + 2*B;
				indy[0] = 0;
				indy[1] = B;
				indy[2] = grid_size.Vy - 1+B;
				indy[3] = grid_size.Vy - 1 + 2*B;
				*/
				for (z = 0; z < 2; z++)for (y = 0; y < cy+1; y++)
				{
					for (x = 0; x < cx+1; x++)
					{
						ind = index(indx[x], indy[y], this->grid_size.Vx + this->B * 2);
						fprintf(fp, "%1.4f ", data[V][ind]);
					}
					fprintf(fp, "\n");
				}
			}
			
		}
	}
	int vtk_Output_Cells(int np2, FILE *fp,int lowres)
	{
		int x, y, z;
		int v0, v1, v2, v3, v4, v5, v6, v7;
		if (lowres == 0)
		{
			if (this->dim == 3)
			{
				for (z = 0; z < this->grid_size.Vz - 1 + this->B * 2; z++)for (y = 0; y < this->grid_size.Vy - 1 + this->B * 2; y++)for (x = 0; x < this->grid_size.Vx - 1 + this->B * 2; x++)
				{
					v0 = np2 + index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v1 = np2 + index(x + 1, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v2 = np2 + index(x, y + 1, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v3 = np2 + index(x + 1, y + 1, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v4 = np2 + index(x, y, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v5 = np2 + index(x + 1, y, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v6 = np2 + index(x, y + 1, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v7 = np2 + index(x + 1, y + 1, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", v0, v1, v2, v3, v4, v5, v6, v7);
				}
				return (this->grid_size.Vx + 2 * this->B)*(this->grid_size.Vy + 2 * this->B)*(this->grid_size.Vz + 2 * this->B);
			}
			else if (this->dim == 2)
			{
				for (z = 0; z < 1; z++)for (y = 0; y < this->grid_size.Vy - 1 + this->B * 2; y++)for (x = 0; x < this->grid_size.Vx - 1 + this->B * 2; x++)
				{
					v0 = np2 + index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v1 = np2 + index(x + 1, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v2 = np2 + index(x, y + 1, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v3 = np2 + index(x + 1, y + 1, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v4 = np2 + index(x, y, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v5 = np2 + index(x + 1, y, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v6 = np2 + index(x, y + 1, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					v7 = np2 + index(x + 1, y + 1, z + 1, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
					fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", v0, v1, v2, v3, v4, v5, v6, v7);
				}
				return (this->grid_size.Vx + this->B * 2)*(this->grid_size.Vy + this->B * 2) * 2;
			}
		}
		else
		{
			if (this->dim == 3)
			{
				int inc, cx, cy,cz;
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				cz = ceil(grid_size.Vz / (pow(2, lowres)));
				for (z = 0; z < cz; z++)for (y = 0; y < cy; y++)for (x = 0; x < cx; x++)
				{
					v0 = np2 + index(x, y, z, cx + 1, cy + 1);
					v1 = np2 + index(x + 1, y, z, cx + 1, cy + 1);
					v2 = np2 + index(x, y + 1, z, cx + 1, cy + 1);
					v3 = np2 + index(x + 1, y + 1, z, cx + 1, cy + 1);
					v4 = np2 + index(x, y, z + 1, cx + 1, cy + 1);
					v5 = np2 + index(x + 1, y, z + 1, cx + 1, cy + 1);
					v6 = np2 + index(x, y + 1, z + 1, cx + 1, cy + 1);
					v7 = np2 + index(x + 1, y + 1, z + 1, cx + 1, cy + 1);
					fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", v0, v1, v2, v3, v4, v5, v6, v7);
				}
				return (cx+1)*(cy+1)*(cz+1);
			}
			else if (this->dim == 2)
			{
				int inc, cx, cy;
				cx = ceil(grid_size.Vx / (pow(2, lowres)));
				cy = ceil(grid_size.Vy / (pow(2, lowres)));
				for (z = 0; z < 1; z++)for (y = 0; y < cy; y++)for (x = 0; x < cx; x++)
				{
					v0 = np2 + index(x, y, z, cx + 1, cy + 1);
					v1 = np2 + index(x + 1, y, z, cx + 1, cy + 1);
					v2 = np2 + index(x, y + 1, z, cx + 1, cy + 1);
					v3 = np2 + index(x + 1, y + 1, z, cx + 1, cy + 1);
					v4 = np2 + index(x, y, z + 1, cx + 1, cy + 1);
					v5 = np2 + index(x + 1, y, z + 1, cx + 1, cy + 1);
					v6 = np2 + index(x, y + 1, z + 1, cx + 1, cy + 1);
					v7 = np2 + index(x + 1, y + 1, z + 1, cx + 1, cy + 1);
					fprintf(fp, "8 %d %d %d %d %d %d %d %d\n",v0, v1, v2, v3, v4, v5, v6, v7);
				}
				return 2*(cx+1)*(cy+1);
			}
		}
	}
	
	// ************** TMP FUNCTIONS ********************
	void Output_To_Screen(int num)
	{
		int x, y, z, f, ind;
		FILE *fp;
		char filename[BUFSIZ];
		sprintf(filename, "out_%d.txt", num);
		fp = fopen(filename, "a");

		if (dim == 3)
		{
			fprintf(fp, "Origin: %d %d %d\n", this->origin.Vx, this->origin.Vy, this->origin.Vz);
			for (f = 0; f < this->nf; f++)
			{
				fprintf(fp, "field = %d\n", f);
				for (z = 0; z < this->grid_size.Vz + 2 * this->B; z++)
				{
					fprintf(fp, "z=%d\n", this->origin.Vz - this->B*this->dx + z*this->dx);
					for (y = this->grid_size.Vy + 2 * this->B - 1; y >= 0; y--)
					{
						for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
						{
							ind = index(x, y, z, this->grid_size.Vx + 2 * this->B, this->grid_size.Vy + 2 * this->B);
							fprintf(fp, "%3.4f ", this->data[f][ind]);
						}
						fprintf(fp, "\n");
					}
					fprintf(fp, "\n");
				}
				fprintf(fp, " ********** \n");
			}
			fclose(fp);
		}
		else if (dim == 2)
		{
			fprintf(fp, "Origin: %d %d\n", this->origin.Vx, this->origin.Vy);
			for (f = 0; f < 1;f++)// this->nf; f++)
			{
				fprintf(fp, "field = %d\n", f);
				for (y = this->grid_size.Vy + 2 * this->B - 1; y >= 0; y--)
				{
					for (x = 0; x < this->grid_size.Vx + 2 * this->B; x++)
					{
						ind = index(x, y, this->grid_size.Vx + 2 * this->B);
						fprintf(fp, "%3.2f ", this->data[f][ind]);
					}
					fprintf(fp, "\n");
				}
				fprintf(fp, " ********** \n");

			}
			fclose(fp);
		}
	}
	void initialize(int i)
	{
		int x, y, z,f, ind;
		int xpos;
		if (this->dim == 3)
		{
			for (f = 0; f < this->nf; f++)
			{
				for (z = 0; z < this->N + 2 * this->B; z++)for (y = 0; y < this->N + 2 * this->B; y++)for (x = 0; x < this->N + 2 * this->B; x++)
				{
					ind = index(x, y, z, this->N + 2 * this->B, this->N + 2 * this->B);
					xpos = this->origin.Vx + (x - this->B)*this->dx;
					data[f][ind] = xpos;
				}
			}
		}
		else if (this->dim == 2)
		{
			for (f = 0; f < this->nf; f++)
			{
				for (y = 0; y < this->N + 2 * this->B; y++)for (x = 0; x < this->N + 2 * this->B; x++)
				{
					ind = index(x, y, this->N + 2 * this->B);
					xpos = this->origin.Vx + (x - this->B)*this->dx;
					data[f][ind] = xpos;
				}
			}
		}
	}
};
#endif
