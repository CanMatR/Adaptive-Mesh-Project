/* Copyright (c) 2015 Government of Canada  */
#ifndef controlDEF
#define controlDEF
#include <time.h>
#include "region.h"
#include "utility.h"
#include "geometry.h"
#include "mpi.h"
#include "mesh.h"

class regionControl
{
public:
	int region_id;
	double **send_val;
	double **recv_val;
	int **send_buffer;//buffer of coordinates and domains to be sent
	int **recv_buffer;
	int *recv_buffer_size;
	int *send_buffer_size;
	int *send_val_size;
	int *recv_val_size;
	int recv_val_index;
	int num_comms;

	regionControl()
	{
		region_id = -1;
		send_val = NULL;
		send_buffer = NULL;
		recv_buffer = NULL;
		recv_buffer_size = 0;
		send_buffer_size = 0;
	}
	~regionControl()
	{}
	void output_recv_buffer(int comm)
	{
		printf("Output rcv_buffer for Region : %d\n", region_id);
		int i;
		for (i = 0; i<recv_buffer_size[comm]; i++)
		{
			printf("%d\n", recv_buffer[comm][i]);
		}
	}
};
class domainControl
{
public:
	int size;//size of the domain at lowest grid spacing i.e. global uniform mesh
	int domain_weight;//Computational weight of the domain(i.e. number of nodes)
	cart_vector origin;//Origin of the domain

	int region_id;//The id of the region the domain belongs to.
	//Communication
	int **recv_coords;//Coordinates the current cpu requests from this domain
	double **recv_value;
	int recv_index;
	int *recv_num;
	int num_comms;
	char outfile[200];
	//Functions
	domainControl()
	{
		region_id = 0;
		domain_weight = 0;
		size = 0;
		origin.Vx = 0;
		origin.Vy = 0;
		origin.Vz = 0;
		recv_coords = NULL;
		recv_value = NULL;
	}
	~domainControl()
	{

	}
	void set_origin(cart_vector origin)
	{
		this->origin = origin;
	}
	void free_memory(int comm)
	{
		if (recv_coords[comm] != NULL)
		{
			delete recv_coords[comm];
			recv_coords[comm] = NULL;
		}

		if (recv_value[comm] != NULL)
		{
			delete recv_value[comm];
			recv_value[comm] = NULL;
		}
	}
	void allocate_memory(int d, int comm)
	{
		int num;
		num = recv_num[comm];
		if (recv_coords[comm] == NULL)
			recv_coords[comm] = new int[num*(d + 1)];//coords + id
		else
		{
			printf("Error memory allocated to recv_coords when already allocated in allocate_memory\n");
			exit(0);
		}
		if (recv_value[comm] == NULL)
			recv_value[comm] = new double[num];//changed num*d to num
		else
		{
			printf("Error memory allocated to recv_value when already allocated in allocate_memory\n");
			exit(0);
		}
	}
};
class control
{
	region *Reg;// Region associated with this control processor
	//parallel use
	//control area definitions
	int D[3];//Number of domains in each direction
	regionControl *RegionMap;
	domainControl *DomainMap;//uniform mesh of domaincontrol with D[0]*D[1]*D[2] elements

public:
	int myid;
	int ncpus;
	int num_comms;
	int dim;//number of dimensions
	control()
	{
		RegionMap = NULL;
		DomainMap = NULL;
		Reg = NULL;
	}
	control(int myid, int ncpus)
	{
		RegionMap = NULL;
		DomainMap = NULL;
		Reg = NULL;
		this->myid = myid;
		this->ncpus = ncpus;
	}	
	~control()
	{
		if (RegionMap != NULL)
			delete RegionMap;
		if (DomainMap != NULL)
			delete DomainMap;
		if (Reg != NULL)
			delete Reg;
	}
	void Adaption_Cycle(int t)
	{
		int i;
/*		Reg->Unadapt_Domains();
		this->Build_Communication_Requests();
		for (i = 1; i < MESH_PARAMS.comm_calls; i++)
		{
			this->Communicate(i);
		}
		*/
		Reg->Adapt_Domains();
		this->Build_Adaption_Communication_Requests();
		Adaption_Parallel();
		this->Build_Communication_Requests();
	}
	void Adaption_Parallel()
	{
		int i, j;
		int m = 2, n = 0;
		int dx;
		int dom_id;
		cart_vector coord;
		for (i = 0; i<ncpus; i++)
		{
			j = 2;
			m = 2;
			while (j<RegionMap[i].send_buffer_size[0])
			{
				if (j == m)//get domain header information
				{
					dom_id = RegionMap[i].send_buffer[0][m];
					m = m + 2 + RegionMap[i].send_buffer[0][j + 1] * (this->dim + 1);//4 is the number of coordinate information variables, i.e. x,y,z, id
					j += 2;
				}
				else
				{
					coord.Vx = RegionMap[i].send_buffer[0][j];
					coord.Vy = RegionMap[i].send_buffer[0][j + 1];
					if (this->dim == 3)
					{
						coord.Vz = RegionMap[i].send_buffer[0][j + 2];
						dx = RegionMap[i].send_buffer[0][j + 3];
					}
					else
						dx = RegionMap[i].send_buffer[0][j + 2];
					Reg->Adapt_Parallel(dom_id, coord, dx);
					j += (this->dim + 1);
				}
			}
		}
	}
	
	void Create_Initialization()
	{
		int x, y, z,ind;
		FILE *fp;
		int num_procs;
		cart_vector origin;
		char inputline[BUFSIZ];
		char filename[200];
		domain *t_Dom;
		dom_params DP;

		fp = fopen("init_geometry.dat", "r");
		fgets(inputline, BUFSIZ, fp);//num_comms
			sscanf(inputline, "%d", &(this->num_comms));
			MESH_PARAMS.num_comms = this->num_comms;
			MESH_PARAMS.comm_calls = this->num_comms;
		fgets(inputline, BUFSIZ, fp);//Number of dimensions
			sscanf(inputline, "%d", &(this->dim));
			MESH_PARAMS.DIMENSION = this->dim;
		fgets(inputline, BUFSIZ, fp);//Number of processors
			sscanf(inputline, "%d", &(num_procs));
		fgets(inputline, BUFSIZ, fp);//Dx domains
			sscanf(inputline, "%d %d %d", &(this->D[0]), &(this->D[1]), &(this->D[2]));
			MESH_PARAMS.Num_Dom.Vx = this->D[0];
			MESH_PARAMS.Num_Dom.Vy = this->D[1];
			MESH_PARAMS.Num_Dom.Vz = this->D[2];
		fgets(inputline, BUFSIZ, fp);// mesh geometry
			sscanf(inputline, "%d %d %d", &(DP.maxRes),&DP.uniRes,&DP.buffer);
			MESH_PARAMS.maxRes = DP.maxRes;
			MESH_PARAMS.uniRes = DP.uniRes;
			MESH_PARAMS.buffer_size = DP.buffer;
			//update numcomms with max Res;
			this->num_comms = (this->num_comms-1)*MESH_PARAMS.maxRes+1;//comm 0 adaption, all others by level
			MESH_PARAMS.num_comms = this->num_comms;
		fgets(inputline, BUFSIZ, fp);//Boundary Symmetry
			sscanf(inputline, "%d %d %d", &(MESH_PARAMS.periodic.Vx), &MESH_PARAMS.periodic.Vy, &MESH_PARAMS.periodic.Vz);
		fgets(inputline, BUFSIZ, fp);// number of fields
			sscanf(inputline, "%d %d", &(DP.num_fields), &DP.num_aux_fields);
			MESH_PARAMS.num_fields = DP.num_fields;
		this->Create_Domain_Map();//create domain maps with default values
		DP.numComms = this->num_comms;
		DP.dim = this->dim;
		DP.size = (int)pow(2, DP.maxRes + DP.uniRes - 1) + 1;
		
		if (MESH_PARAMS.periodic.Vx)
			DP.global_size.Vx = this->D[0] * (DP.size-1);
		else
			DP.global_size.Vx = this->D[0] * (DP.size - 1) + 1;

		if (MESH_PARAMS.periodic.Vy)
			DP.global_size.Vy = this->D[1] * (DP.size-1);
		else
			DP.global_size.Vy = this->D[1] * (DP.size - 1) + 1;
		if (MESH_PARAMS.periodic.Vz)
			DP.global_size.Vz = this->D[2] * (DP.size-1);
		else
			DP.global_size.Vz = this->D[2] * (DP.size - 1) + 1;
		MESH_PARAMS.global_size.Vx = DP.global_size.Vx;
		MESH_PARAMS.global_size.Vy = DP.global_size.Vy;
		MESH_PARAMS.global_size.Vz = DP.global_size.Vz;
		//printf("%d\n",DP.size);
		//fflush(stdout);
		for (z = 0; z < this->D[2]; z++)for (y = 0; y < this->D[1]; y++)for (x = 0; x < this->D[0]; x++)
		{
			ind = index(x, y, z, this->D[0], this->D[1]);
			origin.Vx = x*(DP.size-1);
			origin.Vy = y*(DP.size-1);
			origin.Vz = z*(DP.size-1);
			sprintf(filename,"domain_%d_%d_%d.dat",origin.Vx,origin.Vy,origin.Vz);
			sprintf(this->DomainMap[ind].outfile,"%s", filename);
			this->DomainMap[ind].set_origin(origin);
			this->DomainMap[ind].size = DP.size;
			this->Init_Domain_Dist();
		}
		//function to set region ids in domainmap
		this->create_region();
		this->create_region_map();
		for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
		{
			ind = index(x, y, z, D[0], D[1]);
			if (DomainMap[ind].region_id == myid)
			{
				DP.O.Vx = x*(DP.size-1);
				DP.O.Vy = y*(DP.size-1);
				DP.O.Vz = z*(DP.size-1);
				DP.dom_id = ind;
				t_Dom = new domain(DP);
				Reg->Dom.push_back(t_Dom);
				Reg->num_domain++;
			}
		}
		// *********** NOW LOAD GEOMETRY FOR INITIAL ADAPTIONS **********
		//Change this after determining communication method
		fgets(inputline, BUFSIZ, fp);//Field Communications
		int i,j,comm_num,fld_num;
		communications comm;
		for (i = 0; i < DP.numComms; i++)
		{
			COMMS.push_back(comm);
		}
		for (i = 0; i < DP.num_fields+DP.num_aux_fields; i++)
		{
			fgets(inputline, BUFSIZ, fp);
			sscanf(inputline, "%d %d", &(fld_num), &(comm_num));
			COMMS[comm_num].push_comm(fld_num);
	//		for (j = 0; j < DP.maxRes; j++)
//			{
//				ind = (comm_num - 1)*DP.maxRes + j;
//				COMMS[ind].push_comm(fld_num);
//		
//			}
		}
	
		int num_geom;
		adapt_geometry **geom;
		fgets(inputline, BUFSIZ, fp);//Initial Geometry
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d", &(num_geom));
		
		geom = Init_Load_Seed_Geom(fp,num_geom);
		Seed_Adaption(geom,num_geom);
		
		//delete geometry
		for (i = 0; i < num_geom; i++)
			if (geom[i] != NULL)
				delete geom[i];
		if (geom != NULL)
			delete geom;

		//Initialize field values for fields
		fclose(fp);
	}
	adapt_geometry **Init_Load_Seed_Geom(FILE *fp,int num_geom)
	{
		int i,type;
		char inputline[BUFSIZ];
		adapt_geometry **ad_geo;
		circle *circ;
		sphere *sphr;
		ad_geo = new adapt_geometry*[num_geom];
		cart_vector_d origin;
		double R;
		for (i = 0; i < num_geom; i++)
		{
			fgets(inputline, BUFSIZ, fp);
			sscanf(inputline, "%d", &(type));
			switch (type)
			{
			case 1://circle(2D) /sphere(3D)
				if (this->dim == 3)
				{
					sphr = new sphere();
					fgets(inputline, BUFSIZ, fp);
					sscanf(inputline, "%lf %lf %lf %lf", &(R), &(origin.Vx), &(origin.Vy), &(origin.Vz));
					sphr->set_geom(origin, R);
					ad_geo[i] = sphr;
				}
				else if (this->dim == 2)
				{
					circ = new circle();
					origin.Vz = 0.0;
					fgets(inputline, BUFSIZ, fp);
					sscanf(inputline, "%lf %lf %lf", &(R), &(origin.Vx), &(origin.Vy));
					circ->set_geom(origin, R);
					ad_geo[i] = circ;
				}
				break;
			}
		}
		return ad_geo;
	}
	void Init_Domain_Dist()
	{
		int ind,x,y,z;
		int xp, yp, xyp;
		int indz,xpz,ypz,xypz;
		int mx,my,mz;
		int idx, idy, idz;
		int Np;
		int proc = 0;
		for (z = 0; z < this->D[2]; z++)
		  {
			proc = 0;
			for (y = 0; y < this->D[1]; y++)for (x = 0; x < this->D[0]; x++)
			{
				ind = index(x, y, z, this->D[0], this->D[1]);
				this->DomainMap[ind].region_id = ind%this->ncpus;
/*				xp = index(x + 1, y, z, this->D[0], this->D[1]);
				yp = index(x, y + 1, z, this->D[0], this->D[1]);
				xyp = index(x + 1, y + 1, z, this->D[0], this->D[1]);

				indz = index(x, y, z+1, this->D[0], this->D[1]);
                                xpz = index(x + 1, y, z+1, this->D[0], this->D[1]);
                                ypz = index(x, y + 1, z+1, this->D[0], this->D[1]);
                                xypz = index(x + 1, y + 1, z+1, this->D[0], this->D[1]);

				//Distribute domains into regions
				//stripe in x
//				this->DomainMap[ind].region_id = ind%this->ncpus;
				this->DomainMap[ind].region_id = proc;
				this->DomainMap[xp].region_id = proc + 1;
				this->DomainMap[yp].region_id = proc + 2;
				this->DomainMap[xyp].region_id = proc + 3;

		    		this->DomainMap[indz].region_id = proc+4;
				this->DomainMap[xpz].region_id = proc+5;
				this->DomainMap[ypz].region_id = proc+6;
				this->DomainMap[xypz].region_id = proc+7;

				//			this->DomainMap[ind].region_id = ind%this->ncpus;
				proc += 8;
				*/
				/*			if (this->dim == 3)
							{
							Np = 2;
							mx = (int)(ceil((float)this->D[0]/(float)Np));
							my = (int)(ceil((float)this->D[0] / (float)Np));
							mz = (int)(ceil((float)this->D[0] / (float)Np));
							idx = (int)floor(x / mx);
							idy = (int)floor(y / my);
							idz = (int)floor(z / mz);
							this->DomainMap[ind].region_id = index(idx,idy,idz,Np,Np);
							}
							else
							{
							this->DomainMap[ind].region_id = x%this->ncpus;
							}
							*/
			}
		}
//		fgets(inputline, BUFSIZ, fp);//region_id
	//	sscanf(inputline, "%d", &reg_id);
		
	}
	void Seed_Adaption(adapt_geometry **ad_geo, int num_geom)
	{
		Reg->Seed_Adapt_Domains(ad_geo,num_geom);
	}
	void Load_Restart(int time)
	{
		Load_Control(time);
		Load_Region(time);
	}
	void Load_Region(int time)
	{
		int x, y, z, ind;
		//loop over domain map, create a domain if it belongs to this region, then load the domain in from the appropriate file
		for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
		{
			ind = index(x, y, z, D[0], D[1]);
			if (DomainMap[ind].region_id == myid)
			{
				Reg->Create_Domain(DomainMap[ind].outfile, time, ind, num_comms);
			}
		}
	}
	void Load_Control(int time)
	{
		int x, y, z, ind;
		int size, reg_id;
		cart_vector origin;
		int num_procs;
		char ctrl_rest_file[BUFSIZ];
		char inputline[BUFSIZ];
		sprintf(ctrl_rest_file, "./restart_%d/control.restart", time);
		FILE *fp;
		fp = fopen(ctrl_rest_file, "r");

		fgets(inputline, BUFSIZ, fp);//num_comms
		sscanf(inputline, "%d", &(this->num_comms));

		fgets(inputline, BUFSIZ, fp);//Number of dimensions
		sscanf(inputline, "%d", &(this->dim));

		fgets(inputline, BUFSIZ, fp);//Number of processors
		sscanf(inputline, "%d", &(num_procs));
		if (num_procs != this->ncpus)
		{
			printf("Number of processors used %d not equal to number of processors in restart %d.  Auto redistribution to be programmed at later date.\n", this->ncpus, num_procs);
			exit(0);
		}
		fgets(inputline, BUFSIZ, fp);//Dx domains
		sscanf(inputline, "%d %d %d", &(this->D[0]), &(this->D[1]), &(this->D[2]));

		fgets(inputline, BUFSIZ, fp);//skip domainmap header
		this->Create_Domain_Map();//create domain maps with default values

		for (z = 0; z < this->D[2]; z++)for (y = 0; y < this->D[1]; y++)for (x = 0; x < this->D[0]; x++)
		{
			ind = index(x, y, z, this->D[0], this->D[1]);
			fgets(inputline, BUFSIZ, fp);//x origin
			sscanf(inputline, "%d %d %d", &(origin.Vx), &(origin.Vy), &(origin.Vz));
			fgets(inputline, BUFSIZ, fp);//size
			sscanf(inputline, "%d", &size);
			fgets(inputline, BUFSIZ, fp);//region_id
			sscanf(inputline, "%d", &reg_id);
			fgets(inputline, BUFSIZ, fp);//region_id
			sscanf(inputline, "%s", &(this->DomainMap[ind].outfile));
			this->DomainMap[ind].set_origin(origin);
			this->DomainMap[ind].size = size;
			this->DomainMap[ind].region_id = reg_id;
		}
		
		this->create_region();
		this->create_region_map();
		fclose(fp);
	}
	void create_region_map()
	{
		RegionMap = new regionControl[ncpus];
		int i, j;
		for (i = 0; i<ncpus; i++)
		{
			//Move to RegionMap init function
			RegionMap[i].region_id = i;
			RegionMap[i].num_comms = num_comms;
			RegionMap[i].recv_buffer = new int*[num_comms];
			RegionMap[i].recv_buffer_size = new int[num_comms];
			RegionMap[i].recv_val = new double*[num_comms];
			RegionMap[i].recv_val_size = new int[num_comms];
			RegionMap[i].send_buffer = new int*[num_comms];
			RegionMap[i].send_buffer_size = new int[num_comms];
			RegionMap[i].send_val = new double*[num_comms];
			RegionMap[i].send_val_size = new int[num_comms];
			for (j = 0; j<num_comms; j++)
			{
				RegionMap[i].recv_buffer[j] = NULL;
				RegionMap[i].recv_val[j] = NULL;
				RegionMap[i].send_buffer[j] = NULL;
				RegionMap[i].send_val[j] = NULL;
			}
		}
	}
	void create_region()
	{
		int i;
		/**** Create region for control****/
		Reg = new region();
		Reg->Regid = myid;
		Reg->global_size.Vx = MESH_PARAMS.global_size.Vx;
		Reg->global_size.Vy = MESH_PARAMS.global_size.Vy;
		Reg->global_size.Vz = MESH_PARAMS.global_size.Vz;

//		Reg->global_size.Vx = D[0] * DomainMap[0].size;
//		Reg->global_size.Vy = D[1] * DomainMap[0].size;
//		Reg->global_size.Vz = D[2] * DomainMap[0].size;
		Reg->dim = dim;
		Reg->recv_val = new double*[num_comms];
		Reg->recv_val_size = new int[num_comms];
		Reg->recv_coords_array = new int*[num_comms];
		for (i = 0; i<num_comms; i++)
		{
			Reg->recv_val[i] = NULL;
			Reg->recv_val_size[i] = 0;
			Reg->recv_coords_array[i] = NULL;
		}
	}
	void Create_Domain_Map()
	{
		int x, y, z, ind, i;
		cart_vector origin;
		DomainMap = new domainControl[D[0] * D[1] * D[2]];
		for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
		{
			ind = index(x, y, z, D[0], D[1]);
			DomainMap[ind].domain_weight = 0;
			DomainMap[ind].size = 0;
			DomainMap[ind].num_comms = num_comms;
			DomainMap[ind].recv_coords = new int*[num_comms];
			DomainMap[ind].recv_num = new int[num_comms];
			DomainMap[ind].recv_value = new double*[num_comms];
			DomainMap[ind].set_origin(origin);
			for (i = 0; i<num_comms; i++)
			{
				DomainMap[ind].recv_coords[i] = NULL;
				DomainMap[ind].recv_value[i] = NULL;
				DomainMap[ind].recv_num[i] = 0;
			}
		}
	}
	void Get_Field_List(list<field*> *fields)
	{
		//	Reg->delete_field_list(fields);
		Reg->Get_Field_List(fields);
	}
	// **************** COMMUNICATION FUNCTIONS **************/
	void Communicate(int comm)
	{
		int c;
		for (int i=0; i < MESH_PARAMS.maxRes; i++)
		{
			c = 1+(comm-1)*MESH_PARAMS.maxRes + i;
			Create_Send_Value(c);
			Communicate_Send_Value(c);
			Load_Send_Value(c);
			Reg->Local_Communicate(c);
		}
	}
	void Communicate_Send_Value(int comm)
	{
		int i;
		int send_to, receive_from;
		int send_tag, receive_tag;
		MPI_Request request;
		MPI_Status status;

		for (i = 0; i<ncpus; i++)
		{
			send_to = (myid - i + ncpus) % ncpus;
			receive_from = (i + myid) % ncpus;
			send_tag = ncpus*(send_to)+i;
			receive_tag = ncpus*(myid) + i;
			fflush(stdout);
			if (RegionMap[receive_from].recv_val_size[comm] > 0)
				MPI_Irecv(RegionMap[receive_from].recv_val[comm], RegionMap[receive_from].recv_val_size[comm], MPI_DOUBLE, receive_from, receive_tag, MPI_COMM_WORLD, &request);
			if (RegionMap[send_to].send_val_size[comm] > 0)
				MPI_Send(RegionMap[send_to].send_val[comm], RegionMap[send_to].send_val_size[comm], MPI_DOUBLE, send_to, send_tag, MPI_COMM_WORLD);
			if (RegionMap[receive_from].recv_val_size[comm] > 0)
				MPI_Wait(&request, &status);
		}
	}
	void Load_Send_Value(int comm)
	{
		int ind, i, x, y, z, r;
		int dom_size = DomainMap[0].size - 1;
		//load the values into the domainmap that the values came from
		for (i = 0; i<ncpus; i++)
			RegionMap[i].recv_val_index = 0;
		switch (this->dim)
		{
		case 2:
			for (y = 0; y<D[1]; y++)
			{
				for (x = 0; x<D[0]; x++)
				{
					ind = index(x, y, D[0]);
					r = DomainMap[ind].region_id;//the region the domain belongs to
					for (i = 0; i<DomainMap[ind].recv_num[comm]; i++)
					{
						DomainMap[ind].recv_value[comm][i] = RegionMap[r].recv_val[comm][RegionMap[r].recv_val_index];
						RegionMap[r].recv_val_index++;
					}
				}
			}
			//load values into the region in the same order they requested coordinates
			for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
					DomainMap[index(x, y, D[0])].recv_index = 0;

			for (i = 0; i<Reg->recv_val_size[comm]; i++)
			{
				x = Reg->recv_coords_array[comm][(this->dim + 1)*i];//replace d with (d+1) d for the coords +1 for the id
				y = Reg->recv_coords_array[comm][(this->dim + 1)*i + 1];
				x = x / dom_size;
				if (x >= D[0])
					x = D[0] - 1;
				y = y / dom_size;
				if (y >= D[1])
					y = D[1] - 1;
				ind = index(x, y, D[0]);
				Reg->recv_val[comm][i] = DomainMap[ind].recv_value[comm][DomainMap[ind].recv_index];

				DomainMap[ind].recv_index++;
			}
			break;
		case 3:
			for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
			{
				ind = index(x, y, z, D[0], D[1]);
				r = DomainMap[ind].region_id;//the region the domain belongs to
				for (i = 0; i<DomainMap[ind].recv_num[comm]; i++)
				{
					DomainMap[ind].recv_value[comm][i] = RegionMap[r].recv_val[comm][RegionMap[r].recv_val_index];
					RegionMap[r].recv_val_index++;
				}
			}

			//load values into the region in the same order they requested coordinates
			for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
				DomainMap[index(x, y, z, D[0], D[1])].recv_index = 0;

			for (i = 0; i<Reg->recv_val_size[comm]; i++)
			{
				x = Reg->recv_coords_array[comm][(this->dim + 1)*i];//replace d with (d+1) d for the coords +1 for the id
				y = Reg->recv_coords_array[comm][(this->dim + 1)*i + 1];
				z = Reg->recv_coords_array[comm][(this->dim + 1)*i + 2];
				x = x / dom_size;
				if (x >= D[0])
					x = D[0] - 1;
				y = y / dom_size;
				if (y >= D[1])
					y = D[1] - 1;
				z = z / dom_size;
				if (z >= D[2])
					z = D[2] - 1;
				ind = index(x, y,z, D[0],D[1]);


//				ind = index(x / DomainMap[0].size, y / DomainMap[0].size, z / DomainMap[0].size, D[0], D[1]);
				Reg->recv_val[comm][i] = DomainMap[ind].recv_value[comm][DomainMap[ind].recv_index];

				DomainMap[ind].recv_index++;
			}
			break;
		}
		Reg->Load_Values_Into_Domains(comm);
	}
	void Build_Adaption_Communication_Requests()
	{
		Build_Receive_Request(0);
		Build_Region_Receive_Requests(0);
		Syncronize_Region_Receive_Requests(0);
		Allocate_Memory_Send_Value(0);
	}
	void Build_Communication_Requests()
	{
		int i;
		for (i = 1; i < this->num_comms; i++)
		{
			Build_Receive_Request(i);
			Build_Region_Receive_Requests(i);
			Syncronize_Region_Receive_Requests(i);
			Allocate_Memory_Send_Value(i);
		}
	}
	void Build_Receive_Request(int comm)
	{
		int ind, i, x, y, z;
		int *recv_coords = NULL;//Coords to be received, if 2D x,y, if 3D x,y,z; i.e. [i]->x [i+1]->y
		int recv_coords_size = 0, dom_size = DomainMap[0].size-1;
		Reg->Get_Comm_Coords(&recv_coords, &recv_coords_size, comm);
		
		if (this->dim == 2)
		{
			
			//Set recv_num and recv_index to zero and free up memory in domain control for this communication
			for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
			{
				ind = index(x, y, D[0]);
				DomainMap[ind].recv_num[comm] = 0;
				DomainMap[ind].recv_index = 0;
				DomainMap[ind].free_memory(comm);
			}
			
			//Count the number of coordinates to be entered into each domain control
			for (i = 0; i<recv_coords_size; i += (this->dim + 1))//d is the dimension i.e. x,y and +1 for the id
			{
				if ((recv_coords[i] < 0) || (recv_coords[i] > D[0] * dom_size+1) || (recv_coords[i + 1] < 0) || (recv_coords[i + 1] > D[1] * dom_size+1))
				{
					fprintf(stderr, "Coordinate request outside global domain\n");
					while (1);
					exit(0);
				}
				x = recv_coords[i] / dom_size;
				if (x >= D[0])
					x = D[0] - 1;
				y = recv_coords[i + 1] / dom_size;
				if (y >= D[1])
					y = D[1] - 1;
				ind = index(x, y, D[0]);
				
				DomainMap[ind].recv_num[comm]++;				
			}
			
			//Allocate the array sizes in domain control
			for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
			{
				ind = index(x, y, D[0]);
				DomainMap[ind].allocate_memory(this->dim, comm);
			}
			
			//Load coordinates into DomainControl
			for (i = 0; i<recv_coords_size; i += (this->dim + 1))//d is the dimension i.e. x,y and +1 for the id
			{//assume domain sizes are the same to find the index in DomainControl of the coordinate
				x = recv_coords[i] / dom_size;
				if (x >= D[0])
					x = D[0] - 1;
				y = recv_coords[i + 1] / dom_size;
				if (y >= D[1])
					y = D[1] - 1;
				ind = index(x, y, D[0]);

				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index] = recv_coords[i];
				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index + 1] = recv_coords[i + 1];
				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index + 2] = recv_coords[i + 2];
				DomainMap[ind].recv_index += (this->dim + 1);//d is the dimension i.e. x,y and +1 for the id
			}
			
		}
		else if (this->dim == 3)
		{
			
			//Set recv_num and recv_index to zero and free up memory in domain control for this communication
			for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
			{
				ind = index(x, y, z, D[0], D[1]);
				DomainMap[ind].recv_num[comm] = 0;
				DomainMap[ind].recv_index = 0;
				DomainMap[ind].free_memory(comm);
			}
			//Count the number of coordinates to be entered into each domain control
			for (i = 0; i<recv_coords_size; i += (this->dim + 1))//d is the dimension i.e. x,y,z and +1 for the id
			{
				if ((recv_coords[i] < 0) || (recv_coords[i] > D[0] * DomainMap[0].size) || (recv_coords[i + 1] < 0) || (recv_coords[i + 1] > D[1] * DomainMap[0].size) || (recv_coords[i + 2] < 0) || (recv_coords[i + 2] > D[2] * DomainMap[0].size))
				{
					fprintf(stderr, "Coordinate request outside global domain\n");
					while (1);
					exit(0);
				}
				x = recv_coords[i] / dom_size;
				if (x >= D[0])
					x = D[0] - 1;
				y = recv_coords[i + 1] / dom_size;
				if (y >= D[1])
					y = D[1] - 1;
				z = recv_coords[i + 2] / dom_size;
				if (z >= D[2])
					z = D[2] - 1;
				ind = index(x, y, z, D[0], D[1]);

//				ind = index(recv_coords[i] / dom_size, recv_coords[i + 1] / dom_size, recv_coords[i + 2] / dom_size, D[0], D[1]);
				DomainMap[ind].recv_num[comm]++;
			}
			//Allocate the array sizes in domain control
			for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
			{
				ind = index(x, y, z, D[0], D[1]);
				DomainMap[ind].allocate_memory(this->dim, comm);
			}
			//Load coordinates into DomainControl
			for (i = 0; i<recv_coords_size; i += (this->dim + 1))//d is the dimension i.e. x,y and +1 for the id
			{//assume domain sizes are the same to find the index in DomainControl of the coordinate

				x = recv_coords[i] / dom_size;
				if (x >= D[0])
					x = D[0] - 1;
				y = recv_coords[i + 1] / dom_size;
				if (y >= D[1])
					y = D[1] - 1;
				z = recv_coords[i + 2] / dom_size;
				if (z >= D[2])
					z = D[2] - 1;
				ind = index(x, y, z, D[0], D[1]);
				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index] = recv_coords[i];
				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index + 1] = recv_coords[i + 1];
				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index + 2] = recv_coords[i + 2];
				DomainMap[ind].recv_coords[comm][DomainMap[ind].recv_index + 3] = recv_coords[i + 3];
				DomainMap[ind].recv_index += (this->dim + 1);//d is the dimension i.e. x,y and +1 for the id
			}
			
		}
		
		if (recv_coords != NULL)
			delete recv_coords;
			
	}
	/****************************
	Function: Build_Region_Recieve_Requests
	Purpose: To assemble the send buffer for each region and calculating its size.  Created from the domain control recv_coords.
	The send buffer is of type integer and has the following formatted layout:
	R_id	#The region that made the request(should be the current CPU constructing the array)
	D_num	#The number of domains this region requires coordinates of from the target region
	D_ind_1	#The index of the domain in the domain control
	C_num_1	#The number of coordinates from the domain
	X_1
	Y_1
	Z_1
	X_2
	...
	D_ind_2	#The index of the domain in the domain control
	C_num_2	#The number of coordinates from the domain
	X_1
	Y_1
	Z_1
	X_2
	...
	*****************************/
	void Build_Region_Receive_Requests(int comm)
	{
		int i, r, x, y, z, ind;
		int coord_count;
		int Dom_count;
		int rcv_ind;
		if (this->dim == 2)
		{
			for (r = 0; r<ncpus; r++)//loop over regions
			{
				coord_count = 0;
				Dom_count = 0;
				for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
				{
					ind = index(x, y, D[0]);
					if (DomainMap[ind].region_id == RegionMap[r].region_id)
					{
						if (DomainMap[ind].recv_num[comm] > 0)
						{
							coord_count += DomainMap[ind].recv_num[comm];
							Dom_count++;//temporarily use as a domain count to allocate memory for recv buffer
						}
					}
				}
				RegionMap[r].recv_val_size[comm] = coord_count;
				if (RegionMap[r].recv_val[comm] != NULL)
				{
					delete RegionMap[r].recv_val[comm];
				}
				if (coord_count > 0)
				{
					RegionMap[r].recv_val[comm] = new double[coord_count];//allocate memory to receive incoming values
				}
				else
					RegionMap[r].recv_val[comm] = NULL;
				RegionMap[r].recv_buffer_size[comm] = coord_count*(this->dim + 1) + Dom_count * 2 + 2;//numcoordinates*(2(coords)+1(id)) + number of domains *2 + region id + number domains

				if (RegionMap[r].recv_buffer[comm] != NULL)//need to initialize recv_buffer to NULl
					free(RegionMap[r].recv_buffer[comm]);

				if (RegionMap[r].recv_buffer_size[comm] >2)
				{
					RegionMap[r].recv_buffer[comm] = new int[RegionMap[r].recv_buffer_size[comm]];
					RegionMap[r].recv_buffer[comm][0] = myid;
					RegionMap[r].recv_buffer[comm][1] = Dom_count;
					rcv_ind = 2;
					for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
					{
						ind = index(x, y, D[0]);
						if (DomainMap[ind].region_id == RegionMap[r].region_id)
						{
							if (DomainMap[ind].recv_num[comm] > 0)
							{
								RegionMap[r].recv_buffer[comm][rcv_ind] = ind;
								RegionMap[r].recv_buffer[comm][rcv_ind + 1] = DomainMap[ind].recv_num[comm];
								rcv_ind += 2;
								for (i = 0; i<DomainMap[ind].recv_num[comm]; i++)
								{
									RegionMap[r].recv_buffer[comm][rcv_ind] = DomainMap[ind].recv_coords[comm][3 * i];		//x
									RegionMap[r].recv_buffer[comm][rcv_ind + 1] = DomainMap[ind].recv_coords[comm][3 * i + 1];	//y
									RegionMap[r].recv_buffer[comm][rcv_ind + 2] = DomainMap[ind].recv_coords[comm][3 * i + 2];	//id
									rcv_ind += 3;
								}
							}
						}
					}
				}
				else
				{
					RegionMap[r].recv_buffer_size[comm] = 0;
					RegionMap[r].recv_buffer[comm] = NULL;
				}
			}
		}
		/*
		if (this->dim == 2)
		{
		for (r = 0; r<ncpus; r++)//loop over regions
		{
		coord_count = 0;
		Dom_count = 0;
		for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
		{

		ind = x + y*D[0];
		if (DomainMap[ind].region_id == RegionMap[r].region_id)
		{
		if (DomainMap[ind].recv_num[comm] > 0)
		{
		coord_count += DomainMap[ind].recv_num[comm];
		Dom_count++;//temporarily use as a domain count to allocate memory for recv buffer
		}
		}
		}
		RegionMap[r].recv_val_size[comm] = coord_count;
		RegionMap[r].recv_val[comm] = new double[coord_count];//allocate memory to receive incoming values
		RegionMap[r].recv_buffer_size[comm] = coord_count*(this->dim + 1) + Dom_count * 2 + 2;//numcoordinates*(2(coords)+1(id)) + number of domains *2 + region id + number domains

		if (RegionMap[r].recv_buffer[comm] != NULL)//need to initialize recv_buffer to NULl
		free(RegionMap[r].recv_buffer[comm]);

		if (RegionMap[r].recv_buffer_size[comm] >2)
		{
		RegionMap[r].recv_buffer[comm] = new int[RegionMap[r].recv_buffer_size[comm]];
		RegionMap[r].recv_buffer[comm][0] = myid;
		RegionMap[r].recv_buffer[comm][1] = Dom_count;
		rcv_ind = 2;
		for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
		{
		ind = x + y*D[0];
		if (DomainMap[ind].region_id == RegionMap[r].region_id)
		{
		if (DomainMap[ind].recv_num[comm] > 0)
		{
		RegionMap[r].recv_buffer[comm][rcv_ind] = ind;
		RegionMap[r].recv_buffer[comm][rcv_ind + 1] = DomainMap[ind].recv_num[comm];
		rcv_ind += 2;
		for (i = 0; i<DomainMap[ind].recv_num[comm]; i++)
		{
		RegionMap[r].recv_buffer[comm][rcv_ind] = DomainMap[ind].recv_coords[comm][3 * i];		//x
		RegionMap[r].recv_buffer[comm][rcv_ind + 1] = DomainMap[ind].recv_coords[comm][3 * i + 1];	//y
		RegionMap[r].recv_buffer[comm][rcv_ind + 2] = DomainMap[ind].recv_coords[comm][3 * i + 2];	//id
		rcv_ind += 3;
		}
		}
		}
		}
		}
		else
		{
		RegionMap[r].recv_buffer_size[comm] = 0;
		RegionMap[r].recv_buffer[comm] = NULL;
		}
		}
		}
		*/
		else if (this->dim ==3)
		{
			for (r = 0; r<ncpus; r++)//loop over regions
			{
				coord_count = 0;
				Dom_count = 0;
				for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
				{
					ind = index(x, y, z, D[0], D[1]);
					if (DomainMap[ind].region_id == RegionMap[r].region_id)
					{
						if (DomainMap[ind].recv_num[comm] > 0)
						{
							coord_count += DomainMap[ind].recv_num[comm];
							Dom_count++;//temporarily use as a domain count to allocate memory for recv buffer
						}
					}
				}
				RegionMap[r].recv_val_size[comm] = coord_count;
				if (RegionMap[r].recv_val[comm] != NULL)
				{
					delete RegionMap[r].recv_val[comm];
				}
				if (coord_count > 0)
				{
					RegionMap[r].recv_val[comm] = new double[coord_count];//allocate memory to receive incoming values
				}
				else
					RegionMap[r].recv_val[comm] = NULL;
				RegionMap[r].recv_buffer_size[comm] = coord_count*(this->dim + 1) + Dom_count * 2 + 2;//numcoordinates*(2(coords)+1(id)) + number of domains *2 + region id + number domains
				if (RegionMap[r].recv_buffer[comm] != NULL)//need to initialize recv_buffer to NULl
					free(RegionMap[r].recv_buffer[comm]);

				if (RegionMap[r].recv_buffer_size[comm] >2)
				{
					RegionMap[r].recv_buffer[comm] = new int[RegionMap[r].recv_buffer_size[comm]];
					RegionMap[r].recv_buffer[comm][0] = myid;
					RegionMap[r].recv_buffer[comm][1] = Dom_count;
					rcv_ind = 2;
					for (z = 0; z<D[2]; z++)for (y = 0; y<D[1]; y++)for (x = 0; x<D[0]; x++)
					{
						ind = index(x, y, z, D[0], D[1]);
						if (DomainMap[ind].region_id == RegionMap[r].region_id)
						{
							if (DomainMap[ind].recv_num[comm] > 0)
							{
								RegionMap[r].recv_buffer[comm][rcv_ind] = ind;
								RegionMap[r].recv_buffer[comm][rcv_ind + 1] = DomainMap[ind].recv_num[comm];
								rcv_ind += 2;
								for (i = 0; i<DomainMap[ind].recv_num[comm]; i++)
								{
									RegionMap[r].recv_buffer[comm][rcv_ind] = DomainMap[ind].recv_coords[comm][(this->dim + 1)*i];		//x
									RegionMap[r].recv_buffer[comm][rcv_ind + 1] = DomainMap[ind].recv_coords[comm][(this->dim + 1)*i + 1];	//y
									RegionMap[r].recv_buffer[comm][rcv_ind + 2] = DomainMap[ind].recv_coords[comm][(this->dim + 1)*i + 2];	//z
									RegionMap[r].recv_buffer[comm][rcv_ind + 3] = DomainMap[ind].recv_coords[comm][(this->dim + 1)*i + 3];	//id

									rcv_ind += 4;
								}
							}
						}
					}
				}
				else
				{
					RegionMap[r].recv_buffer_size[comm] = 0;
					RegionMap[r].recv_buffer[comm] = NULL;
				}

			}
		}
	}
	void Syncronize_Region_Receive_Requests(int comm)
	{
		MPI_Request R1, R2;
		MPI_Status S1, S2;
		int i, j = -1;
		int send_tag, receive_tag;


		int receive_from, send_to;
		for (i = 0; i<ncpus; i++)
		{		
			send_to = (myid - i + ncpus) % ncpus;
			receive_from = (i + myid) % ncpus;
			send_tag = ncpus*(send_to)+i;
			receive_tag = ncpus*(myid)+i;
			MPI_Irecv(&(RegionMap[receive_from].send_buffer_size[comm]), 1, MPI_INT, receive_from, receive_tag, MPI_COMM_WORLD, &R1);
			MPI_Send(&(RegionMap[send_to].recv_buffer_size[comm]), 1, MPI_INT, send_to, send_tag, MPI_COMM_WORLD);
			MPI_Wait(&R1, &S1);
		}
		for (i = 0; i < ncpus; i++)
		{
			if (RegionMap[i].send_buffer[comm] != NULL)
				delete RegionMap[i].send_buffer[comm];
			RegionMap[i].send_buffer[comm] = new int[RegionMap[i].send_buffer_size[comm]];
		}
		for (i = 0; i < ncpus; i++)
		{
			send_to = (myid - i + ncpus) % ncpus;
			receive_from = (i + myid) % ncpus;
			MPI_Barrier(MPI_COMM_WORLD);
			if (RegionMap[receive_from].send_buffer_size[comm] > 0)
				MPI_Irecv(RegionMap[receive_from].send_buffer[comm], RegionMap[receive_from].send_buffer_size[comm], MPI_INT, receive_from, myid, MPI_COMM_WORLD, &R2);
			if (RegionMap[send_to].recv_buffer_size[comm] > 0)
				MPI_Send(RegionMap[send_to].recv_buffer[comm], RegionMap[send_to].recv_buffer_size[comm], MPI_INT, send_to, send_to, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			if (RegionMap[receive_from].send_buffer_size[comm] > 0)
				MPI_Wait(&R2, &S2);
		}
	}
	void Create_Send_Value(int comm)
	{
		int i;
		for (i = 0; i<ncpus; i++)
		{
			if (RegionMap[i].send_buffer_size[comm] > 2)
			{
				Reg->Fill_Buffer(RegionMap[i].send_val[comm], RegionMap[i].send_val_size[comm], RegionMap[i].send_buffer[comm], RegionMap[i].send_buffer_size[comm], comm);
			}
		}

	}
	void Allocate_Memory_Send_Value(int comm)
	{
		int r;
		for (r = 0; r<ncpus; r++)
		{
			if (RegionMap[r].send_buffer_size[comm] > 2)
			{
				RegionMap[r].send_val_size[comm] = (RegionMap[r].send_buffer_size[comm] - 2 - 2 * RegionMap[r].send_buffer[comm][1]) / (this->dim + 1);//(d+1) because coords (x,y) + id (1)
				if (RegionMap[r].send_val[comm] != NULL)
				{
					delete RegionMap[r].send_val[comm];
					RegionMap[r].send_val[comm] = NULL;
				}
				RegionMap[r].send_val[comm] = new double[RegionMap[r].send_val_size[comm]];
			}
			else
			{
				RegionMap[r].send_val_size[comm] = 0;
			}
		}
	}
	// ****** OUTPUT FUNCTIONS ********/
	void Output_Load(int time)
	{
		FILE *fp;
		for (int i = 0; i < MESH_PARAMS.num_procs; i++)
		{
			if (this->myid == i)
			{
				fp = fopen("load_data.out", "a");
				Reg->Output_Load(time, fp,this->myid);
				fclose(fp);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	void Output_VTK(int time,int lowres)
	{
		this->Reg->Output_VTK(time, lowres);
	}
	void Output_Restart(int time)
	{
		int x, y, z, ind;
		char filename[BUFSIZ];
		char command[BUFSIZ];
		FILE *fp;
		#ifdef _WIN32
			sprintf(command, "mkdir restart_%d", time);
		#elif __linux
			sprintf(command, "mkdir ./restart_%d", time);
		#elif __unix
			sprintf(command, "mkdir ./restart_%d", time);
		#endif
		system(command);
		/*
		if (this->myid == 0)
		{
		  
			sprintf(filename, "./restart_%d/control.restart", time);
			fp = fopen(filename, "w");
			fprintf(fp, "%d	: Number of Communications\n",this->num_comms);
			fprintf(fp, "%d	: Number of Dimensions\n",this->dim);
			fprintf(fp, "%d	: Number of Processors\n",this->ncpus);
			fprintf(fp, "%d %d %d	: Number of domains in x,y,z\n",this->D[0],this->D[1],this->D[2]);
			fprintf(fp, "#DomainMap Restart (x,y,z,size, region id, domain filename\n");
			for (z = 0; z < D[2]; z++)for (y = 0; y < D[1]; y++)for (x = 0; x < D[0]; x++)
			{
				ind = index(x, y, z, D[0], D[1]);
				fprintf(fp, "%d %d %d	:	Domain Origin\n", this->DomainMap[ind].origin.Vx, this->DomainMap[ind].origin.Vy, this->DomainMap[ind].origin.Vz);
				fprintf(fp, "%d	:	size (number of grid points)\n", this->DomainMap[ind].size);
				fprintf(fp, "%d	:	region id\n", this->DomainMap[ind].region_id);
				fprintf(fp, "domain_%d_%d_%d.dat\n", this->DomainMap[ind].origin.Vx, this->DomainMap[ind].origin.Vy, this->DomainMap[ind].origin.Vz);
			}
			fclose(fp);
		}
		  */
		#ifdef _WIN32
			sprintf(command, "mkdir restart_%d\\region_%d", time,this->myid);
		#elif __linux
			sprintf(command, "mkdir ./restart_%d/region_%d", time,this->myid);
		#elif __unix
			sprintf(command, "mkdir ./restart_%d/region_%d", time,this->myid);
		#endif
		  
		system(command);
		//		this->Reg->Output_Restart(time);
	}
};
#endif
