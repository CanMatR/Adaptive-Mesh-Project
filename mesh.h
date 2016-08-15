#ifndef mesh_paramsDEF
#define mesh_paramsDEF

class mesh_params
{
public:
	int num_procs;
	int DIMENSION;
	cart_vector Num_Dom;
	int maxRes;
	int uniRes;
	int buffer_size;
	int num_fields;
	int num_comms;
	int comm_calls;
	cart_vector global_size;
	cart_vector periodic;
	mesh_params()
	{
		DIMENSION = 3;
		maxRes = 1;
		uniRes = 1;
		buffer_size = 1;
		num_fields = 1;
		num_comms = 1;
		comm_calls = 1;
	}
	~mesh_params()
	{}
	void output_params()
	{
		printf("num_procs  %d\n",this->num_procs);
		printf("Dim %d\n",DIMENSION);
		printf("maxRes %d\n",maxRes);
		printf("uniRes%d\n",uniRes);
		printf("buffer_size %d\n",buffer_size);
		printf("num_fields %d\n",num_fields);
		printf("num_comms %d\n",num_comms);
		printf("Num_Dom %d %d %d\n",Num_Dom.Vx, Num_Dom.Vy,Num_Dom.Vz);
		printf("global_size %d %d %d\n",global_size.Vx, global_size.Vy,global_size.Vz);
		printf("periodic %d %d %d\n",periodic.Vx,periodic.Vy,periodic.Vz);
	}
};

mesh_params MESH_PARAMS;
#endif