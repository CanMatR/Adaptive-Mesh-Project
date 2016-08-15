/* Copyright (c) 2015 Government of Canada  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "utility.h"
#include "phaseFieldAirBinaryModule.h"
#include <mpi.h>
#include "mesh.h"

int num_procs = 1;
int myid = 0;
int main(int argc, char **argv)
{
	int err;
	err = MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	
		phaseFieldAirBinaryModule *module;
		module = new phaseFieldAirBinaryModule(myid, num_procs);
		MESH_PARAMS.num_procs = num_procs;
		module->Create_Initialization();
//		module->Load_Restart(1000);
		fflush(stdout);
		module->Time_Loop(0);
		MPI_Barrier(MPI_COMM_WORLD);

		printf("Finished run\n");
		fflush(stdout);
	fflush(stdout);
	err = MPI_Finalize();

	return 0;
}
