#ifndef phaseFieldAirBinaryModuleDEF
#define phaseFieldAirBinaryModuleDEF

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <list>
#include <stack>
#include <cmath>

#include "field.h"
#include "control.h"

using namespace std;
class simul_params
{
public:
	//model parameters
	double a1, a2, at;
	double obs;	//solid-solid obstable 
	//physical materials parameters
	double c_o, m0, ke;//wt% alloy composition, liquidus slope, partition coefficient
	double Tm, Tl, Ts, DT_o; //Pure Melting Temp, Liquidus Temp, solidus Temp, the freezing range as calculated per the alloy parameters
	double W, E4, Gamma, d_o;	//interface width, anisotropy, Gibb-Thomson coefficient, capillary length
	
	//physical kinetic parameters
	double Dl, Qe, Epsilon_0;// liquid diffusion coefficient, activation barrier for solid diffusion, Ratio of solid-liquid diffusion
	//Noise Parameters
	double Fu;
	double Fu_dt;			//coefficient preceeding uncorrelated noise
	//simulation parameters
	double lbTemp;			//Lower bound temperature, is it used?
	double Q_t;				//cooling rate
	double binary_omega;	//initial undercooling
	double T_init;			//initial temperature
	//dimensionless variables
	double Lambda, tau, Ddimless;//scaling parameter, dimensionless time, dimensionless diffusion

	double dx0, dt;			//mesh spacing and timestep

	//constants
	double Pi, Rg;//Pi and gas constant

	//rotations
	double psi_zs[4];
	double psi_ys[4];
	double czs[4];
	double szs[4];
	double cys[4];
	double sys[4];
	double sx, sy, cx, cy;

	int nsym;

	//Adaption parameters
	double AdaptThreshold;
	double *AdaptWeights;
	
	// time parameters
	int printFreq, AdaptFreq, endt,endt_mult;
	
	// 0 = phase field, 1 = original U field, 2 = flux-calculated U field, where indC0 refers to
	int nf;  //number of fields
	int ng;			//number of grains
	int indP0;		// order parameter 0 index
	int indC1;		// first concentration index
	int ind_dPdt;	// index for dPdt

	//simulate flags
	bool Q_simple;	// instead of 
	

	simul_params()
	{
		//dimensionless variables
		a1 = 0.8839;
		a2 = 0.6267;
		at = 0.35355;
	}
	~simul_params()
	{}
	void Load_Simul_Params(char *filename)
	{
		char inputline[BUFSIZ];
		FILE *fp;
		fp = fopen(filename, "r");
		//constants
		obs = 154.147;
		Rg = 8.3144621;						// Gas Constant in J /( mol K )
		Pi = 3.14159265;

		double alpha, beta;


		//load parameters
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf %lf %lf", &(this->c_o),&(this->m0),&(this->ke));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf", &(Tm));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf %lf %lf", &(this->W), &(this->E4), &(this->Gamma));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf %lf %lf", &(this->Dl), &(this->Qe), &(this->Epsilon_0));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf", &(this->Fu));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf %lf %lf %lf", &(this->T_init), &(this->binary_omega), &(this->Q_t), &(this->lbTemp));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf %lf", &(this->dx0), &(this->dt));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%d %d %d", &(this->printFreq), &(this->AdaptFreq), &(this->endt_mult));
		fgets(inputline, BUFSIZ, fp);
		sscanf(inputline, "%lf %lf", &alpha, &beta);
		fclose(fp);

		this->sx = sin(alpha);
		this->cx = cos(alpha);
		this->sy = sin(beta);
		this->cy = cos(beta);
		//will add to input file when we decide how we want to do it
		AdaptThreshold = 1E-4;
		nf = 5;
		AdaptWeights = new double[nf]; //0.001; // adaption threshold for splitting, in last function of field.h
		for (int i = 0; i < nf; i++)
			AdaptWeights[i] = 0.001;
		Q_simple = false; // instead of 

		//Calculate the rest from the input
		endt = endt_mult * printFreq;

		Tl = Tm - m0*c_o;		//liquidus temperature 
		Ts = Tm - m0*c_o / ke;  //solidus temperature
		DT_o = m0*c_o*(1 - ke); //the freezing range as calculated per the alloy parameters
		d_o = Gamma / DT_o;     //the capillary as described by the material parameters

		//conditions for simulations
		Lambda = a1*W / d_o;    //lambda, inverse of the nucleation barrier and convergence parameter
		Ddimless = a2*Lambda;   //dimensionless diffusion coefficient, in units of tau and W
		if (dt == 0.0)	
			dt = 0.8*dx0*dx0 / (6.*Ddimless); // dimensionless D (in units of tau and W. dt should be given in taus)
		if (T_init == 0.0)
			T_init = Tl;

		tau = a2*Lambda*W*W / Dl;			//time constant, scales atomic attachment kinetics
		
		nsym = 0;
		ng =1;
		indP0 = 0;// concentration in index(x,y,z,indC,N,N,N)
		// 0 = phase field, 1 = original U field, 2 = flux-calculated U field, where indC0 refers to
		indC1 = ng;//4; phase-fields in index(C,x,y,z,indP0...indP0+ng,N,N,N)
		ind_dPdt = indC1 + 1; //indP0+indC0+3+1; // currently it's at index 5 (6th field) 



		Fu_dt = Fu / dt;      //coefficient preceeding uncorrelated noise

		//Rotations
		psi_zs[0] = Pi / 4.;
		psi_zs[1] = 0.;
		psi_zs[2] = Pi/8.;
		psi_zs[3] = Pi / 8.;
		
		psi_ys[0] = Pi / 4.; 
		psi_ys[1] = 0.; 
		psi_ys[2] = Pi / 4.; 
		psi_ys[3] = Pi / 8.;

		czs[0] = cos(psi_zs[0]);
		czs[1] = cos(psi_zs[1]);
		czs[2] = cos(psi_zs[2]);
		czs[3] = cos(psi_zs[3]);

		szs[0] = sin(psi_zs[0]);
		szs[1] = sin(psi_zs[1]);
		szs[2] = sin(psi_zs[2]);
		szs[3] = sin(psi_zs[3]);

		cys[0] = cos(psi_ys[0]);
		cys[1] = cos(psi_ys[1]);
		cys[2] = cos(psi_ys[2]);
		cys[3] = cos(psi_ys[3]);

		sys[0] = sin(psi_ys[0]);
		sys[1] = sin(psi_ys[1]);
		sys[2] = sin(psi_ys[2]);
		sys[3] = sin(psi_ys[3]);


	}
};


class phaseFieldAirBinaryModule
{
	list <field *> *fields;
	list <field *>::iterator itField;
	simul_params prms;
public:
	control *Con;
	int myid;
	int ncpus;
	int flag;

	phaseFieldAirBinaryModule()
	{
		flag = 0;
		Con = NULL;
		fields = new list < field * > ;
	}
	phaseFieldAirBinaryModule(int myid, int ncpus)
	{
		Con = new control(myid, ncpus);
		this->myid = myid;
		this->ncpus = ncpus;
		fields = new list < field * >;
	}
	~phaseFieldAirBinaryModule()
	{
		if (Con != NULL)
			delete Con;
	}
	void Load_Simulation_Parameters(char *filename)
	{
		prms.Load_Simul_Params(filename);
	}
	void Load_Restart(int time)
	{
		Con->Load_Restart(time);
	}
	void Create_Initialization()
	{
	  char filename[BUFSIZ];
	  sprintf(filename,"Simul_Params.dat");
	  this->Load_Simulation_Parameters(filename);
	  //		this->Load_Simulation_Parameters("Simul_Params.dat");
		Con->Create_Initialization();
		Con->Get_Field_List(fields);
		Initialize_Fields();
		
		Con->Output_Restart(0);
		Con->Output_VTK(0,0);
	}
	void Initialize_Fields()
	{
		list <field*>::iterator it;
		it = fields->begin();
		while (it != fields->end())
		{
			Init(*it);
			it++;
		}
	}
	void Init(field *fld)
	{

		double pos, min_pos;
		int nf = fld->nf;
		int B = fld->B;
		int i;
		int x, y, z;
		int dx = fld->dx;
		int ind, im;
		double **data;
		data = fld->data;
		int xo, yo, zo;
		int xp, yp, zp;
		xo = fld->origin.Vx;
		yo = fld->origin.Vy;
		zo = fld->origin.Vz;

		double eu;
		eu = 1 - (1 - prms.ke)*prms.binary_omega;
		double c,c_eq;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;
		
		double R;
		double *xc, *yc,*zc,*Rad;
		if (this->Con->dim == 2)
		{
			FILE *fp;
			int numflds;
			char inputline[BUFSIZ];
			fp = fopen("init_geometry.dat", "r");
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//
			sscanf(inputline, "%d", &(numflds));
			fgets(inputline, BUFSIZ, fp);//num_comms
			for (i = 0; i < numflds; i++)
			{
				fgets(inputline, BUFSIZ, fp);
			}
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			int numseeds;
			sscanf(inputline, "%d", &(numseeds));
			xc = new double[numseeds];
			yc = new double[numseeds];
			zc = new double[numseeds];
			Rad = new double[numseeds];
			for (i = 0; i < numseeds; i++)
			{
				fgets(inputline, BUFSIZ, fp);//num_comms
				fgets(inputline, BUFSIZ, fp);//num_comms
				sscanf(inputline, "%lf %lf %lf %lf", &(Rad[i]), &(xc[i]), &(yc[i]), &(zc[i]));
				Rad[i] = Rad[i] * 0.75;
			}
			fclose(fp);
			for (y = 0; y < fld->grid_size.Vy + 2 * B; y++)for (x = 0; x < fld->grid_size.Vx + 2 * B; x++)
			{
				xp = (x - B)*dx + xo;
				yp = (y - B)*dx + yo;
				ind = index(x, y, fld->grid_size.Vx + 2 * B);
				min_pos = 10000000.0;
				R = 0.0;
				for (i = 0; i < numseeds; i++)
				{
					pos = sqrt((xp - xc[i])*(xp - xc[i]) + (yp - yc[i])*(yp - yc[i]));
					if (pos < min_pos)
					{
						min_pos = pos;
						R = Rad[i];
					}
				}
				data[0][ind] = -tanh((min_pos - R) / sqrt(2));
				data[1][ind] = 0.5*(1 + data[0][ind])*prms.ke*eu + 0.5*(1 - data[0][ind])*eu + (1 - eu)*exp(-(min_pos - (R + 1.))*(min_pos - (R + 1.))/20.);
				data[2][ind] = 0.0;
			

			}
		}
		else
		{
			FILE *fp;
			int numflds;
			char inputline[BUFSIZ];
			fp = fopen("init_geometry.dat", "r");
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//
			sscanf(inputline, "%d", &(numflds));
			fgets(inputline, BUFSIZ, fp);//num_comms
			for (i = 0; i < numflds; i++)
			{
				fgets(inputline, BUFSIZ, fp);
			}
			fgets(inputline, BUFSIZ, fp);//num_comms
			fgets(inputline, BUFSIZ, fp);//num_comms
			int numseeds;
			sscanf(inputline, "%d", &(numseeds));
			xc = new double[numseeds];
			yc = new double[numseeds];
			zc = new double[numseeds];
			Rad = new double[numseeds];
			for (i = 0; i < numseeds; i++)
			{
				fgets(inputline, BUFSIZ, fp);//num_comms
				fgets(inputline, BUFSIZ, fp);//num_comms
				sscanf(inputline, "%lf %lf %lf %lf", &(Rad[i]), &(xc[i]), &(yc[i]), &(zc[i]));
				Rad[i] = Rad[i] * 0.75;
			}
			fclose(fp);
			for (z = 0; z < fld->grid_size.Vz + 2 * B; z++)for (y = 0; y < fld->grid_size.Vy + 2 * B; y++)for (x = 0; x < fld->grid_size.Vx + 2 * B; x++)
			{
				xp = (x - B)*dx + xo;
				yp = (y - B)*dx + yo;
				zp = (z - B)*dx + zo;
				ind = index(x, y, z, fld->grid_size.Vx + 2 * B, fld->grid_size.Vy + 2 * B);
				min_pos = 10000000.0;
				R = 0.0;
				for (i = 0; i < numseeds; i++)
				{
					pos = sqrt((xp - xc[i])*(xp - xc[i]) + (yp - yc[i])*(yp - yc[i]) + (zp - zc[i])*(zp - zc[i]));
					if (pos < min_pos)
					{
						min_pos = pos;
						R = Rad[i];
					}
				}
				data[0][ind] = -tanh((min_pos - R) / sqrt(2));
				data[1][ind] = 0.5*(1 + data[0][ind])*prms.ke*eu + 0.5*(1 - data[0][ind])*eu + (1 - eu)*exp(-(min_pos - (R + 1.))*(min_pos - (R + 1.)) / 20.);
				data[2][ind] = 0.0;
			}

		}
		delete xc;
		delete yc;
		delete zc;
		delete Rad;
	}
	void Time_Loop(int time)
	{
		double *dPdt, *dUdt,*dCdt;
		double T;
		double xMax_interface; 
		int i;
		list <field*>::iterator it;
		Con->Get_Field_List(this->fields);
		this->Initialize_Fields();
		Con->Build_Communication_Requests();
		
		for (i = 1; i<MESH_PARAMS.comm_calls; i++)
		{
			Con->Communicate(i);
		}
		Con->Adaption_Cycle(1);
		Con->Build_Communication_Requests();
		Con->Get_Field_List(this->fields);
		Con->Output_Restart(1);
//		Con->Output_VTK(1, 0);
//		Con->Output_VTK(1, 1);
		Con->Output_VTK(1, 1);
//		Con->Output_VTK(1, 3);
		
		printf("Begin Time Loop\n");
		it = this->fields->begin();
		int B = (*it)->B;
		int N = (*it)->grid_size.Vx*3;
		
		int arraysize = (N + 2 * B)*(N + 2 * B)*(N + 2 * B);
		dPdt = new double[arraysize];
		dUdt = new double[arraysize];
		dCdt = new double[arraysize];

		int t;
		
		for (t = time + 1; t <= prms.endt; t++)
		{
			T = prms.T_init - prms.Q_t*t*prms.dt*prms.tau;
			for (i = 1; i < MESH_PARAMS.comm_calls; i++)
			{
				Con->Communicate(i);
			}
			if (t%prms.printFreq == 0)
			{
				printf("%d\n", t);
				fflush(stdout);
				Con->Output_Restart(t);
				Con->Output_VTK(t, 1);
			}
			if (t % prms.AdaptFreq == 0)
			{

				Con->Adaption_Cycle(t);
				Con->Get_Field_List(this->fields);
			}
			else
			{

			}
			Con->Communicate(1); // communicate phase field and concentration
			it = this->fields->begin();
			while (it != this->fields->end())
			{
				if (this->Con->dim == 3)
					calc_binaryP_3D(*it, dPdt, T); // calculate dPdt, C fluxes for pseudo finite volume															// update Phi if proper finite volume formulation is not used
				if (this->Con->dim == 2)
					calc_binaryP_2D(*it, dPdt, T); // calculate dPdt, C fluxes for pseudo finite volume
				it++;
			}
			
			Con->Communicate(2); // communicate dPdt, Cfluxes, 
			
			it = this->fields->begin();
			while (it != this->fields->end())
			{
				//xMax_interface = find_interface(*it, xMax_interface); 
				if (this->Con->dim == 3)
				{
					calc_binaryC_finiteVolume_3D(*it, dCdt, T); // 
					updateP_3D(*it); // update phi after upadting C, if proper finite volume formulation is used 
				}
				if (this->Con->dim == 2)
				{
					calc_binaryC_finiteVolume_2D(*it, dCdt, T); // 
					updateP_2D(*it); // update phi after upadting C, if proper finite volume formulation is used 
				}
				it++;
			}			
		}
	}
	double find_interface(field *fld, double xMax_interface)
	{
		int B = fld->B;
		double **data;
		data = fld->data;
		int x, y, z;
		int ind,xp;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;
		double xMax_interface0 = xMax_interface; // node with the largest x coordinate, which is at the interface
		for (z = B; z < Nz + B; z++) for (y = B; y < Ny + B; y++) for (x = B; x < Nx + B; x++)
		{
			ind = index(x,     y, z, Nx + 2 * B, Ny + 2 * B);
			//xp  = index(x + 1, y, z, Nx + 2 * B, Ny + 2 * B);
			xp  = index(x, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			if ( data[0][ind]*data[0][xp] < 0 ) // phi(x,y,z) and phi(x+1,y,z) have a different sign --> interface must be between these two nodes
			{
				if (y > xMax_interface0)
					xMax_interface0 = y; 
			}
		}
		return xMax_interface0; 
	}

	void calc_binaryP_3D(field *fld, double *dPdt, double T)
	{
		double **data;
		data = fld->data;
		int B = fld->B;
		//	int N = fld->N;
		int x, y, z;
		int xp, xm, yp, ym, zp, zm, ind;
		int xypp, xypm, xymp, xymm, xzpp, xzpm, xzmp, xzmm, yzpp, yzpm, yzmp, yzmm;
		//seb start
		double lap, omps;
		double a_n, da_n, px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz, px2, py2, pz2, xnorm, xnorm2, xnorm3;
		double px4, px6, py4, py6, pz4, pz6, px2y2, px4y2, px2z2, px4z2, py2z2, py4x2, py4z2, pz4x2, pz4y2, px2y2z2;
		double pxpy, pxpz, pypz, termx, termy, termz, t1x, t1y, t1z, t2x, t2y, t2z, t1xy, t1xz, t1yz, t2xy, t2xz, t2yz;

		double dx = fld->dx*prms.dx0;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;
		//seb end

		int indP = prms.indP0;
		int indC1 = prms.indC1;
		int ind_dPdt = prms.ind_dPdt;

		double sumphi, sumphi_xp, sumphi_xm, sumphi_yp, sumphi_ym, sumphi_zp, sumphi_zm;
		double chemTerm, repulsion, eU, q, ceq, c;

		double DsPerDl = prms.Epsilon_0*exp(-prms.Qe / (prms.Rg*T));

		double t3, A, AyAx;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;
		//Boundary Condition zero flux on Z
		if (fld->origin.Vz == 0)
		{
			int indp;
			z = 0;
			for (y = 0; y < Ny + 2 * B; y++)for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
				indp = index(x, y, z + 2, Nx + 2 * B, Ny + 2 * B);
				fld->data[0][ind] = fld->data[0][indp];
			}
		}
		if (fld->origin.Vz + (fld->grid_size.Vz - 1)*fld->dx == MESH_PARAMS.global_size.Vz - 1)
		{
			int indm;
			z = fld->grid_size.Vz + fld->B * 2 - 1;
			for (y = 0; y < Ny + 2 * B; y++)for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
				indm = index(x, y, z - 2, Nx + 2 * B, Ny + 2 * B);
				fld->data[0][ind] = fld->data[0][indm];
			}
		}
		for (z = B; z<Nz + B; z++) for (y = B; y<Ny + B; y++) for (x = B; x<Nx + B; x++)
		{
			//phi simulation
			ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
			xp = index(x + 1, y, z, Nx + 2 * B, Ny + 2 * B);
			xm = index(x - 1, y, z, Nx + 2 * B, Ny + 2 * B);
			yp = index(x, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			ym = index(x, y - 1, z, Nx + 2 * B, Ny + 2 * B);
			zp = index(x, y, z + 1, Nx + 2 * B, Ny + 2 * B);
			zm = index(x, y, z - 1, Nx + 2 * B, Ny + 2 * B);

			sumphi = data[0][ind];
			sumphi_xp = data[0][xp];
			sumphi_xm = data[0][xm];
			sumphi_yp = data[0][yp];
			sumphi_ym = data[0][ym];
			sumphi_zp = data[0][zp];
			sumphi_zm = data[0][zm];

			xypp = index(x + 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			xypm = index(x + 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
			xymp = index(x - 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			xymm = index(x - 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
			xzpp = index(x + 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
			xzpm = index(x + 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
			xzmp = index(x - 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
			xzmm = index(x - 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
			yzpp = index(x, y + 1, z + 1, Nx + 2 * B, Ny + 2 * B);
			yzpm = index(x, y + 1, z - 1, Nx + 2 * B, Ny + 2 * B);
			yzmp = index(x, y - 1, z + 1, Nx + 2 * B, Ny + 2 * B);
			yzmm = index(x, y - 1, z - 1, Nx + 2 * B, Ny + 2 * B);

			// anisotropy calc:
			// px,py,pz and xnorm are used also in antitrapping current
			px = (data[0][xp] - data[0][xm]) / dx2;
			py = (data[0][yp] - data[0][ym]) / dx2;
			pz = (data[0][zp] - data[0][zm]) / dx2;

			px2 = px*px;
			py2 = py*py;
			pz2 = pz*pz;
			xnorm = px2 + py2 + pz2;

			a_n = 1.;
			da_n = 0.;
//			calcAnisotropy_3D(fld, x, y, z, &a_n, &da_n); // calculate a_n and da_n 
			calcRot_Anisotropy_3D(fld, x, y, z, &a_n, &da_n);
			// end anisotropy calc

			lap = data[0][ym] + data[0][yp]
				+ data[0][xp] + data[0][xm]
				+ data[0][zp] + data[0][zm]
				- 6 * data[0][ind];

			omps = 1. - data[0][ind] * data[0][ind];

			ceq = 0.5*(1. + prms.ke - (1. - prms.ke)*sumphi);

			//c = data[indC0][ind];  
			c = data[indC1][ind];
			eU = c / ceq; // if c > ceq, eU > 0
			// solid-liquid diffusion coefficient interpolation function
			repulsion = 0.0;
			chemTerm = -prms.Lambda / (1. - prms.ke)*(eU - 1. + (T - prms.Tl) / (prms.m0*prms.c_o))*omps*omps;

			dPdt[ind] = lap / dxs + (omps*data[indP][ind] + chemTerm + da_n + repulsion) / (a_n*a_n);
			data[ind_dPdt][ind] = dPdt[ind];
		}
	}
	void calc_binaryC_finiteVolume_3D(field *fld, double *dCdt, double T)
	{
		double **data;
		data = fld->data;
		int B = fld->B;
		//	int N = fld->N;
		int x, y, z;
		int xp, xm, yp, ym, zp, zm, ind;
		int xypp, xypm, xymp, xymm, xzpp, xzpm, xzmp, xzmm, yzpp, yzpm, yzmp, yzmm;

		double dx = fld->dx*prms.dx0;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;
		//seb end

		int indP = prms.indP0;
		int indC1 = prms.indC1;
		int ind_dPdt = prms.ind_dPdt;

		double sumphi, sumphi_xp, sumphi_xm, sumphi_yp, sumphi_ym, sumphi_zp, sumphi_zm;
		double q, ceq, c;

		double DsPerDl = prms.Epsilon_0*exp(-prms.Qe / (prms.Rg*T));

		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;

		double Jright, Jleft, Jtop, Jbot, Jfront, Jback;
		double PhiAve, cAve, c1, ceq1, Q, eu, eu1, du1, eu_dpdt, dPhi_dx, dPhi_dy, dPhi_dz, normal, dPhi_abs_sq;
		int ind1;
		double c_air, c_air1;
		//Boundary Condition zero flux on Z
		if (fld->origin.Vz == 0)
		{
			int indp;
			z = 0;
			for (y = 0; y < Ny + 2 * B; y++)for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
				indp = index(x, y, z + 2, Nx + 2 * B, Ny + 2 * B);
				fld->data[1][ind] = fld->data[1][indp];
			}
		}
		if (fld->origin.Vz + (fld->grid_size.Vz-1)*fld->dx == MESH_PARAMS.global_size.Vz-1)
		{
			int indm;
			z = fld->grid_size.Vz + fld->B * 2 - 1;
			for (y = 0; y < Ny + 2 * B; y++)for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
				indm = index(x, y, z - 2, Nx + 2 * B, Ny + 2 * B);
				fld->data[1][ind] = fld->data[1][indm];
			}
		}
		for (z = B; z<Nz + B; z++) for (y = B; y<Ny + B; y++) for (x = B; x<Nx + B; x++)
		{
			//phi simulation
			ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
			xp = index(x + 1, y, z, Nx + 2 * B, Ny + 2 * B);
			xm = index(x - 1, y, z, Nx + 2 * B, Ny + 2 * B);
			yp = index(x, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			ym = index(x, y - 1, z, Nx + 2 * B, Ny + 2 * B);
			zp = index(x, y, z + 1, Nx + 2 * B, Ny + 2 * B);
			zm = index(x, y, z - 1, Nx + 2 * B, Ny + 2 * B);

			xypp = index(x + 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			xypm = index(x + 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
			xymp = index(x - 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
			xymm = index(x - 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
			xzpp = index(x + 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
			xzpm = index(x + 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
			xzmp = index(x - 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
			xzmm = index(x - 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
			yzpp = index(x, y + 1, z + 1, Nx + 2 * B, Ny + 2 * B);
			yzpm = index(x, y + 1, z - 1, Nx + 2 * B, Ny + 2 * B);
			yzmp = index(x, y - 1, z + 1, Nx + 2 * B, Ny + 2 * B);
			yzmm = index(x, y - 1, z - 1, Nx + 2 * B, Ny + 2 * B);

			Jright = 0;
			Jleft = 0;
			Jtop = 0;
			Jbot = 0;
			Jfront = 0;
			Jback = 0;

			sumphi = data[0][ind];

			ceq = 0.5*(1. + prms.ke - (1. - prms.ke)*sumphi);
			c = data[indC1][ind];
			eu = c / ceq;

			// the "1" in variable name means that the quantity is evaluated at node on the other side of the finite volume boundary, 

			// all values should be evaluated directly at this position. so for any functions of phi or c, one should take the average of phi or c.
			// (not average of the function). In other words, f( average of phi), not average of f(phi)
			// position (i+1/2, j, k)
			ind1 = xp;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve); // NO SOLID DIFFUSION 
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			dPhi_dx = (data[0][ind1] - data[0][ind]) / dx;
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i+1,j-1,k)
			dPhi_dy = (data[0][yp] - data[0][ym] + data[0][xypp] - data[0][xypm]) / (2 * 2 * dx);
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i+1,j,k-1)
			dPhi_dz = (data[0][zp] - data[0][zm] + data[0][xzpp] - data[0][xzpm]) / (2 * 2 * dx);

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz;

			normal = 0.;
			if (dPhi_abs_sq > 1e-6)
				normal = dPhi_dx / sqrt(dPhi_abs_sq);

			Jright = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// position (i-1/2, j, k)
			ind1 = xm;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			dPhi_dx = (data[0][ind1] - data[0][ind]) / dx;
			// average of y-directional derivatives of Phi at (i,j+/-1,k) and (i+1,j+/-1,k)
			dPhi_dy = (data[0][yp] - data[0][ym] + data[0][xymp] - data[0][xymm]) / (2 * 2 * dx);
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i-1,j,k-1)
			dPhi_dz = (data[0][zp] - data[0][zm] + data[0][xzmp] - data[0][xzmm]) / (2 * 2 * dx);

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz;

			normal = 0.;
			if (dPhi_abs_sq > 1e-6)
				normal = dPhi_dx / sqrt(dPhi_abs_sq);

			Jleft = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// position (i, j+1/2, k)
			ind1 = yp;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			//dPhi_dx = ( data[0][ind1]-data[0][ind] )/dx; 
			dPhi_dx = (data[0][xp] - data[0][xm] + data[0][xypp] - data[0][xymp]) / (2 * 2 * dx);
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i-1,j-1,k)
			dPhi_dy = (data[0][ind1] - data[0][ind]) / dx;
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i-1,j,k-1)
			dPhi_dz = (data[0][zp] - data[0][zm] + data[0][yzpp] - data[0][yzpm]) / (2 * 2 * dx);

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz;

			normal = 0.;
			if (dPhi_abs_sq > 1e-6)
				normal = dPhi_dy / sqrt(dPhi_abs_sq);

			Jtop = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// position (i, j-1/2, k)
			ind1 = ym;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			//dPhi_dx = ( data[0][ind1]-data[0][ind] )/dx; 
			dPhi_dx = (data[0][xp] - data[0][xm] + data[0][xypm] - data[0][xymm]) / (2 * 2 * dx);
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i-1,j-1,k)
			dPhi_dy = (data[0][ind1] - data[0][ind]) / dx;
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i-1,j,k-1)
			dPhi_dz = (data[0][zp] - data[0][zm] + data[0][yzmp] - data[0][yzmm]) / (2 * 2 * dx);

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz;

			normal = 0.;
			if (dPhi_abs_sq > 1e-6)
				normal = dPhi_dy / sqrt(dPhi_abs_sq);

			Jbot = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;


			// position (i, j, k+1/2)
			ind1 = zp;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			//dPhi_dx = ( data[0][ind1]-data[0][ind] )/dx; 
			dPhi_dx = (data[0][xp] - data[0][xm] + data[0][xzpp] - data[0][xzmp]) / (2 * 2 * dx);
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i-1,j-1,k)
			dPhi_dy = (data[0][yp] - data[0][ym] + data[0][yzpp] - data[0][yzmp]) / (2 * 2 * dx);
			//dPhi_dz = ( data[0][zp]-data[0][zm] + data[0][yzpp]-data[0][yzpm] )/(2*2*dx); 
			dPhi_dz = (data[0][ind1] - data[0][ind]) / dx;

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz;

			normal = 0.;
			if (dPhi_abs_sq > 1e-6)
				normal = dPhi_dz / sqrt(dPhi_abs_sq);

			Jfront = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;


			// position (i, j, k-1/2)
			ind1 = zm;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			//dPhi_dx = ( data[0][ind1]-data[0][ind] )/dx; 
			dPhi_dx = (data[0][xp] - data[0][xm] + data[0][xzpm] - data[0][xzmm]) / (2 * 2 * dx);
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i-1,j-1,k)
			dPhi_dy = (data[0][yp] - data[0][ym] + data[0][yzpm] - data[0][yzmm]) / (2 * 2 * dx);
			//dPhi_dz = ( data[0][zp]-data[0][zm] + data[0][yzpp]-data[0][yzpm] )/(2*2*dx); 
			dPhi_dz = (data[0][ind1] - data[0][ind]) / dx;

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz;

			normal = 0.;
			if (dPhi_abs_sq > 1e-6)
				normal = dPhi_dz / sqrt(dPhi_abs_sq);

			Jback = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// the left, bot, and back fluxes are already pointing away from the finite volume center,
			// so they have the same sign as right, top, front
			dCdt[ind] = -(Jright + Jleft + Jtop + Jbot + Jfront + Jback) / dx; // Eq. 6.80 in the book  
		}

		// should the update be done after all the concentrations have been looped through?
		// no: the other domains (*fields) do not know about each other before communication is done
		for (z = B; z<Nz + B; z++)for (y = B; y<Ny + B; y++)for (x = B; x<Nx + B; x++)
		{
			ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
			data[indC1][ind] += prms.dt*dCdt[ind];
		}

	}
	void updateP_3D(field *fld)
	{
		double dx = fld->dx*prms.dx0;
		int B = fld->B;
		double **data;
		data = fld->data;
		int x, y, z, ind;
		double dx2 = 2.*dx;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;
		for (z = B; z < Nz + B; z++) for (y = B; y < Ny + B; y++) for (x = B; x < Nx + B; x++)
		{
			ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);

			data[0][ind] += prms.dt*data[prms.ind_dPdt][ind];
		}
	}
	void calcRot_Anisotropy_3D(field *fld, int x, int y, int z, double *a_n, double *da_n)
	{
		//double a_n, da_n;
		double **data = fld->data;
		double dx = fld->dx*prms.dx0;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;
		int B = fld->B;

		int xp, xm, yp, ym, zp, zm, ind;
		double px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz, px2, py2, pz2, xnorm, xnorm2, xnorm3;
		double px4, px6, py4, py6, pz4, pz6, px2y2, px4y2, px2z2, px4z2, py2z2, py4x2, py4z2, pz4x2, pz4y2, px2y2z2;
		double pxpy, pxpz, pypz, termx, termy, termz, t1x, t1y, t1z, t2x, t2y, t2z, t1xy, t1xz, t1yz, t2xy, t2xz, t2yz;
		int xypp, xypm, xymp, xymm, xzpp, xzpm, xzmp, xzmm, yzpp, yzpm, yzmp, yzmm;

		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;

		//phi simulation
		ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
		xp = index(x + 1, y, z, Nx + 2 * B, Ny + 2 * B);
		xm = index(x - 1, y, z, Nx + 2 * B, Ny + 2 * B);
		yp = index(x, y + 1, z, Nx + 2 * B, Ny + 2 * B);
		ym = index(x, y - 1, z, Nx + 2 * B, Ny + 2 * B);
		zp = index(x, y, z + 1, Nx + 2 * B, Ny + 2 * B);
		zm = index(x, y, z - 1, Nx + 2 * B, Ny + 2 * B);

		//seb start
		xypp = index(x + 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
		xypm = index(x + 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
		xymp = index(x - 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
		xymm = index(x - 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
		xzpp = index(x + 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
		xzpm = index(x + 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
		xzmp = index(x - 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
		xzmm = index(x - 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
		yzpp = index(x, y + 1, z + 1, Nx + 2 * B, Ny + 2 * B);
		yzpm = index(x, y + 1, z - 1, Nx + 2 * B, Ny + 2 * B);
		yzmp = index(x, y - 1, z + 1, Nx + 2 * B, Ny + 2 * B);
		yzmm = index(x, y - 1, z - 1, Nx + 2 * B, Ny + 2 * B);

		px = (data[0][xp] - data[0][xm]) / dx2;
		py = (data[0][yp] - data[0][ym]) / dx2;
		pz = (data[0][zp] - data[0][zm]) / dx2;

		px2 = px*px;
		py2 = py*py;
		pz2 = pz*pz;
		xnorm = px2 + py2 + pz2;

		if (xnorm<1.e-6) {
			*a_n = 1;
			*da_n = 0;
			return;
		}

		pxx = (data[0][xp] - 2.*data[0][ind] + data[0][xm]) / dxs;
		pyy = (data[0][yp] - 2.*data[0][ind] + data[0][ym]) / dxs;
		pzz = (data[0][zp] - 2.*data[0][ind] + data[0][zm]) / dxs;

		pxy = (data[0][xypp] - data[0][xymp] - data[0][xypm] + data[0][xymm]) / dxs4;
		pxz = (data[0][xzpp] - data[0][xzmp] - data[0][xzpm] + data[0][xzmm]) / dxs4;
		pyz = (data[0][yzpp] - data[0][yzmp] - data[0][yzpm] + data[0][yzmm]) / dxs4;

		double pxx_c, pyy_c, pzz_c, pxy_c, pxz_c, pyz_c;
		double cx = prms.cx;
		double cy = prms.cy;
		double sx = prms.sx;
		double sy = prms.sy;
		double pn = sqrt(xnorm);
		double ipn = 1.0 / pn;
		double nx_r = (px*cy + pz*sy)*ipn;
		double ny_r = (px*sx*sy - pz*cy*sx + py*cx)*ipn;
		double nz_r = (-px*cx*sy + pz*cx*cy + py*sx)*ipn;
		double nr_4 = nx_r*nx_r*nx_r*nx_r + ny_r*ny_r*ny_r*ny_r + nz_r*nz_r*nz_r*nz_r;

		pxx_c = 0.384e3 * prms.E4 * (((((nz_r * nz_r * cx * cx / 0.2e1 + sx * sx * ny_r * ny_r / 0.2e1) * nr_4 + 0.2e1 / 0.3e1 * pow(nz_r, 0.6e1) * cx * cx + 0.2e1 / 0.3e1 * sx * sx * pow(ny_r, 0.6e1) - 0.3e1 / 0.8e1 * nz_r * nz_r * cx * cx - 0.3e1 / 0.8e1 * sx * sx * ny_r * ny_r - 0.4e1 / 0.3e1 * pow(nz_r, 0.3e1) * cx * sx * pow(ny_r, 0.3e1)) * sy * sy - 0.4e1 / 0.3e1 * cy * pow(nx_r, 0.3e1) * (cx * pow(nz_r, 0.3e1) - pow(ny_r, 0.3e1) * sx) * sy - nr_4 * nr_4 / 0.6e1 + (cy * cy * nx_r * nx_r / 0.2e1 + 0.1e1 / 0.8e1) * nr_4 + 0.2e1 / 0.3e1 * cy * cy * pow(nx_r, 0.6e1) - 0.3e1 / 0.8e1 * cy * cy * nx_r * nx_r) * prms.E4 + (nz_r * nz_r * cx * cx / 0.8e1 + sx * sx * ny_r * ny_r / 0.8e1) * sy * sy - nr_4 / 0.24e2 + cy * cy * nx_r * nx_r / 0.8e1) * pn * pn + 0.2e1 * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * px * ((cx * pow(nz_r, 0.3e1) - pow(ny_r, 0.3e1) * sx) * sy - cy * pow(nx_r, 0.3e1)) * pn + (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * nr_4 * px * px) * pow(pn, -0.2e1);
		pyy_c = 0.384e3 * prms.E4 * (((-nr_4 * nr_4 / 0.6e1 + (cx * cx * ny_r * ny_r / 0.2e1 + sx * sx * nz_r * nz_r / 0.2e1 + 0.1e1 / 0.8e1) * nr_4 + 0.2e1 / 0.3e1 * cx * cx * pow(ny_r, 0.6e1) + 0.4e1 / 0.3e1 * pow(nz_r, 0.3e1) * cx * sx * pow(ny_r, 0.3e1) + 0.2e1 / 0.3e1 * sx * sx * pow(nz_r, 0.6e1) - 0.3e1 / 0.8e1 * cx * cx * ny_r * ny_r - 0.3e1 / 0.8e1 * sx * sx * nz_r * nz_r) * prms.E4 + cx * cx * ny_r * ny_r / 0.8e1 + sx * sx * nz_r * nz_r / 0.8e1 - nr_4 / 0.24e2) * pn * pn - 0.2e1 * py * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * (cx * pow(ny_r, 0.3e1) + pow(nz_r, 0.3e1) * sx) * pn + py * py * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * nr_4) * pow(pn, -0.2e1);
		pzz_c = 0.384e3 * (((((nz_r * nz_r * cx * cx / 0.2e1 + sx * sx * ny_r * ny_r / 0.2e1) * nr_4 + 0.2e1 / 0.3e1 * pow(nz_r, 0.6e1) * cx * cx + 0.2e1 / 0.3e1 * sx * sx * pow(ny_r, 0.6e1) - 0.3e1 / 0.8e1 * nz_r * nz_r * cx * cx - 0.3e1 / 0.8e1 * sx * sx * ny_r * ny_r - 0.4e1 / 0.3e1 * pow(nz_r, 0.3e1) * cx * sx * pow(ny_r, 0.3e1)) * cy * cy + 0.4e1 / 0.3e1 * cy * pow(nx_r, 0.3e1) * (cx * pow(nz_r, 0.3e1) - pow(ny_r, 0.3e1) * sx) * sy - nr_4 * nr_4 / 0.6e1 + (sy * sy * nx_r * nx_r / 0.2e1 + 0.1e1 / 0.8e1) * nr_4 + 0.2e1 / 0.3e1 * sy * sy * pow(nx_r, 0.6e1) - 0.3e1 / 0.8e1 * sy * sy * nx_r * nx_r) * prms.E4 + (nz_r * nz_r * cx * cx / 0.8e1 + sx * sx * ny_r * ny_r / 0.8e1) * cy * cy - nr_4 / 0.24e2 + sy * sy * nx_r * nx_r / 0.8e1) * pn * pn - 0.2e1 * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * ((cx * pow(nz_r, 0.3e1) - pow(ny_r, 0.3e1) * sx) * cy + sy * pow(nx_r, 0.3e1)) * pz * pn + (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * nr_4 * pz * pz) * prms.E4 * pow(pn, -0.2e1);
		pxy_c = -0.768e3 * prms.E4 * (((0.2e1 / 0.3e1 * pow(ny_r, 0.3e1) * pow(nz_r, 0.3e1) * cx * cx * sy + (-0.2e1 / 0.3e1 * (pow(ny_r, 0.4e1) + ny_r * ny_r * nz_r * nz_r + pow(nz_r, 0.4e1) + 0.3e1 / 0.4e1 * nr_4 - 0.9e1 / 0.16e2) * (ny_r - nz_r) * (ny_r + nz_r) * sy * sx - 0.2e1 / 0.3e1 * pow(ny_r, 0.3e1) * cy * pow(nx_r, 0.3e1)) * cx - 0.2e1 / 0.3e1 * pow(nz_r, 0.3e1) * sx * (pow(ny_r, 0.3e1) * sx * sy + cy * pow(nx_r, 0.3e1))) * prms.E4 - cx * sx * sy * (ny_r - nz_r) * (ny_r + nz_r) / 0.8e1) * pn * pn + ((-py * pow(nz_r, 0.3e1) * sy + px * pow(ny_r, 0.3e1)) * cx + (py * pow(ny_r, 0.3e1) * sy + px * pow(nz_r, 0.3e1)) * sx + py * cy * pow(nx_r, 0.3e1)) * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * pn - nr_4 * py * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * px) * pow(pn, -0.2e1);
		pxz_c = 0.768e3 * prms.E4 * (((-0.2e1 / 0.3e1 * pow(nx_r, 0.3e1) * (-cx * pow(nz_r, 0.3e1) + pow(ny_r, 0.3e1) * sx) * cy * cy - 0.2e1 / 0.3e1 * sy * (-pow(nx_r, 0.6e1) + (-0.3e1 / 0.4e1 * nr_4 + 0.9e1 / 0.16e2) * nx_r * nx_r + sx * sx * pow(ny_r, 0.6e1) - 0.2e1 * pow(nz_r, 0.3e1) * cx * sx * pow(ny_r, 0.3e1) + 0.3e1 / 0.4e1 * sx * sx * (-0.3e1 / 0.4e1 + nr_4) * ny_r * ny_r + nz_r * nz_r * cx * cx * (pow(nz_r, 0.4e1) + 0.3e1 / 0.4e1 * nr_4 - 0.9e1 / 0.16e2)) * cy + 0.2e1 / 0.3e1 * sy * sy * pow(nx_r, 0.3e1) * (-cx * pow(nz_r, 0.3e1) + pow(ny_r, 0.3e1) * sx)) * prms.E4 - cy * sy * (nz_r * nz_r * cx * cx + sx * sx * ny_r * ny_r - nx_r * nx_r) / 0.8e1) * pn * pn + ((-pz * pow(nx_r, 0.3e1) + px * (-cx * pow(nz_r, 0.3e1) + pow(ny_r, 0.3e1) * sx)) * cy - sy * (px * pow(nx_r, 0.3e1) + pz * (-cx * pow(nz_r, 0.3e1) + pow(ny_r, 0.3e1) * sx))) * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * pn + nr_4 * pz * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * px) * pow(pn, -0.2e1);
		pyz_c = 0.768e3 * (((0.2e1 / 0.3e1 * pow(ny_r, 0.3e1) * pow(nz_r, 0.3e1) * cx * cx * cy + (-0.2e1 / 0.3e1 * (pow(ny_r, 0.4e1) + ny_r * ny_r * nz_r * nz_r + pow(nz_r, 0.4e1) + 0.3e1 / 0.4e1 * nr_4 - 0.9e1 / 0.16e2) * sx * (ny_r - nz_r) * (ny_r + nz_r) * cy + 0.2e1 / 0.3e1 * pow(ny_r, 0.3e1) * sy * pow(nx_r, 0.3e1)) * cx - 0.2e1 / 0.3e1 * pow(nz_r, 0.3e1) * sx * (cy * pow(ny_r, 0.3e1) * sx - sy * pow(nx_r, 0.3e1))) * prms.E4 - cx * cy * sx * (ny_r - nz_r) * (ny_r + nz_r) / 0.8e1) * pn * pn + ((-py * cy * pow(nz_r, 0.3e1) - pz * pow(ny_r, 0.3e1)) * cx + py * pow(ny_r, 0.3e1) * cy * sx - py * sy * pow(nx_r, 0.3e1) - pow(nz_r, 0.3e1) * pz * sx) * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * pn + nr_4 * pz * (0.1e1 / 0.12e2 + (nr_4 - 0.1e1 / 0.4e1) * prms.E4) * py) * prms.E4 * pow(pn, -0.2e1);

		*a_n = 1. - 3.*prms.E4 + 4.*prms.E4 *(nr_4);
		*da_n = pxx_c*pxx + pyy_c*pyy + pzz_c*pzz + pxy_c*pxy + pxz_c*pxz + pyz_c*pyz;

		//Old anisotropy
		/*
		pxpy = px*py;
		pxpz = px*pz;
		pypz = py*pz;

		px4 = px2*px2;
		px6 = px2*px4;
		py4 = py2*py2;
		py6 = py2*py4;
		pz4 = pz2*pz2;
		pz6 = pz2*pz4;

		px2y2 = px2*py2;
		px4y2 = px4*py2;
		px2z2 = px2*pz2;
		px4z2 = px4*pz2;
		py2z2 = py2*pz2;
		py4x2 = py4*px2;
		py4z2 = py4*pz2;
		pz4x2 = pz4*px2;
		pz4y2 = pz4*py2;
		px2y2z2 = px2*py2*pz2;

		termx = px2y2 - py4 - pz4 + px2z2;
		termy = px2y2 - px4 - pz4 + py2z2;
		termz = py2z2 - px4 - py4 + px2z2;

		t1x = pxx*(px4y2 + 4.*py4x2 - py6 + px4z2 + 6.*px2y2z2 - pz4y2 - py4z2 + 4.*pz4x2 - pz6);
		t1y = pyy*(py4z2 + 4.*pz4y2 - pz6 + py4x2 + 6.*px2y2z2 - px4z2 - pz4x2 + 4.*px4y2 - px6);
		t1z = pzz*(pz4x2 + 4.*px4z2 - px6 + pz4y2 + 6.*px2y2z2 - py4x2 - px4y2 + 4.*py4z2 - py6);

		t2x = pxx*px2*termx*termx;
		t2y = pyy*py2*termy*termy;
		t2z = pzz*pz2*termz*termz;

		t1xy = pxy*pxpy*(-2.*px2y2 - px2z2 - py2z2 + pz4);
		t1xz = pxz*pxpz*(-2.*px2z2 - px2y2 - py2z2 + py4);
		t1yz = pyz*pypz*(-2.*py2z2 - px2y2 - px2z2 + px4);

		t2xy = 2.*pxy*pxpy*termx*termy;
		t2xz = 2.*pxz*pxpz*termx*termz;
		t2yz = 2.*pyz*pypz*termy*termz;

		xnorm2 = xnorm*xnorm;
		xnorm3 = xnorm2*xnorm;

		*a_n = 1. - 3.*prms.E4 + 4.*prms.E4 *(px4 + py4 + pz4) / xnorm2;
		*da_n = 16.*prms.E4*(*a_n)*(t1x + t1y + t1z) / xnorm3 + 64.*prms.E4*(*a_n) *(t1xy + t1xz + t1yz) / xnorm3 + 256.*prms.E4*prms.E4*((t2x + t2y + t2z) / xnorm2 + (t2xy + t2xz + t2yz) / xnorm2) / xnorm3;
		*/
		return;
	}
	void calcAnisotropy_3D(field *fld, int x, int y, int z, double *a_n, double *da_n)
	{
		//double a_n, da_n;
		double **data = fld->data;
		double dx = fld->dx*prms.dx0;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;
		int B = fld->B;

		int xp, xm, yp, ym, zp, zm, ind;
		double px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz, px2, py2, pz2, xnorm, xnorm2, xnorm3;
		double px4, px6, py4, py6, pz4, pz6, px2y2, px4y2, px2z2, px4z2, py2z2, py4x2, py4z2, pz4x2, pz4y2, px2y2z2;
		double pxpy, pxpz, pypz, termx, termy, termz, t1x, t1y, t1z, t2x, t2y, t2z, t1xy, t1xz, t1yz, t2xy, t2xz, t2yz;
		int xypp, xypm, xymp, xymm, xzpp, xzpm, xzmp, xzmm, yzpp, yzpm, yzmp, yzmm;

		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		int Nz = fld->grid_size.Vz;

		//phi simulation
		ind = index(x, y, z, Nx + 2 * B, Ny + 2 * B);
		xp = index(x + 1, y, z, Nx + 2 * B, Ny + 2 * B);
		xm = index(x - 1, y, z, Nx + 2 * B, Ny + 2 * B);
		yp = index(x, y + 1, z, Nx + 2 * B, Ny + 2 * B);
		ym = index(x, y - 1, z, Nx + 2 * B, Ny + 2 * B);
		zp = index(x, y, z + 1, Nx + 2 * B, Ny + 2 * B);
		zm = index(x, y, z - 1, Nx + 2 * B, Ny + 2 * B);

		//seb start
		xypp = index(x + 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
		xypm = index(x + 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
		xymp = index(x - 1, y + 1, z, Nx + 2 * B, Ny + 2 * B);
		xymm = index(x - 1, y - 1, z, Nx + 2 * B, Ny + 2 * B);
		xzpp = index(x + 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
		xzpm = index(x + 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
		xzmp = index(x - 1, y, z + 1, Nx + 2 * B, Ny + 2 * B);
		xzmm = index(x - 1, y, z - 1, Nx + 2 * B, Ny + 2 * B);
		yzpp = index(x, y + 1, z + 1, Nx + 2 * B, Ny + 2 * B);
		yzpm = index(x, y + 1, z - 1, Nx + 2 * B, Ny + 2 * B);
		yzmp = index(x, y - 1, z + 1, Nx + 2 * B, Ny + 2 * B);
		yzmm = index(x, y - 1, z - 1, Nx + 2 * B, Ny + 2 * B);

		px = (data[0][xp] - data[0][xm]) / dx2;
		py = (data[0][yp] - data[0][ym]) / dx2;
		pz = (data[0][zp] - data[0][zm]) / dx2;

		px2 = px*px;
		py2 = py*py;
		pz2 = pz*pz;
		xnorm = px2 + py2 + pz2;

		if (xnorm<1.e-6){
			*a_n = 1;
			*da_n = 0;
			return;
		}

		pxx = (data[0][xp] - 2.*data[0][ind] + data[0][xm]) / dxs;
		pyy = (data[0][yp] - 2.*data[0][ind] + data[0][ym]) / dxs;
		pzz = (data[0][zp] - 2.*data[0][ind] + data[0][zm]) / dxs;

		pxy = (data[0][xypp] - data[0][xymp] - data[0][xypm] + data[0][xymm]) / dxs4;
		pxz = (data[0][xzpp] - data[0][xzmp] - data[0][xzpm] + data[0][xzmm]) / dxs4;
		pyz = (data[0][yzpp] - data[0][yzmp] - data[0][yzpm] + data[0][yzmm]) / dxs4;

		pxpy = px*py;
		pxpz = px*pz;
		pypz = py*pz;

		px4 = px2*px2;
		px6 = px2*px4;
		py4 = py2*py2;
		py6 = py2*py4;
		pz4 = pz2*pz2;
		pz6 = pz2*pz4;

		px2y2 = px2*py2;
		px4y2 = px4*py2;
		px2z2 = px2*pz2;
		px4z2 = px4*pz2;
		py2z2 = py2*pz2;
		py4x2 = py4*px2;
		py4z2 = py4*pz2;
		pz4x2 = pz4*px2;
		pz4y2 = pz4*py2;
		px2y2z2 = px2*py2*pz2;

		termx = px2y2 - py4 - pz4 + px2z2;
		termy = px2y2 - px4 - pz4 + py2z2;
		termz = py2z2 - px4 - py4 + px2z2;

		t1x = pxx*(px4y2 + 4.*py4x2 - py6 + px4z2 + 6.*px2y2z2 - pz4y2 - py4z2 + 4.*pz4x2 - pz6);
		t1y = pyy*(py4z2 + 4.*pz4y2 - pz6 + py4x2 + 6.*px2y2z2 - px4z2 - pz4x2 + 4.*px4y2 - px6);
		t1z = pzz*(pz4x2 + 4.*px4z2 - px6 + pz4y2 + 6.*px2y2z2 - py4x2 - px4y2 + 4.*py4z2 - py6);

		t2x = pxx*px2*termx*termx;
		t2y = pyy*py2*termy*termy;
		t2z = pzz*pz2*termz*termz;

		t1xy = pxy*pxpy*(-2.*px2y2 - px2z2 - py2z2 + pz4);
		t1xz = pxz*pxpz*(-2.*px2z2 - px2y2 - py2z2 + py4);
		t1yz = pyz*pypz*(-2.*py2z2 - px2y2 - px2z2 + px4);

		t2xy = 2.*pxy*pxpy*termx*termy;
		t2xz = 2.*pxz*pxpz*termx*termz;
		t2yz = 2.*pyz*pypz*termy*termz;

		xnorm2 = xnorm*xnorm;
		xnorm3 = xnorm2*xnorm;

		*a_n = 1. - 3.*prms.E4 + 4.*prms.E4 *(px4 + py4 + pz4) / xnorm2;
		*da_n = 16.*prms.E4*(*a_n)*(t1x + t1y + t1z) / xnorm3 + 64.*prms.E4*(*a_n) *(t1xy + t1xz + t1yz) / xnorm3 + 256.*prms.E4*prms.E4*((t2x + t2y + t2z) / xnorm2 + (t2xy + t2xz + t2yz) / xnorm2) / xnorm3;

		return;
	}

	void calc_binaryP_2D(field *fld, double *dPdt, double T)
	{
		double **data;
		data = fld->data;
		int B = fld->B;

		double xo = fld->origin.Vx;
		double yo = fld->origin.Vy;
		double zo = fld->origin.Vz;
		int x, y;
		int xp, xm, yp, ym, zp, zm, ind;
		int xypp, xypm, xymp, xymm, xzpp, xzpm, xzmp, xzmm, yzpp, yzpm, yzmp, yzmm;
		//seb start
		double lap, omps;
		double a_n, da_n, px, py, pxx, pyy, pxy, px2, py2, xnorm, xnorm2, xnorm3;
		double px4, px6, py4, py6, px2y2, px4y2, py4x2;
		double pxpy, pxpz, pypz, termx, termy, termz, t1x, t1y, t2x, t2y, t1xy, t2xy;

		double dx = fld->dx*prms.dx0;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;
		//seb end

		int indP = prms.indP0;
		int indC1 = prms.indC1;
		int ind_dPdt = prms.ind_dPdt;

		double sumphi, sumphi_xp, sumphi_xm, sumphi_yp, sumphi_ym; //, sumphi_zp, sumphi_zm;
		double chemTerm, repulsion, eu, q, c_eq, c;

		double DsPerDl = prms.Epsilon_0*exp(-prms.Qe / (prms.Rg*T));

		double t3, A, AyAx;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		//int Nz = fld->grid_size.Vz;

		double airRepulsion = 0;
		double c_air; // air concentration 
		double dist_x, dist_y, dist;


		//Set Boundary condition
		if (fld->origin.Vy == 0)
		{
			int indp;
			y = 0;
			for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, Nx + 2 * B);
				indp = index(x, y+2, Nx + 2 * B);
				fld->data[0][ind] = fld->data[0][indp];
			}
		}
		if (fld->origin.Vy + fld->grid_size.Vy*fld->dx == MESH_PARAMS.global_size.Vy)
		{
			int indm;
			y = fld->grid_size.Vy + fld->B * 2 - 1;
			for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, Nx + 2 * B);
				indm = index(x, y-2, Nx + 2 * B);
				fld->data[0][ind] = fld->data[0][indm];
			}
		}

		//for (z = B; z<Nz +B; z++) for (y = B; y<Ny +B; y++) for (x = B; x<Nx +B; x++)
		for (y = B; y<Ny + B; y++) for (x = B; x<Nx + B; x++)
		{
			//phi simulation
			ind = index(x, y, Nx + 2 * B);
			xp = index(x + 1, y, Nx + 2 * B);
			xm = index(x - 1, y, Nx + 2 * B);
			yp = index(x, y + 1, Nx + 2 * B);
			ym = index(x, y - 1, Nx + 2 * B);

			sumphi = data[0][ind];

			xypp = index(x + 1, y + 1, Nx + 2 * B);
			xypm = index(x + 1, y - 1, Nx + 2 * B);
			xymp = index(x - 1, y + 1, Nx + 2 * B);
			xymm = index(x - 1, y - 1, Nx + 2 * B);

			// anisotropy calc:
			// px,py,pz and xnorm are used also in antitrapping current
			px = (data[0][xp] - data[0][xm]) / dx2;
			py = (data[0][yp] - data[0][ym]) / dx2;

			px2 = px*px;
			py2 = py*py;
			xnorm = px2 + py2; // + pz2;

			a_n = 1.;
			da_n = 0.;
			calcAnisotropy_2D(fld, x, y, &a_n, &da_n); // calculate a_n and da_n 

			lap = data[0][ym] + data[0][yp]
				+ data[0][xp] + data[0][xm]
				- 4 * data[0][ind];

			omps = 1. - data[0][ind] * data[0][ind];

			c_eq = 0.5*(1. + prms.ke - (1. - prms.ke)*sumphi);

			c = data[indC1][ind];
			eu = c / c_eq; // if c > ceq, eU > 0
			// solid-liquid diffusion coefficient interpolation function
			repulsion = 0.0;
			chemTerm = -prms.Lambda / (1. - prms.ke)*(eu - 1. + (T - prms.Tl) / (prms.m0*prms.c_o))*omps*omps;

			dPdt[ind] = lap / dxs + (omps*data[indP][ind] + chemTerm + da_n + repulsion + airRepulsion) / (a_n*a_n);
			data[ind_dPdt][ind] = dPdt[ind];
		}
	}
	void calc_binaryC_finiteVolume_2D(field *fld, double *dCdt, double T)
	{
		double **data;
		data = fld->data;
		int B = fld->B;
		//	int N = fld->N;
		int x, y; //, z;
		int xp, xm, yp, ym, ind;
		int xypp, xypm, xymp, xymm; //xzpp, xzpm, xzmp, xzmm, yzpp, yzpm, yzmp, yzmm;

		double dx = fld->dx*prms.dx0;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;
		//seb end

		int indP = prms.indP0;
		int indC1 = prms.indC1;
		int ind_dPdt = prms.ind_dPdt;

		double sumphi, sumphi_xp, sumphi_xm, sumphi_yp, sumphi_ym;
		//double chemTerm,repulsion,eU,q,ceq,c; 
		double q, ceq, c;

		double DsPerDl = prms.Epsilon_0*exp(-prms.Qe / (prms.Rg*T));

		//double t3, A, AyAx;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		//int Nz = fld->grid_size.Vz;

		double Jright, Jleft, Jtop, Jbot;
		double PhiAve, cAve, c1, ceq1, Q, eu, eu1, du1, eu_dpdt, dPhi_dx, dPhi_dy, normal, dPhi_abs_sq;
		int ind1;
		double Phi_abs_sq_threshold = 1e-8; // threshold below which normal is just assumed to be zero

		double u_air, u_air1, c_air, c_air1; // air's contribution to dimensionless chemical potential
		// air's contribution to chemical potential: d/dc 0.5*obs_air_C*c_air^2*c^2 = obs_air_C*c_air^2*c
		//for (z = B; z<Nz +B; z++) for (y = B; y<Ny +B; y++) for (x = B; x<Nx +B; x++)
		// 
		//Set Boundary condition
		if (fld->origin.Vy == 0)
		{
			int indp;
			y = 0;
			for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, Nx + 2 * B);
				indp = index(x, y + 2, Nx + 2 * B);
				fld->data[1][ind] = fld->data[1][indp];
			}
		}
		if (fld->origin.Vy + fld->grid_size.Vy*fld->dx == MESH_PARAMS.global_size.Vy)
		{
			int indm;
			y = fld->grid_size.Vy + fld->B * 2 - 1;
			for (x = 0; x < Nx + 2 * B; x++)
			{
				ind = index(x, y, Nx + 2 * B);
				indm = index(x, y - 2, Nx + 2 * B);
				fld->data[1][ind] = fld->data[1][indm];
			}
		}

		for (y = B; y<Ny + B; y++) for (x = B; x<Nx + B; x++)
		{
			//phi simulation
			ind = index(x, y, Nx + 2 * B);
			xp = index(x + 1, y, Nx + 2 * B);
			xm = index(x - 1, y, Nx + 2 * B);
			yp = index(x, y + 1, Nx + 2 * B);
			ym = index(x, y - 1, Nx + 2 * B);

			xypp = index(x + 1, y + 1, Nx + 2 * B);
			xypm = index(x + 1, y - 1, Nx + 2 * B);
			xymp = index(x - 1, y + 1, Nx + 2 * B);
			xymm = index(x - 1, y - 1, Nx + 2 * B);


			// u = molVol/kB*T * (mu-mu_inf) = ln(c/ceq) is non-dimensional chemical potential
			// if c > ceq --> u = ln(c/ceq) > 0 --> negative driving force dPdt ~ -lambda/(1-k) * (e^u - 1) 


			Jright = 0;
			Jleft = 0;
			Jtop = 0;
			Jbot = 0;
			//		Jfront = 0;
			//		Jback  = 0;

			sumphi = data[0][ind];

			ceq = 0.5*(1. + prms.ke - (1. - prms.ke)*sumphi);
			c = data[indC1][ind];
			eu = c / ceq;

			// the "1" in variable name means that the quantity is evaluated at node on the other side of the finite volume boundary, 

			// all values should be evaluated directly at this position. so for any functions of phi or c, one should take the average of phi or c.
			// (not average of the function). In other words, f( average of phi), not average of f(phi)
			// position (i+1/2, j, k)
			ind1 = xp;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);

//			Q = calcQ(PhiAve);
			if (prms.Q_simple)
				Q= 0.5*(1 - PhiAve);
			else
				Q= (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);
			
			
			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			dPhi_dx = (data[0][ind1] - data[0][ind]) / dx;
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i+1,j-1,k)
			dPhi_dy = (data[0][yp] - data[0][ym] + data[0][xypp] - data[0][xypm]) / (2 * 2 * dx);
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i+1,j,k-1)
			//dPhi_dz = ( data[0][zp]-data[0][zm] + data[0][xzpp]-data[0][xzpm] )/(2*2*dx); 

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy; // + dPhi_dz*dPhi_dz; 

			normal = 0.;
			if (dPhi_abs_sq > Phi_abs_sq_threshold)
				normal = dPhi_dx / sqrt(dPhi_abs_sq);

			Jright = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// position (i-1/2, j, k)
			ind1 = xm;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);

//			Q = calcQ(PhiAve);
			if (prms.Q_simple)
				Q = 0.5*(1 - PhiAve);
			else
				Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);

			eu1 = c1 / ceq1;

			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			dPhi_dx = (data[0][ind1] - data[0][ind]) / dx;
			// average of y-directional derivatives of Phi at (i,j+/-1,k) and (i+1,j+/-1,k)
			dPhi_dy = (data[0][yp] - data[0][ym] + data[0][xymp] - data[0][xymm]) / (2 * 2 * dx);
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i-1,j,k-1)
			//dPhi_dz = ( data[0][zp]-data[0][zm] + data[0][xzmp]-data[0][xzmm] )/(2*2*dx); 

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy; // + dPhi_dz*dPhi_dz; 

			normal = 0.;
			if (dPhi_abs_sq > Phi_abs_sq_threshold)
				normal = dPhi_dx / sqrt(dPhi_abs_sq);

			Jleft = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// position (i, j+1/2, k)
			ind1 = yp;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
//			Q = calcQ(PhiAve);
			if (prms.Q_simple)
				Q = 0.5*(1 - PhiAve);
			else
				Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);

			//Q = (1-PhiAve)/( 1+ke - (1-ke)*PhiAve ); 
			eu1 = c1 / ceq1;


			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary
			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			//dPhi_dx = ( data[0][ind1]-data[0][ind] )/dx; 
			dPhi_dx = (data[0][xp] - data[0][xm] + data[0][xypp] - data[0][xymp]) / (2 * 2 * dx);
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i-1,j-1,k)
			dPhi_dy = (data[0][ind1] - data[0][ind]) / dx;
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i-1,j,k-1)
			//dPhi_dz = ( data[0][zp]-data[0][zm] + data[0][yzpp]-data[0][yzpm] )/(2*2*dx); 

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy; // + dPhi_dz*dPhi_dz; 
			//dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz; 

			normal = 0.;
			if (dPhi_abs_sq > Phi_abs_sq_threshold)
				normal = dPhi_dy / sqrt(dPhi_abs_sq);

			Jtop = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// position (i, j-1/2, k)
			ind1 = ym;
			c1 = data[indC1][ind1];
			ceq1 = 0.5*(1 + prms.ke - (1 - prms.ke)*data[0][ind1]);
			PhiAve = 0.5*(data[0][ind] + data[0][ind1]);
			cAve = 0.5*(c + c1);
			//Q = (1-PhiAve)/( 1+ke - (1-ke)*PhiAve ); 
//			Q = calcQ(PhiAve);
			if (prms.Q_simple)
				Q = 0.5*(1 - PhiAve);
			else
				Q = (1 - PhiAve) / (1 + prms.ke - (1 - prms.ke)*PhiAve);

			eu1 = c1 / ceq1;
			// d/dx u = d/dx ln(e^u) = 1/e^u * d/dx e^u. Left-side multiplier is average e^u at finite volume boundary


			du1 = 2 * (eu1 - eu) / (eu1 + eu) / dx;
			eu_dpdt = 0.5*(eu1 + eu)*0.5*(data[ind_dPdt][ind1] + data[ind_dPdt][ind]); // dPdt NEEDS TO BE LOOPED THROUGH AND COMMUNICATED 
			// IN THE WHOLE SYSTEM BEFORE THIS 
			//dPhi_dx = ( data[0][ind1]-data[0][ind] )/dx; 
			dPhi_dx = (data[0][xp] - data[0][xm] + data[0][xypm] - data[0][xymm]) / (2 * 2 * dx);
			// average of y-directional derivatives of Phi at (i,j+1,k) and (i-1,j-1,k)
			dPhi_dy = (data[0][ind1] - data[0][ind]) / dx;
			// average of z-directional derivatives of Phi at (i,j,k+1) and (i-1,j,k-1)
			//dPhi_dz = ( data[0][zp]-data[0][zm] + data[0][yzmp]-data[0][yzmm] )/(2*2*dx); 

			dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy; // + dPhi_dz*dPhi_dz; 
			//dPhi_abs_sq = dPhi_dx*dPhi_dx + dPhi_dy*dPhi_dy + dPhi_dz*dPhi_dz; 

			normal = 0.;
			if (dPhi_abs_sq > Phi_abs_sq_threshold)
				normal = dPhi_dy / sqrt(dPhi_abs_sq);

			Jbot = -prms.Ddimless*Q*cAve*du1 - prms.at*(1 - prms.ke)*eu_dpdt*normal;

			// the left, bot, and back fluxes are already pointing away from the finite volume center,
			// so they have the same sign as right, top, front
			dCdt[ind] = -(Jright + Jleft + Jtop + Jbot) / dx; // Eq. 6.80 in the book  
		}

		// should the update be done after all the concentrations have been looped through?
		// no: the other domains (*fields) do not know about each other before communication is done
		for (y = B; y<Ny + B; y++)for (x = B; x<Nx + B; x++)
		{
			ind = index(x, y, Nx + 2 * B);
			data[indC1][ind] += prms.dt*dCdt[ind];
		}

	}
	void updateP_2D(field *fld)
	{
		double dx = fld->dx*prms.dx0;
		int B = fld->B;
		double **data;
		data = fld->data;
		int x, y, ind;
		double dx2 = 2.*dx;
		int ind_dPdt = prms.ind_dPdt;
		int Nx = fld->grid_size.Vx;
		int Ny = fld->grid_size.Vy;
		for (y = B; y < Ny + B; y++) for (x = B; x < Nx + B; x++)
		{
			ind = index(x, y, Nx + 2 * B);
			data[0][ind] += prms.dt*data[ind_dPdt][ind];
		}
	}
	void calcAnisotropy_2D(field *fld, int x, int y, double *a_n, double *da_n)
	{
		double **data;
		data = fld->data;
		int B = fld->B;
		int xp, xm, yp, ym, ind;
		int xypp, xypm, xymp, xymm;

		int Nx = fld->grid_size.Vx;

		double dx = fld->dx*prms.dx0;
		double px, py, pz, pxx, pyy, pxy, px2, py2, xnorm, xnorm2, xnorm3;
		double px4, px6, py4, py6, px2y2, px4y2, py4x2;
		double pxpy, pxpz, termx, termy, t1x, t1y, t2x, t2y, t1xy, t2xy;
		double dx2 = 2.*dx;
		double dxs = dx*dx;
		double dxs4 = 4.*dx*dx;

		//phi simulation
		ind = index(x, y, Nx + 2 * B);
		xp = index(x + 1, y, Nx + 2 * B);
		xm = index(x - 1, y, Nx + 2 * B);
		yp = index(x, y + 1, Nx + 2 * B);
		ym = index(x, y - 1, Nx + 2 * B);
		xypp = index(x + 1, y + 1, Nx + 2 * B);
		xypm = index(x + 1, y - 1, Nx + 2 * B);
		xymm = index(x - 1, y - 1, Nx + 2 * B);
		xymp = index(x - 1, y + 1, Nx + 2 * B);

		// anisotropy calc:

		px = (data[0][xp] - data[0][xm]) / dx2;
		py = (data[0][yp] - data[0][ym]) / dx2;

		px2 = px*px;
		py2 = py*py;
		xnorm = px2 + py2;

		*a_n = 1.;
		*da_n = 0.;

		if (xnorm>1.e-6) {

			pxx = (data[0][xp] - 2.*data[0][ind] + data[0][xm]) / dxs;
			pyy = (data[0][yp] - 2.*data[0][ind] + data[0][ym]) / dxs;
			pxy = (data[0][xypp] - data[0][xymp] - data[0][xypm] + data[0][xymm]) / dxs4;

			pxpy = px*py;

			px4 = px2*px2;
			px6 = px2*px4;
			py4 = py2*py2;
			py6 = py2*py4;

			px2y2 = px2*py2;
			px4y2 = px4*py2;
			py4x2 = py4*px2;

			termx = px2y2 - py4;
			termy = px2y2 - px4;

			t1x = pxx*(px4y2 + 4.*py4x2 - py6);
			t1y = pyy*(py4x2 + 4.*px4y2 - px6);

			t2x = pxx*px2*termx*termx;
			t2y = pyy*py2*termy*termy;

			t1xy = pxy*pxpy*(-2.*px2y2);

			t2xy = 2.*pxy*pxpy*termx*termy;

			xnorm2 = xnorm*xnorm;
			xnorm3 = xnorm2*xnorm;

			*a_n = 1. - 3.*prms.E4 + 4.*prms.E4 *(px4 + py4) / xnorm2;
			*da_n = 16.*prms.E4*(*a_n)*(t1x + t1y) / xnorm3 + 64.*prms.E4*(*a_n) *(t1xy) / xnorm3 + 256.*prms.E4*prms.E4*((t2x + t2y) / xnorm2 + (t2xy) / xnorm2) / xnorm3;
		}
		return;
	}
};
#endif
