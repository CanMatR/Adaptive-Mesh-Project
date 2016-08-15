/* Copyright (c) 2013 Government of Canada  */
#ifndef utilityDEF
#define utilityDEF

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <list>
#include <cmath>

static inline int index(int a, int b, int Na)
{
	return (a + (b)*Na);
}
static inline int index(int a, int b, int c, int Na, int Nb)
{
	return (a + (b + c*Nb)*Na);
}
static inline int index(int a, int b, int c, int d, int Na, int Nb, int Nc)
{
	return a + (b + (c + d*Nc)*Nb)*Na;
}


class cart_vector
{
public:
	int Vx, Vy, Vz;
	cart_vector(){ Vx = 0; Vy = 0; Vz = 0; }
	cart_vector(int Valx, int Valy, int Valz)
	{
		Vx = Valx;
		Vy = Valy;
		Vz = Valz;
	}
	~cart_vector(){}
	cart_vector& operator= (const cart_vector &rhs)
	{
		this->Vx = rhs.Vx;
		this->Vy = rhs.Vy;
		this->Vz = rhs.Vz;
		return *this;
	}
	void output()
	{
		printf("origin : x: %d  y: %d  z:  %d | ", Vx, Vy, Vz);
	}

};
class cart_vector_d
{
public:
	double Vx, Vy, Vz;
	cart_vector_d(){}
	cart_vector_d(double Valx, double Valy, double Valz)
	{
		Vx = Valx;
		Vy = Valy;
		Vz = Valz;
	}
	~cart_vector_d(){}
	cart_vector_d& operator= (const cart_vector_d &rhs)
	{
		Vx = rhs.Vx;
		Vy = rhs.Vy;
		Vz = rhs.Vz;
		return *this;
	}
	void output()
	{
		printf("origin : x: %f  y: %f  z:  %f | ", Vx, Vy, Vz);
	}

};

#endif