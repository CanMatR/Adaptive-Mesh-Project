/* Copyright (c) 2015 Government of Canada  */
#ifndef geometryDEF
#define geometryDEF
#include "utility.h"


class adapt_geometry
{

public:
	adapt_geometry()
	{
	}
	~adapt_geometry()
	{}
	virtual bool check_inside(double x, double y, double z)
	{
		return 0;
	}
	virtual bool check_inside(double x, double y)
	{
		return 0;
	}
	virtual bool check_contained_in(double x, double y, double z, double dx, double dy, double dz)
	{
		return 0;
	}
	virtual bool check_contained_in(double x, double y, double dx, double dy)
	{
		return 0;
	}
	virtual void print_params()
	{

	}
};
class sphere : public adapt_geometry
{
	cart_vector_d origin;
	double radius;
	
public:
	void set_geom(cart_vector_d origin, double radius)
	{
		this->origin = origin;
		this->radius = radius;
	}
	void print_params()
	{
		printf("R:%f\n",this->radius);
		this->origin.output();
	}
	bool check_inside(double x, double y, double z)
	{
		if ((x - origin.Vx)*(x - origin.Vx) + (y - origin.Vy)*(y - origin.Vy) + (z - origin.Vz)*(z - origin.Vz) <= this->radius*this->radius)
		{
			return 1;
		}
		return 0;
	}
	bool check_contained_in(double x, double y, double z, double dx, double dy, double dz)
	{
		if ((origin.Vx >= x  && origin.Vx <= x + dx) && (origin.Vy >= y  && origin.Vy <= y + dy) && (origin.Vz >= z  && origin.Vz <= z + dz))
			return 1;
		return 0;
	}
};
class circle : public adapt_geometry
{
	cart_vector_d origin;
	double radius;

public:
	void set_geom(cart_vector_d origin, double radius)
	{
		this->origin = origin;
		this->radius = radius;
	}
	void print_params()
	{
		printf("R:%f\n", this->radius);
		this->origin.output();
	}
	bool check_contained_in(double x, double y, double dx, double dy)
	{
		if ((origin.Vx >= x  && origin.Vx <= x + dx) && (origin.Vy >= y  && origin.Vy <= y + dy))
			return 1;
		return 0;
	}
	bool check_inside(double x, double y)
	{
		if ((x - origin.Vx)*(x - origin.Vx) + (y - origin.Vy)*(y - origin.Vy) <= this->radius*this->radius)
		{
			return 1;
		}
		return 0;
	}
};

#endif