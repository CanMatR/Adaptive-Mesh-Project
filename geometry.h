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


class plane : public adapt_geometry
{
        cart_vector_d origin;
        cart_vector_d width;

public:
        void set_geom(cart_vector_d origin, cart_vector_d width)
        {
                this->origin = origin;
                this->width = width;
        }
        void print_params()
        {
                printf("Wx Wy :%f %f \n", width.Vx , width.Vy);
                this->origin.output();
        }

        bool check_contained_in(double x, double y, double dx)
        {
                if ((((origin.Vx >= x) && (origin.Vx <= x + dx))) && (((origin.Vy >= y)  && (origin.Vy <= y + dx))))
                {
                        return 1;
                }
                return 0;
        }

        bool check_inside(double x, double y)
        {

          double xp;
          double yp;
          xp = fabs(x-origin.Vx);
          yp = fabs(y-origin.Vy);

	  if ((((x - origin.Vx)*(x - origin.Vx)) <= (width.Vx * width.Vx)) && (((y - origin.Vy) * (y - origin.Vy)) <= (width.Vy * width.Vy)))
            {
            return 1;
            }

        return 0;
        }
};

class cubid : public adapt_geometry
{
        cart_vector_d origin;
        cart_vector_d width;

public:
        void set_geom(cart_vector_d origin, cart_vector_d width)
        {
                this->origin = origin;
                this->width = width;
        }
        void print_params()
        {
                printf("Wx Wy Wz :%f %f %f \n", width.Vx , width.Vy,width.Vz);
                this->origin.output();
        }

        bool check_contained_in(double x, double y,double z, double dx)
        {

                if ((((origin.Vx >= x) && (origin.Vx <= x + dx))) && (((origin.Vy >= y)  && (origin.Vy <= y + dx))) && (((origin.Vz >= z)  && (origin.Vz <= z + dx))))
              {
                        return 1;
              }
                return 0;
        }

        bool check_inside(double x, double y,double z)
        {

          double xp;
          double yp;
          double zp;
          xp = fabs(x-origin.Vx);
          yp = fabs(y-origin.Vy);
          zp = fabs(z-origin.Vz);

          if ((((x - origin.Vx)*(x-origin.Vx)) <= (width.Vx * width.Vx)) && (((y- origin.Vy) * (y-origin.Vy)) <= (width.Vy * width.Vy)) && (((z- origin.Vz) * (z-origin.Vz)) <= (width.Vz * width.Vz)))
            {
            return 1;
            }
                return 0;
        }
};

#endif
