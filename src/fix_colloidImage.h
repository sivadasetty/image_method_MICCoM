#ifdef FIX_CLASS

FixStyle(colloid/image, FixColloidImage)

#else

#ifndef LMP_FIX_COLLOID_IMAGE_H
#define LMP_FIX_COLLOID_IMAGE_H

#include "fix.h"
#include "stdio.h"
#include <vector>

namespace LAMMPS_NS {

class FixColloidImage : public Fix 
{
private:
	double xlo = 0.0;
	double xhi = 1.0;
	const int ngauss = 5;
	const double _xg0[5] = {-0.9061798459386640,-0.5384693101056831,0.00000000000000000,0.5384693101056831,0.9061798459386640};
	const double _wg0[5] = {0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891};

public:
  	FixColloidImage(class LAMMPS *, int, char **);
  	~FixColloidImage();
 
 	int itype0;
  	
	double eouter, einner;
  	
	double **sigma;

        double **delta; // for lj/expand (added by Siva Dasetty)

	int nlocal;

	int *tag;
	
	double *q;
	
	double **f;

	double **x;

	int *type;

	double qqrd2e;
	
	// helper variables for image method kernel function 
	double _e, _g;
	double Rxkj, Rykj, Rzkj, Rkj2, rkj;
	double Rxij, Ryij, Rzij, Rij2, rij;
        double ukj, vkj, wkj;
	double aux1, aux2;
        double auxv_x_integ, auxv_y_integ, auxv_z_integ, aux3_integ, aux3Sqrt_integ;
        double auxv_x_delta, auxv_y_delta, auxv_z_delta, aux3_delta, aux3Sqrt_delta;
        double aa;
	double prefactor, etotal;

	//
  	int setmask();
  	void init();
  	void setup(int);
  	void min_setup(int);
  	void post_force(int);
  	void min_post_force(int);
  	double compute_scalar();
  	double compute_vector(int);
  	void initial_integrate(int);

  /*
    function called to calculate the polarization correction between {i;j;k}, which means
    k-th charge (either colloid or ion) polarizes j-th colloid, and the induced charges act
    on i-th object (either colloid or ion). so the force is F_i.  
  */
 	void force_pol (int, int, int);
  
	double energy_pol (int, int, int);
   
};
  
}


#endif
#endif

/* ERROR/WARNING messages:

*/
