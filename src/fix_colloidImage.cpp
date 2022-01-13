/*----------------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Original author: Jiyuan Li (dePablo group, UChicago)
   Contributing author: Siva Dasetty (Ferguson lab, UChicago)
------------------------------------------------------------------------- */


#include <math.h>
#include <stdlib.h>
#include "float.h"
#include <string.h>
#include "fix_colloidImage.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "pair.h"
#include "group.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10
//#define __DEBUG
//#define __DEBUG_GAUSS
//#define _TEST_SCALING
enum{TETHER,COUPLE};

/* ---------------------------------------------------------------------- */

FixColloidImage::FixColloidImage(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  
  itype0 = 0;
  
  for(int i=3; i<narg-1; i++)
  {
    if (strcmp(arg[i],"ion-type-start") == 0) itype0 = atoi(arg[i+1]);
    else if (strcmp(arg[i],"einner") == 0) einner = atof(arg[i+1]);
  }
  
  eouter = force->dielectric;
  
#ifdef __DEBUG
  fprintf(screen, "[Fix/Colloid/Image] Start running ...\n");
  fprintf(screen, "    ion types starting ID:     %d\n", itype0);
  fprintf(screen, "    outer dielectric constant: %lf\n", eouter);
  fprintf(screen, "    inner dielectric constant: %lf\n", einner);
#endif

}

/* ---------------------------------------------------------------------- */

FixColloidImage::~FixColloidImage()
{

}

/* ---------------------------------------------------------------------- */

int FixColloidImage::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= MIN_POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixColloidImage::init()
{
  _e = (1.0 - einner/eouter) / (1.0 + einner/eouter);
  _g = 1.0 / (1.0 + einner/eouter);

  int dim;
  sigma = (double**)(force->pair->extract("sigma", dim));
  delta = (double**)(force->pair->extract("delta", dim));
}


/* ---------------------------------------------------------------------- */

void FixColloidImage::post_force(int vflag)
{
  
  type  = (int*)(atom->type);
  etotal = 0;
  nlocal = atom->nlocal;
 

 // polarization force from {i;j;i}, calculate force on k-th object and total energy;
  //j-th object is a colloid. 
  std::vector<double> fijk(3),fkji(3); 

 
  for(int j=0; j<nlocal; j++) if(type[j]<itype0)
  { 
    for(int i=0; i<nlocal; i++) if(i != j) 
    { 
      for(int k=0; k<nlocal; k++) if (k != j)
      {      
	
          #ifdef __DEBUG
          printf("  I-J-K %d-%d-%d\n", tag[i], tag[j], tag[k]);
          #endif
          // Interaction: [I] <--- [J] <--- [K]
          // F[i] = f{i;j;k)
          // F[j] = - (f{k;j;i} + f{i;j;k} )
          // F[k] = f{k;j;i}

	  		  
          if (i == k)
          {
	   force_pol(i,j,k);
          // f[i][0] += 2 * fijk[0];
	  //f[i][1] += 2 * fijk[1];
	   //f[i][2] += 2 * fijk[2];

	   //f[j][0] -= 2 * fijk[0];
	  // f[j][1] -= 2 * fijk[1];
	  // f[j][2] -= 2 * fijk[2];
          }
          else if(i < k) 
          {
            	force_pol(i,j,k);
		force_pol(k,j,i);
         	//f[i][0] += 2 * fijk[0]; 
          	//f[i][1] += 2 * fijk[1];
         	//f[i][2] += 2 * fijk[2];

          	//f[k][0] += 2 * fkji[0];
          	//f[k][1] += 2 * fkji[1];
          	//f[k][2] += 2 * fkji[2];

          	//f[j][0] -= 2 * (fijk[0] + fkji[0]);
          	//f[j][1] -= 2 * (fijk[1] + fkji[1]);
          	//f[j][2] -= 2 * (fijk[2] + fkji[2]);
	 }
      }//end k
    }// end i  
   }// end j         
}
 
/* ---------------------------------------------------------------------- */

void FixColloidImage::initial_integrate(int vflag)
{

}

/* ---------------------------------------------------------------------- */

void FixColloidImage::setup(int vflag)
{
 post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixColloidImage::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixColloidImage::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixColloidImage::compute_scalar()
{
  return etotal;
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixColloidImage::compute_vector(int n)
{
 return 0.0;
}

/*----------------------------------------------------------------------
    return polarization correction to the forces
-------------------------------------------------------------------------*/

void FixColloidImage::force_pol (int ith, int jth, int kth)
  {
    #ifdef _TEST_SCALING
	printf("compute (%d, %d, %d)\n", ith, jth, kth);
    #endif
    nlocal = atom->nlocal;
    tag   = (int*)(atom->tag);
    f = (double**)(atom->f);
    x = (double**)(atom->x);
    q = (double*)(atom->q);
    aa = ( sigma[type[jth]][type[jth]] + delta[type[jth]][type[jth]] ) / 2;
    //std::cout << "aa: " << aa << "\t" << "delta: " << delta << "\n" <<  std::endl;
    qqrd2e = force->qqrd2e;  
    //==================kernel function starts========================
    // Rxkj, Rykj, Rzkj: vector points from j-th to k-th, in x, y, z direction respectively.
    Rxkj = x[kth][0] - x[jth][0];
    Rykj = x[kth][1] - x[jth][1];
    Rzkj = x[kth][2] - x[jth][2];
    // rkj: module of the vector between jth and kth particle
    Rkj2 = Rxkj*Rxkj + Rykj*Rykj + Rzkj*Rzkj;
    rkj = sqrt (Rkj2);
    // ukj, vkj, wkj: unit vector points from jth to kth, in x, y, z, direction respectively.
    ukj = Rxkj / rkj;
    vkj = Rykj / rkj;
    wkj = Rzkj / rkj;
    // Rxij, Ryij, Rzij: vector points from ith to jth particle, in x, y, z direction respectively
    Rxij = x[ith][0] - x[jth][0];
    Ryij = x[ith][1] - x[jth][1];
    Rzij = x[ith][2] - x[jth][2];
    Rij2 = Rxij * Rxij + Ryij * Ryij + Rzij * Rzij;
    rij = sqrt(Rij2);

    //Auxiliary variables
    aux1 = _e * aa / rkj; 
    aux2 = aa * aa / rkj;
    // Tmp storage of polarization force
    std::vector<double> force_pol(3);
    
    // 1)delta_s,1 term of equation.29 in method paper
    auxv_x_delta = Rxij - aux2 * ukj;
    auxv_y_delta = Ryij - aux2 * vkj;
    auxv_z_delta = Rzij - aux2 * wkj;
    aux3Sqrt_delta = sqrt(auxv_x_delta * auxv_x_delta + auxv_y_delta * auxv_y_delta + auxv_z_delta * auxv_z_delta);
    aux3_delta = aux3Sqrt_delta * aux3Sqrt_delta * aux3Sqrt_delta;
       
    force_pol[0] = auxv_x_delta / aux3_delta;
    force_pol[1] = auxv_y_delta / aux3_delta;
    force_pol[2] = auxv_z_delta / aux3_delta;
    // 2) integral term of equation.29 in method paper 
    std::vector<double> _xg(5);
    for (int ig = 0; ig < ngauss; ++ig)
    {
      _xg[ig] = (xhi - xlo) / 2.0 * _xg0[ig] + (xhi + xlo) / 2.0; 
      auxv_x_integ = Rxij - ( pow(_xg[ig] , 1.0/_g) * aux2 * ukj );
      auxv_y_integ = Ryij - ( pow(_xg[ig] , 1.0/_g) * aux2 * vkj );
      auxv_z_integ = Rzij - ( pow(_xg[ig] , 1.0/_g) * aux2 * wkj );
      aux3Sqrt_integ = sqrt(auxv_x_integ * auxv_x_integ + auxv_y_integ * auxv_y_integ + auxv_z_integ * auxv_z_integ);
      aux3_integ = aux3Sqrt_integ * aux3Sqrt_integ * aux3Sqrt_integ;
      force_pol[0] += (xhi-xlo)/2 * (- auxv_x_integ / aux3_integ * _wg0[ig]); 
      force_pol[1] += (xhi-xlo)/2 * (- auxv_y_integ /aux3_integ * _wg0[ig]);
      force_pol[2] += (xhi-xlo)/2 * (- auxv_z_integ /aux3_integ * _wg0[ig]);
    }

    // final polarization force for term {i;j;k}
    force_pol[0] *= 1.0 / 2.0 * qqrd2e * aux1 * q[ith] * q[kth];
    force_pol[1] *= 1.0 / 2.0 * qqrd2e * aux1 * q[ith] * q[kth];
    force_pol[2] *= 1.0 / 2.0 * qqrd2e * aux1 * q[ith] * q[kth];

    //std::cout << "force_pol 0: " << force_pol[0] << "\t" << "force_pol 1: " << force_pol[1] << "\t" << "force_pol 2: " << force_pol[2] << "\n" << std::endl;
    
    // update total force
    f[ith][0] += 2 * force_pol[0];
    f[ith][1] += 2 * force_pol[1];
    f[ith][2] += 2 * force_pol[2];
    f[jth][0] -= 2 * force_pol[0];
    f[jth][1] -= 2 * force_pol[1];
    f[jth][2] -= 2 * force_pol[2];
}


  /*------------------------------------------------------------------------ 
      return the polarization correction to the total energy of the system 
  ----------------------------------------------------------------------------*/

  double FixColloidImage::energy_pol (int i, int j, int k)
  {
    int *type  = atom->type;
    double aa = ( sigma[type[i]][type[i]] + delta[type[i]][type[i]] ) / 2 ;
    int ngauss = 5;
    double *q  = atom->q;
    double **x = atom->x;
  
//    Rxkj, Rykj, Rzkj: vector points from j-th to k-th, in x, y, z direction respectively.
  //  rkj: module of the vector between j-th and k-th
   // ukj, vkj, wkj: unit vector points from j-th to k-th, in x, y, z direction repectively

    double Rxkj = x[k][0] - x[j][0];
    double Rykj = x[k][1] - x[j][1];
    double Rzkj = x[k][2] - x[j][2];
    double Rkj2 = Rxkj*Rxkj + Rykj*Rykj + Rzkj*Rzkj;
    double rkj = sqrt (Rkj2);
    double ukj = Rxkj / rkj;
    double vkj = Rykj / rkj;
    double wkj = Rzkj / rkj;


//    vector points from i-th to j-th

    double Rxij = x[i][0] - x[j][0];
    double Ryij = x[i][1] - x[j][1];
    double Rzij = x[i][2] - x[j][2];
    double Rij2 = Rxij * Rxij + Ryij * Ryij + Rzij * Rzij;
    double rij = sqrt(Rij2);

    double aux1, aux2, aux3, auxv_x, auxv_y, auxv_z, energy_pol;

    aux1 = _e * aa / rkj; 
    aux2 = aa * aa / rkj;

    //delta_s,1 term of equation.29 in method paper
    auxv_x = Rxij - aux2 * ukj;
    auxv_y = Ryij - aux2 * vkj;
    auxv_z = Rzij - aux2 * wkj;
    double aux3Sqrt = sqrt(auxv_x*auxv_x + auxv_y*auxv_y + auxv_z*auxv_z);

    //energy_pol = 1.0 / aux3Sqrt;

    
    //  Integration term of equation.14 in method paper
    

    //Transform Gauss points for integral of [0,1] from that of [-1,1]
    double xlo = 0.0;
    double xhi = 1.0;
    double _xg[] = {-0.9061798459386640,
                -0.5384693101056831,
                0.00000000000000000,
                0.5384693101056831,
                0.9061798459386640};
    double _wg[] = { 0.2369268850561891,
              0.4786286704993665,
              0.5688888888888889,
              0.4786286704993665,
              0.2369268850561891};

    //start integration, the integral is transformed by substituting s by u, where u = s^g, to
    //remove the numerical instability. 
    for (int ig = 0; ig < ngauss; ++ig)
    {
      _xg[ig] = (xhi - xlo) / 2.0 * _xg[ig] + (xhi + xlo) / 2.0; 
      auxv_x = Rxij - ( pow(_xg[ig] , 1.0/_g) * aux2 * ukj );
      auxv_y = Ryij - ( pow(_xg[ig] , 1.0/_g) * aux2 * vkj );
      auxv_z = Rzij - ( pow(_xg[ig] , 1.0/_g) * aux2 * wkj );
      aux3Sqrt = sqrt(auxv_x*auxv_x + auxv_y*auxv_y + auxv_z*auxv_z);
     // energy_pol = energy_pol - 1.0 / aux3Sqrt * _wg[ig];

//      #ifdef __DEBUG_GAUSS
  //    fprintf(screen, "%d-th gauss point is: %f; weight is: %f\n", ig, _xg[ig], _wg[ig]);
    //  #endif

  }

    // end integration

    //multiple by pre-factor in eq.14)
    energy_pol = 1.0 / 2.0 * energy_pol * aux1 * q[i] * q[k];
    return energy_pol;
  }

