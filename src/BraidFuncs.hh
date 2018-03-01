#ifndef _BRAIDFUNCS_H_
#define _BRAIDFUNCS_H_

/*-------- System --------*/

/*-------- Third Party --------*/
#include <deal.II/numerics/vector_tools.h>

#include "braid.h"
#include "braid_test.h"

/*-------- Project --------*/
#include "HeatEquation.hh"


typedef struct _braid_Vector_struct
{
  dealii::Vector<double> data;
} my_Vector;

// Wrap the Heat Equation in a struct
typedef struct _braid_App_struct
{
  HeatEquation<2> eq;
} my_App;

int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status);

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr);

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr);

int
my_Free(braid_App    app,
        braid_Vector u);

int
my_Sum(braid_App app,
       double alpha,
       braid_Vector x,
       double beta,
       braid_Vector y);

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr);

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus);

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus);

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus);

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus);

#endif // _BRAIDFUNCS_H_
