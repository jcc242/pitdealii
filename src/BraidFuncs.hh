#ifndef _BRAIDFUNCS_H_
#define _BRAIDFUNCS_H_

/*-------- System --------*/

/*-------- Third Party --------*/
#include <deal.II/numerics/vector_tools.h>

#include "braid.h"
#include "braid_test.h"

/*-------- Project --------*/
#include "HeatEquation.hh"

/**
 * \brief Struct that contains the data 
 */
typedef struct _braid_Vector_struct
{
  dealii::Vector<double> data;
} my_Vector;

// Wrap the Heat Equation in a struct
typedef struct _braid_App_struct
{
  HeatEquation<2> eq;
  int final_step;
} my_App;


/**
 * @brief my_Step - Takes a step in time, advancing the u vector
 *
 * @param app - The braid app struct
 * @param ustop - The solution data at the end of this time step
 * @param fstop - RHS data (such as forcing function?)
 * @param u - The solution data at the beginning of this time step
 * @param status - Status structure that contains various info of this time
 *
 * @return Success (0) or failure (1)
 **/
int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status);


/**
 * @brief my_Init - Initializes a solution data at the given time
 * For now, initializes the solution to zero no matter what time we are at
 *
 * @param app - The braid app struct containing user data
 * @param t - Time at which the solution is initialized
 * @param u_ptr - The solution data that needs to be filled
 *
 * @return Success (0) or failure (1)
 **/
int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr);


/**
 * @brief my_Clone - Clones a vector into a new vector
 *
 * @param app - The braid app struct containing user data
 * @param u - The existing vector containing data
 * @param v_ptr - The empty vector that needs to be filled
 *
 * @return Success (0) or failure (1)
 **/
int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr);


/**
 * @brief my_Free - Deletes a vector
 *
 * @param app - The braid app struct containing user data
 * @param u - The vector that needs to be deleted
 *
 * @return Success (0) or failure (1)
 **/
int
my_Free(braid_App    app,
        braid_Vector u);


/**
 * @brief my_Sum - Sums two vectors in an AXPY operation
 * The operation is y = alpha*x + beta*y
 *
 * @param app - The braid app struct containing user data
 * @param alpha - The coefficient in front of x
 * @param x - A vector that is multiplied by alpha then added to y
 * @param beta - The coefficient of y
 * @param y - A vector that is multiplied by beta then summed with x
 *
 * @return Success (0) or failure (1)
 **/
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
