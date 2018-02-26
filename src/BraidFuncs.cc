/*-------- Project --------*/
#include "BraidFuncs.hh"

int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int level, i;
   double deltaT;
   
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;

   HeatEquation<2>& heateq = app->eq;

   heateq.step(deltaT);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
  my_Vector *u;

  *u_ptr = u;

  return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
  my_Vector *v;
  v->data = u->data;
  *v_ptr = v;

  return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
  free(u);

  return 0;
}

int my_Sum(braid_App app,
           double alpha,
           braid_Vector x,
           double beta,
           braid_Vector y)
{
  Vector<double>& vec = y->data;
  vec.sadd(beta, alpha, x->data);
  return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
  int    i;
  double dot = 0.0;
  dot = u->data.l2_norm();
  *norm_ptr = sqrt(dot);

  return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{

   return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
  return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
  return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
 
  return 0;
}
