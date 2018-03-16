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
  int level;
  double deltaT;

  int index;
  braid_StepStatusGetLevel(status, &level);
  braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
  braid_StepStatusGetTIndex(status, &index);

  deltaT = tstop - tstart;

  dealii::Vector<double>& solution = u->data;

  HeatEquation<2>& heateq = app->eq;

  heateq.step(solution, deltaT, tstart, index);

  return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
  my_Vector *u = new(my_Vector);
  int size = app->eq.size();
  u->data.reinit(size);

  if(t < 0)
    {
      double halfsize = u->data.size()/2.;
      for (int i = 0; i != u->data.size(); ++i)
        {
          u->data[i] = (i-halfsize)*(i-halfsize)*t;
        }
    }

  *u_ptr = u;

  return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
  my_Vector *v = new(my_Vector);
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
  double dot = 0.0;
  dot = u->data.l2_norm();
  *norm_ptr = dot;
  return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
  // int        index, rank, level, done;
  // double     t, error;

  // std::cout << "Access hit here: " << std::endl;
  // braid_AccessStatusGetT(astatus, &t);
  // braid_AccessStatusGetTIndex(astatus, &index);
  // braid_AccessStatusGetLevel(astatus, &level);
  // braid_AccessStatusGetDone(astatus, &done);

  // std::cout << "Result: " << std::endl
  //           << "\ttime: " << t << "\n"
  //           << "\tindex: " << index << "\n"
  //           << "\tlevel: " << level << "\n"
  //           << "\tdone?: " << done << "\n"
  //           << "NumElements: " << u->data.size() << std::endl;

  // std::cout << "Vec: ";
  // for(int i = 0; i != u->data.size(); ++i)
  //   {
  //     std::cout << " " << u->data[i];
  //   }
  // std::cout << std::endl;

  return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
  int size = app->eq.size();
  *size_ptr = (size+1)*sizeof(double);
  return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
  double *dbuffer = (double*)buffer;
  int size = u->data.size();
  dbuffer[0] = size;
  for(int i=0; i != size; ++i)
    {
      dbuffer[i+1] = (u->data)[i];
    }
  braid_BufferStatusSetSize(bstatus, (size+1)*sizeof(double));
  return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
  my_Vector *u = NULL;
  double *dbuffer = (double*)buffer;
  int size = dbuffer[0];
  u = new(my_Vector);
  u->data.reinit(size);

  for(int i = 0; i != size; ++i)
    {
      (u->data)[i] = dbuffer[i+1];
    }
  *u_ptr = u;
  return 0;
}
