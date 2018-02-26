#ifndef _HEATEQUATION_H_
#define _HEATEQUATION_H_

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

using namespace dealii;

template <int dim>
class HeatEquation
{
public:
    HeatEquation();
    void define();
    void step(double deltaT);

private:
    void setup_system();
    void solve_time_step();
    void output_results() const;
    void refine_mesh (const unsigned int min_grid_level,
                      const unsigned int max_grid_level);

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       old_solution;
    Vector<double>       system_rhs;

    double               time;
    double               time_step;
    unsigned int         timestep_number;

    const double         theta;

    // These were originally in the run() function but because
    // I am splitting the run() function up into define and step
    // they need to become member data
    Vector<double> tmp;
    Vector<double> forcing_terms;
};

template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide ()
        :
        Function<dim>(),
        period (0.2)
    {}

    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;

private:
    const double period;
};

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
    virtual double value (const Point<dim>  &p,
                          const unsigned int component = 0) const;
};

#include "HeatEquationImplem.hh"

#endif // _HEATEQUATION_H_
