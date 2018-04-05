#include <iomanip>
#include <math.h>
#include "Utilities.hh"

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int component) const
{
  (void) component;
  Assert (component == 0, ExcIndexRange(component, 0, 1));
  Assert (dim == 2, ExcNotImplemented());

  double time = this->get_time();

  if ((p[0] > 0.5) && (p[1] > -0.5))
    {
      return std::exp(-0.5*(time-0.125)*(time-0.125)/(0.005));
    }
  else if ((p[0] > -0.5) && (p[1] > 0.5))
    {
      return std::exp(-0.5*(time-0.375)*(time-0.375)/(0.005));
    }
  else
    {
      return 0;
    }
}

template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &/*p*/,
                                   const unsigned int component) const
{
  (void) component;
  Assert (component == 0, ExcIndexRange(component, 0, 1));
  return 0;
}

template <int dim>
HeatEquation<dim>::HeatEquation ()
  :
  fe(1),
  dof_handler(triangulation),
  theta(0.5)
{
}

template <int dim>
void HeatEquation<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  pout() << std::endl
         << "==========================================="
         << std::endl
         << "Number of active cells: " << triangulation.n_active_cells()
         << std::endl
         << "Number of degrees of freedom: " << dof_handler.n_dofs()
         << std::endl
         << std::endl;

  constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ true);
  sparsity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparsity_pattern);
  laplace_matrix.reinit(sparsity_pattern);
  system_matrix.reinit(sparsity_pattern);

  MatrixCreator::create_mass_matrix(dof_handler,
                                    QGauss<dim>(fe.degree+1),
                                    mass_matrix);
  MatrixCreator::create_laplace_matrix(dof_handler,
                                       QGauss<dim>(fe.degree+1),
                                       laplace_matrix);

  solution.reinit(dof_handler.n_dofs());

  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void HeatEquation<dim>::solve_time_step()
{
  SolverControl solver_control(1000, 1.e-10 * system_rhs.l2_norm());
  SolverCG<> cg(solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);

  cg.solve(system_matrix, solution, system_rhs,
           preconditioner);

  constraints.distribute(solution);

  pout() << "     " << solver_control.last_step()
         << " CG iterations." << std::endl;
}



template <int dim>
void HeatEquation<dim>::output_results(int a_time_idx,
                                       double a_time,
                                       Vector<double>& a_solution) const
{

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.time = a_time;
  vtk_flags.cycle = a_time_idx;

  DataOut<dim> data_out;
  data_out.set_flags(vtk_flags);

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(a_solution, "U");

  data_out.build_patches();

  const std::string filename = "solution-parallel-"
    + Utilities::int_to_string(a_time_idx, 3) +
    ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}

// This function won't make much sense in real parallel in time...
template <int dim>
void HeatEquation<dim>::define()
{
  const unsigned int initial_global_refinement = 4;

  GridGenerator::hyper_L (triangulation);
  triangulation.refine_global (initial_global_refinement);

  setup_system();

  tmp.reinit (solution.size());
  forcing_terms.reinit (solution.size());

  VectorTools::interpolate(dof_handler,
                           Functions::ZeroFunction<dim>(),
                           old_solution);
  solution = old_solution;

  int time_step = 0;
  double time = 0.;

  output_results(time_step, time, solution);
}

template<int dim>
void HeatEquation<dim>::step(Vector<double>& braid_data,
                             double deltaT,
                             double a_time,
                             int a_time_idx)
{
  // Set old solution to the braid data
  old_solution = braid_data;
  a_time += deltaT;
  ++a_time_idx;

  mass_matrix.vmult(system_rhs, old_solution);

  laplace_matrix.vmult(tmp, old_solution);

  system_rhs.add(-(1 - theta) * deltaT, tmp);

  RightHandSide<dim> rhs_function;
  rhs_function.set_time(a_time);
  VectorTools::create_right_hand_side(dof_handler,
                                      QGauss<dim>(fe.degree+1),
                                      rhs_function,
                                      tmp);

  forcing_terms = tmp;
  forcing_terms *= deltaT * theta;

  rhs_function.set_time(a_time - deltaT);
  VectorTools::create_right_hand_side(dof_handler,
                                      QGauss<dim>(fe.degree+1),
                                      rhs_function,
                                      tmp);

  forcing_terms.add(deltaT * (1 - theta), tmp);

  system_rhs += forcing_terms;

  system_matrix.copy_from(mass_matrix);
  system_matrix.add(theta * deltaT, laplace_matrix);

  constraints.condense (system_matrix, system_rhs);

  {
    BoundaryValues<dim> boundary_values_function;
    boundary_values_function.set_time(a_time);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             boundary_values_function,
                                             boundary_values);

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
    // system_rhs different here
  }

  solve_time_step();


  old_solution = solution;
  // Also set braid solution to the new solution
  braid_data = solution;
}

template<int dim>
int HeatEquation<dim>::size() const
{
  return solution.size();
}

template<int dim> void
HeatEquation<dim>::dump_vec(std::ofstream& file,
                            const Vector<double>& vector) const
{
  file << "Dumping vec:";
  for(int i=0; i != vector.size(); ++i)
    {
      file << "\n" << vector[i];
    }
  file << std::endl;
}
