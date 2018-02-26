#include <fstream>

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int component) const
{
  (void) component;
  Assert (component == 0, ExcIndexRange(component, 0, 1));
  Assert (dim == 2, ExcNotImplemented());

  const double time = this->get_time();
  const double point_within_period = (time/period - std::floor(time/period));

  if ((point_within_period >= 0.0) && (point_within_period <= 0.2))
    {
      if ((p[0] > 0.5) && (p[1] > -0.5))
        return 1;
      else
        return 0;
    }
  else if ((point_within_period >= 0.5) && (point_within_period <= 0.7))
    {
      if ((p[0] > -0.5) && (p[1] > 0.5))
        return 1;
      else
        return 0;
    }
  else
    return 0;
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
  time (0.0),
  time_step(1. / 500),
  timestep_number (0),
  theta(0.5)
{}

template <int dim>
void HeatEquation<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  std::cout << std::endl
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
  SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
  SolverCG<> cg(solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);

  cg.solve(system_matrix, solution, system_rhs,
           preconditioner);

  constraints.distribute(solution);

  std::cout << "     " << solver_control.last_step()
            << " CG iterations." << std::endl;
}



template <int dim>
void HeatEquation<dim>::output_results() const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "U");

  data_out.build_patches();

  const std::string filename = "solution-"
    + Utilities::int_to_string(timestep_number, 3) +
    ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}


template <int dim>
void HeatEquation<dim>::refine_mesh (const unsigned int min_grid_level,
                                     const unsigned int max_grid_level)
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(fe.degree+1),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                     estimated_error_per_cell,
                                                     0.6, 0.4);

  if (triangulation.n_levels() > max_grid_level)
    for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
         cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();
  for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active(min_grid_level);
       cell != triangulation.end_active(min_grid_level); ++cell)
    cell->clear_coarsen_flag ();

  SolutionTransfer<dim> solution_trans(dof_handler);

  Vector<double> previous_solution;
  previous_solution = solution;
  triangulation.prepare_coarsening_and_refinement();
  solution_trans.prepare_for_coarsening_and_refinement(previous_solution);

  triangulation.execute_coarsening_and_refinement ();
  setup_system ();

  solution_trans.interpolate(previous_solution, solution);
  constraints.distribute (solution);
}

template <int dim>
void HeatEquation<dim>::define()
{
  const unsigned int initial_global_refinement = 2;
  const unsigned int n_adaptive_pre_refinement_steps = 4;

  GridGenerator::hyper_L (triangulation);
  triangulation.refine_global (initial_global_refinement);

  setup_system();

  unsigned int pre_refinement_step = 0;

  Vector<double> tmp;
  Vector<double> forcing_terms;

  tmp.reinit (solution.size());
  forcing_terms.reinit (solution.size());


  VectorTools::interpolate(dof_handler,
                           Functions::ZeroFunction<dim>(),
                           old_solution);
  solution = old_solution;

  output_results();

  // I think we can safely set the time here
  time = 0.;
}

template<int dim>
void HeatEquation<dim>::step(double deltaT)
{

   time += deltaT;
   ++timestep_number;

   std::cout << "Time step " << timestep_number << " at t=" << time
             << std::endl;

   mass_matrix.vmult(system_rhs, old_solution);

   laplace_matrix.vmult(tmp, old_solution);
   system_rhs.add(-(1 - theta) * time_step, tmp);

   RightHandSide<dim> rhs_function;
   rhs_function.set_time(time);
   VectorTools::create_right_hand_side(dof_handler,
                                       QGauss<dim>(fe.degree+1),
                                       rhs_function,
                                       tmp);
   forcing_terms = tmp;
   forcing_terms *= time_step * theta;

   rhs_function.set_time(time - time_step);
   VectorTools::create_right_hand_side(dof_handler,
                                       QGauss<dim>(fe.degree+1),
                                       rhs_function,
                                       tmp);

   forcing_terms.add(time_step * (1 - theta), tmp);

   system_rhs += forcing_terms;

   system_matrix.copy_from(mass_matrix);
   system_matrix.add(theta * time_step, laplace_matrix);

   constraints.condense (system_matrix, system_rhs);

   {
     BoundaryValues<dim> boundary_values_function;
     boundary_values_function.set_time(time);

     std::map<types::global_dof_index, double> boundary_values;
     VectorTools::interpolate_boundary_values(dof_handler,
                                              0,
                                              boundary_values_function,
                                              boundary_values);

     MatrixTools::apply_boundary_values(boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
   }

   solve_time_step();

   output_results();

   old_solution = solution;
}