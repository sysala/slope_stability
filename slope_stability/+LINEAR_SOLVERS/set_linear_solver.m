function [linear_system_solver] = set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q)


is_agmg_present = LINEAR_SOLVERS.check_agmg_present(agmg_folder); % Check for AGMG in specified folder

if solver_type == "DFGMRES_AGMG" && is_agmg_present == 0
    warning('AGMG not found. Switching to DIRECT solver (backslash).');
    solver_type = 'DFGMRES_ICHOL';
end

switch solver_type
    case 'DIRECT'
        linear_system_solver = LINEAR_SOLVERS.DIRECT_BACKSLASH();
    case 'DFGMRES_ICHOL'
        preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_ICHOL(A, Q);
        linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
        linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing);

    case 'DFGMRES_AGMG'
        preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_AGMG(A, Q);
        linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
        linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing);
    otherwise
        error('Bad choice of the solver type');
end