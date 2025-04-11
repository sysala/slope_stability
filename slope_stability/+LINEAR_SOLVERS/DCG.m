classdef DCG < LINEAR_SOLVERS.DFGMRES
    %--------------------------------------------------------------------------
    % DIRECT_BACKSLASH implements a direct solver using MATLAB's backslash operator.
    %
    % This class inherits from LINEAR_SOLVERS.DFGMRES to provide a consistent
    % interface with other iterative solvers. It does not require any preconditioning,
    % deflation basis expansion, or orthogonalization, so those methods are overridden
    % with empty implementations.
    %
    % The constructor takes no input arguments.
    %--------------------------------------------------------------------------
    
    methods
        function obj = DCG(preconditioner_builder, tolerance, max_iterations, tolerance_deflation_basis, verbose)
            %--------------------------------------------------------------------------
            % DCG Constructor.
            %
            %   obj = DCG(preconditioner_builder, tolerance, max_iterations, tolerance_deflation_basis, verbose)
            %
            % Inputs:
            %   preconditioner_builder   - Function or handle to build the preconditioner.
            %   tolerance                - Convergence tolerance for GMRES.
            %   max_iterations           - Maximum number of iterations allowed.
            %   tolerance_deflation_basis- Tolerance for deflation basis orthogonalization.
            %   verbose                  - Logical flag for verbose output.
            %
            % The constructor calls the superclass (DFGMRES) constructor.
            %--------------------------------------------------------------------------
            obj@LINEAR_SOLVERS.DFGMRES(preconditioner_builder, tolerance, max_iterations, tolerance_deflation_basis, verbose);
        end
        
        function [u, nit] = solve_core(obj, A, b)
            %--------------------------------------------------------------------------
            % solve_core is the core implementation of the flexible GMRES solver.
            %
            %   [u, nit] = obj.solve_core(A, b)
            %
            % It calls the function flexible_GMRES_deflate from the LINEAR_SOLVERS
            % package with the provided parameters.
            %--------------------------------------------------------------------------
            [u, nit, residuals] = LINEAR_SOLVERS.dcg(A, b, ...
                obj.preconditioner, obj.deflation_basis, ...
                obj.max_iterations, obj.tolerance, []);
        end
    end
end
