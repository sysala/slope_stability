classdef DIRECT_BACKSLASH < LINEAR_SOLVERS.DFGMRES
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
        function obj = DIRECT_BACKSLASH(~, ~, ~, ~, ~)
            %--------------------------------------------------------------------------
            % DIRECT_BACKSLASH Constructor.
            %
            %   obj = DIRECT_BACKSLASH()
            %
            % This constructor initializes the DIRECT_BACKSLASH solver by calling the
            % superclass constructor with dummy parameters.
            %--------------------------------------------------------------------------
            dummy_prec_builder = @(A) [];  % Dummy preconditioner builder (not used).
            dummy_tolerance = 1e-6;
            dummy_max_iterations = 1;
            dummy_tol_defl = 1e-6;
            dummy_verbose = false;
            obj@LINEAR_SOLVERS.DFGMRES(dummy_prec_builder, dummy_tolerance, dummy_max_iterations, dummy_tol_defl, dummy_verbose);
        end
        
        function preconditioner = setup_preconditioner_core(~, ~)
            %--------------------------------------------------------------------------
            % setup_preconditioner_core: No preconditioning is needed for the direct solver.
            %
            %   preconditioner = obj.setup_preconditioner_core(A)
            %
            % This method returns an empty preconditioner since the direct solver does not
            % require any preconditioning.
            %--------------------------------------------------------------------------
            preconditioner = [];
        end
        

        function [u, nit] = solve_core(~, A, b)
            %--------------------------------------------------------------------------
            % solve_core computes the solution using MATLAB's backslash operator.
            %
            %   [u, nit] = obj.solve_core(A, b)
            %
            % This method directly computes the solution u = A \ b and returns a
            % zero iteration count (nit = 0) since no iterations are needed.
            %--------------------------------------------------------------------------
            u = A \ b;
            nit = 0;
        end
        
        function obj = expand_deflation_basis(obj, ~)
            %--------------------------------------------------------------------------
            % expand_deflation_basis: The direct solver does not require deflation.
            %
            %   obj = obj.expand_deflation_basis(additional_vectors)
            %
            % This method is overridden with an empty implementation as the direct
            % solver does not utilize a deflation basis.
            %--------------------------------------------------------------------------
        end
    end
end
