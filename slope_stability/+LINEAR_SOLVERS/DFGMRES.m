classdef DFGMRES < handle
    %--------------------------------------------------------------------------
    % DFGMRES implements a Deflated Flexible GMRES solver with preconditioning.
    %
    % This class solves linear systems using a flexible GMRES method that
    % incorporates deflation of the Krylov subspace. It uses a provided
    % preconditioner builder to construct the preconditioner, and collects
    % iteration and timing data through an IterationCollector.
    %
    % PROPERTIES:
    %   deflation_basis         - Matrix whose columns form the deflation basis.
    %   preconditioner          - The current preconditioner.
    %   preconditioner_builder  - Function handle or object to build a preconditioner.
    %   tolerance               - Convergence tolerance for the GMRES solver.
    %   max_iterations          - Maximum number of iterations allowed.
    %   tolerance_deflation_basis - Tolerance for orthogonalizing the deflation basis.
    %   verbose                 - Flag for verbose output.
    %   iteration_collector     - Shared IterationCollector object for storing solver metrics.
    %   instance_id             - Unique instance identifier registered in the iteration_collector.
    %
    % METHODS:
    %   DFGMRES                 - Constructor.
    %   setup_preconditioner     - Sets up the preconditioner using the preconditioner_builder.
    %   A_orthogonalize         - Orthogonalizes the deflation basis with respect to A.
    %   solve                   - Solves the linear system A*u = b using FGMRES with deflation.
    %   solve_core              - Core flexible GMRES solver function.
    %   expand_deflation_basis  - Expands the deflation basis with additional vectors.
    %   copy                    - Creates a deep copy of the solver object with a new instance ID.
    %--------------------------------------------------------------------------
    
    properties
        deflation_basis
        preconditioner
        preconditioner_builder
        tolerance
        max_iterations
        tolerance_deflation_basis
        verbose
        iteration_collector  % Shared among all copies.
        instance_id          % Unique ID in iteration_collector.
    end
    
    methods
        function obj = DFGMRES(preconditioner_builder, tolerance, max_iterations, tolerance_deflation_basis, verbose)
            %--------------------------------------------------------------------------
            % DFGMRES Constructor.
            %
            %   obj = DFGMRES(preconditioner_builder, tolerance, max_iterations, ...
            %                 tolerance_deflation_basis, verbose)
            %
            % Inputs:
            %   preconditioner_builder  - Function or object handle to build the preconditioner.
            %   tolerance               - Convergence tolerance for GMRES.
            %   max_iterations          - Maximum number of iterations.
            %   tolerance_deflation_basis - Tolerance used when orthogonalizing the deflation basis.
            %   verbose                 - Logical flag for verbose output.
            %
            %--------------------------------------------------------------------------
            obj.deflation_basis = [];
            obj.preconditioner = [];
            obj.preconditioner_builder = preconditioner_builder;
            obj.tolerance = tolerance;
            obj.max_iterations = max_iterations;
            obj.tolerance_deflation_basis = tolerance_deflation_basis;
            obj.verbose = verbose;
            
            % Initialize the iteration collector and register this instance.
            obj.iteration_collector = LINEAR_SOLVERS.IterationCollector();
            obj.instance_id = obj.iteration_collector.register_instance();
        end
        
        function obj = setup_preconditioner(obj, A)
            %--------------------------------------------------------------------------
            % setup_preconditioner constructs and stores the preconditioner for matrix A.
            %
            %   obj = obj.setup_preconditioner(A)
            %
            % This method calls the internal function setup_preconditioner_core to
            % build a preconditioner for the provided matrix A using the preconditioner_builder.
            % It then records the time required for preconditioner setup in the iteration
            % collector.
            %
            % INPUT:
            %   A - Matrix for which the preconditioner is to be constructed.
            %
            % OUTPUT:
            %   The preconditioner is stored in obj.preconditioner.
            %--------------------------------------------------------------------------
            t_start = tic;
            obj.preconditioner = obj.setup_preconditioner_core(A);
            elapsed_time = toc(t_start);
            obj.iteration_collector.store_preconditioner_time(obj.instance_id, elapsed_time);
        end
        
        function preconditioner = setup_preconditioner_core(obj, A)
            %--------------------------------------------------------------------------
            % setup_preconditioner_core builds a preconditioner for matrix A.
            %
            %   preconditioner = obj.setup_preconditioner_core(A)
            %
            % This internal method applies the preconditioner_builder to the given
            % matrix A to construct a preconditioner. This method encapsulates the
            % preconditioner construction logic.
            %
            % INPUT:
            %   A - Matrix for which the preconditioner is to be constructed.
            %
            % OUTPUT:
            %   preconditioner - The preconditioner for matrix A.
            %--------------------------------------------------------------------------
            preconditioner = obj.preconditioner_builder(A);
        end


        function obj = A_orthogonalize(obj, A)
            %--------------------------------------------------------------------------
            % A_orthogonalize orthogonalizes the deflation basis with respect to A.
            %
            %   obj = obj.A_orthogonalize(A)
            %
            % This method updates the deflation basis by orthogonalizing it against
            % A using a specified tolerance. The time taken is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            obj.deflation_basis = LINEAR_SOLVERS.A_orthogonalize(obj.deflation_basis, A, obj.tolerance_deflation_basis);
            elapsed_time = toc(t_start);
            obj.iteration_collector.store_orthogonalization_time(obj.instance_id, elapsed_time);
        end

        function [u] = solve(obj, A, b)
            %--------------------------------------------------------------------------
            % solve computes the solution u of the linear system A*u = b.
            %
            %   u = obj.solve(A, b)
            %
            % This method wraps the core flexible GMRES solver, records the number
            % of iterations and time taken, and optionally prints verbose output.
            %--------------------------------------------------------------------------
            t_start = tic;
            [u, nit] = obj.solve_core(A, b);
            elapsed_time = toc(t_start);
            obj.iteration_collector.store_iteration(obj.instance_id, nit, elapsed_time);
            
            if obj.verbose
                fprintf('%d|', nit);
            end
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
            [u, nit, ~] = LINEAR_SOLVERS.dfgmres(A, b, ...
                obj.preconditioner, obj.deflation_basis, ...
                obj.max_iterations, obj.tolerance, []);
        end

        function obj = expand_deflation_basis(obj, additional_vectors)
            %--------------------------------------------------------------------------
            % expand_deflation_basis adds new vectors to the deflation basis.
            %
            %   obj = obj.expand_deflation_basis(additional_vectors)
            %
            % This method appends the columns in additional_vectors to the existing
            % deflation basis.
            %--------------------------------------------------------------------------
            obj.deflation_basis = [obj.deflation_basis, additional_vectors];
        end

        function new_obj = copy(obj)
            %--------------------------------------------------------------------------
            % copy creates a copy of the DFGMRES object with a new instance ID.
            %
            %   new_obj = obj.copy()
            %
            % The new copy shares the same iteration collector and preconditioner,
            % but is registered as a separate instance for data collection.
            %--------------------------------------------------------------------------
            new_obj = feval(class(obj), obj.preconditioner_builder, obj.tolerance, ...
                            obj.max_iterations, obj.tolerance_deflation_basis, obj.verbose);
            new_obj.deflation_basis = obj.deflation_basis;  
            new_obj.preconditioner = obj.preconditioner; 
            new_obj.iteration_collector = obj.iteration_collector;
            new_obj.instance_id = obj.iteration_collector.register_instance();
        end
    end
end
