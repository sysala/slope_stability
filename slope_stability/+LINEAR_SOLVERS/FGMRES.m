classdef FGMRES < LINEAR_SOLVERS.DFGMRES
    %--------------------------------------------------------------------------
    % FGMRES implements a Flexible GMRES solver without deflation basis expansion.
    %
    % This class inherits from LINEAR_SOLVERS.DFGMRES and overrides the method
    % for expanding the deflation basis, effectively disabling deflation. All
    % other functionalities (preconditioning, solving, timing collection) are
    % inherited from the DFGMRES class.
    %
    %--------------------------------------------------------------------------
    
    methods
        function obj = FGMRES(preconditioner_builder, tolerance, max_iterations, tolerance_deflation_basis, verbose)
            %--------------------------------------------------------------------------
            % FGMRES Constructor.
            %
            %   obj = FGMRES(preconditioner_builder, tolerance, max_iterations, tolerance_deflation_basis, verbose)
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
        
        function expand_deflation_basis(~, ~)
            %--------------------------------------------------------------------------
            % expand_deflation_basis Disables the expansion of the deflation basis.
            %
            % This method overrides the superclass implementation to prevent any
            % changes to the deflation basis.
            %
            % Inputs:
            %   (none)
            %
            % Outputs:
            %   (none)
            %--------------------------------------------------------------------------
            % No action is taken.
        end
    end
end
