classdef CONSTITUTIVE < handle & matlab.mixin.Copyable
    %--------------------------------------------------------------------------
    % CONSTITUTIVE is a class for handling the constitutive model for
    % elastic-perfectly plastic materials based on the Mohr-Coulomb yield
    % criterion. It computes stress, consistent tangent operators, and 
    % residual forces from strain fields at integration points.
    %
    % The class stores material parameters, quadrature weights, and the 
    % strain-displacement matrix. It also measures runtimes for various 
    % operations such as reduction, stress evaluation, and tangent assembly.
    %
    % Properties:
    %   c0          - Effective cohesion at integration points.
    %   phi         - Effective friction angle (in degrees or radians).
    %   psi         - Dilatancy angle.
    %   Davis_type  - Character flag indicating the Davis' reduction approach.
    %   c_bar       - Reduced cohesion parameter computed by reduction.
    %   sin_phi     - Sine of the reduced friction angle.
    %   dim         - Spatial dimension (2 or 3).
    %   shear       - Shear modulus at integration points.
    %   bulk        - Bulk modulus at integration points.
    %   lame        - Lame's coefficient at integration points.
    %   B           - Strain-displacement matrix.
    %   S           - Stress tensor at integration points.
    %   DS          - Consistent tangent operator.
    %   WEIGHT      - Quadrature weight coefficients.
    %   n_strain    - Number of strain components (6 for 3D, 3 for 2D).
    %   n_int       - Number of integration points.
    %   AUX, iD, jD, vD_pre - Auxiliary arrays for assembling sparse matrices.
    %
    %   Runtime measurements (vectors):
    %     time_reduction, time_stress, time_stress_tangent,
    %     time_build_F, time_build_F_K_tangent, time_potential.
    %
    % Methods:
    %   CONSTITUTIVE   - Constructor.
    %   reduction      - Perform material reduction with a given factor.
    %   constitutive_problem_stress - Compute stress from displacement U.
    %   constitutive_problem_stress_tangent - Compute stress and tangent DS.
    %   build_F        - Assemble the residual force vector.
    %   build_F_K_tangent - Assemble both the residual vector and the tangent matrix.
    %   build_F_all, build_F_K_tangent_all - Combined operations with reduction.
    %   build_F_reduced, build_F_K_tangent_reduced - Operations without reduction.
    %   potential      - Compute the integrated potential.
    %   get_total_time and similar - Retrieve accumulated runtimes.
    %   get_*_time_vector - Return runtime vectors.
    %--------------------------------------------------------------------------
    
    properties
        c0
        phi
        psi
        Davis_type
        c_bar
        sin_phi
        dim
        shear
        bulk
        lame
        B
        S
        DS
        WEIGHT
        n_strain
        n_int
        AUX
        iD
        jD
        vD_pre

        % Runtime measurements (vectors)
        time_reduction
        time_stress
        time_stress_tangent
        time_build_F
        time_build_F_K_tangent
        time_potential
    end
    
    methods
        function obj = CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim)
            %--------------------------------------------------------------------------
            % CONSTITUTIVE Construct an instance of the CONSTITUTIVE class.
            %
            %   obj = CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame,
            %                      WEIGHT, n_strain, n_int, dim)
            %
            % Inputs:
            %   B          - Strain-displacement matrix.
            %   c0         - Effective cohesion.
            %   phi        - Effective friction angle.
            %   psi        - Dilatancy angle.
            %   Davis_type - Davis' reduction approach ('A', 'B', or 'C').
            %   shear      - Shear modulus.
            %   bulk       - Bulk modulus.
            %   lame       - Lame's coefficient.
            %   WEIGHT     - Quadrature weight coefficients.
            %   n_strain   - Number of strain components.
            %   n_int      - Number of integration points.
            %   dim        - Dimension (2 or 3).
            %
            %--------------------------------------------------------------------------
            obj.B = B;
            obj.c0 = c0;
            obj.phi = phi;
            obj.psi = psi;
            obj.Davis_type = Davis_type;
            obj.shear = shear;
            obj.bulk = bulk;
            obj.lame = lame;
            obj.WEIGHT = WEIGHT;
            obj.n_strain = n_strain;
            obj.n_int = n_int;
            obj.dim = dim;

            obj.AUX = reshape(1:obj.n_strain * obj.n_int, obj.n_strain, obj.n_int);
            obj.iD = repmat(obj.AUX, obj.n_strain, 1);
            obj.jD = kron(obj.AUX, ones(obj.n_strain, 1));
            obj.vD_pre = repmat(obj.WEIGHT, obj.n_strain^2, 1);

            % Initialize runtime measurement vectors as empty.
            obj.time_reduction = [];
            obj.time_stress = [];
            obj.time_stress_tangent = [];
            obj.time_build_F = [];
            obj.time_build_F_K_tangent = [];
            obj.time_potential = [];
        end
        
        function obj = reduction(obj, lambda)
            %--------------------------------------------------------------------------
            % reduction Reduces material strength parameters.
            %
            %   obj = obj.reduction(lambda)
            %
            % This method computes the reduced cohesion c_bar and the sine of the
            % reduced friction angle using a specified strength reduction factor.
            %
            % Inputs:
            %   lambda - Strength reduction factor.
            %
            % Outputs:
            %   c_bar   - Reduced cohesion.
            %   sin_phi - Sine of the reduced friction angle.
            %
            % The runtime is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            [obj.c_bar, obj.sin_phi] = CONSTITUTIVE_PROBLEM.reduction(obj.c0, obj.phi, obj.psi, lambda, obj.Davis_type);
            elapsed_time = toc(t_start);
            obj.time_reduction(end + 1) = elapsed_time;
        end
        
        function obj = constitutive_problem_stress(obj, U)
            %--------------------------------------------------------------------------
            % constitutive_problem_stress Computes the stress tensor from displacement U.
            %
            %   obj = obj.constitutive_problem_stress(U)
            %
            % Input:
            %   U - Displacement vector.
            %
            % This method evaluates the strain at integration points via the
            % strain-displacement matrix B and computes the stress S using the
            % appropriate 2D or 3D constitutive problem.
            %
            % The runtime is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            E = obj.B * U(:);  % Strain at integration points.
            E = reshape(E, obj.n_strain, []);
            if obj.dim == 2
                [obj.S] = CONSTITUTIVE_PROBLEM.constitutive_problem_2D(E, obj.c_bar, obj.sin_phi, obj.shear, obj.bulk, obj.lame);
            elseif obj.dim == 3
                [obj.S] = CONSTITUTIVE_PROBLEM.constitutive_problem_3D(E, obj.c_bar, obj.sin_phi, obj.shear, obj.bulk, obj.lame);
            else
                error("wrong dimension");
            end
            elapsed_time = toc(t_start);
            obj.time_stress(end + 1) = elapsed_time;
        end

        function obj = constitutive_problem_stress_tangent(obj, U)
            %--------------------------------------------------------------------------
            % constitutive_problem_stress_tangent Computes stress and tangent operator.
            %
            %   obj = obj.constitutive_problem_stress_tangent(U)
            %
            % Input:
            %   U - Displacement vector.
            %
            % This method computes both the stress tensor S and the consistent tangent
            % operator DS by evaluating the strain field and applying the appropriate
            % constitutive model (2D or 3D).
            %
            % The runtime is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            E = obj.B * U(:);  % Strain at integration points.
            E = reshape(E, obj.n_strain, []);
            if obj.dim == 2
                [obj.S, obj.DS] = CONSTITUTIVE_PROBLEM.constitutive_problem_2D(E, obj.c_bar, obj.sin_phi, obj.shear, obj.bulk, obj.lame);
            elseif obj.dim == 3
                [obj.S, obj.DS] = CONSTITUTIVE_PROBLEM.constitutive_problem_3D(E, obj.c_bar, obj.sin_phi, obj.shear, obj.bulk, obj.lame);
            else
                error("wrong dimension");
            end
            elapsed_time = toc(t_start);
            obj.time_stress_tangent(end + 1) = elapsed_time;
        end
        
        function F = build_F(obj)
            %--------------------------------------------------------------------------
            % build_F Assembles the residual force vector.
            %
            %   F = obj.build_F()
            %
            % This method builds the force vector F from the stress tensor S and the
            % strain-displacement matrix B. The result is reshaped to match the 
            % spatial dimension.
            %
            % The runtime is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            F = obj.B' * reshape(repmat(obj.WEIGHT, obj.n_strain, 1) .* obj.S(1:obj.n_strain, :), [], 1);
            F = reshape(F, obj.dim, []);
            elapsed_time = toc(t_start);
            obj.time_build_F(end + 1) = elapsed_time;
        end

        function [F, K_tangent] = build_F_K_tangent(obj)
            %--------------------------------------------------------------------------
            % build_F_K_tangent Assembles the residual force and tangent stiffness.
            %
            %   [F, K_tangent] = obj.build_F_K_tangent()
            %
            % This method computes the force vector F and the tangent stiffness matrix
            % K_tangent by assembling the sparse material tangent D_p and applying the
            % strain-displacement matrix.
            %
            % The tangent matrix is symmetrized.
            %
            % The runtime is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            F = obj.build_F();
            vD = obj.vD_pre .* obj.DS;
            D_p = sparse(obj.iD(:), obj.jD(:), vD(:), obj.n_strain * obj.n_int, obj.n_strain * obj.n_int);
            K_tangent = obj.B' * D_p * obj.B;
            K_tangent = (K_tangent + K_tangent') / 2;
            elapsed_time = toc(t_start);
            obj.time_build_F_K_tangent(end + 1) = elapsed_time;
        end

        function F = build_F_all(obj, lambda, U)
            %--------------------------------------------------------------------------
            % build_F_all Performs reduction, stress evaluation, and assembles F.
            %
            %   F = obj.build_F_all(lambda, U)
            %
            % This method first reduces material parameters, then computes the stress,
            % and finally assembles the residual force vector F.
            %--------------------------------------------------------------------------
            obj.reduction(lambda);
            obj.constitutive_problem_stress(U);
            F = obj.build_F();
        end

        function [F, K_tangent] = build_F_K_tangent_all(obj, lambda, U)
            %--------------------------------------------------------------------------
            % build_F_K_tangent_all Performs reduction, computes stress/tangent, and builds F and K.
            %
            %   [F, K_tangent] = obj.build_F_K_tangent_all(lambda, U)
            %
            % This method combines material reduction, stress and tangent evaluation, and
            % assembly of both the residual force vector F and the tangent stiffness matrix.
            %--------------------------------------------------------------------------
            obj.reduction(lambda);
            obj.constitutive_problem_stress_tangent(U);
            [F, K_tangent] = obj.build_F_K_tangent();
        end

        function F = build_F_reduced(obj, U)
            %--------------------------------------------------------------------------
            % build_F_reduced Builds the residual force vector without material reduction.
            %
            %   F = obj.build_F_reduced(U)
            %
            % This method computes the stress and assembles F, assuming that reduction
            % has already been performed.
            %--------------------------------------------------------------------------
            obj.constitutive_problem_stress(U);
            F = obj.build_F();
        end

        function [F, K_tangent] = build_F_K_tangent_reduced(obj, U)
            %--------------------------------------------------------------------------
            % build_F_K_tangent_reduced Builds F and tangent matrix without reduction.
            %
            %   [F, K_tangent] = obj.build_F_K_tangent_reduced(U)
            %
            % This method computes the stress and consistent tangent operator and
            % assembles F and K_tangent, assuming reduction is already done.
            %--------------------------------------------------------------------------
            obj.constitutive_problem_stress_tangent(U);
            [F, K_tangent] = obj.build_F_K_tangent();
        end

        function Psi_integrated = potential(obj, U)
            %--------------------------------------------------------------------------
            % potential Computes the integrated potential.
            %
            %   Psi_integrated = obj.potential(U)
            %
            % This method calculates the potential function Psi at integration
            % points from the strain field and integrates it using the quadrature
            % weights.
            %
            % The runtime is recorded.
            %--------------------------------------------------------------------------
            t_start = tic;
            E = obj.B * U(:);   % Strain at integration points.
            E = reshape(E, obj.n_strain, []);
            if obj.dim == 2
                Psi = CONSTITUTIVE_PROBLEM.potential_2D(E, obj.c_bar, obj.sin_phi, obj.shear, obj.bulk, obj.lame);
            elseif obj.dim == 3
                Psi = CONSTITUTIVE_PROBLEM.potential_3D(E, obj.c_bar, obj.sin_phi, obj.shear, obj.bulk, obj.lame);
            else
                error("wrong dimension");
            end
            Psi_integrated = obj.WEIGHT * Psi';
            elapsed_time = toc(t_start);
            obj.time_potential(end + 1) = elapsed_time;
        end
        
        % Methods to get total runtime measurements.
        function total_time = get_total_time(obj)
            %--------------------------------------------------------------------------
            % get_total_time Returns the sum of all recorded runtimes.
            %--------------------------------------------------------------------------
            total_time = sum(obj.time_reduction) + sum(obj.time_stress) + ...
                         sum(obj.time_stress_tangent) + sum(obj.time_build_F) + ...
                         sum(obj.time_build_F_K_tangent) + sum(obj.time_potential);
        end

        function total_time = get_total_reduction_time(obj)
            total_time = sum(obj.time_reduction);
        end

        function total_time = get_total_stress_time(obj)
            total_time = sum(obj.time_stress);
        end

        function total_time = get_total_stress_tangent_time(obj)
            total_time = sum(obj.time_stress_tangent);
        end

        function total_time = get_total_build_F_time(obj)
            total_time = sum(obj.time_build_F);
        end

        function total_time = get_total_build_F_K_tangent_time(obj)
            total_time = sum(obj.time_build_F_K_tangent);
        end

        function total_time = get_total_potential_time(obj)
            total_time = sum(obj.time_potential);
        end

        % Methods to get runtime vectors.
        function time_vector = get_reduction_time_vector(obj)
            time_vector = obj.time_reduction;
        end

        function time_vector = get_stress_time_vector(obj)
            time_vector = obj.time_stress;
        end

        function time_vector = get_stress_tangent_time_vector(obj)
            time_vector = obj.time_stress_tangent;
        end

        function time_vector = get_build_F_time_vector(obj)
            time_vector = obj.time_build_F;
        end

        function time_vector = get_build_F_K_tangent_time_vector(obj)
            time_vector = obj.time_build_F_K_tangent;
        end

        function time_vector = get_potential_time_vector(obj)
            time_vector = obj.time_potential;
        end
    end
end
