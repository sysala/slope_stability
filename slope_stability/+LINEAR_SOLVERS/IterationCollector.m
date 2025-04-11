classdef IterationCollector < handle
    %--------------------------------------------------------------------------
    % IterationCollector collects and stores iteration and timing data from
    % iterative solvers.
    %
    % This class maintains cell arrays for iteration counts, total solve times,
    % preconditioner setup times, and orthogonalization times for different
    % solver instances. It provides methods to register new instances, store
    % timing data, and retrieve aggregated statistics.
    %
    %--------------------------------------------------------------------------
    
    properties
        iterations              % Cell array storing iteration counts per instance.
        times                   % Cell array storing solve times per instance.
        preconditioner_times    % Cell array storing preconditioner setup times.
        orthogonalization_times % Cell array storing orthogonalization times.
    end
    
    methods
        function obj = IterationCollector()
            %--------------------------------------------------------------------------
            % IterationCollector Construct an IterationCollector instance.
            %
            %   obj = IterationCollector()
            %
            % Initializes empty cell arrays for iterations, times, preconditioner
            % times, and orthogonalization times.
            %--------------------------------------------------------------------------
            obj.iterations = {}; 
            obj.times = {};
            obj.preconditioner_times = {};
            obj.orthogonalization_times = {};
        end
        
        function instance_id = register_instance(obj)
            %--------------------------------------------------------------------------
            % register_instance Registers a new solver instance.
            %
            %   instance_id = obj.register_instance()
            %
            % This method allocates new cells for storing iteration counts, solve
            % times, preconditioner times, and orthogonalization times for the new
            % instance and returns its unique instance ID.
            %--------------------------------------------------------------------------
            instance_id = length(obj.iterations) + 1;
            obj.iterations{instance_id} = {};
            obj.times{instance_id} = {};
            obj.preconditioner_times{instance_id} = {};
            obj.orthogonalization_times{instance_id} = {};
        end
        
        function store_iteration(obj, instance_id, nit, time)
            %--------------------------------------------------------------------------
            % store_iteration Stores an iteration count and its corresponding time.
            %
            %   obj.store_iteration(instance_id, nit, time)
            %
            % Inputs:
            %   instance_id - Identifier for the solver instance.
            %   nit         - Number of iterations for the current step.
            %   time        - Solve time for the current step.
            %--------------------------------------------------------------------------
            obj.iterations{instance_id}{end+1} = nit;
            obj.times{instance_id}{end+1} = time;
        end

        function store_preconditioner_time(obj, instance_id, time)
            %--------------------------------------------------------------------------
            % store_preconditioner_time Stores preconditioner setup time.
            %
            %   obj.store_preconditioner_time(instance_id, time)
            %
            % Inputs:
            %   instance_id - Identifier for the solver instance.
            %   time        - Time taken for preconditioner setup.
            %--------------------------------------------------------------------------
            obj.preconditioner_times{instance_id}{end+1} = time;
        end

        function store_orthogonalization_time(obj, instance_id, time)
            %--------------------------------------------------------------------------
            % store_orthogonalization_time Stores orthogonalization time.
            %
            %   obj.store_orthogonalization_time(instance_id, time)
            %
            % Inputs:
            %   instance_id - Identifier for the solver instance.
            %   time        - Time taken for orthogonalization.
            %--------------------------------------------------------------------------
            obj.orthogonalization_times{instance_id}{end+1} = time;
        end

        % Get the sum of all iterations across all instances.
        function total_iterations = get_total_iterations(obj)
            %--------------------------------------------------------------------------
            % get_total_iterations Returns the total number of iterations.
            %
            %   total_iterations = obj.get_total_iterations()
            %
            % Sums the iteration counts stored across all registered instances.
            %--------------------------------------------------------------------------
            total_iterations = sum(cellfun(@(inst) sum(cell2mat(inst)), obj.iterations));
        end

        % Get the sum of all solve times across all instances.
        function total_time = get_total_time(obj)
            %--------------------------------------------------------------------------
            % get_total_time Returns the total time (solve + preconditioner + orthogonalization).
            %
            %   total_time = obj.get_total_time()
            %
            % Aggregates solve times, preconditioner setup times, and 
            % orthogonalization times.
            %--------------------------------------------------------------------------
            total_time = obj.get_total_solve_time() + obj.get_total_preconditioner_time() + obj.get_total_orthogonalization_time();
        end
        
        % Get the sum of preconditioner setup times across all instances.
        function total_preconditioner_time = get_total_preconditioner_time(obj)
            %--------------------------------------------------------------------------
            % get_total_preconditioner_time Returns the total preconditioner setup time.
            %
            %   total_preconditioner_time = obj.get_total_preconditioner_time()
            %--------------------------------------------------------------------------
            total_preconditioner_time = sum(cellfun(@(inst) sum(cell2mat(inst)), obj.preconditioner_times));
        end

        % Get the sum of orthogonalization times across all instances.
        function total_orthogonalization_time = get_total_orthogonalization_time(obj)
            %--------------------------------------------------------------------------
            % get_total_orthogonalization_time Returns the total orthogonalization time.
            %
            %   total_orthogonalization_time = obj.get_total_orthogonalization_time()
            %--------------------------------------------------------------------------
            total_orthogonalization_time = sum(cellfun(@(inst) sum(cell2mat(inst)), obj.orthogonalization_times));
        end
        
        % Get the sum of solve times across all instances.
        function total_solve_time = get_total_solve_time(obj)
            %--------------------------------------------------------------------------
            % get_total_solve_time Returns the total solve time.
            %
            %   total_solve_time = obj.get_total_solve_time()
            %--------------------------------------------------------------------------
            total_solve_time = sum(cellfun(@(inst) sum(cell2mat(inst)), obj.times));
        end
        
        % Get all iteration counts as a single vector.
        function iteration_vector = get_iterations_vector(obj)
            %--------------------------------------------------------------------------
            % get_iterations_vector Returns all iteration counts in one vector.
            %
            %   iteration_vector = obj.get_iterations_vector()
            %--------------------------------------------------------------------------
            iteration_vector = cell2mat(cellfun(@cell2mat, obj.iterations, 'UniformOutput', false));
        end

        % Get all preconditioner setup times as a single vector.
        function preconditioner_time_vector = get_preconditioner_time_vector(obj)
            %--------------------------------------------------------------------------
            % get_preconditioner_time_vector Returns all preconditioner setup times.
            %
            %   preconditioner_time_vector = obj.get_preconditioner_time_vector()
            %--------------------------------------------------------------------------
            preconditioner_time_vector = cell2mat(cellfun(@cell2mat, obj.preconditioner_times, 'UniformOutput', false));
        end

        % Get all orthogonalization times as a single vector.
        function orthogonalization_time_vector = get_orthogonalization_time_vector(obj)
            %--------------------------------------------------------------------------
            % get_orthogonalization_time_vector Returns all orthogonalization times.
            %
            %   orthogonalization_time_vector = obj.get_orthogonalization_time_vector()
            %--------------------------------------------------------------------------
            orthogonalization_time_vector = cell2mat(cellfun(@cell2mat, obj.orthogonalization_times, 'UniformOutput', false));
        end
        
        % Get all solve times as a single vector.
        function solve_time_vector = get_solve_time_vector(obj)
            %--------------------------------------------------------------------------
            % get_solve_time_vector Returns all solve times in one vector.
            %
            %   solve_time_vector = obj.get_solve_time_vector()
            %--------------------------------------------------------------------------
            solve_time_vector = cell2mat(cellfun(@cell2mat, obj.times, 'UniformOutput', false));
        end
    end
end
