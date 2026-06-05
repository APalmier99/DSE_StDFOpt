%% StoMADS algorithm (https://doi.org/10.1007/s10589-020-00249-0 or https://rdcu.be/cS8fR)
% (Unconstrained) Stochastic variant of MADS, StoMADS is developed by
% K. Joseph Dzahini, PhD. (https://scholar.google.com/citations?user=qVz8p6kAAAAJ&hl=fr&oi=ao)
% with the joint collaboration of Profs. C. Audet, M. Kokkolaras and S. Le Digabel.

% See also (https://doi.org/10.1007/s10589-021-00329-9 or https://rdcu.be/cS8fS)
% for an expected complexity analysis of StoMADS and variants, and
% (https://doi.org/10.1007/s10107-022-01787-7 or https://rdcu.be/cS8gs) for StoMADS-PB
% a constraint version (with stochastically noisy constraints) of StoMADS.

%%
%  Argonne National Laboratory (USA) / Polytechnique Montreal (Canada)

%  Kwassi Joseph Dzahini, PhD. September 2022.
%%
function z = stomads_algorithm(myfun, x0, stomads_option, probspecs)   %
%%
if nargin == 3
    probspecs = struct('n', 'Dimension');
    probspecs.Dimension = length(x0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_vec = inf(stomads_option.MaxFuncEval, 1);
myfun_true = probspecs.myfun_true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_epsilon = stomads_option.GammaEpsilon;
sampl_sze = stomads_option.SampleSize;
lw_bounds = stomads_option.LowerBounds;
if isempty(lw_bounds)
    stomads_option.LowerBounds = -Inf * ones(length(x0), 1);
    lw_bounds = stomads_option.LowerBounds;
end
if isrow(lw_bounds)
    lw_bounds = lw_bounds';
end
up_bounds = stomads_option.UpperBounds;
if isempty(up_bounds)
    stomads_option.UpperBounds = Inf * ones(length(x0), 1);
    up_bounds = stomads_option.UpperBounds;
end
if isrow(up_bounds)
    up_bounds = up_bounds';
end
prev_sampl_sze = sampl_sze;
max_iter = stomads_option.MaxNumberIters;
poll_size = stomads_option.InitPollSize;
eta = 1;
tau = stomads_option.Tau;
cur_sol = x0;
if stomads_option.experiments == 1i
    probspecs.Dimension = probspecs.n;
end
if isempty(probspecs.Dimension)
    probspecs.Dimension = length(cur_sol);
end
if isrow(cur_sol)
    cur_sol = cur_sol';
end
%% Warnings
stomads_warnings;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stomads_option.warning == 1
    z = [cur_sol, f_vec, myfun_true(cur_sol)];
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if stomads_option.MeshType == -1
    poll_size_coeffs = init_poll_size_coeffs(x0, stomads_option);
    b_0 = poll_size_coeffs(:, 2);
    poll_size_vector = poll_size_coeffs(:, 1) .* 10.^(poll_size_coeffs(:, 2));
    stomads_option.InitMeshSizeVector = 10.^(poll_size_coeffs(:, 2) - ...
        abs(poll_size_coeffs(:, 2) - b_0));
elseif stomads_option.MeshType == 1
    poll_size_vector = init_poll_size_coeffs(x0, stomads_option);
    stomads_option.InitMeshSizeVector = min([poll_size_vector, poll_size_vector.^2], [], 2);
    if isrow(poll_size_vector)
        poll_size_vector = poll_size_vector';
    end
else
    stomads_option.Anisotropy = 0;
end
success = 0;
acceptance_flag = 0;
success_flag = 0;
failure_flag = 0;
if (stomads_option.FixSeed == 1) && (isempty(stomads_option.SeedValue))
    stomads_option.SeedValue = floor(1 + pi);
end
% To make sure rng is set for reproducibility
if stomads_option.FixSeed == 1
    rand('state', stomads_option.SeedValue);
    randn('state', stomads_option.SeedValue);
end
if stomads_option.MaxFuncEval < Inf
    max_iter = 10^14;
end
global nfEval nfEvalExceeded fEval_History fEval_Stats
nfEval = 1;
nfEvalExceeded = 0;
fEval_History = [];
fEval_Stats = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfEval_old = nfEval;
f_vec(nfEval_old) = myfun_true(cur_sol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:max_iter
    %% Search
    % (Optional, no strategies proposed yet. In case you propose one, please make sure
    % the trial points are always rounded to the mesh)
    % Search strategies may include, random model search (possibly in random subspaces
    % (see e.g., https://arxiv.org/abs/2207.06452)),
    % stochastic Nelder-Mead search (see https://doi.org/10.1007/s10589-018-0016-0 and
    % https://doi.org/10.1016/j.ejor.2012.02.028)
    %% Poll set
    if stomads_option.MeshType == -1
        mesh_size_vector = 10.^(poll_size_coeffs(:, 2) - abs(poll_size_coeffs(:, 2) - b_0));
        basis = gmesh_poll_basis(poll_size_coeffs, mesh_size_vector);
    elseif stomads_option.MeshType == 1
        mesh_size_vector = min([poll_size_vector, poll_size_vector.^2], [], 2);
        basis = xmesh_poll_basis(poll_size_vector, mesh_size_vector);
    else
        mesh_size = min(poll_size, poll_size^2);
        basis = basic_poll_basis(mesh_size, poll_size, probspecs);
    end
    pos_span_set = [basis, -basis];
    if success == 0
        if (stomads_option.MeshType == -1) || (stomads_option.MeshType == 1)
            poll_set = cur_sol + diag(mesh_size_vector) * pos_span_set;
        else
            poll_set = cur_sol + mesh_size * pos_span_set;
        end
    else
        if (stomads_option.MeshType == -1) || (stomads_option.MeshType == 1)
            if stomads_option.OrthoN_PlusOne == 1
                poll_set = cur_sol + diag(mesh_size_vector) * ...
                    ortho_n_plus_one(pos_span_set, last_success_direction);
            else
                poll_set = cur_sol + diag(mesh_size_vector) * ...
                    order_last(pos_span_set, last_success_direction);
            end
        else
            if stomads_option.OrthoN_PlusOne == 1
                poll_set = cur_sol + mesh_size * ...
                    ortho_n_plus_one(pos_span_set, last_success_direction);
            else
                poll_set = cur_sol + mesh_size * ...
                    order_last(pos_span_set, last_success_direction);
            end
        end
    end
    % Deleting poll points violating bounds constraints
    violat_vec = sum(poll_set < lw_bounds, 1) + sum(poll_set > up_bounds, 1);
    violat_index = find(violat_vec > 0);
    if isempty(violat_index) == 0
        poll_set(:, violat_index) = [];
    end
    if isempty(poll_set)
        if stomads_option.MeshType == -1
            poll_size_coeffs = decrease_poll_size_coeffs(poll_size_coeffs);
            poll_size_vector = poll_size_coeffs(:, 1) .* 10.^(poll_size_coeffs(:, 2));
        elseif stomads_option.MeshType == 1
            poll_size_vector = tau^2 * poll_size_vector;
        else
            poll_size = tau^2 * poll_size;
        end
        continue
    end
    %% Estimates computation
    % Estimate at the current solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if stomads_option.noise_coeff == 0
        sampl_sze = ceil(poll_size.^(-4));
    else
        sampl_sze = ceil(stomads_option.noise_coeff^2 * poll_size.^(-4));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    funct_estm1 = MCestimate(cur_sol, sampl_sze, myfun, probspecs, stomads_option);
    if funct_estm1 ~= 1i
        if success_flag == 1
            cur_obj = (sampl_sze * funct_estm1 + PrevSamplSze_Succ * funct_estm2) / ...
                (sampl_sze + PrevSamplSze_Succ);
            prev_sampl_sze = sampl_sze + PrevSamplSze_Succ;
            success_flag = 0;
            PrevSamplSze_Succ = 0;
        elseif failure_flag == 1
            cur_obj = (prev_sampl_sze * cur_obj + sampl_sze * funct_estm1) / ...
                (prev_sampl_sze + sampl_sze);
            prev_sampl_sze = prev_sampl_sze + sampl_sze;
            failure_flag = 0;
        else
            cur_obj = funct_estm1;
        end
    else
        break
    end
    % Estimates at the trial points
    for j = 1:size(poll_set, 2)
        funct_estm2 = MCestimate(poll_set(:, j), sampl_sze, myfun, probspecs, stomads_option);
        if funct_estm2 ~= 1i
            if stomads_option.MeshType == -1
                lambda = max(eta^2, min(mesh_size_vector));
            elseif stomads_option.MeshType == 1
                lambda = max(eta^2, min(poll_size_vector.^2));
            else
                lambda = poll_size^2;
            end
            if funct_estm2 <= cur_obj - gamma_epsilon * lambda
                success = 1;
                last_success_direction = pos_span_set(:, j);   % Last direction of success
                success_point = poll_set(:, j);
                acceptance_flag = 1;
                break
            end
        else
            nfEvalExceeded = 1;
            break
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    true_f = myfun_true(cur_sol);
    for j = nfEval_old+1: min(nfEval, stomads_option.MaxFuncEval)
        f_vec(j) = true_f;
    end
    nfEval_old = nfEval;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nfEvalExceeded == 1
        break
    end
    %% Checking success and failure
    if acceptance_flag == 1      %  Success
        cur_sol = success_point;
        if stomads_option.MeshType == -1
            poll_size_coeffs = increase_poll_size_coeffs(poll_size_coeffs, ...
                last_success_direction, mesh_size_vector, poll_size_vector, stomads_option);
            poll_size_vector = poll_size_coeffs(:, 1) .* 10.^(poll_size_coeffs(:, 2));
            eta = tau^(-2) * eta;
        elseif stomads_option.MeshType == 1
            poll_size_vector = incr_poll_size_vec(poll_size_vector, ...
                mesh_size_vector, last_success_direction, stomads_option);
            eta = tau^(-2) * eta;
        else
            poll_size = tau^(-2) * poll_size;
        end
        acceptance_flag = 0;
        if stomads_option.UsePreviousSamples == 1
            success_flag = 1;
            PrevSamplSze_Succ = sampl_sze;
        end
    else         % Failure
        if stomads_option.MeshType == -1
            poll_size_coeffs = decrease_poll_size_coeffs(poll_size_coeffs);
            poll_size_vector = poll_size_coeffs(:, 1) .* 10.^(poll_size_coeffs(:, 2));
            eta = tau^2 * eta;
        elseif stomads_option.MeshType == 1
            poll_size_vector = tau^2 * poll_size_vector;
            eta = tau^2 * eta;
        else
            poll_size = tau^2 * poll_size;
        end
        if stomads_option.UsePreviousSamples == 1
            failure_flag = 1;
        end
    end
    %% Display (or not) true function values and mesh/poll size parameters
    %% during optimization process
    if stomads_option.DisplayOutputs == 1
        if stomads_option.experiments == 1i
            smooth_fun = @(x)calfun_sample(x, probspecs, 'smooth');
            if (stomads_option.MeshType == -1) || (stomads_option.MeshType == 1)
                disp([smooth_fun(cur_sol), max(mesh_size_vector)]);
            else
                disp([smooth_fun(cur_sol), mesh_size]);
            end
        else
            if (stomads_option.MeshType == -1) || (stomads_option.MeshType == 1)
                disp([myfun(cur_sol), max(mesh_size_vector)]);
            else
                disp([myfun(cur_sol), mesh_size]);
            end
        end
    end
end
%%
if stomads_option.experiments == 1i  % When using 'stomads_experiments.m'
    smooth_fun = @(x)calfun_sample(x, probspecs, 'smooth');
    if (stomads_option.MeshType == -1) || (stomads_option.MeshType == 1)
        z = [smooth_fun(cur_sol), max(mesh_size_vector)];
    else
        z = [smooth_fun(cur_sol), mesh_size];
    end
else  % When using 'stomads_applications.m'
    z = [];  % Can be modified depending on the preferred output
end

%% Display or not the best solution found

if stomads_option.DisplaySolution == 1
    if stomads_option.MaxFuncEval < Inf
        tot_evals = stomads_option.MaxFuncEval;
        Sol_message = [' Best solution found after ', num2str(tot_evals), ' function evaluations: '];
    else
        tot_iters = stomads_option.MaxNumberIters;
        Sol_message = [' Best solution found after ', num2str(tot_iters), ' iterations: '];
    end
    disp(Sol_message);
    disp(num2str(cur_sol));
    fprintf(' \n ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   Do not comment/delete/modify
% txt_files_generator;

fx_star = myfun_true(cur_sol);

for j = nfEval + 1: stomads_option.MaxFuncEval
    f_vec(j) = fx_star;
end

z = struct();
z.x_star = cur_sol;
z.f_star = fx_star;
z.f_vec = f_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
