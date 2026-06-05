% RUN ALL ALGORITHMS ON ONE PROBLEM (CLEAN RELEASE)
function run_all_algorithms_on_problem(problem)
% problem: struct with fields
%   - name    (string)
%   - dim     (integer)
%   - fhandle (function handle)
%   - startp  (column vector)

    % General algorithm options
    algs_option.r = 100;
    algs_option.q = 2;
    algs_option.corr = false;
    algs_option.eps = 0;
    algs_option.alpha_start = 1.0;
    algs_option.function_budget = 10000 * (problem.dim + 1);
    algs_option.noise_coeff = 0.1;

    % SDS params
    sds_option.theta_sds   = 0.1;
    sds_option.theta_c_sds = 0.999;
    sds_option.theta_exp_sds = 1.001;

    % Sequential SDS options
    seq_sds_option.gamma     = sds_option.theta_sds;
    seq_sds_option.theta_c   = sds_option.theta_c_sds;
    seq_sds_option.theta_exp = sds_option.theta_exp_sds;
    seq_sds_option.correlation = algs_option.corr;
    seq_sds_option.seq_variant  = 1;

    % SDS v2 extras
    seq_sds_option.alpha_threshold = 0.5;
    seq_sds_option.theta_c1 = 0.99;
    seq_sds_option.theta_exp1 = 1 +  (seq_sds_option.theta_exp - 1) * problem.dim;

    % DSE-DD params (tested set)
    dse_option.gamma_dse_c  = 0.999;
    dse_option.gamma_dse_e  = 1/0.999;
    dse_option.theta_dse    = 0.1;
    dse_option.alpha_threshold = 0.5;
    dse_option.gamma_dse_c1 = 0.99;
    dse_option.gamma_dse_e1 = 1 +  (dse_option.gamma_dse_e - 1) * problem.dim;
    dse_option.ls_threshold = 1e-2;

    % Run SDS
    fprintf('Running SDS...\n');
    tic;
    [x_sds, f_sds, fvec_sds] = direct_search_ns_noise_v2(problem, algs_option.function_budget, mergeStructs(algs_option, seq_sds_option));
    time_sds = toc;

    % Run DSE-DD
    fprintf('Running DSE-DD...\n');
    tic;
    [x_dse, f_dse, fvec_dse] = sds_sls_v2(problem, mergeStructs(algs_option, dse_option));
    time_dse = toc;

    % Run StoMADS (optional, requires third_party files)
    fprintf('Running StoMADS...\n');
    try
        s_options.theta_c = 0.999;
        s_options.alpha_start = 2;
        s_options.gamma = 0.1;
        s_options.noise_coeff = algs_option.noise_coeff;
        s_options.sample_size = 1;
        tic;
        [x_sto, f_sto, fvec_sto] = stomads_wrapper(problem, algs_option.function_budget,  s_options);
        time_sto = toc;
    catch ME
        warning(ME.identifier, 'StoMADS unavailable or failed: %s. Skipping StoMADS.', ME.message);
        x_sto = NaN; f_sto = NaN; time_sto = 0;
        fvec_sto = NaN(algs_option.function_budget,1);
    end

    % Random nonconvex search baseline
    fprintf('Running Random-Search (Sect. 7)...\n');
    n = problem.dim;
    f = problem.fhandle;
    x0 = problem.startp(:);

    rnd_option.B = eye(n);
    mu = 1e-4; rnd_option.mu = mu;
    rnd_option.h = 1e-3;
    rnd_option.N = algs_option.function_budget;

    tic;
    [x_rnd, f_rnd, fvec_rnd, ~] = random_nonconvex_search_noisy(f, x0, mergeStructs(algs_option, rnd_option));
    time_rnd = toc;

    % Shamir Algorithm 1 + 2 
    fprintf('Running Shamir Algorithm 1 + 2...\n');
    shamir_option.N = algs_option.function_budget;

    shamir_option.mu = 1e-3;
    shamir_option.h  = 1e-3;
    shamir_option.nu = 1e-2;
    shamir_option.D  = 1e-2;
    shamir_option.batch_size = 1;

    % Use this for performance/data profiles.
    shamir_option.output_rule = 'paper';

    tic;
    [x_sham, f_sham, fvec_sham, ~] = shamir(f, x0, mergeStructs(algs_option, shamir_option));
    time_sham = toc;

    % Save consolidated results (repo-level results folder)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    results_dir = fullfile(repo_root,'results');
    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    save(fullfile(results_dir,'results_placeholder.mat'), ...
         'x_sds','f_sds','fvec_sds','time_sds', ...
         'x_dse','f_dse','fvec_dse','time_dse', ...
         'x_sto','f_sto','fvec_sto','time_sto', ...
         'x_rnd','f_rnd','fvec_rnd','time_rnd', ...
         'x_sham','f_sham','fvec_sham',"time_sham");
end

function S = mergeStructs(varargin)
    S = struct();
    for i = 1:numel(varargin)
        Si = varargin{i};
        fn = fieldnames(Si);
        for k = 1:numel(fn)
            S.(fn{k}) = Si.(fn{k});
        end
    end
end
