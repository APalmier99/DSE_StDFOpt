function make_profiles(num_runs)
% MAKE_PROFILES   Load solver histories and plot data & performance profiles.
%
% Usage:
%   make_profiles(num_runs)
%
% Reads results from results/results_num_runs_<num_runs>.mat and writes
% EPS figures under results/plots/make_profiles_num_runs_<num_runs>/

setup_paths();
repo_root = fileparts(fileparts(mfilename('fullpath')));
folder_save = fullfile(repo_root, 'results', 'plots', sprintf('make_profiles_num_runs_%d', num_runs));
if ~exist(folder_save,'dir'), mkdir(folder_save); end

solver_order = ["dse", "sds", "sto", "rnd", "sham"];
legend_names = ["DSE", "SDS", "STOMADS", "GS", "OSNNOA"];
filename = fullfile(repo_root, 'results', sprintf('results_num_runs_%d.mat', num_runs));
function_budget = 10000;
tolset = [1, 2, 4];

for itol=1:length(tolset)
    tol = tolset(itol);

    % 1) load results
    data = load(filename, 'resultcell','dimvec','solverNames');
    C = data.resultcell;             % P × S cell array
    dims = data.dimvec;              % P×1
    if isfield(data,'solverNames')
        solver_names = string(data.solverNames);
    else
        solver_names = ["sds", "dse", "sto", "rnd", "sham"];
    end
    [found, solver_indices] = ismember(solver_order, solver_names);
    if ~all(found)
        error('Missing solver results for: %s', strjoin(solver_order(~found), ', '));
    end
    C = C(:, solver_indices);
    [P, S] = size(C);

    % 2) build H(e,p,s) and N(p) for data & performance profiles
    N = dims;
    maxEvals = max(N*function_budget);
    H = inf(maxEvals, P, S);

    for p = 1:P
        for s = 1:S
            fvec = C{p,s}; fvec = fvec(:);
            L = min(length(fvec), maxEvals);
            % best = inf;
            for e = 1:L
                % best = min(best, fvec(e));  % sort by best so far, not by last fvec(e)
                % H(e,p,s) = best;
                H(e,p,s) = fvec(e);  % use last fvec(e), not by best so far
            end
        end
    end

    % 3) data profile
    gate = 10^(-tol);
    figure; hl1 = data_profile_original(H, N, gate);
    title(sprintf('Data Profile (\\tau = %g)', gate));
    legend(hl1, legend_names, 'Location','southeast');
    xlabel('Groups of (n+1) evaluations of f');
    ylabel('Proportion of problems solved');
    grid on;
    saveas(gcf, fullfile(folder_save, sprintf('dp_tol%d.eps', tol)), 'epsc');
    close;

    % 4) performance profile (log-scale)
    figure; hl2 = perf_profile_original(H, gate, 1);
    title(sprintf('Performance Profile (\\tau = %g)', gate));
    legend(hl2, legend_names, 'Location','southeast');
    xlabel('Ratio of evaluations of f');
    ylabel('Proportion of problems solved');
    grid on;
    saveas(gcf, fullfile(folder_save, sprintf('pp_tol%d.eps', tol)), 'epsc');
    close;
end
end
