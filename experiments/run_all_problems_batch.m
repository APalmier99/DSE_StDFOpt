function run_all_problems_batch(num_runs)
    if nargin<1, num_runs = 1; end

    % Add this repo to the MATLAB path (one level up)
    setup_paths();
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    output_folder = fullfile(repo_root,'results');
    if ~exist(output_folder,'dir'), mkdir(output_folder); end

    % Load problems list (needs a problems.mat snapshot)
    problems = [];
    candidateMats = { ...
        fullfile(repo_root,'problems','problems.mat'), ...
        fullfile(repo_root,'problems.mat'), ...
        'problems.mat' ...
    };
    for c = 1:numel(candidateMats)
        if exist(candidateMats{c}, 'file') == 2
            S = load(candidateMats{c});
            if isfield(S,'problems')
                problems = S.problems;
            end
            break;
        end
    end
    if isempty(problems)
        fprintf('No problems.mat found. Please provide a problems.mat file under the problems/ folder or the repo root.\n');
        return;
    end
    problems = normalize_problems(problems); % ensure cell array with required fields
    P = numel(problems);

    solverNames = {'sds','dse','sto','rnd', 'sham'}; ns = numel(solverNames);

    totalRows  = P * num_runs;
    resultcell = cell(totalRows, ns);
    dimvec     = zeros(totalRows,1);

    for k = 1:P
        prob = problems{k};
        for run_idx = 1:num_runs
            fprintf('\n--- Problem #%d (%s), run %d/%d ---\n', k, prob.name, run_idx, num_runs);

            run_all_algorithms_on_problem(prob);
            S = load(fullfile(output_folder,'results_placeholder.mat'));

            r = (k-1)*num_runs + run_idx;
            resultcell{r,1} = S.fvec_sds;
            resultcell{r,2} = S.fvec_dse;
            resultcell{r,3} = S.fvec_sto;
            resultcell{r,4} = S.fvec_rnd;
            resultcell{r,5} = S.fvec_sham;
            dimvec(r) = prob.dim + 1;
        end
    end

    filename = fullfile(output_folder, sprintf('results_num_runs_%d.mat', num_runs));
    save(filename, 'resultcell','dimvec','solverNames','-v7.3');
    fprintf('Saved aggregated results to: %s\n', filename);
end

function problems = normalize_problems(problems)
    % Ensure problems is a cell array of structs with fields:
    % name (char), dim (scalar), fhandle (function handle), startp (dimx1)

    % Convert struct array to cell if needed
    if isstruct(problems)
        problems = arrayfun(@(s) s, problems, 'UniformOutput', false);
    end
    if ~iscell(problems)
        error('normalize_problems:badType','Expected problems to be cell or struct array.');
    end

    for i = 1:numel(problems)
        pr = problems{i};
        % If provided as a struct with fname string, resolve to fhandle
        if ~isfield(pr,'fhandle') && isfield(pr,'fname')
            try
                pr.fhandle = str2func(pr.fname);
            catch
                error('normalize_problems:badFname','Invalid function name in problems{%d}.fname', i);
            end
        end
        % Default start point if missing and dim available
        if ~isfield(pr,'startp') && isfield(pr,'dim') && isnumeric(pr.dim)
            pr.startp = zeros(pr.dim,1);
        end
        % Basic validation
        req = {'name','dim','fhandle','startp'};
        missing = req(~isfield(pr, req));
        if ~isempty(missing)
            error('normalize_problems:missingFields','Problem %d missing fields: %s', i, strjoin(missing, ', '));
        end
        problems{i} = pr;
    end
end
