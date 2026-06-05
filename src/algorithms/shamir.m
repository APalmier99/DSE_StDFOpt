function [x, fx, f_hist, grad_norm_hist, info] = shamir(fhandle, x0, options)
%KORNOWSKI_SHAMIR_NOISY
% Implementation of Kornowski-Shamir Algorithm 1 + 2.
%
%   num_eval = 1;
%   f_hist(1) = fhandle(x0);
%
% Then each noisy function evaluation used inside the algorithm increments
% num_eval. The true value fhandle(x) used for profiling is NOT counted.
%
% Main options:
%   options.function_budget  total budget for profile indexing
%   options.N                length of f_hist
%   options.mu or rho         smoothing radius rho
%   options.h  or eta         stepsize eta
%   options.nu               Goldstein neighborhood parameter
%   options.D                clipping radius
%   options.batch_size       number of two-point estimators averaged
%   options.noise_coeff      additive Gaussian noise std, as in baseline
%   options.output_rule      'last' or 'paper'
% 
% Exact paper output:
%   options.output_rule = 'paper';
% Alternatively, to use the last iterate:
%   options.output_rule = 'last' 
% 
% The paper output samples uniformly among block averages
% x_k = (1/M) sum_{m=1}^M z_{(k-1)M+m},
% where M = floor(nu/D).

if nargin < 3 || isempty(options)
    options = struct();
end

x0 = x0(:);
d = numel(x0);

% ------------------------------------------------------------
% Options
% ------------------------------------------------------------
function_budget = get_option(options, 'function_budget', get_option(options, 'N', 1000));
N               = get_option(options, 'N', function_budget);

rho = get_option(options, 'rho', get_option(options, 'mu', 1e-4));
eta = get_option(options, 'eta', get_option(options, 'h',  1e-3));

nu = get_option(options, 'nu', rho);
D  = get_option(options, 'D',  nu);

batch_size = get_option(options, 'batch_size', get_option(options, 'k', 1));
output_rule = get_option(options, 'output_rule', 'last');

function_budget = max(1, floor(function_budget));
N               = max(1, floor(N));
batch_size      = max(1, floor(batch_size));

if rho <= 0
    error('kornowski_shamir_noisy:badRho', 'rho/mu must be positive.');
end
if eta <= 0
    error('kornowski_shamir_noisy:badEta', 'eta/h must be positive.');
end
if nu <= 0
    error('kornowski_shamir_noisy:badNu', 'nu must be positive.');
end
if D <= 0
    error('kornowski_shamir_noisy:badD', 'D must be positive.');
end
if D > nu
    error('kornowski_shamir_noisy:badD', ...
        'Need D <= nu so that M = floor(nu/D) >= 1.');
end

evals_per_iter = 2 * batch_size;

% Same convention as baseline:
% f_hist(1) corresponds to the initial point.
num_eval = 1;
last_recorded_eval = 1;

max_budget_iters = floor((function_budget - num_eval) / evals_per_iter);
max_budget_iters = max(0, max_budget_iters);

T = get_option(options, 'T', get_option(options, 'max_iter', max_budget_iters));
T = floor(T);
T = max(0, min(T, max_budget_iters));

% Paper block size
M = floor(nu / D);
M = max(1, M);

K_plan = floor(T / M);
x_candidates = zeros(d, max(1, K_plan));
num_candidates = 0;
z_block_sum = zeros(d, 1);
z_block_count = 0;

% ------------------------------------------------------------
% Initialize history
% ------------------------------------------------------------
f_hist = inf(N, 1);
grad_norm_hist = inf(N, 1);

x_prev = x0;
Delta = zeros(d, 1);

fx0 = fhandle(x_prev);
f_hist(1) = fx0;
grad_norm_hist(1) = NaN;

fx_current = fx0;
last_grad_norm = NaN;

actual_T = 0;

% ------------------------------------------------------------
% Main loop
% ------------------------------------------------------------
for t = 1:T

    if num_eval + evals_per_iter > function_budget
        break;
    end

    % Algorithm 1:
    % x_t = x_{t-1} + Delta_t
    % z_t = x_{t-1} + s_t Delta_t
    s = rand();

    x_curr = x_prev + Delta;
    z = x_prev + s * Delta;

    % Algorithm 2: averaged two-point spherical estimator
    [g, used_evals] = ks_grad_estimator(fhandle, z, rho, batch_size, options);

    num_eval = num_eval + used_evals;

    % Clipped displacement update:
    % Delta_{t+1} = min(1, D / ||Delta_t - eta g_t||) ...
    %               * (Delta_t - eta g_t)
    Delta_trial = Delta - eta * g;
    Delta_trial_norm = norm(Delta_trial);

    if Delta_trial_norm > D
        Delta = (D / Delta_trial_norm) * Delta_trial;
    else
        Delta = Delta_trial;
    end

    % Move to x_t
    x_prev = x_curr;
    actual_T = actual_T + 1;

    % --------------------------------------------------------
    % Baseline-compatible profiling
    % --------------------------------------------------------
    % This true function call is for plotting only and is not counted.
    fx_current = fhandle(x_prev);
    last_grad_norm = norm(g);

    profile_idx = min(num_eval, N);

    if profile_idx > last_recorded_eval
        f_hist(last_recorded_eval + 1 : profile_idx) = fx_current;
        grad_norm_hist(last_recorded_eval + 1 : profile_idx) = last_grad_norm;
    end

    last_recorded_eval = profile_idx;

    % --------------------------------------------------------
    % Paper candidate construction
    % --------------------------------------------------------
    z_block_sum = z_block_sum + z;
    z_block_count = z_block_count + 1;

    if z_block_count == M
        num_candidates = num_candidates + 1;

        if num_candidates <= size(x_candidates, 2)
            x_candidates(:, num_candidates) = z_block_sum / M;
        end

        z_block_sum(:) = 0;
        z_block_count = 0;
    end
end

% ------------------------------------------------------------
% Output rule
% ------------------------------------------------------------
% For exact paper output, use a uniformly sampled block average.
if strcmpi(output_rule, 'paper') && num_candidates >= 1
    selected_candidate = randi(num_candidates);
    x = x_candidates(:, selected_candidate);
else
    % Alternatively, use the last iterate.
    selected_candidate = NaN;
    x = x_prev;
end

fx = fhandle(x);

% Fill remaining history as in baseline.
if last_recorded_eval < N
    f_hist(last_recorded_eval + 1 : N) = fx;
    grad_norm_hist(last_recorded_eval + 1 : N) = last_grad_norm;
end

% ------------------------------------------------------------
% Info
% ------------------------------------------------------------
info = struct();
info.rho = rho;
info.eta = eta;
info.nu = nu;
info.D = D;
info.M = M;
info.K_plan = K_plan;
info.num_candidates = num_candidates;
info.selected_candidate = selected_candidate;
info.output_rule = output_rule;
info.batch_size = batch_size;
info.evals_per_iter = evals_per_iter;
info.function_budget = function_budget;
info.N = N;
info.num_eval = num_eval;
info.T_requested = T;
info.T_actual = actual_T;
info.last_iterate = x_prev;
info.last_displacement = Delta;
info.last_grad_norm = last_grad_norm;
info.candidates = x_candidates(:, 1:num_candidates);

end

% ========================================================================
% Two-point spherical gradient estimator
% ========================================================================
function [g, used_evals] = ks_grad_estimator(fhandle, z, rho, batch_size, options)

d = numel(z);
g = zeros(d, 1);

for j = 1:batch_size

    w = randn(d, 1);
    w_norm = norm(w);

    while w_norm == 0
        w = randn(d, 1);
        w_norm = norm(w);
    end

    w = w / w_norm;

    f_plus  = noisy_value(fhandle, z + rho * w, options);
    f_minus = noisy_value(fhandle, z - rho * w, options);

    g = g + (d / (2 * rho)) * (f_plus - f_minus) * w;
end

g = g / batch_size;
used_evals = 2 * batch_size;

end

% ========================================================================
% Noisy oracle, matching baseline convention
% ========================================================================
function y = noisy_value(fhandle, x, options)

y = fhandle(x);

noise_coeff = get_option(options, 'noise_coeff', 0);

if noise_coeff ~= 0
    y = y + noise_coeff * randn();
end

end

% ========================================================================
% Option getter
% ========================================================================
function value = get_option(options, name, default_value)

if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = default_value;
end

end