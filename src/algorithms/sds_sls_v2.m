function [x_k, f_k, f_k_vec] = sds_sls_v2(problem, options)
% DSE-DD style algorithm with stochastic line search and Sobol directions.

% Unpack problem and options
function_budget   = options.function_budget;
alpha_start       = options.alpha_start;
alpha_threshold   = options.alpha_threshold;
dim               = problem.dim;

theta_dse         = options.theta_dse;
gamma_c           = options.gamma_dse_c;
gamma_e           = options.gamma_dse_e;
gamma_c1          = options.gamma_dse_c1;
gamma_e1          = options.gamma_dse_e1;

q                 = options.q;
correlation       = options.corr;
eps               = options.eps;
noise_coeff       = options.noise_coeff;
sample_num_exp    = sample_computer(q, correlation, eps);

alpha_k           = alpha_threshold;
alpha_k_vec       = ones(1, 2*dim) * alpha_start;

% Build basis directions for Phase 1
mat_basis = zeros(dim, 2*dim);
for i = 1:dim
    mat_basis(i, 2*i - 1) = 1;
    mat_basis(i, 2*i)     = -1;
end

% Initialize
f_k_vec      = inf(function_budget, 1);
x_k          = problem.startp(:);
f_obj        = problem.fhandle;
f_k          = f_obj(x_k);

% budget accounting
[~, m0] = noisegenerator(alpha_start, sample_num_exp, correlation, noise_coeff);
tot_f_eval_old = 1;
tot_f_eval = max(1, ceil(m0));
for j = 1: min(tot_f_eval, function_budget)
    f_k_vec(j) = f_k;
end

i_basis_curr = 1;
phase_two    = false;

% Sobol sequence for phase 2 (fallback to Gaussian if unavailable)
try
    p = sobolset(dim,'Skip',1e3,'Leap',1e2);
    p = scramble(p,'MatousekAffineOwen');
    s = qrandstream(p);
catch
    s = [];
end

ls_threshold = options.ls_threshold;

while tot_f_eval < function_budget

    for j = tot_f_eval_old+1: min(tot_f_eval, function_budget)
        f_k_vec(j) = f_k;
    end

    if max(alpha_k_vec) < alpha_threshold
        phase_two = true;
    end

    if phase_two
        if ~isempty(s)
            u = rand(s,1,dim)';
            z = sqrt(2) * erfinv(2*u - 1);
            g_k = z / norm(z);
        else
            z = randn(dim,1);
            g_k = z / max(1e-15, norm(z));
        end
        alpha_k_curr = alpha_k;
    else
        g_k          = mat_basis(:, i_basis_curr);
        alpha_k_curr = alpha_k_vec(i_basis_curr);
    end

    if alpha_k_curr <= ls_threshold
        [delta_out, delta_next, num_sample] = stochastic_linesearch_sampled( ...
            x_k, alpha_k_curr, g_k, theta_dse, [gamma_c, gamma_e], ...
            f_obj, sample_num_exp, correlation, noise_coeff);

        if delta_out > 0
            x_k = x_k + delta_out * g_k;
            f_k = f_obj(x_k);
        end
        alpha_next = delta_next;

    else
        [noise, num_sample] = noisegenerator(alpha_k_curr, sample_num_exp, correlation, noise_coeff);
        num_sample = 2*num_sample;
        x_kp = x_k + alpha_k_curr*g_k;
        f_kp = f_obj(x_kp);
        if f_k - f_kp + noise > theta_dse * alpha_k^q
            x_k = x_kp;
            f_k = f_kp;
            if phase_two
                alpha_next = alpha_k_curr * gamma_e;
            else
                alpha_next = alpha_k_curr * gamma_e1;
            end
        else
            if phase_two
                alpha_next = alpha_k_curr * gamma_c;
            else
                alpha_next = alpha_k_curr * gamma_c1;
            end
        end
    end

    if phase_two
        alpha_k = alpha_next;
    else
        alpha_k_vec(i_basis_curr) = alpha_next;
        % advance basis index
        prev_idx = i_basis_curr;
        i_basis_curr = i_basis_curr + 1;
        if i_basis_curr > 2*dim, i_basis_curr = 1; end
        while alpha_k_vec(i_basis_curr) < alpha_threshold && i_basis_curr ~= prev_idx
            i_basis_curr = i_basis_curr + 1;
            if i_basis_curr > 2*dim, i_basis_curr = 1; end
        end
    end

    tot_f_eval_old = tot_f_eval;
    tot_f_eval = tot_f_eval + num_sample;
end

for j = tot_f_eval_old + 1: function_budget
    f_k_vec(j) = f_k;
end

end

function [delta_out, delta_next, total_samples] = stochastic_linesearch_sampled( ...
    x, delta, d, theta, gamma, f_obj, exp_samples, corr, noise_coeff)

gamma_c = gamma(1);
gamma_e = gamma(2);

[noise1, m1] = noisegenerator(delta, exp_samples, corr, noise_coeff);
fx_hat       = f_obj(x) + noise1;
[noise2, m2] = noisegenerator(delta, exp_samples, corr, noise_coeff);
fxd_hat      = f_obj(x + delta * d) + noise2;
total_samples = m1 + m2;

if fxd_hat > fx_hat - theta * delta^2
    delta_out  = 0;
    delta_next = gamma_c * delta;
    return;
end

while fxd_hat <= fx_hat - theta * delta^2
    delta = delta * gamma_e;
    [noise2, m2] = noisegenerator(delta, exp_samples, corr, noise_coeff);
    fxd_hat      = f_obj(x + delta * d) + noise2;
    total_samples = total_samples + m2;
end

delta_out  = delta / gamma_e;
delta_next = delta_out;
end
