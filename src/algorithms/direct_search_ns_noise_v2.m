function [x_k, f_k, f_k_vec] = direct_search_ns_noise_v2(problem, function_budget, s_options)
alpha_start = s_options.alpha_start;
gamma = s_options.gamma;
theta_c = s_options.theta_c;
theta_exp = s_options.theta_exp;
q = s_options.q;
correlation = s_options.correlation;
eps = s_options.eps;
noise_coeff = s_options.noise_coeff;
alpha_threshold = s_options.alpha_threshold;
theta_c1 = s_options.theta_c1;
theta_exp1 = s_options.theta_exp1;

sample_num_exp = sample_computer(q, correlation, eps);

f_k_vec = inf(function_budget, 1);
x_k = problem.startp;
dim = problem.dim;
alpha_k = alpha_threshold;
alpha_k_vec = ones(1, 2*dim) * alpha_start;
mat_basis = zeros(dim, 2*dim);
for i = 1:dim
    mat_basis(i, 2*i - 1) = 1;
    mat_basis(i, 2*i) = -1;
end
f_obj = problem.fhandle;
f_k = f_obj(x_k);

% track budget using the sampling rule
[~, m0] = noisegenerator(alpha_start, sample_num_exp, correlation, noise_coeff);
tot_f_eval_old = 1;
tot_f_eval = max(1, ceil(m0));

i = 1;
i_basis_curr = 1;
phase_two = false;
f_k_vec(1) = f_k;
while tot_f_eval < function_budget
    for j = tot_f_eval_old+1: min(tot_f_eval, function_budget)
        f_k_vec(j) = f_k;
    end
    if max(alpha_k_vec) < alpha_threshold
        phase_two = true;
    end
    if phase_two
        g_k = randn(dim, 1); g_k = g_k/norm(g_k);
        alpha_k_curr = alpha_k;
    else
        g_k = mat_basis(:, i_basis_curr);
        alpha_k_curr = alpha_k_vec(i_basis_curr);
    end
    [noise, num_sample] = noisegenerator(alpha_k_curr, sample_num_exp, correlation, noise_coeff);
    x_kp = x_k + alpha_k_curr*g_k;
    f_kp = f_obj(x_kp);
    i_basis_prev = i_basis_curr;
    if f_k - f_kp + noise > gamma * alpha_k^q
        x_k = x_kp;
        f_k = f_kp;
        if phase_two
            alpha_k_curr = alpha_k_curr * theta_exp;
        else
            alpha_k_curr = alpha_k_curr * theta_exp1;
        end
    else
        if phase_two
            alpha_k_curr = alpha_k_curr * theta_c;
        else
            alpha_k_curr = alpha_k_curr * theta_c1;
        end
        i_basis_curr = mod((i_basis_curr + 1), 2*dim);
        if i_basis_curr == 0, i_basis_curr = 2*dim; end
        while alpha_k_vec(i_basis_curr) < alpha_threshold
            i_basis_curr = mod((i_basis_curr + 1), 2*dim);
            if i_basis_curr == 0, i_basis_curr = 2*dim; end
            if i_basis_curr == i_basis_prev
                break
            end
        end
    end
    if phase_two
        alpha_k = alpha_k_curr;
    else
        alpha_k_vec(i_basis_prev) = alpha_k_curr;
    end
    tot_f_eval_old = tot_f_eval;
    tot_f_eval = tot_f_eval + 2 * num_sample;
    i = i + 1;
end
for j = tot_f_eval_old + 1: function_budget
    f_k_vec(j) = f_k;
end
