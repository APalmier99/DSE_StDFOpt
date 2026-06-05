function [x, fx, f_hist, grad_norm_hist] = random_nonconvex_search_noisy(fhandle, x0, options)
% Random nonconvex search baseline with noisy finite differences.

q = options.q; 
corr = options.corr;
eps = options.eps;
noise_coeff = options.noise_coeff;
function_budget = options.function_budget;

B = options.B; mu = options.mu; h = options.h; N = options.N;

x = x0(:);
n = numel(x);
f_hist         = inf(N,1);
grad_norm_hist = inf(N,1);

num_eval = 1; num_eval_old = num_eval;

% record true f for profiling
f_hist(1) = fhandle(x);
grad_norm_hist(1) = NaN;

L    = chol(B,'lower');
Linv = L \ eye(n);

% exp_samp = sample_computer(q, corr, eps);

for k = 1:N
    u = Linv * randn(n,1);

    % Using simple Gaussian noise here for speed; budget still enforced
    noise1 = noise_coeff*randn(1,1); m1 = 1;
    if num_eval + m1 > function_budget, break; end
    fx_noisy = fhandle(x) + noise1; num_eval = num_eval + m1;

    noise2 = noise_coeff*randn(1,1); m2 = 1;
    if num_eval + m2 > function_budget, break; end
    fpu_noisy = fhandle(x + mu*u) + noise2; num_eval  = num_eval + m2;

    df  = (fpu_noisy - fx_noisy)/mu;
    gmu = df*(B*u);
    x   = x - h*(B\gmu);

    fx = fhandle(x);
    f_hist(num_eval_old+1: min(num_eval, N)) = fx;
    grad_norm_hist(num_eval_old+1: num_eval) = sqrt(gmu'*(B\gmu));
    num_eval_old = num_eval;
end

for j = num_eval + 1: N
    f_hist(j) = fx;
end

end
