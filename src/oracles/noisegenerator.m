function [noise, num_sample] = noisegenerator(alpha, exp_samples, corr, noise_coeff)
% Unified noise model: average of num_sample noisy evals; here we return
% an equivalent Gaussian perturbation and the implied sample count.
% num_sample scales like ceil(C^2 * alpha^{exp_samples}).

if noise_coeff == 0
    noise = 0;
    num_sample = max(1, ceil(alpha^(exp_samples)));
    return
end

num_sample = max(1, ceil(noise_coeff^2 * alpha^(exp_samples)));

if corr==false
    std = (2*num_sample)^(-0.5) * noise_coeff;
else
    std = num_sample^(-0.5) * alpha * noise_coeff;
end

noise = std * randn();

end
