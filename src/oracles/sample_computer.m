function sample_num_exp = sample_computer(q, corr, eps)
% Returns exponent for the sampling rule m(h) ~ h^{exp}
% Default regime used in our experiments yields exp = -4.
if q > 2
    sample_num_exp = -q^2 - eps;
elseif corr == true
    sample_num_exp = 2 - 2*q - eps;
elseif q == 1.5
    sample_num_exp = - 2*q;
else
    sample_num_exp = -4 - eps;
end
