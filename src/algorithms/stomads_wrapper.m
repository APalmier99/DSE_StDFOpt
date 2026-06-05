function [x_star, f_star, f_vec] = stomads_wrapper(problem, function_budget, s_options)
% Thin wrapper around StoMADS to match our interface and noise model.

x_start = problem.startp;
f_obj_true = @(x)problem.fhandle(x);
f_obj_noisy = @(x) f_obj_true(x);

stomads_option.experiments = 0;

stomads_option.Anisotropy = 0;
stomads_option.AnisoParam = 0.1;
stomads_option.DisplayOutputs = 0;
stomads_option.DisplaySolution = 0;
stomads_option.FixSeed = 0;
stomads_option.GammaEpsilon = s_options.gamma;
stomads_option.HistoryFile = 0;
stomads_option.InitPollSize = s_options.alpha_start;
stomads_option.LowerBounds = [];
stomads_option.MaxFuncEval = function_budget;
stomads_option.MaxNumberIters = 1000;
stomads_option.MeshType = 0;
stomads_option.OrthoN_PlusOne = 1;
stomads_option.SampleSize = 1;
stomads_option.SeedValue = [];
stomads_option.SolutionFile = 0;
stomads_option.StatsFile = 0;
stomads_option.Tau = s_options.theta_c;
stomads_option.UpperBounds = [];
stomads_option.UsePreviousSamples = 0;
stomads_option.warning = 0;

% Some StoMADS versions expect an output folder path in SaveLocation.
% Provide a default under this repo's results directory to avoid errors.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
stomads_option.SaveLocation = fullfile(repo_root, 'results', 'stomads_txt');
if ~exist(stomads_option.SaveLocation,'dir')
	try mkdir(stomads_option.SaveLocation); catch, end
end

stomads_option.MaxFuncEval = function_budget;
stomads_option.SampleSize = 1;
stomads_option.UsePreviousSamples = 0;
stomads_option.noise_coeff = s_options.noise_coeff;

probspecs = struct();
probspecs.Dimension = length(x_start);
probspecs.myfun_true = f_obj_true;

% StoMADS expects its own internal functions on path
sto_root = fullfile(fileparts(mfilename('fullpath')),'..','..','third_party','StoMADS');
addpath(genpath(sto_root));

X = stomads_algorithm(f_obj_noisy, x_start, stomads_option, probspecs);
x_star = X.x_star;
f_star = X.f_star;
f_vec = X.f_vec;
end
