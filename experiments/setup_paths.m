function setup_paths()
% Add root, src, problems, third_party to path for a fresh MATLAB session.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

% Put experiments folder on the path so run_all_problems_batch resolves here
addpath(fullfile(root,'experiments'));

% Add src but exclude legacy 'old_scripts' to avoid name clashes and wrong outputs
srcPath = genpath(fullfile(root,'src'));
% Remove any subpaths containing 'old_scripts'
if ~isempty(srcPath)
	parts = strsplit(srcPath, pathsep);
	for i = 1:numel(parts)
		p = parts{i};
		if ~isempty(p)
			if contains(p, ['old_scripts' filesep]) || endsWith(p, 'old_scripts')
				% skip adding this legacy folder
			else
				addpath(p);
			end
		end
	end
end

addpath(genpath(fullfile(root,'problems')));
addpath(genpath(fullfile(root,'third_party')));
end
