How to integrate the sampling rule and comparable output with StoMADS

Goal
- Make StoMADS use a sample size m(h) ≍ ceil(C^2 · h^{-4}) when computing noisy estimates, matching the unified oracle protocol used in this repo.
- Make StoMADS return the same output shape as the other algorithms, so the harness can consume it directly (struct with fields `x_star`, `f_star`, `f_vec`).

What you replace
- Replace StoMADS’ original `stomads_algorithm.m` with the comparable version provided in this repository. This modified file enforces the sampling rule and returns the expected struct.

Where to find the modified file here
- `third_party/StoMADS/StoMADS_Main_Files/stomads_algorithm.m`

Steps (in your local StoMADS copy)
1) Overwrite the StoMADS driver
   - Copy this repo’s `third_party/StoMADS/StoMADS_Main_Files/stomads_algorithm.m` over the file with the same name in your StoMADS installation (usually under `.../stomads/StoMADS_Main_Files/`).
   - This version does two things:
     - Uses m(h) = ceil(noise_coeff^2 · h^{-4}) (or ceil(h^{-4}) if `noise_coeff=0`) when computing Monte Carlo estimates.
     - Returns a struct with fields `x_star`, `f_star`, `f_vec` (best point, best value, and true f-history), matching the rest of this codebase.

   Path precedence tip: if you prefer not to modify your StoMADS install, ensure MATLAB finds the modified file first by putting this repo’s `third_party/StoMADS/StoMADS_Main_Files/` ahead of any installed StoMADS paths (e.g., call `addpath(genpath(...))` for this folder before others).

2) Update the Monte Carlo estimator
   - Ensure `StoMADS_Main_Files/MCestimate.m` adds Gaussian noise with variance 1/m and supports a `noise_coeff` multiplier. Minimal change:

        z = my_fun(x) + randn() * (sample_size^(-0.5)) * stomads_option.noise_coeff;

   - Alternatively, replace the function with this version:

        function z = MCestimate(x, sample_size, my_fun, probspecs, stomads_option)

            global nfEval fEval_History fEval_Stats nfEvalExceeded

            y = zeros(sample_size, 1);
            nfEval = min([nfEval + sample_size, stomads_option.MaxFuncEval + 1]);

            if stomads_option.noise_coeff == 0
                z = my_fun(x);
            else
                z = my_fun(x) + normrnd(0, sample_size^(-0.5)) * stomads_option.noise_coeff;
            end

            if nfEval > stomads_option.MaxFuncEval
                nfEvalExceeded = 1;
            end

            if nfEvalExceeded == 1
                z = 1i;
            end

        end

3) Put files where the wrapper expects them
   - Copy your updated `StoMADS_Main_Files/` into this repo at:
     `third_party/StoMADS/StoMADS_Main_Files/`.
   - The wrapper `src/algorithms/stomads_wrapper.m` adds this folder to the path automatically during runs.

4) Optional: Outputs and logs
   - Some StoMADS revisions expect `stomads_option.SaveLocation` for text outputs. The wrapper sets a default under `results/stomads_txt/`.

Licensing note
- Keep the StoMADS license with the third-party code. This repository does not redistribute StoMADS; you point the wrapper at your local copy with the modified files.
