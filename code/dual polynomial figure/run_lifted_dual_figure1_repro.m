function result = run_lifted_dual_figure1_repro()
%RUN_LIFTED_DUAL_FIGURE1_REPRO Reproduce the lifted dual-polynomial test.
%%
% Paper: "A Mathematical Theory of Near-field Super-Resolution."
% Author: Sajad Daei, KTH Royal Institute of Technology.
% Email: sajado@kth.se
%
% This code was written by Sajad Daei and is provided as supplementary
% reproducibility material for:
%
%   S. Daei, G. Fodor, and M. Skoglund,
%   "A Mathematical Theory of Near-field Super-Resolution."
%
% Permission is granted to download, inspect, and run the code for academic,
% non-commercial peer-review, verification, and reproducibility purposes.
%
% Any publication, preprint, report, presentation, or software work that uses
% this code, adapts this code, or uses numerical output derived from this code
% must cite the paper above.
%
% Redistribution, modification, incorporation into other software, commercial
% use, or public release of modified versions requires prior written permission
% from Sajad Daei.
%
% This code is provided "as is", without warranty of any kind.
%
%
% This function reproduces the finite Jacobi--Anger lifted SDP experiment
% reported in the paper.  It generates Fresnel-model observations, solves
% the robust lifted dual SDP, extracts the two range--angle support points,
% refits their complex amplitudes, and exports machine-readable diagnostics
% together with the paper figure.
%
% Requirements:
%   - MATLAB R2019b or later
%   - CVX 2.2
%   - SDPT3 4.0 (selected explicitly below)
%
% Output:
%   result : structure containing parameters, diagnostics, support estimates,
%            solver status, and paths to generated artifacts.
%
% Run from any directory with
%   result = run_cpam_lifted_dual_figure1_repro();
this_file = mfilename('fullpath');
if isempty(this_file)
    repro_root = pwd;
else
    repro_root = fileparts(this_file);
end
results_dir = fullfile(repro_root, 'results');
figures_dir = fullfile(repro_root, 'figures');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

diary(fullfile(results_dir, 'cpam_repro_log.txt'));
cleanupObj = onCleanup(@() diary('off'));

% Fixed random seed.  It controls the optional synthetic noise realization.
random_seed = 12;
rng(random_seed, 'twister');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) BASIC USER PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_r = 16;               % number of array elements
fc  = 1e9;              % carrier frequency (Hz)
c0  = 3e8;              % value used in the reported experiment (m/s)

lambda = c0/fc;
k_fac  = 2*pi/lambda;

% Jacobi-Anger truncation orders.
% I1 controls angular Bessel expansion.
% I2 controls range/curvature Bessel expansion.
I1 = 20;
I2 = 8;

% Discrete range grid and the paper-reported two-source scene.  These values
% are explicit so that the archive cannot silently drift to a different run.
N_d = 8;
scene.target_bins = [3; 6];
scene.theta_rad = [0.30*pi; 0.75*pi];
scene.coefficients = [ ...
    -1.00146807223807 - 0.343404412876242*1i; ...
     1.15019479172557 + 0.133020150213613*1i];
nspikes = numel(scene.target_bins);

% Noise level.
sigma_noise = 0;
delta_prob  = 1e-1;

% Use false to validate the lifted JA recovery model.
% Use true to test Fresnel-generated data with model mismatch included in eta.
% The lifted SDP is a computational surrogate and is not part of the QPAC
% exact-recovery theorem.  Set true only for the Fresnel-reference experiment,
% where eta includes an explicit analytic truncation/model-mismatch bound
% evaluated in ordinary double precision.
USE_FRESNEL_REFERENCE = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) DESIGN d, r_min, r_max, AND r_list FROM I1 AND I2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

range_cfg.gamma_z1 = 0.50;
% z1_max = k D_ap <= gamma_z1 I1

range_cfg.gamma_z2 = 0.12;
% z2(r_min) = k D_ap^2/(4r_min) <= gamma_z2 I2

range_cfg.z2_min = 0.08;
% z2(r_max) >= z2_min
% Smaller z2_min gives larger r_max but weaker range identifiability.

range_cfg.reactive_safety = 1.05;

range_cfg.fraunhofer_multiple = 10.0;
% Allows r_max beyond classical Fraunhofer distance if curvature is still
% above z2_min.

% Good choices:
%   'log_curvature'    : recommended for wide range intervals.
%   'linear_curvature' : uniform curvature resolution.
%   'log_range'        : log-spaced range.
%   'linear_range'     : uniform range.
range_cfg.grid_type = 'log_curvature';

[d, D_ap, r_min, r_max, r_list, range_info] = ...
    designRangeGridFromI1I2(N_r, lambda, I1, I2, N_d, range_cfg);

n_idx = (0:N_r-1).';

fprintf('\nArray and truncation-safe range diagnostics:\n');
disp(struct2table(range_info));

fprintf('\nDerived parameters:\n');
fprintf('lambda = %.6e m\n', lambda);
fprintf('d      = %.6e m\n', d);
fprintf('D_ap   = %.6e m\n', D_ap);
fprintf('r_min  = %.6e m\n', r_min);
fprintf('r_max  = %.6e m\n', r_max);
fprintf('span   = %.6e m\n', r_max-r_min);
fprintf('kD     = %.6f\n', k_fac*D_ap);
fprintf('classical radiating Fresnel interval = [%.6e, %.6e] m\n', ...
    range_info.R_reactive, range_info.R_fraunhofer);
fprintf('complete range grid is inside that interval: %d\n', ...
    range_info.full_grid_in_radiating_fresnel);

range_grid_table = table((1:N_d).', r_list(:), ...
    (k_fac*D_ap^2./(4*r_list(:))), ...
    'VariableNames', {'RangeBin','RangeM','CurvatureArgumentZ2'});
writetable(range_grid_table, fullfile(results_dir, 'range_grid.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) JA INDEX SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_vals = -I1:I1;
q_vals = -I2:I2;

L0 = length(l_vals);
Q0 = length(q_vals);
P  = L0*Q0;

fprintf('\nChosen JA truncation:\n');
fprintf('I1 = %d, I2 = %d, P = %d\n', I1, I2, P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) TRUE TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_bins = scene.target_bins(:);
theta_true = scene.theta_rad(:);
coef = scene.coefficients(:);

if any(target_bins < 1) || any(target_bins > N_d) || ...
        any(target_bins ~= round(target_bins))
    error('scene.target_bins must contain valid integer indices in 1:N_d.');
end
if length(theta_true) ~= nspikes || length(coef) ~= nspikes
    error('The target-bin, angle, and coefficient vectors must have equal length.');
end
if any(theta_true <= 0) || any(theta_true >= pi)
    error('All source angles must lie strictly inside (0,pi).');
end

r_true = r_list(target_bins).';
r_true = r_true(:);

fprintf('\nPaper-reported source scene:\n');
scene_table = table((1:nspikes).', target_bins, r_true, theta_true, ...
    real(coef), imag(coef), ...
    'VariableNames', {'Source','RangeBin','RangeM','AngleRad', ...
    'CoefficientReal','CoefficientImag'});
disp(scene_table);
writetable(scene_table, fullfile(results_dir, 'scene.csv'));

true_sources_in_radiating_fresnel = ...
    all(r_true >= range_info.R_reactive) && ...
    all(r_true <= range_info.R_fraunhofer);
true_sources_beyond_fraunhofer = all(r_true > range_info.R_fraunhofer);

fprintf('all true sources in classical radiating Fresnel interval: %d\n', ...
    true_sources_in_radiating_fresnel);
fprintf('all true sources beyond the Fraunhofer boundary: %d\n', ...
    true_sources_beyond_fraunhofer);

if ~true_sources_in_radiating_fresnel
    fprintf(['Interpretation: "Fresnel" denotes the quadratic Fresnel ', ...
        'forward model here; it does not claim that the active ranges lie ', ...
        'inside the classical radiating near-field interval.\n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) NATURAL TRUNCATION DIAGNOSTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I1_main_angle    = floor(k_fac*D_ap);
I2_main_distance = floor(k_fac*D_ap^2 ./ (4*r_true(:)));

fprintf('\nNatural truncation diagnostics:\n');
fprintf('I1_main_angle ~= floor(kD) = %d\n', I1_main_angle);

source_idx = (1:nspikes).';

diag_table = table(source_idx, r_true(:), theta_true(:), ...
                   I2_main_distance(:), ...
    'VariableNames', {'source','r_true','theta_true','I2_main_distance'});

disp(diag_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6) ANGULAR LIFTING: v_theta(theta) = S v_normal(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[v_theta_func, v_normal_func, S, offset] = vThetaFunction1(l_vals, q_vals); %#ok<NASGU>

M = size(S,2);

v_normal_at_t = cell(nspikes,1);
v_theta_at_t  = cell(nspikes,1);

for i = 1:nspikes
    v_normal_at_t{i} = v_normal_func(theta_true(i));
    v_theta_at_t{i}  = v_theta_func(theta_true(i));
end

fprintf('\nLifted angular dimension:\n');
fprintf('M = %d\n', M);
fprintf('Full dual SDP block size = %d\n', (M+1)*N_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7) FRESNEL REFERENCE CHANNEL h_fresnel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_fresnel = cell(nspikes,1);
h_fresnel = zeros(N_r,1);

for i = 1:nspikes
    a_fresnel{i} = fresnelSteeringVector(N_r, lambda, d, n_idx, ...
                                          r_true(i), theta_true(i));

    h_fresnel = h_fresnel + coef(i)*a_fresnel{i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8) DIRECT TRUNCATED JA CHANNEL h_JA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_JA = cell(nspikes,1);
h_JA = zeros(N_r,1);

for i = 1:nspikes
    a_JA{i} = jacobiAngerSteeringVector(N_r, lambda, d, n_idx, ...
                                        l_vals, q_vals, ...
                                        r_true(i), theta_true(i));

    h_JA = h_JA + coef(i)*a_JA{i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9) FACTORIZED ATOMS a(r,theta) = C(r)v(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_r = cell(nspikes,1);
a_fromFunctions = cell(nspikes,1);

fprintf('\nSanity check C(r)v(theta) versus direct JA:\n');

for i = 1:nspikes
    C_r{i} = buildCr(N_r, lambda, d, n_idx, l_vals, q_vals, r_true(i));
    a_fromFunctions{i} = C_r{i} * v_theta_at_t{i};

    fprintf('source %d: ||C(r)v(theta)-a_JA||_2 = %.3e\n', ...
        i, norm(a_fromFunctions{i}-a_JA{i}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10) BUILD DISTANCE DICTIONARY D = [vec(C(r_1)),...,vec(C(r_Nd))]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = zeros(N_r*P, N_d);

for i = 1:N_d
    C_ri   = buildCr(N_r, lambda, d, n_idx, l_vals, q_vals, r_list(i));
    D(:,i) = C_ri(:);
end

idx_min = zeros(nspikes,1);

for i = 1:nspikes
    [~, idx_min(i)] = min(abs(r_list-r_true(i)));
end

alpha_sparse = cell(nspikes,1);

for i = 1:nspikes
    alpha_sparse{i} = zeros(N_d,1);
    alpha_sparse{i}(idx_min(i)) = 1;
end

d_vec = D(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11) LIFTED JA MEASUREMENT y_me(n) = <X, Phi_n S>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(N_d, M);

for ii = 1:nspikes
    X = X + coef(ii) * alpha_sparse{ii} * (v_normal_at_t{ii}).';
end

Phi_n = cell(N_r,1);

for n = 1:N_r
    Phi_n{n} = zeros(N_d, P);
end

for i = 1:N_d
    blockStart = (i-1)*(P*N_r) + 1;

    for j = 1:P
        subStart = blockStart + (j-1)*N_r;
        subEnd   = subStart + N_r - 1;

        subVec = d_vec(subStart:subEnd);

        for n = 1:N_r
            Phi_n{n}(i,j) = subVec(n);
        end
    end
end

Phiten = zeros(N_d, P, N_r);
PhiS1  = zeros(N_d, M, N_r);

for n = 1:N_r
    Phiten(:,:,n) = Phi_n{n};
end

y_me = zeros(N_r,1);

for n = 1:N_r
    PhiS = Phiten(:,:,n) * S;
    PhiS1(:,:,n) = PhiS;

    y_me(n) = sum(sum(X .* PhiS));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 12) MODEL ERROR DIAGNOSTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_fresnel_JA_direct = norm(h_JA - h_fresnel);
err_lift_vs_JA        = norm(y_me - h_JA);
err_lift_vs_fresnel   = norm(y_me - h_fresnel);

rel_err_fresnel_JA_direct = err_fresnel_JA_direct / max(norm(h_fresnel), eps);
rel_err_lift_vs_JA        = err_lift_vs_JA / max(norm(h_JA), eps);
rel_err_lift_vs_fresnel   = err_lift_vs_fresnel / max(norm(h_fresnel), eps);

fprintf('\nActual errors:\n');
fprintf('||h_JA - h_fresnel||_2       = %.6e\n', err_fresnel_JA_direct);
fprintf('relative direct JA error     = %.6e\n', rel_err_fresnel_JA_direct);
fprintf('||y_me - h_JA||_2            = %.6e\n', err_lift_vs_JA);
fprintf('relative internal lift error = %.6e\n', rel_err_lift_vs_JA);
fprintf('||y_me - h_fresnel||_2       = %.6e\n', err_lift_vs_fresnel);
fprintf('relative model error         = %.6e\n', rel_err_lift_vs_fresnel);

% Distance dictionary coherence.
G = abs(D' * D);
G = G ./ sqrt(diag(G) * diag(G).');
G(1:N_d+1:end) = 0;

fprintf('\nDistance dictionary coherence:\n');
fprintf('max off-diagonal coherence = %.6f\n', max(G(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 13) THEORETICAL JA TRUNCATION ERROR BOUND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eta_model_th, model_parts] = predictFresnelJATruncationBound( ...
    N_r, lambda, d, n_idx, r_true, coef, I1, I2);

fprintf('\nTheoretical Fresnel-to-JA truncation bound:\n');
fprintf('analytic JA truncation/model bound (double precision) = %.6e\n', eta_model_th);
fprintf('relative eta_model_th = %.6e\n', eta_model_th/max(norm(h_fresnel),eps));

disp(struct2table(model_parts));

lift_consistency_tolerance = 1e-11*max(norm(h_JA),1);
lift_consistency_ok = err_lift_vs_JA <= lift_consistency_tolerance;
model_bound_contains_error = ...
    err_fresnel_JA_direct <= eta_model_th*(1+1e-10);

if ~lift_consistency_ok
    error(['The factorized lift does not agree with the direct truncated ', ...
           'Jacobi--Anger channel within the stated tolerance.']);
end
if ~model_bound_contains_error
    error(['The realized Fresnel-to-Jacobi--Anger mismatch exceeds the ', ...
           'computed analytic truncation bound.']);
end

if eta_model_th > 0.01*norm(h_fresnel)
    warning(['eta_model_th is large relative to ||h_fresnel||. ', ...
             'Recovery from Fresnel data may fail or become non-informative. ', ...
             'Increase I2, reduce gamma_z2, reduce range span, or use JA data ', ...
             'when validating the convex machinery.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14) ADD COMPLEX GAUSSIAN NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_noise_th = predictNoiseRadiusCN(sigma_noise, N_r, delta_prob);

w_raw = (sigma_noise/sqrt(2)) * (randn(N_r,1) + 1j*randn(N_r,1));

w = w_raw;
nw = norm(w);

if nw > eta_noise_th && eta_noise_th > 0
    w = (eta_noise_th / nw) * w_raw;
end

y_JA_noisy      = y_me + w;
y_fresnel_noisy = h_fresnel + w;

eta_exact_total = eta_model_th + eta_noise_th;

fprintf('\nNoise statistics:\n');
fprintf('actual ||w||_2        = %.6e\n', norm(w));
fprintf('theoretical eta_noise = %.6e\n', eta_noise_th);
fprintf('eta_model + eta_noise = %.6e\n', eta_exact_total);

if USE_FRESNEL_REFERENCE
    y_obs   = y_fresnel_noisy;
    eta_sdp = eta_exact_total;
    fprintf('\nSDP uses y_obs = h_fresnel + w.\n');
    fprintf('This tests Fresnel data with JA model mismatch included in eta.\n');
else
    y_obs   = y_JA_noisy;
    eta_sdp = eta_noise_th;

    fprintf('\nSDP uses y_obs = y_me + w.\n');
    fprintf('This validates the finite lifted JA recovery model itself.\n');
end

% Numerical feasibility floor for the conic solver.  For the paper-reported
% Fresnel experiment the analytic model bound is larger, so this floor does
% not change the reported robustness radius.
eta_numerical_floor = 1e-6*norm(y_obs);
eta_sdp = max(eta_sdp, eta_numerical_floor);

fprintf('numerical eta floor = %.6e\n', eta_numerical_floor);
fprintf('eta_sdp             = %.6e\n', eta_sdp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 15) NOISY DUAL SDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('cvx_begin','file') ~= 2
    error(['CVX is required to run this reproducibility script. ', ...
           'Install and initialize CVX with an SDP solver such as SDPT3, ', ...
           'then rerun this file.']);
end

solver_name = "SDPT3";
solver_precision = "best";

cvx_clear

cvx_begin sdp
    cvx_solver sdpt3
    % Keep the strict setting used for the paper run.  SDPT3 may legitimately
    % return Inaccurate/Solved at residuals of a few parts in 10^6.
    cvx_precision best

    variable q_me(N_r,1) complex
    variable t_dual nonnegative
    variable X_total((M+1)*N_d, (M+1)*N_d) hermitian

    maximize( real(y_obs' * q_me) - eta_sdp*t_dual )

    subject to
        norm(q_me) <= t_dual;

        % Explicit block diagonality.
        for i_blk = 1:N_d
            idx_i = (i_blk-1)*(M+1) + (1:(M+1));

            for j_blk = i_blk+1:N_d
                idx_j = (j_blk-1)*(M+1) + (1:(M+1));

                % X_total is Hermitian, so the conjugate lower block follows.
                X_total(idx_i, idx_j) == 0;
            end
        end

        % Per-distance-bin dual certificate.
        for i = 1:N_d
            idx = (i-1)*(M+1) + (1:(M+1));

            X_total(idx, idx) >= 0;
            X_total(idx(end), idx(end)) == 1;
            trace(X_total(idx, idx)) == 2;

            % Correct adjoint:
            %
            % If y_n = sum_j X_j Phi_{j,n}, then
            %
            %   A^*(q)_j = sum_n q_n conj(Phi_{j,n})
            %
            % because the dual objective is real(y^H q).
            dual_col_i = sum( ...
                repmat(q_me,1,M) .* conj(squeeze(PhiS1(i,:,:)).'), ...
                1).';

            X_total(idx(1:M), idx(end)) == dual_col_i;

            % Positive trigonometric lags are sufficient because the moment
            % block is Hermitian; negative-lag equations are conjugates.
            for lag = 1:(M-1)
                MQ = diag(ones(M-lag,1), lag);
                trace(MQ * X_total(idx(1:M), idx(1:M))) == 0;
            end
        end
cvx_end

fprintf('\nCVX status: %s\n', cvx_status);
fprintf('CVX optimal value: %.6e\n', cvx_optval);

solver_accepted = strcmp(cvx_status,'Solved') || ...
                  strcmp(cvx_status,'Inaccurate/Solved');

solver_table = table(solver_name, solver_precision, string(cvx_status), ...
    cvx_optval, solver_accepted, eta_numerical_floor, eta_sdp, ...
    USE_FRESNEL_REFERENCE, ...
    'VariableNames', {'Solver','RequestedPrecision','CVXStatus', ...
    'CVXOptimalValue','SolverAccepted','EtaNumericalFloor','EtaSDP', ...
    'UseFresnelReference'});
writetable(solver_table, fullfile(results_dir, 'solver_report.csv'));

if ~solver_accepted || ~isfinite(cvx_optval) || any(~isfinite(q_me))
    error('CVX failed. Do not use q_me or peak picking from this run.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16) DUAL POLYNOMIAL SURFACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dualpoly1 = @(t,i_d) abs( ...
    v_normal_func(t)' * ...
    sum( ...
        repmat(q_me,1,M) .* conj(squeeze(PhiS1(i_d,:,:)).'), ...
        1).' );

res_angle = 0.001;
t_ind = 0:res_angle:pi;

fun11 = zeros(length(t_ind), N_d);

for j = 1:N_d
    for i = 1:length(t_ind)
        fun11(i,j) = dualpoly1(t_ind(i), j);
    end
end

sampled_dual_max = max(fun11(:));
dual_bound_tolerance = 5e-3;
sampled_dual_bound_ok = sampled_dual_max <= 1+dual_bound_tolerance;

fprintf('\nSampled dual-certificate check:\n');
fprintf('max sampled |Q(r,theta)| = %.12f\n', sampled_dual_max);
fprintf('max sampled value <= 1 + %.1e: %d\n', ...
    dual_bound_tolerance, sampled_dual_bound_ok);

if ~sampled_dual_bound_ok
    warning(['The sampled dual polynomial exceeds one beyond the stated ', ...
             'numerical tolerance.  Treat the certificate as unresolved.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 17) PEAK PICKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fun_temp = fun11;

estimated_angles_raw        = zeros(nspikes,1);
estimated_distance_bins_raw = zeros(nspikes,1);
peak_amplitudes_raw         = zeros(nspikes,1);

suppression_radius_angle = max(1, round(0.01/(t_ind(2)-t_ind(1))));
suppression_radius_dist  = 1;

for k = 1:nspikes
    [max_val, linIdx] = max(fun_temp(:));

    [row_idx, col_idx] = ind2sub(size(fun_temp), linIdx);

    estimated_angles_raw(k)        = t_ind(row_idx);
    estimated_distance_bins_raw(k) = col_idx;
    peak_amplitudes_raw(k)         = max_val;

    row_min = max(1, row_idx - suppression_radius_angle);
    row_max = min(length(t_ind), row_idx + suppression_radius_angle);

    col_min = max(1, col_idx - suppression_radius_dist);
    col_max = min(length(r_list), col_idx + suppression_radius_dist);

    fun_temp(row_min:row_max, col_min:col_max) = 0;
end

estimated_distances_raw = r_list(estimated_distance_bins_raw(:));
estimated_distances_raw = estimated_distances_raw(:);

% Match range--angle pairs jointly.  Independently sorting the two coordinate
% lists can combine an angle from one source with the range of another.
angle_match_scale = max(max(theta_true)-min(theta_true), res_angle);
range_match_scale = max(max(r_true)-min(r_true), eps);
[match_order, normalized_joint_matching_cost] = ...
    bestJointSupportAssignment(theta_true, r_true, ...
    estimated_angles_raw, estimated_distances_raw, ...
    angle_match_scale, range_match_scale);

estimated_angles        = estimated_angles_raw(match_order);
estimated_distance_bins = estimated_distance_bins_raw(match_order);
estimated_distances     = estimated_distances_raw(match_order);
peak_amplitudes         = peak_amplitudes_raw(match_order);

results_table = table(source_idx, estimated_angles, ...
    estimated_distance_bins, estimated_distances, peak_amplitudes, ...
    'VariableNames', {'Source','AngleEstimate','DistanceBinEstimate', ...
    'DistanceEstimate','DualPolyValue'});

truth_table = table(source_idx, theta_true, target_bins, r_true, ...
    'VariableNames', {'Source','TrueAngle','TrueDistanceBin','TrueDistance'});

fprintf('\nEstimated support:\n');
disp(results_table);

fprintf('Ground truth:\n');
disp(truth_table);

angle_error = max(abs(theta_true(:)-estimated_angles(:)));
range_error = max(abs(r_true(:)-estimated_distances(:)));
peak_height_ok = all(peak_amplitudes >= 0.95);

fprintf('\nSupport errors:\n');
fprintf('max paired angle error = %.6e rad\n', angle_error);
fprintf('max paired range error = %.6e m\n', range_error);
fprintf('normalized joint matching cost = %.6e\n', ...
    normalized_joint_matching_cost);
fprintf('all recovered peak values >= 0.95: %d\n', peak_height_ok);

% Persist support results before amplitude estimation.
writetable(results_table, fullfile(results_dir, 'support_estimates.csv'));
writetable(truth_table, fullfile(results_dir, 'support_truth.csv'));
support_error_table = table(angle_error, range_error, ...
    normalized_joint_matching_cost, peak_height_ok, ...
    'VariableNames', {'MaxPairedAngleErrorRad','MaxPairedRangeErrorM', ...
    'NormalizedJointMatchingCost','PeakHeightOK'});
writetable(support_error_table, fullfile(results_dir, 'support_errors.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 18) AMPLITUDE ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th_est = estimated_angles;
r_est  = estimated_distances(:);
s_hat  = length(th_est);

Ame = zeros(N_r, s_hat);

for ell = 1:s_hat
    C_rl = buildCr(N_r, lambda, d, n_idx, l_vals, q_vals, r_est(ell));
    v_theta_l = v_theta_func(th_est(ell));

    Ame(:,ell) = C_rl * v_theta_l;
end

c_hat = Ame \ y_obs;
amplitude_fit_residual = norm(Ame*c_hat-y_obs);
relative_amplitude_fit_residual = ...
    amplitude_fit_residual/max(norm(y_obs),eps);
amplitude_matrix_condition = cond(Ame);

fprintf('\nEstimated channel amplitudes:\n');
disp(c_hat);

fprintf('True channel amplitudes:\n');
disp(coef);
fprintf('amplitude-refit residual          = %.6e\n', amplitude_fit_residual);
fprintf('relative amplitude-refit residual = %.6e\n', ...
    relative_amplitude_fit_residual);
fprintf('condition number of refit matrix  = %.6e\n', ...
    amplitude_matrix_condition);

amplitude_table = table((1:length(c_hat)).', real(c_hat(:)), imag(c_hat(:)), ...
    real(coef(:)), imag(coef(:)), abs(c_hat(:)-coef(:)), ...
    'VariableNames', {'Source','RealEstimated','ImagEstimated','RealTrue', ...
    'ImagTrue','AbsoluteCoefficientError'});
writetable(amplitude_table, fullfile(results_dir, 'amplitudes.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 19) PAPER FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Paper Figure 1: dual-polynomial localization surface.  A surface viewed
% from above preserves the nonuniform physical range coordinates; imagesc
% would incorrectly display the log-spaced range rows as linearly spaced.
fig_dual = figure('Color','w', 'Name','Lifted dual polynomial');
surf(t_ind/pi, r_list, fun11.', ...
    'EdgeColor','none', 'FaceColor','interp', 'HandleVisibility','off');
view(2);
axis tight;
set(gca, 'YScale','log', 'Layer','top', 'FontSize',10);
caxis([0, max(1.01, 1.01*sampled_dual_max)]);

xlabel('\theta/\pi');
ylabel('range $r$ (m)', 'Interpreter','latex');
title('|Dual polynomial| over (r,\theta)');
colormap parula;
colorbar;
hold on;
marker_height = max(fun11(:)) + 0.02;

plot3(theta_true(:)/pi, r_true(:), ...
    marker_height*ones(nspikes,1), ...
    'ko', 'MarkerFaceColor','w', 'MarkerSize', 8, 'LineWidth', 1.2, ...
    'DisplayName','True support');

plot3(th_est(:)/pi, r_est(:), ...
    marker_height*ones(nspikes,1), ...
    'r^', 'MarkerFaceColor','r', 'MarkerSize', 8, 'LineWidth', 1.0, ...
    'DisplayName','Estimated support');

legend('Location','southoutside', ...
    'Orientation','horizontal');

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 20) SAVE REPRODUCIBILITY ARTIFACTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save all numerical arrays needed to regenerate paper tables and plots.
metadata = struct();
metadata.schema_version = '1.0.0';
metadata.random_seed = random_seed;
metadata.matlab_version = version;
metadata.use_fresnel_reference = USE_FRESNEL_REFERENCE;
metadata.solver = char(solver_name);
metadata.requested_cvx_precision = char(solver_precision);
metadata.generated_at = datestr(now, 30);
metadata.note = ['Generated by run_cpam_lifted_dual_figure1_repro.m. ', ...
    'Compare numerical values with stated tolerances, not bitwise.'];

save(fullfile(results_dir, 'cpam_results.mat'), ...
    'metadata', 'N_r', 'fc', 'c0', 'lambda', 'k_fac', 'I1', 'I2', ...
    'N_d', 'nspikes', 'scene', 'scene_table', 'range_grid_table', ...
    'sigma_noise', 'delta_prob', ...
    'd', 'D_ap', 'r_min', 'r_max', 'r_list', 'range_info', ...
    'l_vals', 'q_vals', 'M', 'target_bins', 'r_true', 'theta_true', 'coef', ...
    'true_sources_in_radiating_fresnel', 'true_sources_beyond_fraunhofer', ...
    'h_fresnel', 'h_JA', 'y_me', 'y_obs', 'eta_model_th', 'eta_noise_th', ...
    'eta_exact_total', 'eta_numerical_floor', 'eta_sdp', ...
    'err_fresnel_JA_direct', 'err_lift_vs_JA', ...
    'err_lift_vs_fresnel', 'rel_err_fresnel_JA_direct', 'rel_err_lift_vs_JA', ...
    'rel_err_lift_vs_fresnel', 'lift_consistency_ok', ...
    'model_bound_contains_error', 'q_me', 'cvx_status', 'cvx_optval', ...
    'solver_accepted', 't_ind', 'fun11', 'sampled_dual_max', ...
    'sampled_dual_bound_ok', 'estimated_angles_raw', ...
    'estimated_distance_bins_raw', 'estimated_distances_raw', ...
    'peak_amplitudes_raw', 'estimated_angles', 'estimated_distance_bins', ...
    'estimated_distances', 'peak_amplitudes', 'angle_error', 'range_error', ...
    'normalized_joint_matching_cost', 'peak_height_ok', 'c_hat', ...
    'amplitude_fit_residual', 'relative_amplitude_fit_residual', ...
    'amplitude_matrix_condition');

% Write compact diagnostics table used by the paper.
diagnostics_table = table( ...
    err_fresnel_JA_direct, rel_err_fresnel_JA_direct, ...
    err_lift_vs_JA, rel_err_lift_vs_JA, ...
    err_lift_vs_fresnel, rel_err_lift_vs_fresnel, ...
    eta_model_th, eta_noise_th, eta_exact_total, ...
    eta_numerical_floor, eta_sdp, lift_consistency_ok, ...
    model_bound_contains_error, sampled_dual_max, sampled_dual_bound_ok, ...
    'VariableNames', {'ErrFresnelJADirect','RelErrFresnelJADirect', ...
    'ErrLiftVsJA','RelErrLiftVsJA','ErrLiftVsFresnel','RelErrLiftVsFresnel', ...
    'EtaModelTheoretical','EtaNoiseTheoretical','EtaModelPlusNoise', ...
    'EtaNumericalFloor','EtaSDP','LiftConsistencyOK', ...
    'ModelBoundContainsError','SampledDualMaximum','SampledDualBoundOK'});
writetable(diagnostics_table, fullfile(results_dir, 'diagnostics.csv'));

% Save only the paper-facing figure.
saveFigureForPaper(fig_dual, figures_dir, 'dualpol1');

fprintf('\nSaved reproducibility outputs to:\n');
fprintf('  %s\n', results_dir);
fprintf('  %s\n', figures_dir);

result = struct();
result.metadata = metadata;
result.results_dir = results_dir;
result.figures_dir = figures_dir;
result.results_table = results_table;
result.truth_table = truth_table;
result.support_error_table = support_error_table;
result.amplitude_table = amplitude_table;
result.diagnostics_table = diagnostics_table;
result.solver_table = solver_table;
result.scene = scene;
result.scene_table = scene_table;
result.range_grid_table = range_grid_table;
result.range_info = range_info;
result.cvx_status = cvx_status;
result.cvx_optval = cvx_optval;
result.solver_accepted = solver_accepted;
result.eta_numerical_floor = eta_numerical_floor;
result.eta_sdp = eta_sdp;
result.angle_error = angle_error;
result.range_error = range_error;
result.normalized_joint_matching_cost = normalized_joint_matching_cost;
result.peak_height_ok = peak_height_ok;
result.sampled_dual_max = sampled_dual_max;
result.sampled_dual_bound_ok = sampled_dual_bound_ok;
result.amplitude_fit_residual = amplitude_fit_residual;
result.relative_amplitude_fit_residual = relative_amplitude_fit_residual;
result.amplitude_matrix_condition = amplitude_matrix_condition;
result.figure_pdf = fullfile(figures_dir, 'dualpol1.pdf');
result.figure_png = fullfile(figures_dir, 'dualpol1.png');
result.figure_matlab = fullfile(figures_dir, 'dualpol1.fig');
result.experimentType = 'finite-harmonic full-circle SDP numerical illustration';
result.isQPACRecoveryExperiment = false;
result.optimizerConnection = ['The same-vector QPAC perturbation estimate does not ', ...
    'establish optimizer identity between the lifted SDP and the exact Fresnel TV problem.'];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveFigureForPaper(fh, figures_dir, basename)
%SAVEFIGUREFORPAPER Save the paper-facing figure as vector PDF, PNG, and FIG.
pdf_name = fullfile(figures_dir, [basename '.pdf']);
png_name = fullfile(figures_dir, [basename '.png']);
fig_name = fullfile(figures_dir, [basename '.fig']);

if exist('exportgraphics','file') == 2 || exist('exportgraphics','builtin') == 5
    exportgraphics(fh, pdf_name, 'ContentType','vector');
    exportgraphics(fh, png_name, 'Resolution', 300);
else
    print(fh, pdf_name, '-dpdf', '-painters');
    print(fh, png_name, '-dpng', '-r300');
end
savefig(fh, fig_name);
end

function [d, D_ap, r_min, r_max, r_list, info] = ...
    designRangeGridFromI1I2(N_r, lambda, I1, I2, N_d, cfg)
% DESIGNRANGEGRIDFROMI1I2
%
% Designs aperture spacing and range grid from JA truncation orders I1 and I2.
%
% Angular argument:
%
%   z1_max = k D_ap.
%
% Enforce:
%
%   z1_max <= gamma_z1 I1.
%
% Quadratic/range argument:
%
%   z2(r) = k D_ap^2/(4r).
%
% Enforce:
%
%   z2(r_min) <= gamma_z2 I2,
%   z2(r_max) >= z2_min.

if nargin < 6
    cfg = struct();
end

if ~isfield(cfg, 'gamma_z1')
    cfg.gamma_z1 = 0.75;
end

if ~isfield(cfg, 'gamma_z2')
    cfg.gamma_z2 = 0.30;
end

if ~isfield(cfg, 'z2_min')
    cfg.z2_min = 0.05;
end

if ~isfield(cfg, 'reactive_safety')
    cfg.reactive_safety = 1.05;
end

if ~isfield(cfg, 'fraunhofer_multiple')
    cfg.fraunhofer_multiple = 30.0;
end

if ~isfield(cfg, 'grid_type')
    cfg.grid_type = 'log_curvature';
end

if I1 <= 0
    error('I1 must be positive.');
end

if I2 <= 0
    error('I2 must be positive. I2 = 0 removes range-sensitive q-expansion.');
end

if N_r < 2
    error('N_r must be at least 2.');
end

if N_d < 2
    error('N_d must be at least 2.');
end

if cfg.gamma_z1 <= 0 || cfg.gamma_z1 >= 1
    error('cfg.gamma_z1 must be in (0,1).');
end

if cfg.gamma_z2 <= 0 || cfg.gamma_z2 >= 1
    error('cfg.gamma_z2 must be in (0,1).');
end

if cfg.z2_min <= 0
    error('cfg.z2_min must be positive.');
end

if cfg.fraunhofer_multiple <= 0
    error('cfg.fraunhofer_multiple must be positive.');
end

k = 2*pi/lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aperture design from I1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_ap = cfg.gamma_z1*I1/k;

d = D_ap/(N_r-1);

z1_max = k*D_ap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical scales.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_reactive = 0.62*sqrt(D_ap^3/lambda);
R_F        = 2*D_ap^2/lambda;

r_reactive_min = cfg.reactive_safety*R_reactive;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range design from I2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z2_max_allowed = cfg.gamma_z2*I2;

r_min_JA = k*D_ap^2/(4*z2_max_allowed);

r_min = max(r_reactive_min, r_min_JA);

r_max_curvature = k*D_ap^2/(4*cfg.z2_min);

r_max_fraunhofer = cfg.fraunhofer_multiple*R_F;

r_max = min(r_max_curvature, r_max_fraunhofer);

if r_min >= r_max
    error(['No valid range interval: r_min >= r_max. ', ...
           'Increase I1, increase I2, decrease z2_min, ', ...
           'or increase fraunhofer_multiple.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curvature-compatible grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z2_at_r_min = k*D_ap^2/(4*r_min);
z2_at_r_max = k*D_ap^2/(4*r_max);

switch lower(cfg.grid_type)

    case 'log_curvature'
        z2_grid = logspace(log10(z2_at_r_min), ...
                           log10(z2_at_r_max), N_d);
        r_list = k*D_ap^2 ./ (4*z2_grid);
        r_list = sort(r_list, 'ascend');

    case 'linear_curvature'
        z2_grid = linspace(z2_at_r_min, z2_at_r_max, N_d);
        r_list = k*D_ap^2 ./ (4*z2_grid);
        r_list = sort(r_list, 'ascend');

    case 'log_range'
        r_list = logspace(log10(r_min), log10(r_max), N_d);

    case 'linear_range'
        r_list = linspace(r_min, r_max, N_d);

    otherwise
        error(['cfg.grid_type must be one of: ', ...
               'log_curvature, linear_curvature, log_range, linear_range.']);
end

r_list = r_list(:).';

z2_grid_check = k*D_ap^2 ./ (4*r_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = struct();

info.lambda = lambda;
info.k = k;

info.I1 = I1;
info.I2 = I2;

info.gamma_z1 = cfg.gamma_z1;
info.gamma_z2 = cfg.gamma_z2;
info.z2_min = cfg.z2_min;

info.d = d;
info.d_over_lambda = d/lambda;
info.D_aperture = D_ap;

info.z1_max = z1_max;
info.z1_max_over_I1 = z1_max/I1;

info.R_reactive = R_reactive;
info.R_fraunhofer = R_F;

info.r_reactive_min = r_reactive_min;
info.r_min_JA = r_min_JA;
info.r_min = r_min;

info.r_max_curvature = r_max_curvature;
info.r_max_fraunhofer = r_max_fraunhofer;
info.r_max = r_max;

info.range_span = r_max-r_min;
info.range_ratio = r_max/r_min;

info.z2_at_r_min = z2_at_r_min;
info.z2_at_r_max = z2_at_r_max;

info.max_z2_on_grid = max(z2_grid_check);
info.min_z2_on_grid = min(z2_grid_check);
info.max_z2_over_I2 = max(z2_grid_check)/I2;

info.fraunhofer_multiple = cfg.fraunhofer_multiple;
info.full_grid_in_radiating_fresnel = ...
    (r_min >= R_reactive) && (r_max <= R_F);
info.lower_fresnel_margin = r_min-R_reactive;
info.upper_fresnel_margin = R_F-r_max;

end

function a = fresnelSteeringVector(N_r, lambda, d, n_idx, r, theta)
% FRESNELSTEERINGVECTOR
%
% Fresnel path:
%
%   Delta_n = -x_n cos(theta) + x_n^2 sin^2(theta)/(2r)
%
% Steering:
%
%   a_n = exp(-j k Delta_n)

k = 2*pi/lambda;

x = n_idx(:)*d;

delta_fresnel = -x*cos(theta) + (x.^2/(2*r))*sin(theta)^2;

a = exp(-1j*k*delta_fresnel);

end

function a = jacobiAngerSteeringVector(N_r, lambda, d, n_idx, ...
                                       l_vals, q_vals, r, theta)
% JACOBIANGERSTEERINGVECTOR
%
% Truncated JA approximation of the Fresnel steering vector.
%
% Fresnel factorization:
%
%   a_n =
%       exp(-j k x_n^2/(4r))
%       exp( j k x_n cos(theta))
%       exp( j k x_n^2 cos(2theta)/(4r))
%
% JA expansions:
%
%   exp(j z cos(theta))
%       = sum_l j^l J_l(z) exp(jl theta)
%
%   exp(j beta cos(2theta))
%       = sum_q j^q J_q(beta) exp(j2q theta)

k = 2*pi/lambda;

a = zeros(N_r,1);

for nn = 1:N_r
    n   = n_idx(nn);
    x_n = n*d;

    E_n = exp(-1j*k*(x_n^2/(4*r)));

    total_sum = 0;

    for q = q_vals
        for l = l_vals
            jphase = (1j)^(l+q);

            Jl_ = besselj(l, k*x_n);
            Jq_ = besselj(q, k*(x_n^2/(4*r)));

            ang_ = exp(1j*(l+2*q)*theta);

            total_sum = total_sum + jphase*Jl_*Jq_*ang_;
        end
    end

    a(nn) = E_n*total_sum;
end

end

function C_r = buildCr(N_r, lambda, d, n_idx, l_vals, q_vals, r)
% BUILDCR
%
% Constructs C(r) in C^{N_r x (|l_vals||q_vals|)}.

L0 = length(l_vals);
Q0 = length(q_vals);
P  = L0*Q0;

C_r = zeros(N_r, P);

k = 2*pi/lambda;

idx = 1;

for q_i = 1:Q0
    q_val = q_vals(q_i);

    for l_i = 1:L0
        l_val = l_vals(l_i);

        for nn = 1:N_r
            n   = n_idx(nn);
            x_n = n*d;

            E_n = exp(-1j*k*(x_n^2/(4*r)));

            jpow = (1j)^(l_val+q_val);

            Jl_ = besselj(l_val, k*x_n);

            Jq_ = besselj(q_val, k*(x_n^2/(4*r)));

            C_r(nn,idx) = E_n*jpow*Jl_*Jq_;
        end

        idx = idx + 1;
    end
end

end

function [v_theta_func, v_normal_func, S, offset] = vThetaFunction1(l_vals, q_vals)
% VTHETAFUNCTION1
%
% Builds S so that:
%
%   v_theta(theta) = S v_normal(theta),
%
% where v_theta contains:
%
%   exp(j(l+2q)theta).

P = length(l_vals)*length(q_vals);

exponents = zeros(P,1);

idx = 1;

for q = q_vals
    for l = l_vals
        exponents(idx) = l+2*q;
        idx = idx+1;
    end
end

offset = -min(exponents);

shifted_exponents = exponents+offset;

M = max(shifted_exponents)+1;

S = zeros(P,M);

for k = 1:P
    S(k, shifted_exponents(k)+1) = 1;
end

v_normal_func = @(t) exp(1j*(((0:M-1).')-offset)*t);

v_theta_func = @(t) S*v_normal_func(t);

end

function [eta_model_th, parts] = predictFresnelJATruncationBound( ...
    N_r, lambda, d, n_idx, r_true, coef, I1, I2)
% PREDICTFRESNELJATRUNCATIONBOUND
%
% Bound for:
%
%   ||h_fresnel - y_JA||_2.
%
% This only accounts for JA truncation error, because the reference model is
% already Fresnel.

k = 2*pi/lambda;

nspikes = length(r_true);

eta_trunc_each = zeros(nspikes,1);
eta_total_each = zeros(nspikes,1);

for ell = 1:nspikes
    r = r_true(ell);

    e_trunc_vec = zeros(N_r,1);

    for nn = 1:N_r
        n   = n_idx(nn);
        x_n = n*d;

        z1 = k*x_n;
        z2 = k*(x_n^2/(4*r));

        eps_l = besselTailBound(z1, I1);
        eps_q = besselTailBound(z2, I2);

        e_trunc_vec(nn) = eps_l + eps_q + eps_l*eps_q;
    end

    eta_trunc_each(ell) = norm(e_trunc_vec, 2);

    eta_total_each(ell) = abs(coef(ell))*eta_trunc_each(ell);
end

eta_model_th = sum(eta_total_each);

parts.eta_trunc_each = eta_trunc_each;
parts.eta_total_each = eta_total_each;

end

function tail = besselTailBound(z, I)
%BESSELTAILBOUND Two-sided JA tail at the fixed argument z.
%
% Computes
%     2*sum_{m=I+1}^\infty |J_m(z)|
% using direct summation plus an analytic remainder bound.


X = abs(z);

if X == 0
    tail = 0;
    return;
end

m0 = I + 1;
K  = max(m0 + 20, ceil(X) + 30);

orders  = m0:(K-1);
partial = sum(abs(besselj(orders, X)));

% For m >= K:
% |J_m(X)| <= (X/2)^m/m! * exp(X^2/(4(m+1))).
x = X/2;

logUK = K*log(x) - gammaln(K+1) + X^2/(4*(K+1));
UK    = exp(logUK);

rho = x/(K+1);
remainder = UK/(1-rho);

tail = 2*(partial + remainder)*(1 + 100*eps);
end

function eta_noise = predictNoiseRadiusCN(sigma_noise, N, delta_prob)
% PREDICTNOISERADIUSCN
%
% For:
%
%   w ~ CN(0, sigma_noise^2 I_N)
%
% Since:
%
%   ||w||_2^2 = (sigma_noise^2/2) chi^2_{2N},
%
% Laurent-Massart gives:
%
%   ||w||_2^2 <= sigma_noise^2 ...
%       * (N + sqrt(2N log(1/delta)) + log(1/delta)).

x = log(1/delta_prob);

eta_noise = sigma_noise * sqrt(N + sqrt(2*N*x) + x);

end

function [best_order, best_cost] = bestJointSupportAssignment( ...
    theta_true, range_true, theta_est, range_est, angle_scale, range_scale)
%BESTJOINTSUPPORTASSIGNMENT Pair estimates jointly in angle and range.
%
% The true support is used only to report numerical recovery errors after
% optimization.  It does not affect the SDP or the extracted peak locations.

theta_true = theta_true(:);
range_true = range_true(:);
theta_est  = theta_est(:);
range_est  = range_est(:);

n_sources = length(theta_true);
if length(range_true) ~= n_sources || ...
        length(theta_est) ~= n_sources || length(range_est) ~= n_sources
    error('True and estimated support arrays must have equal lengths.');
end
if n_sources > 8
    error('Brute-force support matching is restricted to at most 8 sources.');
end

candidate_orders = perms(1:n_sources);
best_order = (1:n_sources).';
best_cost = Inf;

for i_order = 1:size(candidate_orders,1)
    order_now = candidate_orders(i_order,:).';
    normalized_angle_error = ...
        (theta_est(order_now)-theta_true)/angle_scale;
    normalized_range_error = ...
        (range_est(order_now)-range_true)/range_scale;
    cost_now = mean(normalized_angle_error.^2 + ...
                    normalized_range_error.^2);

    if cost_now < best_cost
        best_cost = cost_now;
        best_order = order_now;
    end
end

end
