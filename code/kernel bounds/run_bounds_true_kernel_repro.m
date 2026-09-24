% Copyright (c) 2026 Sajad Daei.
% Email: sajado@kth.se
%
% Self-contained code to reproduce the pointwise bound, support-uniform bound, and true kernel.
%
% This single MATLAB file evaluates and plots the pointwise and
% support-uniform QPAC bounds against the true kernel magnitude for the
% eight Hermite and derivative channels
%
%     K, H, dK, dH, d2K, d3K, d2H, d3H.
%
% No companion MATLAB source files, data files, or precomputed tables are
% required. All helper routines used by the bound computation are included
% below in this file.
%
% This code was written by Sajad Daei and is provided as supplementary
% reproducibility material for:
%
%   S. Daei, G. Fodor, and M. Skoglund,
%   "A Mathematical Theory of Near-Field Super-Resolution."
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
% RUN_BOUNDS_TRUE_KERNEL_REPRO
%
% Reviewer-facing pairwise-bound experiment for the CPAM/QPAC paper.
%
% The code returns and plots only the three paper-relevant quantities:
%
%   1) true pointwise kernel magnitude;
%   2) pointwise bound;
%   3) support-uniform bound.
%
% Previous comparison-only curves and branch-only diagnostics are intentionally
% absent from the returned data and from the plots. This file is therefore
% suitable for the final paper figure/table story where the full-phase
% residue-correlation machinery is the reported method.
%
% Mathematical choices:
%   - theta is allowed in (0,pi), so tau=d*cos(theta)/r may be negative;
%   - the pointwise bound applies the scalar Bbest to the exact sampled
%     coefficient vector of each channel;
%   - the support-uniform bound uses the full-phase residue-correlation
%     support-class machinery, derivative and residue-linear branches, and the
%     tangent-arc/dictionary envelopes for normalized and higher channels;
%   - no artificial unit cap is applied to pointwise higher-derivative channels;
%   - normalized support channels K,H,dK,dH still receive the Hilbert-space
%     unit cap, as in the theorem-level envelope construction.
%
% Run in MATLAB:
%     out = run_bounds_true_kernel_repro();
%
% Optional example:
%     opts.sliceType = 'theta';      % 'theta', 'range', or 'both'
%     opts.makePlots = true;
%     out = run_bounds_true_kernel_repro(opts);
%
% Returned fields:
%     out.thetaSlice.data{j}.Ktrue,       out.thetaSlice.data{j}.Kactual,
%     out.thetaSlice.data{j}.Ksupport,   and similarly for all channels.
%
% Here "actual" means the pointwise bound, and
% "support" means the new support-uniform bound.
%
function out = run_bounds_true_kernel_repro(opts)
if nargin < 1
    opts = struct();
end
if ~isfield(opts, 'sliceType')
    opts.sliceType = 'both';
end
if ~isfield(opts, 'class') || ~isstruct(opts.class)
    opts.class = struct();
end
if ~isfield(opts.class, 'subcells') || ~isstruct(opts.class.subcells)
    opts.class.subcells = struct();
end
if ~isfield(opts.class.subcells, 'theta')
    opts.class.subcells.theta = 1;
end
if ~isfield(opts.class.subcells, 'range')
    opts.class.subcells.range = 1;
end
close all; clc;

%% ============================================================
% User parameters
%% ============================================================
fc_list = get_opt(opts, 'fc_list', 10e9);
Nr = get_opt(opts, 'Nr', 64);
% numPlot = get_opt(opts, 'numPlot', 100);

mShape = get_opt(opts, 'mShape', 6);
weightFamily = get_opt(opts, 'weightFamily', 'binomial');

r_min = get_opt(opts, 'r_min',.1);
r_max = get_opt(opts, 'r_max', 10);
theta_min = get_opt(opts, 'theta_min', 0.15);
theta_max = get_opt(opts, 'theta_max', pi - 0.05);

r_support = get_opt(opts, 'r_support', 2);
theta_support = get_opt(opts, 'theta_support', pi/7);
r_eval_fixed = get_opt(opts, 'r_eval_fixed', 5);
theta_eval_fixed = get_opt(opts, 'theta_eval_fixed', pi/1.5);
% theta1 = linspace(0.15, 0.75, 35);
% theta2 = linspace(0.75, 1.35, 90);
% theta3 = linspace(1.35, pi-0.05, 45);
% 
% theta_eval_grid = unique([theta1 theta2 theta3]);
% % theta_eval_grid = get_opt(opts, 'theta_eval_grid', linspace(theta_min, theta_max, numPlot));
% numPlot = numel(theta_eval_grid);


c_min = cos(pi - 0.05);
c_max = cos(0.15);

c_grid = linspace(c_max, c_min, 50);
theta_grid = acos(c_grid);

theta_eval_grid = theta_grid;
% numPlot = numel(theta_grid);

r1 = linspace(r_min, 8, 50);
r2 = linspace(8,r_max,20 );
r_eval_grid = unique([r1 r2]);
% r_eval_grid = get_opt(opts, 'r_eval_grid', linspace(r_min, r_max, numPlot));

sliceType = lower(char(get_opt(opts, 'sliceType', 'both')));

useHalfLambda = get_opt(opts, 'useHalfLambda', true);
d_fixed = get_opt(opts, 'd_fixed', []);

makePlots = get_opt(opts, 'makePlots', true);
fontSize = get_opt(opts, 'fontSize', 24);
quiet = get_opt(opts, 'quiet', false);
ratioTruthFloor = get_opt(opts, 'ratioTruthFloor', 1e-8);
assertSupport = get_opt(opts, 'assertSupport', true);
assertActual = get_opt(opts, 'assertActual', true);

channels = {'K','H','dK','dH','d2K','d3K','d2H','d3H'};
normalizedChannels = {'K','H','dK','dH'};
higherChannels = {'d2K','d3K','d2H','d3H'};

%% ============================================================
% Geometry / interval options
%% ============================================================
tauDefault = struct();
tauDefault.guard = 1e-12;
tauDefault.tightBilin2D = true;
tauDefault.bilinMaxDepth2D = 2;
tauDefault.bilinMaxCells2D = 64;
tauDefault.bilinTolAbs2D = 1e-10;
tauDefault.bilinTolRel2D = 1e-7;
tauDefault.bilinSampleOrder = 10;
tauDefault.maxDepth2D = 2;
tauDefault.maxCells2D = 32;
tauDefault.tolAbs2D = 1e-10;
tauDefault.tolRel2D = 1e-7;
tauDefault.sampleOrder2D = 10;

tauOpts = merge_structs(tauDefault, get_opt(opts, 'tau', struct()));

%% ============================================================
% Actual pointwise B_best options
%% ============================================================
actualDefault = struct();
actualDefault.enable = true;
actualDefault.useL1Fallback = true;
actualDefault.useDerivativeBranch = true;
actualDefault.useResidueCorrelationBranch = true;
actualDefault.useFullPhaseResidueCorrelationBranch = true;
actualDefault.useResidueLinearBranch = true;
actualDefault.residueQMax = 32;
% Q=1 makes the full-phase residue-correlation branch equal to the exact
% pointwise kernel.  Keep it off by default so the branch remains a genuine
% residue-split envelope.  Set true only when exact singleton evaluation is
% intended.
actualDefault.fullPhaseCorrIncludeQ1 = false;
actualDefault.residueCorrIncludeQ1 = true;
actualDefault.fullPhaseCorrSafety = 1e-12;
actualDefault.residueUseNuCap = false;
actualDefault.residueNuMax = 0.15;
actualDefault.residueMinDen = 1e-12;
actualDefault.pCorrOrderList = 1:8;

actualOpts = merge_structs(actualDefault, get_opt(opts, 'actual', struct()));
actualOpts.pCorrOrderList = actualOpts.pCorrOrderList(:).';


%% ============================================================
% Support class options
%% ============================================================
classDefault = struct();
classDefault.enable = true;
classDefault.theoremMode = true;
classDefault.uniformMode = 'samplecell';
classDefault.sampleCellPad = 1e-12;
classDefault.useOneFallback = true;
classDefault.useL1Fallback = true;

classDefault.derivative = struct();
classDefault.derivative.enable = true;
classDefault.derivative.d0List = [
    0.15 0.20 0.30 0.40 0.55 0.70 ...
    0.85 1.00 1.20 1.40 1.70 2.00 2.40 2.80 3.05];

classDefault.derivative.rhoBarList = [
    1e-8 3e-8 1e-7 3e-7 ...
    1e-6 3e-6 1e-5 3e-5 ...
    1e-4 3e-4 1e-3 3e-3 ...
    1e-2 3e-2 0.07 0.15 0.3 0.6 1.0];
classDefault.derivative.curveDNCells = 32;

classDefault.corr = struct();
classDefault.corr.enable = true;
classDefault.corr.QList = [1 2 8 16 32];
classDefault.corr.pList = actualOpts.pCorrOrderList;
% Full-phase residue-correlation is the default support-uniform correlation
% envelope.  It keeps the lag phase h*beta_s+h^2*gamma and combines a
% phase-aware cosine majorant with the Abel denominator magnitude bound.
classDefault.corr.useFullPhaseResidueCorrelationBranch = true;
classDefault.corr.fullPhaseCorrIncludeQ1 = false;
classDefault.corr.fullPhaseCorrSafety = 1e-12;
classDefault.corr.autoDeltaSafety = 0.98;
classDefault.corr.autoDeltaMin = 1e-12;

classDefault.lin = struct();
classDefault.lin.enable = true;
classDefault.lin.qList = [1 2 3 4 5 6 8 10 12 16 20 24 32];
classDefault.lin.AMode = 'phase';
classDefault.lin.AWindow = 0.12;
classDefault.lin.autoGammaSafety = 0.98;
classDefault.lin.autoGammaMin = 1e-12;

classDefault.higher = struct();
classDefault.higher.method = 'dictionary';
classDefault.higher.useFastCauchyCap = true;

classDefault.subcells = struct();
classDefault.subcells.theta = 1;
classDefault.subcells.range = 1;

classDefault.rSuppRange = [NaN NaN];
classDefault.thetaSuppRange = [NaN NaN];
classDefault.rEvalRange = [NaN NaN];
classDefault.thetaEvalRange = [NaN NaN];

classOpts = merge_structs(classDefault, get_opt(opts, 'class', struct()));
classOpts = normalize_class_options(classOpts);

if any(isnan(classOpts.rSuppRange))
    classOpts.rSuppRange = [r_support r_support];
end
if any(isnan(classOpts.thetaSuppRange))
    classOpts.thetaSuppRange = [theta_support theta_support];
end

%% ============================================================
% Validation and preprocessing
%% ============================================================
validate_inputs(Nr, useHalfLambda, d_fixed, r_min, r_max, theta_min, theta_max, ...
    r_support, theta_support, r_eval_fixed, theta_eval_fixed, theta_eval_grid, ...
    r_eval_grid, sliceType, actualOpts, classOpts, tauOpts);

w = build_flat_end_weight_r4(Nr, mShape, weightFamily);
W0 = sum(w);
b = w(:)/W0;
n = (0:(Nr-1)).';

[xTangent, yTangent, Gtan, E] = build_tangent_geometry_from_b(b);

pMax = max([actualOpts.pCorrOrderList(:); classOpts.corr.pList(:); 1]);
Qall = unique([1:actualOpts.residueQMax, classOpts.corr.QList(:).'], 'stable');

coeffK = precompute_scalar_coeff_data(b, pMax, Qall);
geom = build_tau_geometry_cache(b, xTangent, yTangent, Gtan, E, pMax, Qall, coeffK, tauOpts);

if ~quiet
    fprintf('============================================================\n');
    fprintf('bound all-channel comparison against true kernel\n');
    fprintf('Nr=%d, fc=[ ', Nr); fprintf('%.4g ', fc_list/1e9); fprintf('] GHz\n');
    fprintf('theta domain: (0,pi), tau may be negative\n');
    fprintf('channels: K H dK dH d2K d3K d2H d3H\n');
    fprintf('uniformMode=%s, actual=%d, derivative=%d, corr=%d, linear=%d\n', ...
        char(classOpts.uniformMode), logical(actualOpts.enable), ...
        logical(classOpts.derivative.enable), logical(classOpts.corr.enable), logical(classOpts.lin.enable));
    fprintf('higher method=%s; subcells theta=%d, range=%d\n', ...
        char(classOpts.higher.method), classOpts.subcells.theta, classOpts.subcells.range);
    fprintf('No local Taylor and no expensive higher-correlation polynomial branch.\n');
    fprintf('============================================================\n\n');
end

out = struct();
out.Nr = Nr;
out.fc_list = fc_list;
out.w = w;
out.W0 = W0;
out.b = b;
out.n = n;
out.xTangent = xTangent;
out.yTangent = yTangent;
out.Gtan = Gtan;
out.E = E;
out.coeffK = coeffK;
out.geom = geom;
out.classOpts = classOpts;
out.actualOpts = actualOpts;
out.tauOpts = tauOpts;
out.channels = channels;
out.normalizedChannels = normalizedChannels;
out.higherChannels = higherChannels;
out.ratioTruthFloor = ratioTruthFloor;

if ismember(sliceType, {'theta','both'})
    out.thetaSlice = run_one_slice('theta');
end
if ismember(sliceType, {'range','both'})
    out.rangeSlice = run_one_slice('range');
end

if ~quiet
    fprintf('Done. Returned true kernel, pointwise bound, and support-uniform bound.\n');
end

%% ============================================================
% Nested slice driver
%% ============================================================
    function sliceOut = run_one_slice(whichSlice)
        switch lower(char(whichSlice))
            case 'theta'
                xGrid = theta_eval_grid(:).';
                xLabel = '\theta';
                fixedText = sprintf('r_{eval}=%.4g', r_eval_fixed);
            case 'range'
                xGrid = r_eval_grid(:).';
                xLabel = 'r';
                fixedText = sprintf('\\theta_{eval}=%.4g', theta_eval_fixed);
            otherwise
                error('Unknown slice type.');
        end

        data = cell(1, numel(fc_list));

        for jf = 1:numel(fc_list)
            fc = fc_list(jf);
            lambda = 3e8/fc;
            if useHalfLambda
                d = lambda/2;
            else
                d = d_fixed;
            end
            k_lambda = 2*pi/lambda;

            if ~quiet
                apertureD = (Nr - 1)*d;
                Rreactive = 0.62*sqrt(apertureD^3/lambda);
                Rfraunhofer = 2*apertureD^2/lambda;
                fprintf('fc %.4g GHz: lambda=%.4g m, d=%.4g m, D=%.4g m, reactive=%.4g m, Fraunhofer=%.4g m\n', ...
                    fc/1e9, lambda, d, apertureD, Rreactive, Rfraunhofer);
            end

            classGlobal = make_slice_class_ranges(classOpts, whichSlice, xGrid, r_eval_fixed, theta_eval_fixed);
            if classOpts.enable && strcmpi(char(classOpts.uniformMode), 'global')
                specGlobal = build_support_spec(Nr, classGlobal, k_lambda, d, whichSlice, r_eval_fixed, theta_eval_fixed);
            else
                specGlobal = [];
            end

            R = init_result_arrays(numel(xGrid), channels);
            R.fc = fc;
            R.lambda = lambda;
            R.d = d;
            R.k_lambda = k_lambda;
            R.xGrid = xGrid;

            for ii = 1:numel(xGrid)
                switch lower(char(whichSlice))
                    case 'theta'
                        r_eval = r_eval_fixed;
                        theta_eval = xGrid(ii);
                    case 'range'
                        r_eval = xGrid(ii);
                        theta_eval = theta_eval_fixed;
                end

                if classOpts.enable && strcmpi(char(classOpts.uniformMode), 'samplecell')
                    classLocal = make_sample_class_ranges(classOpts, whichSlice, xGrid, ii, r_eval_fixed, theta_eval_fixed);
                    classSpec = build_support_spec(Nr, classLocal, k_lambda, d, whichSlice, r_eval_fixed, theta_eval_fixed);
                elseif classOpts.enable
                    classLocal = classGlobal;
                    classSpec = specGlobal;
                else
                    classLocal = classOpts;
                    classSpec = disabled_support_spec();
                end

                [truth, actual, support, aux] = evaluate_one_pair( ...
                    Nr, n, coeffK, geom, actualOpts, classLocal, classSpec, tauOpts, ...
                    k_lambda, d, r_support, theta_support, r_eval, theta_eval, channels, higherChannels);

                for ic = 1:numel(channels)
                    ch = channels{ic};
                    R.([ch 'true'])(ii) = truth.(ch);
                    R.([ch 'actual'])(ii) = actual.(ch).value;
                    R.([ch 'support'])(ii) = support.(ch).value;
                    R.([ch 'actualBranch'])(ii) = actual.(ch).bestBranch;
                    R.([ch 'supportBranch'])(ii) = support.(ch).bestBranch;
                end

                R.omega1(ii) = aux.omega1;
                R.omega2(ii) = aux.omega2;
                R.dNplus(ii) = aux.dNplus;
                R.tauEval(ii) = aux.tauEval;
                R.tauSupp(ii) = aux.tauSupp;
                R.aux{ii} = aux;
            end

            if assertActual
                assert_bound_valid(R, channels, 'actual', 1e-9);
            end
            if assertSupport
                assert_bound_valid(R, channels, 'support', 1e-9);
            end

            data{jf} = R;
        end

        if ~quiet
            fprintf('------------------------------------------------------------\n');
            fprintf('Slice = %s (%s)\n', whichSlice, fixedText);
            for jj = 1:numel(data)
                fprintf('  fc = %.4g GHz\n', data{jj}.fc/1e9);
                for ic = 1:numel(channels)
                    ch = channels{ic};
                    print_new_bound_diag(ch, data{jj}.([ch 'true']), data{jj}.([ch 'actual']), data{jj}.([ch 'support']), ratioTruthFloor);
                end
            end
            fprintf('------------------------------------------------------------\n\n');
        end

        figs = [];
        if makePlots
            figs = plot_new_bounds_slice(data, whichSlice, xLabel, fixedText, fontSize, channels);
        end

        sliceOut = struct();
        sliceOut.whichSlice = whichSlice;
        sliceOut.xGrid = xGrid;
        sliceOut.data = data;
        sliceOut.figures = figs;
    end
end

%% =====================================================================
% Pair evaluation
%% =====================================================================
function [truth, actual, support, aux] = evaluate_one_pair( ...
    Nr, n, coeffK, geom, actualOpts, classOpts, classSpec, tauOpts, ...
    k_lambda, d, r_s, theta_s, r_e, theta_e, channels, higherChannels)

[omega1, omega2] = phase_pair(k_lambda, d, r_s, theta_s, r_e, theta_e);
phase = exp(1i*(n*omega1 + (n.^2)*omega2));

tauSupp = d*cos(theta_s)/r_s;
tauEval = d*cos(theta_e)/r_e;
alphaEval = d/r_e;
kappa = k_lambda*d;

hSupp = physical_h_vec(geom, tauSupp);
hEval = physical_h_vec(geom, tauEval);
[Q2eval, Q3eval] = higher_Q_profiles(geom, kappa, alphaEval, theta_e);

coeff = struct();
coeff.K   = coeffK;
coeff.H   = precompute_scalar_coeff_data(1i*geom.b(:).*hSupp(:), geom.pMax, geom.Qall);
coeff.dK  = precompute_scalar_coeff_data(-1i*geom.b(:).*hEval(:), geom.pMax, geom.Qall);
coeff.dH  = precompute_scalar_coeff_data(geom.b(:).*hEval(:).*hSupp(:), geom.pMax, geom.Qall);
coeff.d2K = precompute_scalar_coeff_data(geom.b(:).*Q2eval(:), geom.pMax, geom.Qall);
coeff.d3K = precompute_scalar_coeff_data(geom.b(:).*Q3eval(:), geom.pMax, geom.Qall);
coeff.d2H = precompute_scalar_coeff_data(1i*geom.b(:).*hSupp(:).*Q2eval(:), geom.pMax, geom.Qall);
coeff.d3H = precompute_scalar_coeff_data(1i*geom.b(:).*hSupp(:).*Q3eval(:), geom.pMax, geom.Qall);

truth = struct();
actual = struct();
for ic = 1:numel(channels)
    ch = channels{ic};
    truth.(ch) = abs(sum(coeff.(ch).a(:).*phase));
    if actualOpts.enable
        actual.(ch) = actual_pointwise_Bbest(coeff.(ch), omega1, omega2, actualOpts);
    else
        actual.(ch) = disabled_bound_info('actual-disabled', inf);
    end
end

if classOpts.enable
    [Bsupport, detail] = cpam_support_class_bound_magic(coeffK, geom, classSpec, classOpts, tauOpts, channels, higherChannels);
    support = struct();
    for ic = 1:numel(channels)
        ch = channels{ic};
        support.(ch) = make_bound_info(Bsupport.(ch), detail.branch.(ch), inf);
    end
else
    support = disabled_bound_struct('support-disabled', channels);
    detail = struct();
end

aux = struct();
aux.omega1 = omega1;
aux.omega2 = omega2;
aux.dNplus = exact_dN_plus(Nr, omega1, omega2);
aux.tauSupp = tauSupp;
aux.tauEval = tauEval;
aux.alphaEval = alphaEval;
aux.Q2eval = Q2eval;
aux.Q3eval = Q3eval;
aux.classSpec = classSpec;
aux.supportDetail = detail;
end

%% =====================================================================
% Magic support class bound
%% =====================================================================
function [B, detail] = cpam_support_class_bound_magic(coeffK, geom, spec, classOpts, tauOpts, channels, higherChannels)
% Optional physical subcell maximum. This is not local Taylor. It simply
% covers the support by smaller support classes and takes the max.
nt = max(1, ceil(get_struct_opt(classOpts.subcells, 'theta', 1)));
nr = max(1, ceil(get_struct_opt(classOpts.subcells, 'range', 1)));

if (nt > 1 || nr > 1) && isfinite_range(spec.rEvalRange) && isfinite_range(spec.thetaEvalRange)
    B = zero_channel_branch(channels);
    branch = branch_names('subcell-max', channels);
    detailCells = cell(nr, nt);
    rEdges = linspace(spec.rEvalRange(1), spec.rEvalRange(2), nr+1);
    thEdges = linspace(spec.thetaEvalRange(1), spec.thetaEvalRange(2), nt+1);
    classLocal = classOpts;
    classLocal.subcells.theta = 1;
    classLocal.subcells.range = 1;
    for ir = 1:nr
        for it = 1:nt
            classLocal.rEvalRange = [rEdges(ir) rEdges(ir+1)];
            classLocal.thetaEvalRange = [thEdges(it) thEdges(it+1)];
            subSpec = build_support_spec(spec.Nr, classLocal, spec.k_lambda, spec.d, ...
                spec.whichSlice, spec.rEvalFixed, spec.thetaEvalFixed);
            [Bs, ds] = cpam_support_class_bound_magic_core(coeffK, geom, subSpec, classLocal, tauOpts, channels, higherChannels);
            detailCells{ir,it} = ds;
            for k = 1:numel(channels)
                f = channels{k};
                if Bs.(f) > B.(f)
                    B.(f) = Bs.(f);
                    branch.(f) = string(sprintf('subcell-max:%s', char(ds.branch.(f))));
                end
            end
        end
    end
    detail = struct();
    detail.branch = branch;
    detail.subcells = detailCells;
    detail.phaseBox = spec.phaseBox;
    detail.tauSupp = spec.tauSupp;
    detail.tauEval = spec.tauEval;
    detail.dNplusLower = spec.dNplusLower;
    return;
end

[B, detail] = cpam_support_class_bound_magic_core(coeffK, geom, spec, classOpts, tauOpts, channels, higherChannels);
end

function [B, detail] = cpam_support_class_bound_magic_core(coeffK, geom, spec, classOpts, tauOpts, channels, higherChannels)

guard = get_struct_opt(tauOpts, 'guard', 1e-12);
B = inf_channel_branch(channels);
branch = branch_names('none', channels);
branchRaw = struct();

% K: fixed normalized channel.
infoK = scalar_support_phase_bound(coeffK, spec, classOpts);
C = inf_channel_branch(channels);
C.K = infoK.value;
[B, branch] = update_best_branch(B, branch, C, sprintf('K:%s', char(infoK.bestBranch)), guard, channels);
branchRaw.K = infoK;

% Arc data for support/evaluation tangent channels.
arcS = tau_arc_data(geom.G, spec.tauSupp);
arcE = tau_arc_data(geom.G, spec.tauEval);

% H: arc hull in tau_s.
BH = -inf;
for js = 1:numel(arcS.tau)
    hs = physical_h_vec(geom, arcS.tau(js));
    coeff = precompute_scalar_coeff_data(1i*geom.b(:).*hs(:), geom.pMax, geom.Qall);
    infj = scalar_support_phase_bound(coeff, spec, classOpts);
    BH = max(BH, infj.value);
end
C = inf_channel_branch(channels);
C.H = arcS.kappa * BH;
[B, branch] = update_best_branch(B, branch, C, 'H:arc-fixed-scalar', guard, channels);

% dK: arc hull in tau_e.
BdK = -inf;
for je = 1:numel(arcE.tau)
    he = physical_h_vec(geom, arcE.tau(je));
    coeff = precompute_scalar_coeff_data(-1i*geom.b(:).*he(:), geom.pMax, geom.Qall);
    infj = scalar_support_phase_bound(coeff, spec, classOpts);
    BdK = max(BdK, infj.value);
end
C = inf_channel_branch(channels);
C.dK = arcE.kappa * BdK;
[B, branch] = update_best_branch(B, branch, C, 'dK:arc-fixed-scalar', guard, channels);

% dH: product of two arc hulls.
BdH = -inf;
for je = 1:numel(arcE.tau)
    he = physical_h_vec(geom, arcE.tau(je));
    for js = 1:numel(arcS.tau)
        hs = physical_h_vec(geom, arcS.tau(js));
        coeff = precompute_scalar_coeff_data(geom.b(:).*he(:).*hs(:), geom.pMax, geom.Qall);
        infij = scalar_support_phase_bound(coeff, spec, classOpts);
        BdH = max(BdH, infij.value);
    end
end
C = inf_channel_branch(channels);
C.dH = arcE.kappa * arcS.kappa * BdH;
[B, branch] = update_best_branch(B, branch, C, 'dH:arc2-fixed-scalar', guard, channels);

% Higher derivative dictionary bounds.
high = higher_dictionary_support_bounds(geom, spec, classOpts, tauOpts);
C = inf_channel_branch(channels);
for k = 1:numel(higherChannels)
    f = higherChannels{k};
    C.(f) = high.(f).value;
end
[B, branch] = update_best_branch(B, branch, C, 'higher-dictionary', guard, channels);

% Optional fast l1/Cauchy caps for higher channels.
if get_struct_opt(classOpts.higher, 'useFastCauchyCap', true)
    fast = higher_fast_l1_cauchy_caps(geom, spec);
    C = inf_channel_branch(channels);
    C.d2K = min(fast.d2K.l1, fast.d2K.cauchy);
    C.d3K = min(fast.d3K.l1, fast.d3K.cauchy);
    C.d2H = min(fast.d2H.l1, fast.d2H.cauchy);
    C.d3H = min(fast.d3H.l1, fast.d3H.cauchy);
    [B, branch] = update_best_branch(B, branch, C, 'higher-fast-cap', guard, channels);
else
    fast = struct();
end

% Normalized channels receive the Hilbert-space unit cap.
normFields = {'K','H','dK','dH'};
for k = 1:numel(normFields)
    f = normFields{k};
    B.(f) = min(1, max(real(B.(f)), 0));
end

% Higher derivative channels receive no unit cap.
for k = 1:numel(higherChannels)
    f = higherChannels{k};
    B.(f) = max(real(B.(f)), 0);
end

detail = struct();
detail.branch = branch;
detail.branchRaw = branchRaw;
detail.higher = high;
detail.fastHigher = fast;
detail.arcSupp = arcS;
detail.arcEval = arcE;
detail.phaseBox = spec.phaseBox;
detail.tauSupp = spec.tauSupp;
detail.tauEval = spec.tauEval;
detail.dNplusLower = spec.dNplusLower;
end

function info = scalar_support_phase_bound(coeff, spec, classOpts)
% Fixed-coefficient support-uniform phase-class bound.
candidates = [];
labels = {};

if get_struct_opt(classOpts, 'useL1Fallback', true)
    candidates(end+1) = coeff.l1; %#ok<AGROW>
    labels{end+1} = 'support-l1'; %#ok<AGROW>
end

if isfield(classOpts, 'derivative') && classOpts.derivative.enable && coeff.flatEndsR4
    w2Int = spec.phaseBox.omega2;
    dNLower = spec.dNplusLower;
    rhoUpper = interval_sin_abs_max(w2Int(1), w2Int(2));
    for id = 1:numel(classOpts.derivative.d0List)
        d0 = classOpts.derivative.d0List(id);
        for ir = 1:numel(classOpts.derivative.rhoBarList)
            rhoBar = classOpts.derivative.rhoBarList(ir);
            if dNLower < d0 - 1e-14 || rhoUpper > rhoBar + 1e-14
                continue;
            end
            beta = derivative_class_beta_sine(d0, rhoBar);
            candidates(end+1) = derivative_class_bound_scalar(beta, coeff); %#ok<AGROW>
            labels{end+1} = sprintf('support-derivative-d%.3g-rho%.3g', d0, rhoBar); %#ok<AGROW>
        end
    end
end

if isfield(classOpts, 'corr') && classOpts.corr.enable
    QList = classOpts.corr.QList(:).';
    pList = classOpts.corr.pList(:).';
    if get_struct_opt(classOpts.corr, 'useFullPhaseResidueCorrelationBranch', true)
        includeQ1 = get_struct_opt(classOpts.corr, 'fullPhaseCorrIncludeQ1', false);
        for Q = QList
            if Q == 1 && ~includeQ1
                continue;
            end
            if Q > coeff.N
                continue;
            end
            BQ = fixedQ_fullphase_residue_correlation_cell_scalar(coeff, Q, spec.phaseBox, pList, classOpts.corr);
            if isfinite(BQ)
                candidates(end+1) = BQ; %#ok<AGROW>
                labels{end+1} = sprintf('support-fprc-Q%d', Q); %#ok<AGROW>
            end
        end
    end
end

if isfield(classOpts, 'lin') && classOpts.lin.enable
    w2Int = spec.phaseBox.omega2;
    qList = classOpts.lin.qList(:).';
    for q = qList
        if q > coeff.N
            continue;
        end
        AList = get_A_list_for_q_support(classOpts.lin, q, spec.phaseBox);
        for Ares = AList
            [epsBar, okEps] = eps_bar_from_omega2_interval(q, Ares, w2Int(1), w2Int(2));
            if ~okEps || epsBar > pi
                continue;
            end
            nuBar = interval_dist_to_lattice_max( ...
                q^2*(w2Int(1) - pi*Ares/q), ...
                q^2*(w2Int(2) - pi*Ares/q), 2*pi);
            if ~isfinite(nuBar) || nuBar > pi
                continue;
            end
            gammaProfile = classOpts.lin.autoGammaSafety * ...
                gamma_profile_lower_from_phase_box(q, Ares, epsBar, spec.phaseBox);
            if any(gammaProfile <= classOpts.lin.autoGammaMin) || any(~isfinite(gammaProfile))
                continue;
            end
            candidates(end+1) = fixed_q_residue_linear_sine_scalar_bound(coeff.a, q, nuBar, gammaProfile); %#ok<AGROW>
            labels{end+1} = sprintf('support-linear-q%d-A%d', q, Ares); %#ok<AGROW>
        end
    end
end

[value, label] = min_finite_candidates(candidates, labels);
info = make_bound_info(value, label, inf);
info.l1 = coeff.l1;
end

%% =====================================================================
% Higher derivative dictionary support
%% =====================================================================
function high = higher_dictionary_support_bounds(geom, spec, classOpts, tauOpts) %#ok<INUSD>
% Bound Q2/Q3 by atomic dictionaries in x,y. Coefficient functions are
% bounded on the physical cell in alpha and cos(theta). H-type channels get
% the tau_s tangent-arc hull.

[atoms2, gamma2] = higher_dictionary_atoms(geom, spec, 2);
[atoms3, gamma3] = higher_dictionary_atoms(geom, spec, 3);
arcS = tau_arc_data(geom.G, spec.tauSupp);

B2K = 0;
for l = 1:numel(gamma2)
    coeff = precompute_scalar_coeff_data(atoms2{l}, geom.pMax, geom.Qall);
    infol = scalar_support_phase_bound(coeff, spec, classOpts);
    B2K = B2K + gamma2(l)*infol.value;
end

B3K = 0;
for l = 1:numel(gamma3)
    coeff = precompute_scalar_coeff_data(atoms3{l}, geom.pMax, geom.Qall);
    infol = scalar_support_phase_bound(coeff, spec, classOpts);
    B3K = B3K + gamma3(l)*infol.value;
end

B2H = 0;
for l = 1:numel(gamma2)
    bestEndpoint = 0;
    for js = 1:numel(arcS.tau)
        hs = physical_h_vec(geom, arcS.tau(js));
        coeff = precompute_scalar_coeff_data(1i*atoms2{l}(:).*hs(:), geom.pMax, geom.Qall);
        infol = scalar_support_phase_bound(coeff, spec, classOpts);
        bestEndpoint = max(bestEndpoint, infol.value);
    end
    B2H = B2H + gamma2(l)*arcS.kappa*bestEndpoint;
end

B3H = 0;
for l = 1:numel(gamma3)
    bestEndpoint = 0;
    for js = 1:numel(arcS.tau)
        hs = physical_h_vec(geom, arcS.tau(js));
        coeff = precompute_scalar_coeff_data(1i*atoms3{l}(:).*hs(:), geom.pMax, geom.Qall);
        infol = scalar_support_phase_bound(coeff, spec, classOpts);
        bestEndpoint = max(bestEndpoint, infol.value);
    end
    B3H = B3H + gamma3(l)*arcS.kappa*bestEndpoint;
end

high = struct();
high.d2K.value = max(real(B2K),0);
high.d3K.value = max(real(B3K),0);
high.d2H.value = max(real(B2H),0);
high.d3H.value = max(real(B3H),0);
high.gamma2 = gamma2;
high.gamma3 = gamma3;
high.arcSupp = arcS;
end

function [atoms, gamma] = higher_dictionary_atoms(geom, spec, order)
% Atoms include the b_n weight. gamma contains nonnegative support-uniform
% coefficient envelopes.
x = geom.x(:);
y = geom.y(:);
b = geom.b(:);
kappa = spec.kappa;
alphaHi = max(spec.alphaEvalRange);

cInt = sort(cos(spec.thetaEvalRange));
[cAbsMin, cMax] = abs_c_range(cInt);
c2Max = cMax^2;
c3Max = cMax^3;
[~, s2Max] = sin2_interval(spec.thetaEvalRange);
sMax = sqrt(max(0,s2Max));
s3Max = sMax^3;

p2m1Max = max_abs_poly_in_c2(cInt, [ -1, 2 ]);  % |2c^2-1|
p3m1Max = max_abs_poly_in_c2(cInt, [ -1, 3 ]);  % |3c^2-1|
cm1c2Max = max_t_times_one_minus_t2(cAbsMin, cMax);
c2m1c2Max = max_t2_times_one_minus_t2(cAbsMin, cMax);

switch order
    case 2
        atoms = cell(5,1);
        atoms{1} = b .* x;
        atoms{2} = b .* y;
        atoms{3} = b .* (x.^2);
        atoms{4} = b .* (x.*y);
        atoms{5} = b .* (y.^2);

        gamma = zeros(5,1);
        gamma(1) = kappa * cMax;
        gamma(2) = kappa * alphaHi * p2m1Max;
        gamma(3) = kappa^2 * s2Max;
        gamma(4) = 2*kappa^2 * alphaHi * cm1c2Max;
        gamma(5) = kappa^2 * alphaHi^2 * c2m1c2Max;

    case 3
        atoms = cell(9,1);
        atoms{1} = b .* x;
        atoms{2} = b .* y;
        atoms{3} = b .* (x.^2);
        atoms{4} = b .* (x.*y);
        atoms{5} = b .* (y.^2);
        atoms{6} = b .* (x.^3);
        atoms{7} = b .* (x.^2.*y);
        atoms{8} = b .* (x.*y.^2);
        atoms{9} = b .* (y.^3);

        gamma = zeros(9,1);
        gamma(1) = kappa * sMax;
        gamma(2) = 4*kappa * alphaHi * cMax * sMax;
        gamma(3) = 3*kappa^2 * cMax * sMax;
        gamma(4) = 3*kappa^2 * alphaHi * p3m1Max * sMax;
        gamma(5) = 3*kappa^2 * alphaHi^2 * cMax * p2m1Max * sMax;
        gamma(6) = kappa^3 * s3Max;
        gamma(7) = 3*kappa^3 * alphaHi * cMax * s3Max;
        gamma(8) = 3*kappa^3 * alphaHi^2 * c2Max * s3Max;
        gamma(9) = kappa^3 * alphaHi^3 * c3Max * s3Max;

    otherwise
        error('order must be 2 or 3.');
end

gamma = max(real(gamma),0);
end

function fast = higher_fast_l1_cauchy_caps(geom, spec)
% Conservative fast amplitude caps. These are used only as a fallback/cap.
kappa = spec.kappa;
N = geom.N;
x = geom.x(:);
y = geom.y(:);
b = geom.b(:);

thInt = sort(spec.thetaEvalRange(:).');
tauEvalInt = sort(spec.tauEval(:).');
tauSuppInt = sort(spec.tauSupp(:).');
alphaInt = sort(spec.alphaEvalRange(:).');

smax = interval_sin_abs_max(thInt(1), thInt(2));
[cLo,cHi] = cos_interval(thInt);
[zLo,zHi] = alpha_cos2_interval(alphaInt, thInt);

Aub = zeros(N,1);
Bub = zeros(N,1);
A4ub = zeros(N,1);
hsub = zeros(N,1);
for nn = 1:N
    Aub(nn) = max_abs_affine_1d(x(nn), y(nn), tauEvalInt);
    Bub(nn) = max_abs_bilinear_box(x(nn), y(nn), [cLo cHi], [zLo zHi]);
    A4ub(nn) = max_abs_affine_1d(x(nn), 4*y(nn), tauEvalInt);
    hsub(nn) = max_abs_linear_over_sqrtq(x(nn), y(nn), geom.G, tauSuppInt);
end

Q2ub = kappa*Bub + (kappa^2)*(smax^2).*(Aub.^2);
Q3ub = kappa*smax*A4ub + 3*(kappa^2)*smax.*Aub.*Bub + (kappa^3)*(smax^3).*(Aub.^3);

fast = struct();
fast.d2K.l1 = sum(b.*Q2ub);
fast.d3K.l1 = sum(b.*Q3ub);
fast.d2H.l1 = sum(b.*hsub.*Q2ub);
fast.d3H.l1 = sum(b.*hsub.*Q3ub);
fast.d2K.cauchy = sqrt(sum(b.*Q2ub.^2));
fast.d3K.cauchy = sqrt(sum(b.*Q3ub.^2));
fast.d2H.cauchy = fast.d2K.cauchy;
fast.d3H.cauchy = fast.d3K.cauchy;

fields = {'d2K','d3K','d2H','d3H'};
for k = 1:numel(fields)
    f = fields{k};
    fast.(f).l1 = max(real(fast.(f).l1),0);
    fast.(f).cauchy = max(real(fast.(f).cauchy),0);
end
end

function arc = tau_arc_data(G, tauInt)
tauInt = sort(tauInt(:).');
if any(~isfinite(tauInt))
    tauInt = [0 0];
end
if abs(tauInt(2)-tauInt(1)) <= 1e-15
    tauList = mean(tauInt);
    kappa = 1;
    u0 = tangent_unit_vector(G, tauList);
    u1 = u0;
else
    tauList = [tauInt(1), tauInt(2)];
    u0 = tangent_unit_vector(G, tauInt(1));
    u1 = tangent_unit_vector(G, tauInt(2));
    dotv = max(-1, min(1, real(u0(:).'*u1(:))));
    dpsi = acos(dotv);
    if dpsi >= pi - 1e-12
        kappa = inf;
    else
        kappa = 1/cos(dpsi/2);
    end
end
arc = struct();
arc.tau = tauList;
arc.kappa = kappa;
arc.u0 = u0;
arc.u1 = u1;
end

function u = tangent_unit_vector(G, tau)
v = sqrtm(G) * [1; tau];
nu = norm(v);
if nu <= 0 || ~isfinite(nu)
    error('Invalid tangent arc vector.');
end
u = v/nu;
end

function tf = isfinite_range(x)
tf = numel(x) == 2 && all(isfinite(x));
end

function [lo, hi] = abs_c_range(cInt)
loC = cInt(1); hiC = cInt(2);
hi = max(abs(loC), abs(hiC));
if loC <= 0 && hiC >= 0
    lo = 0;
else
    lo = min(abs(loC), abs(hiC));
end
end

function M = max_abs_poly_in_c2(cInt, coeff01)
% coeff01 = [a0,a1] for |a0 + a1*c^2|.
[tlo, thi] = c2_range(cInt);
vals = abs(coeff01(1) + coeff01(2)*[tlo thi]);
if coeff01(2) ~= 0
    t0 = -coeff01(1)/coeff01(2);
    if t0 >= tlo && t0 <= thi
        vals(end+1) = 0; %#ok<AGROW>
    end
end
M = max(vals);
end

function [tlo, thi] = c2_range(cInt)
loC = cInt(1); hiC = cInt(2);
thi = max(loC^2, hiC^2);
if loC <= 0 && hiC >= 0
    tlo = 0;
else
    tlo = min(loC^2, hiC^2);
end
end

function M = max_t_times_one_minus_t2(tlo, thi)
thi = min(max(thi,0),1);
tlo = min(max(tlo,0),thi);
cand = [tlo thi];
t = 1/sqrt(3);
if t >= tlo && t <= thi
    cand(end+1) = t; %#ok<AGROW>
end
M = max(cand .* (1 - cand.^2));
end

function M = max_t2_times_one_minus_t2(tlo, thi)
thi = min(max(thi,0),1);
tlo = min(max(tlo,0),thi);
cand = [tlo thi];
t = 1/sqrt(2);
if t >= tlo && t <= thi
    cand(end+1) = t; %#ok<AGROW>
end
M = max(cand.^2 .* (1 - cand.^2));
end

%% =====================================================================
% Actual pointwise B_best
%% =====================================================================
function info = actual_pointwise_Bbest(coeff, omega1, omega2, opts)
candidates = [];
labels = {};

if opts.useL1Fallback
    candidates(end+1) = coeff.l1; %#ok<AGROW>
    labels{end+1} = 'actual-l1'; %#ok<AGROW>
end
if opts.useDerivativeBranch
    candidates(end+1) = derivative_bound_pointwise_scalar(coeff, omega1, omega2); %#ok<AGROW>
    labels{end+1} = 'actual-derivative'; %#ok<AGROW>
end
if opts.useResidueCorrelationBranch
    if get_struct_opt(opts, 'useFullPhaseResidueCorrelationBranch', true)
        candidates(end+1) = fullphase_residue_correlation_pointwise_scalar(coeff, omega1, omega2, opts); %#ok<AGROW>
        labels{end+1} = 'actual-fprc'; %#ok<AGROW>
    end
end
if opts.useResidueLinearBranch
    candidates(end+1) = residue_linear_pointwise_scalar(coeff.a, omega1, omega2, opts); %#ok<AGROW>
    labels{end+1} = 'actual-reslin'; %#ok<AGROW>
end

[value, label] = min_finite_candidates(candidates, labels);
info = make_bound_info(value, label, inf);
info.l1 = coeff.l1;
end

function B = derivative_bound_pointwise_scalar(coeff, omega1, omega2)
N = coeff.N;
if N < 9 || ~coeff.flatEndsR4
    B = inf;
    return;
end

dN = exact_dN_plus(N, omega1, omega2);
if dN <= 0 || ~isfinite(dN)
    B = inf;
    return;
end

s0 = sin(dN/2);
rho = abs(sin(omega2));
if s0 <= 0 || ~isfinite(s0)
    B = inf;
    return;
end

beta = derivative_beta_from_s0_rho(s0, rho);
B = derivative_class_bound_scalar(beta, coeff);
end


function Bbest = fullphase_residue_correlation_pointwise_scalar(coeff, omega1, omega2, opts)
% Full-phase residue-correlation pointwise envelope.
% For Q>=2 this is sum_s |T_s| after the residue split n=s+Qm, hence it
% keeps the full lag phase inside each residue class and only uses the
% triangle inequality across residue classes.  If fullPhaseCorrIncludeQ1 is
% true, Q=1 is allowed and the branch becomes exact pointwise evaluation.
a = coeff.a(:);
N = coeff.N;
l1 = coeff.l1;
Bbest = l1;

includeQ1 = get_struct_opt(opts, 'fullPhaseCorrIncludeQ1', false);
if includeQ1
    Qstart = 1;
else
    Qstart = 2;
end
Qmax = min(N, get_struct_opt(opts, 'residueQMax', N));
safety = get_struct_opt(opts, 'fullPhaseCorrSafety', 0);

for Q = Qstart:Qmax
    BQ = fixedQ_fullphase_residue_correlation_pointwise_scalar(a, Q, omega1, omega2);
    Bbest = min(Bbest, BQ);
end

Bbest = max(real(min(l1, Bbest) + safety), 0);
end

function B = fixedQ_fullphase_residue_correlation_pointwise_scalar(a, Q, omega1, omega2)
% Pointwise full-phase residue split.  The sub-sums are evaluated exactly;
% this is equivalent to the full autocorrelation identity on each residue
% class, but is numerically more stable.
a = a(:);
N = numel(a);
if Q < 1 || Q > N
    B = sum(abs(a));
    return;
end

gamma = Q^2 * omega2;
B = 0;
for s0 = 0:(Q-1)
    sub = a((s0+1):Q:N);
    M = numel(sub);
    if M == 0
        continue;
    end
    beta = Q*omega1 + 2*Q*s0*omega2;
    m = (0:(M-1)).';
    Ts = sum(sub(:).*exp(1i*(beta*m + gamma*(m.^2))));
    B = B + abs(Ts);
end
B = max(real(B),0);
end

function B = fixedQ_fullphase_residue_correlation_cell_scalar(coeff, Q, phaseBox, pList, corrOpts)
% Support-uniform full-phase residue-correlation envelope on a phase box.
% It uses, for every residue class and lag, the hybrid upper bound
%     min{ phase-aware interval-cosine real-part bound,
%          Abel denominator magnitude bound }.
% This is a genuine upper bound over the full phase rectangle.
a = coeff.a(:);
N = coeff.N;
l1 = coeff.l1;

if Q < 1 || Q > N || ~isfinite(Q)
    B = l1;
    return;
end

w1Int = sort(phaseBox.omega1(:).');
w2Int = sort(phaseBox.omega2(:).');
if numel(w1Int) ~= 2 || numel(w2Int) ~= 2 || any(~isfinite([w1Int w2Int]))
    B = l1;
    return;
end

H_Q = ceil(N/Q) - 1;
if H_Q < 1
    B = l1;
    return;
end

autoDeltaSafety = get_struct_opt(corrOpts, 'autoDeltaSafety', 0.98);
deltaProfile = delta_profile_from_omega2_interval(H_Q, Q, w2Int(1), w2Int(2), autoDeltaSafety);

statsCell = [];
if isKey(coeff.Qstats, Q)
    statsCell = coeff.Qstats(Q);
end

Btotal = 0;
for s0 = 0:(Q-1)
    sub = a((s0+1):Q:N);
    M = numel(sub);
    if M == 0
        continue;
    end

    E = sum(abs(sub).^2);
    Fupper = E;

    if ~isempty(statsCell) && numel(statsCell) >= s0+1
        stats = statsCell{s0+1};
    else
        stats = coefficient_stats_scalar(sub, max(pList));
    end

    for h = 1:(M-1)
        Ucos = 0;
        for m0 = 0:(M-1-h)
            dlag = sub(m0+h+1) * conj(sub(m0+1));
            amp = abs(dlag);
            if amp == 0
                continue;
            end

            beta1 = h*Q;
            beta2 = 2*Q*s0*h + Q^2*(h^2 + 2*h*m0);
            J = affine_phase_interval(angle(dlag), beta1, w1Int, beta2, w2Int);
            Ucos = Ucos + amp * interval_cos_max(J(1), J(2));
        end

        Umag = fprc_lag_magnitude_bound_from_stats(stats, h, deltaProfile, pList, corrOpts);
        Uh = min(Ucos, Umag);
        Fupper = Fupper + 2*Uh;
    end

    Btotal = Btotal + sqrt(max(real(Fupper),0));
    if Btotal >= l1
        Btotal = l1;
        break;
    end
end

safety = get_struct_opt(corrOpts, 'fullPhaseCorrSafety', 0);
B = max(real(min(l1, Btotal) + safety),0);
end

function U = fprc_lag_magnitude_bound_from_stats(stats, h, deltaProfile, pList, corrOpts)
% Magnitude control for one lag correlation.  This is the Abel
% denominator estimate, but used only as one half of the new hybrid real-part
% majorant.
if h < 1 || h > numel(stats.C)
    U = 0;
    return;
end

U = stats.C(h);
if h <= numel(deltaProfile)
    dh = deltaProfile(h);
    autoDeltaMin = get_struct_opt(corrOpts, 'autoDeltaMin', 1e-12);
    if dh > autoDeltaMin && dh <= pi/2 && isfinite(dh)
        den = 2*sin(dh);
        for p = pList(:).'
            if p <= numel(stats.W)
                Whp = stats.W{p}(h);
                if isfinite(Whp)
                    U = min(U, Whp/den^p);
                end
            end
        end
    end
end
U = max(real(U),0);
end

function J = affine_phase_interval(offset, beta1, I1, beta2, I2)
% Exact interval hull for offset + beta1*w1 + beta2*w2 over a rectangle.
vals = [ ...
    offset + beta1*I1(1) + beta2*I2(1), ...
    offset + beta1*I1(1) + beta2*I2(2), ...
    offset + beta1*I1(2) + beta2*I2(1), ...
    offset + beta1*I1(2) + beta2*I2(2)];
J = [min(vals), max(vals)];
end

function cmax = interval_cos_max(lo, hi)
% Supremum of cos(x) over [lo,hi].
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if hi - lo >= 2*pi
    cmax = 1;
    return;
end
vals = [cos(lo), cos(hi)];
kmin = ceil(lo/(2*pi));
kmax = floor(hi/(2*pi));
if kmin <= kmax
    vals(end+1) = 1; %#ok<AGROW>
end
cmax = max(vals);
cmax = max(min(real(cmax),1),-1);
end

function Bbest = residue_correlation_pointwise_scalar(coeff, omega2, pList, opts)
l1 = coeff.l1;
Bbest = l1;

if opts.residueCorrIncludeQ1
    Qstart = 1;
else
    Qstart = 2;
end

for Q = Qstart:opts.residueQMax
    BQ = fixedQ_corr_pointwise_precomputed_scalar(coeff, Q, omega2, pList);
    Bbest = min(Bbest, BQ);
end

Bbest = max(real(min(l1, Bbest)),0);
end

function Bbest = residue_linear_pointwise_scalar(a, omega1, omega2, opts)
a = a(:);
N = numel(a);
l1 = sum(abs(a));
Bbest = l1;

for q = 1:opts.residueQMax
    for Ares = 0:(2*q-1)
        epsSigned = signed_dist_to_lattice(omega2 - pi*Ares/q, 2*pi);
        nu = signed_dist_to_lattice(q^2*epsSigned, 2*pi);

        if opts.residueUseNuCap && abs(nu) > opts.residueNuMax
            continue;
        end

        Bq = 0;
        for s = 0:(q-1)
            sub = a((s+1):q:N);
            if isempty(sub)
                continue;
            end

            Omega = q*omega1 + pi*Ares*q + 2*q*epsSigned*s;
            Bs = linear_abel_bound_scalar_exact_residual(sub, Omega, nu, opts.residueMinDen);
            Bq = Bq + Bs;

            if Bq >= Bbest
                break;
            end
        end

        Bbest = min(Bbest, Bq);
    end
end

Bbest = max(real(min(l1, Bbest)),0);
end

function B = linear_abel_bound_scalar_exact_residual(aSub, Omega, nu, minDen)
aSub = aSub(:);
M = numel(aSub)-1;
l1 = sum(abs(aSub));

if M < 0
    B = 0;
    return;
end
if M == 0
    B = l1;
    return;
end

den = 2*abs(sin(Omega/2));
if den <= minDen || ~isfinite(den)
    B = l1;
    return;
end

m = (0:M).';
bseq = aSub .* exp(1i*nu*(m.^2));
variation = abs(bseq(1)) + abs(bseq(end)) + sum(abs(diff(bseq)));
B = max(real(min(l1, variation/den)),0);
end

%% =====================================================================
% Scalar coefficient preprocessing and scalar phase-class branches
%% =====================================================================
function coeff = precompute_scalar_coeff_data(a, pMax, Qall)
a = a(:);
N = numel(a);
coeff = struct();
coeff.kind = 'scalar';
coeff.a = a;
coeff.N = N;
coeff.l1 = sum(abs(a));
coeff.diffNorms0to4 = diff_norms_0_to_4_scalar(a);
coeff.flatEndsR4 = N >= 9 && all(abs(a(1:4)) <= 1e-13) && all(abs(a(end-3:end)) <= 1e-13);
coeff.pMax = pMax;
coeff.Qall = Qall;
coeff.Qstats = containers.Map('KeyType','double','ValueType','any');

for Q = unique(Qall(:).')
    statsCell = cell(1,Q);
    for s = 1:Q
        sub = a(s:Q:N);
        statsCell{s} = coefficient_stats_scalar(sub, pMax);
    end
    coeff.Qstats(Q) = statsCell;
end
end

function stats = coefficient_stats_scalar(a, pMax)
a = a(:);
N = numel(a);
stats = struct();
stats.N = N;
stats.l1 = sum(abs(a));
stats.E = sum(abs(a).^2);

if N <= 1
    stats.C = zeros(0,1);
    stats.W = cell(pMax,1);
    for p = 1:pMax
        stats.W{p} = zeros(0,1);
    end
    return;
end

stats.C = zeros(N-1,1);
stats.W = cell(pMax,1);
for p = 1:pMax
    stats.W{p} = inf(N-1,1);
end

for h = 1:(N-1)
    c = a(1+h:N).*conj(a(1:N-h));
    [Csum, Wrow] = variation_profile_for_scalar_sequence(c, pMax);
    stats.C(h) = Csum;
    for p = 1:pMax
        stats.W{p}(h) = Wrow(p);
    end
end
end

function B = derivative_class_bound_scalar(beta, coeff)
B = sum(beta(:).'.*coeff.diffNorms0to4(:).');
B = max(real(min(coeff.l1, B)),0);
end

function B = fixedQ_corr_pointwise_precomputed_scalar(coeff, Q, omega2, pList)
if ~isKey(coeff.Qstats, Q)
    B = coeff.l1;
    return;
end

statsCell = coeff.Qstats(Q);
total = 0;

for s = 1:numel(statsCell)
    stats = statsCell{s};
    Bsub = corr_from_stats_pointwise(stats, Q^2*omega2, pList);
    total = total + min(stats.l1, Bsub);
    if total >= coeff.l1
        break;
    end
end

B = max(real(min(coeff.l1, total)),0);
end

function B = fixedQ_corr_delta_precomputed_scalar(coeff, Q, deltaProfile, pList)
if ~isKey(coeff.Qstats, Q)
    B = coeff.l1;
    return;
end

statsCell = coeff.Qstats(Q);
total = 0;

for s = 1:numel(statsCell)
    stats = statsCell{s};
    Bsub = corr_from_stats_delta(stats, deltaProfile, pList);
    total = total + min(stats.l1, Bsub);
    if total >= coeff.l1
        break;
    end
end

B = max(real(min(coeff.l1, total)),0);
end

function B = fixed_q_residue_linear_sine_scalar_bound(a, q, nuBar, gammaProfile)
a = a(:);
N = numel(a);
l1 = sum(abs(a));
nuBar = min(pi, abs(nuBar));
Bq = 0;

for s = 0:(q-1)
    sub = a((s+1):q:N);
    if isempty(sub)
        continue;
    end

    Ls = sum(abs(sub));
    M = numel(sub)-1;
    gamma = gammaProfile(s+1);
    den = 2*gamma;

    if den <= 0 || ~isfinite(den) || M == 0
        Bs = Ls;
    else
        Vs = abs(sub(1)) + abs(sub(end));
        for m = 0:(M-1)
            T = min(pi, nuBar*(2*m+1));
            Vs = Vs + scalar_phase_jump_envelope(sub(m+1), sub(m+2), T);
        end
        Bs = min(Ls, Vs/den);
    end

    Bq = Bq + Bs;
    if Bq >= l1
        break;
    end
end

B = max(real(min(l1, Bq)),0);
end

function B = corr_from_stats_pointwise(stats, omega2Eff, pList)
sq = stats.E;
for h = 1:(stats.N-1)
    term = stats.C(h);
    sh = abs(sin(h*omega2Eff));
    if sh > 1e-15
        den = 2*sh;
        for p = pList(:).'
            Whp = stats.W{p}(h);
            if isfinite(Whp)
                term = min(term, Whp/den^p);
            end
        end
    end
    sq = sq + 2*term;
end
B = sqrt(max(real(sq),0));
end

function B = corr_from_stats_delta(stats, deltaProfile, pList)
sq = stats.E;
for h = 1:(stats.N-1)
    term = stats.C(h);
    if h <= numel(deltaProfile)
        dh = deltaProfile(h);
        if dh > 0 && isfinite(dh)
            den = 2*sin(dh);
            for p = pList(:).'
                Whp = stats.W{p}(h);
                if isfinite(Whp)
                    term = min(term, Whp/den^p);
                end
            end
        end
    end
    sq = sq + 2*term;
end
B = sqrt(max(real(sq),0));
end

function [Csum, Wrow] = variation_profile_for_scalar_sequence(c, pMax)
c = c(:);
Csum = sum(abs(c));
M = numel(c)-1;
Wrow = inf(1,pMax);

if isempty(c)
    Wrow(:) = 0;
    return;
end

diffs = cell(pMax+1,1);
diffs{1} = c;
for p = 1:pMax
    if numel(diffs{p}) >= 2
        diffs{p+1} = diff(diffs{p});
    else
        diffs{p+1} = [];
    end
end

for p = 1:pMax
    if p > M+1
        continue;
    end

    val = 0;
    for r = 0:(p-1)
        dr = diffs{r+1};
        if isempty(dr)
            val = inf;
            break;
        end
        val = val + 2^(p-1-r)*(abs(dr(1)) + abs(dr(end)));
    end

    dp = diffs{p+1};
    if isempty(dp) && (M-p+1 > 0)
        val = inf;
    elseif ~isempty(dp)
        val = val + sum(abs(dp));
    end

    Wrow(p) = val;
end
end

%% =====================================================================
% Geometry-aware normalized fallback helpers kept for compatibility
%% =====================================================================
function S = sum_poly1d_abs_over_q(Cpoly, G, I)
S = 0;
for m = 1:size(Cpoly,1)
    S = S + max_abs_quad_over_q_exact(Cpoly(m,:), G, I);
end
end

function Wrow = variation_profile_for_poly1d_over_q(Cpoly, G, I, pMax)
M = size(Cpoly,1)-1;
Wrow = inf(1,pMax);

if isempty(Cpoly)
    Wrow(:) = 0;
    return;
end

diffs = cell(pMax+1,1);
diffs{1} = Cpoly;
for p = 1:pMax
    if size(diffs{p},1) >= 2
        diffs{p+1} = diff(diffs{p},1,1);
    else
        diffs{p+1} = zeros(0,3);
    end
end

for p = 1:pMax
    if p > M+1
        continue;
    end

    val = 0;
    for r = 0:(p-1)
        dr = diffs{r+1};
        if size(dr,1) == 0
            val = inf;
            break;
        end
        val = val + 2^(p-1-r)*( ...
            max_abs_quad_over_q_exact(dr(1,:),G,I) + ...
            max_abs_quad_over_q_exact(dr(end,:),G,I));
    end

    dp = diffs{p+1};
    if size(dp,1) == 0 && (M-p+1 > 0)
        val = inf;
    elseif size(dp,1) > 0
        val = val + sum_poly1d_abs_over_q(dp,G,I);
    end

    Wrow(p) = val;
end
end

function S = sum_poly2d_abs_over_qprod(Cpoly, G, Ie, Is, tauOpts)
S = 0;
for m = 1:size(Cpoly,1)
    S = S + max_abs_poly22_over_qprod_adaptive(squeeze(Cpoly(m,:,:)), G, Ie, Is, tauOpts);
end
end

function Wrow = variation_profile_for_poly2d_over_qprod(Cpoly, G, Ie, Is, pMax, tauOpts)
M = size(Cpoly,1)-1;
Wrow = inf(1,pMax);

if isempty(Cpoly)
    Wrow(:) = 0;
    return;
end

diffs = cell(pMax+1,1);
diffs{1} = Cpoly;
for p = 1:pMax
    if size(diffs{p},1) >= 2
        diffs{p+1} = diff(diffs{p},1,1);
    else
        diffs{p+1} = zeros(0,3,3);
    end
end

for p = 1:pMax
    if p > M+1
        continue;
    end

    val = 0;
    for r = 0:(p-1)
        dr = diffs{r+1};
        if size(dr,1) == 0
            val = inf;
            break;
        end
        val = val + 2^(p-1-r)*( ...
            max_abs_poly22_over_qprod_adaptive(squeeze(dr(1,:,:)),G,Ie,Is,tauOpts) + ...
            max_abs_poly22_over_qprod_adaptive(squeeze(dr(end,:,:)),G,Ie,Is,tauOpts));
    end

    dp = diffs{p+1};
    if size(dp,1) == 0 && (M-p+1 > 0)
        val = inf;
    elseif size(dp,1) > 0
        val = val + sum_poly2d_abs_over_qprod(dp,G,Ie,Is,tauOpts);
    end

    Wrow(p) = val;
end
end

function M = max_abs_linear_over_sqrtq(a, b, G, I)
lo = min(I); hi = max(I);
cand = [lo hi];

t = stationary_linear_over_sqrtq(a,b,G);
if isfinite(t) && t >= lo && t <= hi
    cand(end+1) = t; %#ok<AGROW>
end

if abs(b) > 1e-14
    z = -a/b;
    if z >= lo && z <= hi
        cand(end+1) = z; %#ok<AGROW>
    end
end

vals = abs(a + b*cand)./sqrt(q_tau(G,cand));
M = max(vals);
M = max(real(M),0);
end

function B = sum_abs_weighted_linear_over_sqrtq_exact(ax, ay, G, I)
ax = ax(:);
ay = ay(:);
I = sort(I(:).');
loI = I(1);
hiI = I(2);

if hiI <= loI + 1e-15
    tau0 = 0.5*(loI + hiI);
    B = sum(abs(ax + tau0*ay))/sqrt(q_tau(G, tau0));
    B = max(real(B),0);
    return;
end

breaks = [loI hiI];

for m = 1:numel(ax)
    if abs(ay(m)) > 1e-14
        z = -ax(m)/ay(m);
        if z > loI && z < hiI
            breaks(end+1) = z; %#ok<AGROW>
        end
    end
end

breaks = unique(sort(breaks));
B = 0;

for r = 1:(numel(breaks)-1)
    lo = breaks(r);
    hi = breaks(r+1);
    if hi <= lo + 1e-15
        continue;
    end

    mid = 0.5*(lo+hi);
    sgn = sign(ax + mid*ay);
    sgn(sgn == 0) = 1;

    A = sum(sgn(:).*ax(:));
    C = sum(sgn(:).*ay(:));
    B = max(B, max_abs_linear_over_sqrtq(A, C, G, [lo hi]));
end

B = max(real(B),0);
end

function M = max_abs_quad_over_q_exact(p, G, I)
lo = min(I); hi = max(I);
q = [G(1,1), 2*G(1,2), G(2,2)];
cand = [lo hi];

r = poly_sub_asc(poly_conv_asc(poly_der_asc(p), q), poly_conv_asc(p, poly_der_asc(q)));
rootsCand = roots_descending_asc(r);

for k = 1:numel(rootsCand)
    t = real(rootsCand(k));
    if abs(imag(rootsCand(k))) <= 1e-10 && t >= lo && t <= hi
        cand(end+1) = t; %#ok<AGROW>
    end
end

vals = abs(poly_eval_asc(p,cand))./q_tau(G,cand);
M = max(vals);
M = max(real(M),0);
end

function M = max_abs_bilin_over_sqrtqprod_tight(P, G, Ie, Is, tauOpts)
if nargin < 5 || ~isfield(tauOpts, 'tightBilin2D') || ~tauOpts.tightBilin2D
    M = max_abs_bilin_over_sqrtqprod_crude(P, G, Ie, Is);
    return;
end

R = poly2d_mul_bilin(P, P);
opts2 = tauOpts;
opts2.maxDepth2D = get_struct_opt(tauOpts, 'bilinMaxDepth2D', get_struct_opt(tauOpts,'maxDepth2D',2));
opts2.maxCells2D = get_struct_opt(tauOpts, 'bilinMaxCells2D', get_struct_opt(tauOpts,'maxCells2D',64));
opts2.tolAbs2D = get_struct_opt(tauOpts, 'bilinTolAbs2D', get_struct_opt(tauOpts,'tolAbs2D',1e-10));
opts2.tolRel2D = get_struct_opt(tauOpts, 'bilinTolRel2D', get_struct_opt(tauOpts,'tolRel2D',1e-7));
opts2.sampleOrder2D = get_struct_opt(tauOpts, 'bilinSampleOrder', 5);

U = max_abs_poly22_over_qprod_adaptive(R, G, Ie, Is, opts2);
M = sqrt(max(real(U),0));
end

function M = max_abs_bilin_over_sqrtqprod_crude(P, G, Ie, Is)
te = [min(Ie), max(Ie)];
ts = [min(Is), max(Is)];
vals = zeros(4,1);
k = 0;

for i = 1:2
    for j = 1:2
        k = k + 1;
        vals(k) = abs(P(1) + P(2)*te(i) + P(3)*ts(j) + P(4)*te(i)*ts(j));
    end
end

numMax = max(vals);
denMin = sqrt(min_q_interval(G,Ie)*min_q_interval(G,Is));
M = numMax/max(denMin, realmin);
M = max(real(M),0);
end

function M = max_abs_poly22_over_qprod_adaptive(P, G, Ie, Is, tauOpts)
maxDepth = get_struct_opt(tauOpts, 'maxDepth2D', 2);
maxCells = get_struct_opt(tauOpts, 'maxCells2D', 128);
tolAbs = get_struct_opt(tauOpts, 'tolAbs2D', 1e-10);
tolRel = get_struct_opt(tauOpts, 'tolRel2D', 1e-7);
sampleOrder = get_struct_opt(tauOpts, 'sampleOrder2D', 5);

bestLower = 0;
leafUpper = 0;
numCells = 0;
stack = {struct('Ie', sort(Ie), 'Is', sort(Is), 'depth', 0)};

while ~isempty(stack)
    item = stack{end};
    stack(end) = [];
    numCells = numCells + 1;

    ub = poly22_rect_upper_bernstein(P, G, item.Ie, item.Is);
    lb = poly22_rect_sample_lower(P, G, item.Ie, item.Is, sampleOrder);
    bestLower = max(bestLower, lb);

    if item.depth >= maxDepth || numCells >= maxCells || ub <= bestLower + tolAbs + tolRel*max(1,bestLower)
        leafUpper = max(leafUpper, ub);
    else
        we = item.Ie(2)-item.Ie(1);
        ws = item.Is(2)-item.Is(1);
        if we >= ws
            mid = 0.5*(item.Ie(1)+item.Ie(2));
            stack{end+1} = struct('Ie',[item.Ie(1) mid], 'Is',item.Is, 'depth',item.depth+1); %#ok<AGROW>
            stack{end+1} = struct('Ie',[mid item.Ie(2)], 'Is',item.Is, 'depth',item.depth+1); %#ok<AGROW>
        else
            mid = 0.5*(item.Is(1)+item.Is(2));
            stack{end+1} = struct('Ie',item.Ie, 'Is',[item.Is(1) mid], 'depth',item.depth+1); %#ok<AGROW>
            stack{end+1} = struct('Ie',item.Ie, 'Is',[mid item.Is(2)], 'depth',item.depth+1); %#ok<AGROW>
        end
    end
end

M = max(leafUpper, bestLower);
M = max(real(M),0);
end

function ub = poly22_rect_upper_bernstein(P, G, Ie, Is)
Bcoef = monomial22_to_bernstein22_on_rect(P, Ie, Is);
numUpper = max(abs(Bcoef(:)));
denLower = min_q_interval(G,Ie)*min_q_interval(G,Is);
ub = numUpper/max(denLower, realmin);
end

function lb = poly22_rect_sample_lower(P, G, Ie, Is, sampleOrder)
sampleOrder = max(3, ceil(sampleOrder));
k = 0:(sampleOrder-1);
z = 0.5*(1 - cos(pi*k/(sampleOrder-1)));
te = unique([Ie(1) + (Ie(2)-Ie(1))*z, 0.5*(Ie(1)+Ie(2))]);
ts = unique([Is(1) + (Is(2)-Is(1))*z, 0.5*(Is(1)+Is(2))]);

lb = 0;
for i = 1:numel(te)
    qe = q_tau(G,te(i));
    for j = 1:numel(ts)
        val = abs(poly22_eval(P,te(i),ts(j)))/(qe*q_tau(G,ts(j)));
        lb = max(lb, val);
    end
end
end

function B = monomial22_to_bernstein22_on_rect(P, Ie, Is)
a = Ie(1); b = Ie(2); we = b-a;
c = Is(1); d = Is(2); ws = d-c;
Cuv = zeros(3,3);

for i = 0:2
    for j = 0:2
        coeff = P(i+1,j+1);
        for p = 0:i
            ce = nchoosek(i,p)*a^(i-p)*we^p;
            for q = 0:j
                cs = nchoosek(j,q)*c^(j-q)*ws^q;
                Cuv(p+1,q+1) = Cuv(p+1,q+1) + coeff*ce*cs;
            end
        end
    end
end

M = bernstein_power_matrix(2);
B = M \ Cuv / M.';
end

function M = bernstein_power_matrix(n)
M = zeros(n+1,n+1);
for i = 0:n
    for k = i:n
        M(k+1,i+1) = nchoosek(n,i)*nchoosek(n-i,k-i)*(-1)^(k-i);
    end
end
end

function v = poly22_eval(P, te, ts)
v = 0;
for i = 0:2
    for j = 0:2
        v = v + P(i+1,j+1)*te^i*ts^j;
    end
end
end

function P = linprod_poly(a0,a1,b0,b1)
P = [a0*b0, a0*b1 + a1*b0, a1*b1];
end

function R = poly2d_mul_bilin(P, Q)
R = zeros(3,3);
termsP = [0 0 P(1); 1 0 P(2); 0 1 P(3); 1 1 P(4)];
termsQ = [0 0 Q(1); 1 0 Q(2); 0 1 Q(3); 1 1 Q(4)];

for i = 1:4
    for j = 1:4
        de = termsP(i,1) + termsQ(j,1);
        ds = termsP(i,2) + termsQ(j,2);
        R(de+1,ds+1) = R(de+1,ds+1) + termsP(i,3)*termsQ(j,3);
    end
end
end

function qmin = min_q_interval(G, I)
lo = min(I); hi = max(I);
cand = [lo hi];
if G(2,2) > 0
    tv = -G(1,2)/G(2,2);
    if tv >= lo && tv <= hi
        cand(end+1) = tv; %#ok<AGROW>
    end
end
qmin = min(q_tau(G,cand));
qmin = max(qmin, realmin);
end

function tstar = stationary_linear_over_sqrtq(A, C, G)
num0 = C*G(1,1) - A*G(1,2);
num1 = C*G(1,2) - A*G(2,2);
if abs(num1) <= 1e-14
    tstar = NaN;
else
    tstar = -num0/num1;
end
end

function q = q_tau(G, tau)
q = G(1,1) + 2*G(1,2)*tau + G(2,2)*tau.^2;
end

function h = physical_h_vec(geom, tau)
q = q_tau(geom.G, tau);
if any(q <= 0) || any(~isfinite(q))
    error('Invalid physical q(tau).');
end
h = (geom.x(:) + tau*geom.y(:))/sqrt(q);
end

%% =====================================================================
% Phase support specification and class profile helpers
%% =====================================================================
function spec = build_support_spec(Nr, classOpts, k_lambda, d, whichSlice, rEvalFixed, thetaEvalFixed)
spec = struct();
spec.enable = logical(classOpts.enable);
spec.phaseBox = phase_box_for_support_eval_geometry( ...
    k_lambda, d, classOpts.rSuppRange, classOpts.thetaSuppRange, ...
    classOpts.rEvalRange, classOpts.thetaEvalRange);
spec.tauSupp = tau_range_from_geometry(d, classOpts.rSuppRange, classOpts.thetaSuppRange);
spec.tauEval = tau_range_from_geometry(d, classOpts.rEvalRange, classOpts.thetaEvalRange);
spec.alphaEvalRange = alpha_range_from_r(d, classOpts.rEvalRange);
spec.thetaEvalRange = sort(classOpts.thetaEvalRange(:).');
spec.rEvalRange = sort(classOpts.rEvalRange(:).');
spec.rSuppRange = sort(classOpts.rSuppRange(:).');
spec.thetaSuppRange = sort(classOpts.thetaSuppRange(:).');
spec.kappa = k_lambda*d;
spec.k_lambda = k_lambda;
spec.d = d;
spec.Nr = Nr;
spec.whichSlice = whichSlice;
spec.rEvalFixed = rEvalFixed;
spec.thetaEvalFixed = thetaEvalFixed;

spec.dNplusLower = support_dN_plus_lower_curveboxed( ...
    Nr, k_lambda, d, classOpts, whichSlice, rEvalFixed, thetaEvalFixed, ...
    classOpts.derivative.curveDNCells);
end

function spec = disabled_support_spec()
spec = struct();
spec.enable = false;
spec.phaseBox = struct('omega1',[NaN NaN], 'omega2',[NaN NaN]);
spec.tauSupp = [NaN NaN];
spec.tauEval = [NaN NaN];
spec.alphaEvalRange = [NaN NaN];
spec.thetaEvalRange = [NaN NaN];
spec.rEvalRange = [NaN NaN];
spec.rSuppRange = [NaN NaN];
spec.thetaSuppRange = [NaN NaN];
spec.kappa = NaN;
spec.k_lambda = NaN;
spec.d = NaN;
spec.Nr = NaN;
spec.whichSlice = 'none';
spec.rEvalFixed = NaN;
spec.thetaEvalFixed = NaN;
spec.dNplusLower = 0;
end

function ar = alpha_range_from_r(d, rRange)
rRange = sort(rRange(:).');
ar = [d/rRange(2), d/rRange(1)];
end

function dlow = support_dN_plus_lower_curveboxed(Nr, k_lambda, d, classOpts, whichSlice, rEvalFixed, thetaEvalFixed, numCells)
numCells = max(1, ceil(numCells));

switch lower(char(whichSlice))
    case 'theta'
        a = classOpts.thetaEvalRange(1);
        b = classOpts.thetaEvalRange(2);
        edges = linspace(a,b,numCells+1);
        dlow = inf;
        for k = 1:numCells
            tmp = classOpts;
            tmp.rEvalRange = [rEvalFixed rEvalFixed];
            tmp.thetaEvalRange = [edges(k) edges(k+1)];
            box = phase_box_for_support_eval_geometry(k_lambda, d, tmp.rSuppRange, tmp.thetaSuppRange, tmp.rEvalRange, tmp.thetaEvalRange);
            dlow = min(dlow, dN_plus_lower_bound_box(Nr, box.omega1, box.omega2));
        end
    case 'range'
        a = classOpts.rEvalRange(1);
        b = classOpts.rEvalRange(2);
        edges = linspace(a,b,numCells+1);
        dlow = inf;
        for k = 1:numCells
            tmp = classOpts;
            tmp.rEvalRange = [edges(k) edges(k+1)];
            tmp.thetaEvalRange = [thetaEvalFixed thetaEvalFixed];
            box = phase_box_for_support_eval_geometry(k_lambda, d, tmp.rSuppRange, tmp.thetaSuppRange, tmp.rEvalRange, tmp.thetaEvalRange);
            dlow = min(dlow, dN_plus_lower_bound_box(Nr, box.omega1, box.omega2));
        end
    otherwise
        box = phase_box_for_support_eval_geometry(k_lambda, d, classOpts.rSuppRange, classOpts.thetaSuppRange, classOpts.rEvalRange, classOpts.thetaEvalRange);
        dlow = dN_plus_lower_bound_box(Nr, box.omega1, box.omega2);
end

if ~isfinite(dlow)
    dlow = 0;
end
end

function AList = get_A_list_for_q_support(lin, q, phaseBox)
mode = lower(char(get_struct_opt(lin, 'AMode', 'phase')));

switch mode
    case 'all'
        AList = 0:(2*q-1);

    case 'phase'
        window = get_struct_opt(lin, 'AWindow', 0.08);
        lo = phaseBox.omega2(1);
        hi = phaseBox.omega2(2);
        AList = [];

        for Ares = 0:(2*q-1)
            dmin = interval_dist_to_lattice_min(lo - pi*Ares/q, hi - pi*Ares/q, 2*pi);
            if dmin <= window + 1e-14
                AList(end+1) = Ares; %#ok<AGROW>
            end
        end

        if isempty(AList)
            center = 0.5*(lo+hi);
            bestA = 0;
            bestD = inf;
            for Ares = 0:(2*q-1)
                dA = dist_to_lattice(center - pi*Ares/q, 2*pi);
                if dA < bestD
                    bestD = dA;
                    bestA = Ares;
                end
            end
            AList = bestA;
        end

    otherwise
        error('Unsupported lin.AMode. Use phase or all.');
end

AList = unique(AList(:).');
end

function dp = delta_profile_from_omega2_interval(H_Q, Q, omega2Lo, omega2Hi, safety)
dp = zeros(1,H_Q);
for h = 1:H_Q
    dp(h) = safety * interval_dist_to_lattice_min(h*Q^2*omega2Lo, h*Q^2*omega2Hi, pi);
end
end

function gammaProfile = gamma_profile_lower_from_phase_box(q, Ares, epsBar, phaseBox)
alphaLower = alpha_profile_lower_from_phase_box(q, Ares, epsBar, phaseBox);
gammaProfile = max(0, min(1, sin(alphaLower/2)));
end

function alphaProfile = alpha_profile_lower_from_phase_box(q, Ares, epsBar, phaseBox)
alphaProfile = zeros(1,q);
w1Lo = phaseBox.omega1(1);
w1Hi = phaseBox.omega1(2);

for s = 0:(q-1)
    lo = q*w1Lo + pi*Ares*q - 2*q*s*epsBar;
    hi = q*w1Hi + pi*Ares*q + 2*q*s*epsBar;
    alphaProfile(s+1) = interval_dist_to_lattice_min(lo, hi, 2*pi);
end
end

function [epsBar, ok] = eps_bar_from_omega2_interval(q, Ares, omega2Lo, omega2Hi)
lo = omega2Lo - pi*Ares/q;
hi = omega2Hi - pi*Ares/q;
epsBar = interval_dist_to_lattice_max(lo, hi, 2*pi);
ok = isfinite(epsBar) && epsBar <= pi;
end

%% =====================================================================
% Phase geometry and intervals
%% =====================================================================
function [omega1, omega2] = phase_pair(k_lambda, d, r_s, theta_s, r_e, theta_e)
muS = k_lambda*d*cos(theta_s);
muE = k_lambda*d*cos(theta_e);
etaS = (k_lambda*d^2/(2*r_s))*sin(theta_s)^2;
etaE = (k_lambda*d^2/(2*r_e))*sin(theta_e)^2;
omega1 = muS - muE;
omega2 = etaE - etaS;
end

function box = phase_box_for_support_eval_geometry(k_lambda, d, rSuppRange, thetaSuppRange, rEvalRange, thetaEvalRange)
[muSlo, muShi] = mu_interval(k_lambda, d, thetaSuppRange);
[muElo, muEhi] = mu_interval(k_lambda, d, thetaEvalRange);
[etaSlo, etaShi] = eta_interval(k_lambda, d, rSuppRange, thetaSuppRange);
[etaElo, etaEhi] = eta_interval(k_lambda, d, rEvalRange, thetaEvalRange);

box = struct();
box.muSupp = [muSlo muShi];
box.muEval = [muElo muEhi];
box.etaSupp = [etaSlo etaShi];
box.etaEval = [etaElo etaEhi];
box.omega1 = [muSlo-muEhi, muShi-muElo];
box.omega2 = [etaElo-etaShi, etaEhi-etaSlo];
end

function tauInt = tau_range_from_geometry(d, rRange, thetaRange)
rRange = sort(rRange(:).');
thetaRange = sort(thetaRange(:).');
[cLo,cHi] = cos_interval(thetaRange);
invR = [1/rRange(2), 1/rRange(1)];
vals = d*[cLo*invR(1), cLo*invR(2), cHi*invR(1), cHi*invR(2)];
tauInt = [min(vals), max(vals)];
end

function [lo, hi] = mu_interval(k_lambda, d, thetaRange)
thetaRange = sort(thetaRange(:).');
[cLo,cHi] = cos_interval(thetaRange);
vals = k_lambda*d*[cLo cHi];
lo = min(vals);
hi = max(vals);
end

function [lo, hi] = eta_interval(k_lambda, d, rRange, thetaRange)
rRange = sort(rRange(:).');
thetaRange = sort(thetaRange(:).');
[s2lo, s2hi] = sin2_interval(thetaRange);
factor = k_lambda*d^2/2;
lo = factor*s2lo/rRange(2);
hi = factor*s2hi/rRange(1);
end

function [lo, hi] = cos_interval(thetaRange)
a = thetaRange(1); b = thetaRange(2);
vals = [cos(a), cos(b)];
ks = floor(a/(2*pi))-1:ceil(b/(2*pi))+1;
for k = ks
    t = 2*pi*k;
    if t >= a && t <= b
        vals(end+1) = 1; %#ok<AGROW>
    end
    t = pi + 2*pi*k;
    if t >= a && t <= b
        vals(end+1) = -1; %#ok<AGROW>
    end
end
lo = min(vals);
hi = max(vals);
end

function [lo, hi] = sin2_interval(thetaRange)
a = thetaRange(1); b = thetaRange(2);
vals = [sin(a)^2, sin(b)^2];
ks = floor((a-pi/2)/pi)-1 : ceil((b-pi/2)/pi)+1;
for k = ks
    t = pi/2 + k*pi;
    if t >= a && t <= b
        vals(end+1) = 1; %#ok<AGROW>
    end
end
ks = floor(a/pi)-1 : ceil(b/pi)+1;
for k = ks
    t = k*pi;
    if t >= a && t <= b
        vals(end+1) = 0; %#ok<AGROW>
    end
end
lo = min(vals);
hi = max(vals);
end

function [lo, hi] = alpha_cos2_interval(alphaInt, thetaInt)
[c2lo,c2hi] = cos2_interval(thetaInt);
vals = [alphaInt(1)*c2lo, alphaInt(1)*c2hi, alphaInt(2)*c2lo, alphaInt(2)*c2hi];
lo = min(vals);
hi = max(vals);
end

function [lo, hi] = cos2_interval(thetaRange)
a = thetaRange(1); b = thetaRange(2);
vals = [cos(2*a), cos(2*b)];
ks = floor(a/pi)-2:ceil(b/pi)+2;
for k = ks
    t = k*pi;
    if t >= a && t <= b
        vals(end+1) = 1; %#ok<AGROW>
    end
    t = pi/2 + k*pi;
    if t >= a && t <= b
        vals(end+1) = -1; %#ok<AGROW>
    end
end
lo = min(vals);
hi = max(vals);
end

function M = max_abs_affine_1d(a, b, I)
I = sort(I(:).');
vals = [a + b*I(1), a + b*I(2)];
M = max(abs(vals));
M = max(real(M),0);
end

function M = max_abs_bilinear_box(xcoef, ycoef, cInt, zInt)
cInt = sort(cInt(:).');
zInt = sort(zInt(:).');
vals = zeros(4,1);
k = 0;
for i = 1:2
    for j = 1:2
        k = k + 1;
        vals(k) = abs(cInt(i)*xcoef + zInt(j)*ycoef);
    end
end
M = max(vals);
M = max(real(M),0);
end

function classOut = make_slice_class_ranges(classIn, whichSlice, xGrid, rEvalFixed, thetaEvalFixed)
classOut = classIn;
switch lower(char(whichSlice))
    case 'theta'
        classOut.rEvalRange = [rEvalFixed rEvalFixed];
        classOut.thetaEvalRange = [min(xGrid) max(xGrid)];
    case 'range'
        classOut.rEvalRange = [min(xGrid) max(xGrid)];
        classOut.thetaEvalRange = [thetaEvalFixed thetaEvalFixed];
    otherwise
        error('Unknown slice type.');
end
end

function classOut = make_sample_class_ranges(classIn, whichSlice, xGrid, ii, rEvalFixed, thetaEvalFixed)
classOut = classIn;
[xLo, xHi] = grid_cell_interval(xGrid, ii, classIn.sampleCellPad);

switch lower(char(whichSlice))
    case 'theta'
        classOut.rEvalRange = [rEvalFixed rEvalFixed];
        classOut.thetaEvalRange = [xLo xHi];
    case 'range'
        classOut.rEvalRange = [xLo xHi];
        classOut.thetaEvalRange = [thetaEvalFixed thetaEvalFixed];
    otherwise
        error('Unknown slice type.');
end
end

function [lo, hi] = grid_cell_interval(xGrid, ii, pad)
xGrid = xGrid(:).';
M = numel(xGrid);

if M == 1
    dx = max(abs(xGrid(ii)),1)*1e-8;
    lo = xGrid(ii) - dx;
    hi = xGrid(ii) + dx;
else
    if ii == 1
        lo = xGrid(1);
    else
        lo = 0.5*(xGrid(ii-1) + xGrid(ii));
    end

    if ii == M
        hi = xGrid(M);
    else
        hi = 0.5*(xGrid(ii) + xGrid(ii+1));
    end
end

lo = lo - pad;
hi = hi + pad;
end

function dN = exact_dN_plus(N, omega1, omega2)
vals = zeros(N,1);
for n = 0:(N-1)
    vals(n+1) = dist_to_lattice(omega1 + (2*n+1)*omega2, 2*pi);
end
dN = min(vals);
end

function dlow = dN_plus_lower_bound_box(N, w1Int, w2Int)
dlow = inf;
for n = 0:(N-1)
    c = 2*n + 1;
    vals = [w1Int(1) + c*w2Int(1), ...
            w1Int(1) + c*w2Int(2), ...
            w1Int(2) + c*w2Int(1), ...
            w1Int(2) + c*w2Int(2)];
    lo = min(vals);
    hi = max(vals);
    dlow = min(dlow, interval_dist_to_lattice_min(lo, hi, 2*pi));
    if dlow <= 0
        dlow = 0;
        return;
    end
end
end

function m = interval_sin_abs_max(lo, hi)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if hi - lo >= pi
    m = 1;
    return;
end

vals = [abs(sin(lo)), abs(sin(hi))];
kmin = floor((lo - pi/2)/pi)-1;
kmax = ceil((hi + pi/2)/pi)+1;
for k = kmin:kmax
    t = pi/2 + k*pi;
    if t >= lo && t <= hi
        vals(end+1) = 1; %#ok<AGROW>
    end
end
m = max(vals);
end

function dmin = interval_dist_to_lattice_min(lo, hi, period)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if interval_hits_lattice(lo, hi, period)
    dmin = 0;
    return;
end
k1 = round(lo/period);
k2 = round(hi/period);
dmin = min(abs(lo-k1*period), abs(hi-k2*period));
dmin = min(dmin, period/2);
end

function dmax = interval_dist_to_lattice_max(lo, hi, period)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if hi - lo >= period
    dmax = period/2;
    return;
end

vals = [dist_to_lattice(lo, period), dist_to_lattice(hi, period)];
kmin = floor((lo-period/2)/period)-1;
kmax = ceil((hi+period/2)/period)+1;
for k = kmin:kmax
    t = (k + 0.5)*period;
    if t >= lo && t <= hi
        vals(end+1) = period/2; %#ok<AGROW>
    end
end
dmax = max(vals);
end

function tf = interval_hits_lattice(lo, hi, period)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
kmin = ceil(lo/period);
kmax = floor(hi/period);
tf = kmin <= kmax;
end

function d = dist_to_lattice(x, period)
d = abs(mod(x + period/2, period) - period/2);
end

function d = signed_dist_to_lattice(x, period)
d = mod(x + period/2, period) - period/2;
end

%% =====================================================================
% Higher derivative exact pointwise profiles
%% =====================================================================
function [Q2, Q3] = higher_Q_profiles(geom, kappa, alpha, theta)
x = geom.x(:);
y = geom.y(:);
A = x + alpha*cos(theta)*y;
B = cos(theta)*x + alpha*cos(2*theta)*y;
u = -kappa*sin(theta).*A;
ut = -kappa.*B;
utt = kappa*sin(theta).*(x + 4*alpha*cos(theta).*y);
Q2 = -1i*ut - u.^2;
Q3 = -1i*utt - 3*u.*ut + 1i*(u.^3);
end

%% =====================================================================
% Derivative, jump, geometry helpers
%% =====================================================================
function beta = derivative_class_beta_sine(d0, rhoBar)
s0 = sin(d0/2);
if s0 <= 0 || ~isfinite(s0)
    error('Invalid derivative d0.');
end
if rhoBar < 0 || rhoBar > 1 || ~isfinite(rhoBar)
    error('Invalid derivative rhoBar.');
end

beta = zeros(1,5);
beta(1) = 105*rhoBar^4/(16*s0^8);
beta(2) = 105*rhoBar^3/(16*s0^7);
beta(3) = 45*rhoBar^2/(16*s0^6);
beta(4) = 5*rhoBar/(8*s0^5);
beta(5) = 1/(16*s0^4);
end

function beta = derivative_beta_from_s0_rho(s0, rho)
beta = zeros(1,5);
beta(1) = 105*rho^4/(16*s0^8);
beta(2) = 105*rho^3/(16*s0^7);
beta(3) = 45*rho^2/(16*s0^6);
beta(4) = 5*rho/(8*s0^5);
beta(5) = 1/(16*s0^4);
end

function v = diff_norms_0_to_4_scalar(a)
a = a(:);
v = zeros(1,5);
for j = 0:4
    da = diff_order_vec(a, j);
    v(j+1) = sum(abs(da));
end
end

function y = diff_order_vec(x, j)
if j == 0
    y = x;
else
    y = diff(x, j, 1);
end
end

function Y = diff_order_rows(X, j)
if j == 0
    Y = X;
else
    Y = diff(X, j, 1);
end
end

function J = scalar_phase_jump_envelope(a, b, T)
ra = abs(a);
rb = abs(b);

if ra == 0
    J = rb;
    return;
end
if rb == 0
    J = ra;
    return;
end

center = angle(b) - angle(a);
d0 = dist_to_lattice(center, 2*pi);
dmax = min(pi, d0 + abs(T));
J = sqrt(max(ra^2 + rb^2 - 2*ra*rb*cos(dmax), 0));
end

function [x, y, G, E] = build_tangent_geometry_from_b(b)
b = b(:);
n = (0:(numel(b)-1)).';
nbar = sum(b.*n);
n2bar = sum(b.*n.^2);
x = n - nbar;
y = n.^2 - n2bar;

Gxx = sum(b.*x.^2);
Gxy = sum(b.*x.*y);
Gyy = sum(b.*y.^2);
G = [Gxx Gxy; Gxy Gyy];

if min(eig(G)) <= 0
    error('The tangent Gram matrix is not positive definite.');
end

Gsqrt = sqrtm(G);
E = real([x y]/Gsqrt);
end

function geom = build_tau_geometry_cache(b, x, y, G, E, pMax, Qall, coeffK, tauOpts)
geom = struct();
geom.N = numel(b);
geom.b = b(:);
geom.x = x(:);
geom.y = y(:);
geom.G = G;
geom.E = E;
geom.pMax = pMax;
geom.Qall = Qall;
geom.coeffK = coeffK;
geom.tauOpts = tauOpts;
geom.Pdh = dh_bilin_coeffs(geom.b, geom.x, geom.y);
end

function P = dh_bilin_coeffs(b, x, y)
P = zeros(numel(b),4);
P(:,1) = b(:).*x(:).^2;
P(:,2) = b(:).*y(:).*x(:);
P(:,3) = b(:).*x(:).*y(:);
P(:,4) = b(:).*y(:).^2;
end

%% =====================================================================
% Small 1D polynomial utilities
%% =====================================================================
function r = poly_der_asc(p)
p = trim_poly_asc(p(:).');
if numel(p) <= 1
    r = 0;
else
    r = (1:(numel(p)-1)).*p(2:end);
end
end

function r = poly_conv_asc(a,b)
r = conv(trim_poly_asc(a), trim_poly_asc(b));
r = trim_poly_asc(r);
end

function r = poly_sub_asc(a,b)
a = trim_poly_asc(a);
b = trim_poly_asc(b);
L = max(numel(a), numel(b));
aa = zeros(1,L);
bb = zeros(1,L);
aa(1:numel(a)) = a;
bb(1:numel(b)) = b;
r = trim_poly_asc(aa-bb);
end

function v = poly_eval_asc(p,x)
p = trim_poly_asc(p);
v = zeros(size(x));
for k = numel(p):-1:1
    v = v.*x + p(k);
end
end

function r = roots_descending_asc(pAsc)
pAsc = trim_poly_asc(pAsc);
if numel(pAsc) <= 1
    r = [];
else
    r = roots(fliplr(pAsc));
end
end

function p = trim_poly_asc(p)
p = p(:).';
while numel(p) > 1 && abs(p(end)) < 1e-14
    p(end) = [];
end
if isempty(p)
    p = 0;
end
end

%% =====================================================================
% Branch and info helpers
%% =====================================================================
function [B, branch] = update_best_branch(B, branch, C, label, guard, channels)
for k = 1:numel(channels)
    f = channels{k};
    if ~isfield(C, f) || ~isfinite(C.(f))
        continue;
    end
    val = C.(f) + guard;
    if isfinite(val) && val < B.(f)
        B.(f) = val;
        branch.(f) = string(label);
    end
end
end

function B = inf_channel_branch(channels)
B = struct();
for k = 1:numel(channels)
    B.(channels{k}) = inf;
end
end

function B = zero_channel_branch(channels)
B = struct();
for k = 1:numel(channels)
    B.(channels{k}) = 0;
end
end

function s = branch_names(name, channels)
s = struct();
for k = 1:numel(channels)
    s.(channels{k}) = string(name);
end
end

function info = make_bound_info(value, label, cap)
if nargin < 3 || isempty(cap)
    cap = inf;
end
info = struct();
info.raw = max(real(value), 0);
info.value = min(cap, info.raw);
info.bestBranch = string(label);
info.branchValues = struct();
info.rawBranchValues = struct();
info.l1 = NaN;
end

function info = disabled_bound_info(label, cap)
if nargin < 2
    cap = inf;
end
info = make_bound_info(inf, label, cap);
end

function S = disabled_bound_struct(label, channels)
S = struct();
for k = 1:numel(channels)
    S.(channels{k}) = make_bound_info(inf, label, inf);
end
end

function [val, bestBranch] = min_finite_candidates(candidates, labels)
finiteMask = isfinite(candidates);
if ~any(finiteMask)
    val = inf;
    bestBranch = 'none';
else
    cf = candidates(finiteMask);
    lf = labels(finiteMask);
    [val, pos] = min(cf);
    bestBranch = lf{pos};
end
val = max(real(val), 0);
end

%% =====================================================================
% Weights, diagnostics, plotting, validation, options
%% =====================================================================
function w = build_flat_end_weight_r4(N, mShape, family)
n = (0:N-1).';
x = n/max(N-1,1);
leftBin = safe_binom_vector(n,4);
rightBin = safe_binom_vector((N-1)-n,4);

switch lower(char(family))
    case 'binomial'
        w = leftBin.*rightBin;
    case 'binomial_beta'
        w = leftBin.*rightBin.*(x.^mShape).*((1-x).^mShape);
    otherwise
        error('Unknown weightFamily. Use binomial or binomial_beta.');
end

w(1:4) = 0;
w(end-3:end) = 0;

if all(w == 0)
    error('Constructed weight is identically zero. Need N >= 9.');
end
end

function c = safe_binom_vector(nVec, r)
c = zeros(size(nVec));
mask = nVec >= r;
nv = nVec(mask);
c(mask) = exp(gammaln(nv+1) - gammaln(r+1) - gammaln(nv-r+1));
c(abs(c) < 1e-14) = 0;
end

function R = init_result_arrays(M, channels)
R = struct();
for k = 1:numel(channels)
    ch = channels{k};
    R.([ch 'true']) = zeros(1,M);
    R.([ch 'actual']) = zeros(1,M);
    R.([ch 'support']) = zeros(1,M);
    R.([ch 'actualBranch']) = strings(1,M);
    R.([ch 'supportBranch']) = strings(1,M);
end
R.omega1 = zeros(1,M);
R.omega2 = zeros(1,M);
R.dNplus = zeros(1,M);
R.tauEval = zeros(1,M);
R.tauSupp = zeros(1,M);
R.aux = cell(1,M);
end

function print_new_bound_diag(name, truth, actual, support, floorVal)
fprintf('    %-4s min(actual-true)    : %.3e\n', name, min(actual-truth));
fprintf('    %-4s min(support-true)   : %.3e\n', name, min(support-truth));
fprintf('    %-4s median actual/true  : %.3e\n', name, safe_median_ratio(actual, truth, floorVal));
fprintf('    %-4s median support/true : %.3e\n', name, safe_median_ratio(support, truth, floorVal));
fprintf('    %-4s support below true  : %d\n', name, nnz(support-truth < -1e-10));
end

function r = safe_median_ratio(vals, truth, floorVal)
mask = abs(truth) > floorVal & isfinite(vals) & isfinite(truth);
if any(mask)
    r = median(vals(mask)./abs(truth(mask)));
else
    r = NaN;
end
end

function figs = plot_new_bounds_slice(data, whichSlice, xLabel, fixedText, fontSize, channels)
% comparison.
% Each figure shows exactly these three curves:
%   1) true pointwise kernel magnitude;
%   2)  Bbest bound;
%   3)  support-uniform bound.
figs = struct();
labels = channel_plot_labels();
for k = 1:numel(channels)
    ch = channels{k};
    figs.(ch) = plot_one_new_bounds_figure(data, whichSlice, xLabel, fixedText, fontSize, ch, labels.(ch));
end
end

function labels = channel_plot_labels()
labels = struct();
labels.K = '|K|';
labels.H = '|H|';
labels.dK = '|dK|';
labels.dH = '|dH|';
labels.d2K = '|d^2K|';
labels.d3K = '|d^3K|';
labels.d2H = '|d^2H|';
labels.d3H = '|d^3H|';
end

function fig = plot_one_new_bounds_figure(data, whichSlice, xLabel, fixedText, fontSize, channelName, yLabelText)
fig = figure('Color','w','Name',sprintf('%s true kernel and bounds, %s', channelName, whichSlice));
hold on;
leg = {};

for jf = 1:numel(data)
    x = data{jf}.xGrid;
    fcGHz = data{jf}.fc/1e9;

    yTrue = data{jf}.([channelName 'true']);
    yPoint = data{jf}.([channelName 'actual']);
    ySupport = data{jf}.([channelName 'support']);

    [yTruePlot, cap] = finite_plot_vector(yTrue, []);
    [yPointPlot, cap] = finite_plot_vector(yPoint, cap);
    [ySupportPlot, ~] = finite_plot_vector(ySupport, cap);

    if numel(data) == 1
        plot(x, yTruePlot,   'kx--', 'LineWidth', 3, 'MarkerSize', 8);
        plot(x, yPointPlot,  'b-', 'LineWidth', 3);
        plot(x, ySupportPlot,'r-',  'LineWidth', 3);
        leg = {'True kernel', 'Pointwise bound', 'Support-uniform bound'};
    else
        plot(x, yTruePlot,   'x--',  'LineWidth', 3, 'MarkerSize', 8);
        plot(x, yPointPlot,  '-', 'LineWidth', 3);
        plot(x, ySupportPlot,'-',  'LineWidth', 3);
        leg{end+1} = sprintf('True %.3g GHz', fcGHz); %#ok<AGROW>
        leg{end+1} = sprintf('Pointwise bound %.3g GHz', fcGHz); %#ok<AGROW>
        leg{end+1} = sprintf('Support-uniform bound %.3g GHz', fcGHz); %#ok<AGROW>
    end
end

grid on;
xlabel(xLabel);
ylabel(yLabelText);
title(sprintf('%s: (%s)', channelName, fixedText), 'FontWeight','bold');
legend(leg, 'Location','best');
allX = [];
for jf = 1:numel(data)
    allX = [allX, data{jf}.xGrid(:).']; %#ok<AGROW>
end
if ~isempty(allX)
    xmin = min(allX);
    xmax = max(allX);
    if xmax > xmin
        xlim([xmin xmax]);
    else
        pad = max(1e-12, 1e-3*max(abs(xmin),1));
        xlim([xmin-pad xmax+pad]);
    end
end
set(gca, 'FontSize', fontSize);
hold off;
end

function [yp, cap] = finite_plot_vector(y, cap)
% Preserve the stored data exactly, but make nonfinite curves visible in plots.
yp = y;
finiteVals = y(isfinite(y));
if isempty(cap)
    if isempty(finiteVals)
        cap = 1;
    else
        cap = max(abs(finiteVals));
        if cap <= 0 || ~isfinite(cap)
            cap = 1;
        end
        cap = 10*cap;
    end
else
    if ~isempty(finiteVals)
        cap = max(cap, 10*max(abs(finiteVals)));
    end
end
yp(~isfinite(yp)) = cap;
end

function assert_bound_valid(R, channels, kind, tol)
for k = 1:numel(channels)
    ch = channels{k};
    truth = R.([ch 'true']);
    vals = R.([ch kind]);
    bad = vals + tol < truth;
    if any(bad)
        idx = find(bad, 1, 'first');
        error('%s %s bound below true kernel at index %d: bound=%.16e, true=%.16e', ...
            ch, kind, idx, vals(idx), truth(idx));
    end
end
end

function classOpts = normalize_class_options(classOpts)
if ~isfield(classOpts, 'rSuppRange') || isempty(classOpts.rSuppRange)
    classOpts.rSuppRange = [NaN NaN];
end
if ~isfield(classOpts, 'thetaSuppRange') || isempty(classOpts.thetaSuppRange)
    classOpts.thetaSuppRange = [NaN NaN];
end
if ~isfield(classOpts, 'rEvalRange') || isempty(classOpts.rEvalRange)
    classOpts.rEvalRange = [NaN NaN];
end
if ~isfield(classOpts, 'thetaEvalRange') || isempty(classOpts.thetaEvalRange)
    classOpts.thetaEvalRange = [NaN NaN];
end
if ~isfield(classOpts, 'higher') || isempty(classOpts.higher)
    classOpts.higher = struct();
end
if ~isfield(classOpts, 'subcells') || isempty(classOpts.subcells)
    classOpts.subcells = struct('theta',1,'range',1);
end

classOpts.derivative.d0List = unique(classOpts.derivative.d0List(:).', 'stable');
classOpts.derivative.rhoBarList = unique(classOpts.derivative.rhoBarList(:).', 'stable');
classOpts.derivative.curveDNCells = max(1, ceil(get_struct_opt(classOpts.derivative, 'curveDNCells', 16)));
classOpts.corr.QList = unique(classOpts.corr.QList(:).', 'stable');
classOpts.corr.pList = unique(classOpts.corr.pList(:).', 'stable');
classOpts.lin.qList = unique(classOpts.lin.qList(:).', 'stable');
classOpts.subcells.theta = max(1, ceil(get_struct_opt(classOpts.subcells, 'theta', 1)));
classOpts.subcells.range = max(1, ceil(get_struct_opt(classOpts.subcells, 'range', 1)));
end

function validate_inputs(Nr, useHalfLambda, d_fixed, r_min, r_max, theta_min, theta_max, ...
    r_support, theta_support, r_eval_fixed, theta_eval_fixed, theta_eval_grid, ...
    r_eval_grid, sliceType, actualOpts, classOpts, tauOpts)

if Nr < 9
    error('Need Nr >= 9 for fourth-order flat-end taper.');
end
if ~useHalfLambda && isempty(d_fixed)
    error('If useHalfLambda=false, provide opts.d_fixed.');
end
if ~(r_min > 0 && r_min < r_max)
    error('Need 0 < r_min < r_max.');
end
if ~(0 < theta_min && theta_min < theta_max && theta_max < pi)
    error('Need 0 < theta_min < theta_max < pi.');
end
if ~(r_support > 0) || ~(r_eval_fixed > 0)
    error('Support/eval ranges must have positive r.');
end
if ~(0 < theta_support && theta_support < pi) || ~(0 < theta_eval_fixed && theta_eval_fixed < pi)
    error('theta_support and theta_eval_fixed must lie in (0,pi).');
end
if any(theta_eval_grid <= 0) || any(theta_eval_grid >= pi)
    error('theta_eval_grid must lie in (0,pi).');
end
if any(r_eval_grid <= 0)
    error('r_eval_grid must be positive.');
end
if ~ismember(sliceType, {'theta','range','both'})
    error('sliceType must be theta, range, or both.');
end
if isempty(actualOpts.pCorrOrderList) || any(actualOpts.pCorrOrderList < 1) || any(actualOpts.pCorrOrderList ~= floor(actualOpts.pCorrOrderList))
    error('actual.pCorrOrderList must contain positive integers.');
end
if actualOpts.residueQMax < 1 || actualOpts.residueQMax ~= floor(actualOpts.residueQMax)
    error('actual.residueQMax must be a positive integer.');
end
if ~ismember(lower(char(classOpts.uniformMode)), {'samplecell','global'})
    error('class.uniformMode must be samplecell or global.');
end
if tauOpts.maxDepth2D < 0 || tauOpts.maxCells2D < 1 || tauOpts.bilinMaxDepth2D < 0 || tauOpts.bilinMaxCells2D < 1
    error('tau depth/cell parameters must be nonnegative/positive.');
end
end

function val = get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name)
    val = s.(name);
else
    val = defaultVal;
end
end

function val = get_struct_opt(s, name, defaultVal)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end

function s = merge_structs(a,b)
s = a;
if isempty(b)
    return;
end
if ~isstruct(b)
    s = b;
    return;
end
fields = fieldnames(b);
for k = 1:numel(fields)
    f = fields{k};
    if isfield(s,f) && isstruct(s.(f)) && isstruct(b.(f))
        s.(f) = merge_structs(s.(f), b.(f));
    else
        s.(f) = b.(f);
    end
end
end

function c = binomial_diff_coeffs(j)
c = zeros(j+1,1);
for ell = 0:j
    c(ell+1) = (-1)^(j-ell) * nchoosek(j,ell);
end
end
