% Self-contained reproducibility code for the LCS QPAC example.
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
% This single MATLAB file contains the complete LCS example implementation.
% No external MATLAB files or input data files are required.
%
% Run in MATLAB:
%     result = run_lcs_example_qpac();

function result = run_lcs_example_qpac()
%RUN_LCS_EXAMPLE_QPAC Verify QPAC using lag-correlation support bounds.
%
%   result = run_lcs_example_qpac();
%
% The support consists of two sources at the same angle and on distinct
% range rows. Lag correlation bounds the four support-interaction channels.
% Analytic cell bounds and Lipschitz padding control the near and far row
% sums over the continuous angular domain. The output records each candidate
% bound, the selected branch, and the resulting QPAC recovery number.
%
% Calculations use MATLAB double precision. MATLAB R2016b or later is
% required for implicit array expansion.
%

%% 1. Model and numerical settings
par = default_config();
validate_config(par);

[geom,rho] = build_geometry(par);

par.aperture = (par.Nr-1)*par.d;
par.reactiveBoundary = 0.62*sqrt(par.aperture^3/par.lambda);
par.fraunhoferBoundary = 2*par.aperture^2/par.lambda;

% Keep fixed relative margins from both conventional Fresnel boundaries.
par.rangeRows = [1.05*par.reactiveBoundary, ...
                 0.92*par.fraunhoferBoundary];
par.fresnelRegionOK = all(par.rangeRows > par.reactiveBoundary) && ...
    all(par.rangeRows < par.fraunhoferBoundary);
assert(par.fresnelRegionOK, ...
    'QPAC:RangeOutsideFresnel', ...
    'Every range row must lie strictly inside the radiating Fresnel region.');

support = struct();
support.range = par.rangeRows;
support.theta = [par.theta0,par.theta0];

flatEndOK = all(abs(geom.b(1:4)) == 0) && ...
    all(abs(geom.b(end-3:end)) == 0);
assert(flatEndOK,'QPAC:TaperNotFlat', ...
    'The first and last four taper entries must vanish.');

%% 2. Support bounds from lag correlation
% Channel order:
%   K,H,dK,dH,d2K,d3K,d2H,d3H.
% dK,dH are normalized tangent channels; d2/d3 are physical angular
% derivatives on the evaluation side.
channelNames = {'K','H','dK','dH','d2K','d3K','d2H','d3H'};
orientations = [1 2;2 1];
nOrient = size(orientations,1);
nChannel = numel(channelNames);

lcsBound = zeros(nOrient,nChannel);
lcsBestQ = zeros(nOrient,nChannel);
lcsRelativeLagResidual = zeros(nOrient,nChannel);
exactSupport = zeros(nOrient,nChannel);
phaseCoordinates = zeros(nOrient,2);

for io = 1:nOrient
    ie = orientations(io,1);
    is = orientations(io,2);
    rEval = support.range(ie);
    rSource = support.range(is);
    thetaEval = support.theta(ie);
    thetaSource = support.theta(is);

    [omega1,omega2] = phase_coordinates_pair( ...
        rSource,thetaSource,rEval,thetaEval,par);
    phaseCoordinates(io,:) = [omega1,omega2];
    coeff = fixed_channel_coefficients( ...
        rSource,thetaSource,rEval,thetaEval,geom,par);
    exactSupport(io,:) = evaluate_source( ...
        rSource,thetaSource,rEval,thetaEval,geom,par).';

    for ic = 1:nChannel
        lcs = lcs_fixed_point_bound( ...
            coeff{ic},omega1,omega2,par.QList,par.roundoffFactor);
        lcsBound(io,ic) = lcs.bound;
        lcsBestQ(io,ic) = lcs.bestQ;
        lcsRelativeLagResidual(io,ic) = ...
            lcs.maxRelativeLagExpansionResidual;
    end
end

dominanceTolerance = par.roundoffFactor*eps(max(1,abs(lcsBound)));
lcsDominatesExact = all( ...
    lcsBound(:)+dominanceTolerance(:) >= exactSupport(:));
assert(lcsDominatesExact,'QPAC:LCSUnderestimate', ...
    'A printed LCS value is below the corresponding exact channel value.');
supportBestQ = lcsBestQ(:,1:4);
supportUsesNontrivialLCS = all(isfinite(supportBestQ(:))) && ...
    all(supportBestQ(:) >= 2);
assert(supportUsesNontrivialLCS,'QPAC:TrivialSupportCap', ...
    'A support Hermite channel used only the trivial triangle cap.');
auditedLagResidual = lcsRelativeLagResidual( ...
    isfinite(lcsRelativeLagResidual));
lagExpansionAuditOK = ~isempty(auditedLagResidual) && ...
    all(auditedLagResidual <= par.lagResidualTolerance);
assert(lagExpansionAuditOK,'QPAC:LagExpansionAudit', ...
    'The direct and correlation-expanded residue energies disagree.');

uSS = max(lcsBound(:,1:4),[],1);
GSS = [uSS(1),uSS(2);uSS(3),uSS(4)];
MSS = (par.L-1)*GSS;
etaSS = max(abs(eig(MSS)));
etaSSOK = isfinite(etaSS) && etaSS < 1;
assert(etaSSOK,'QPAC:HermiteBlockNotContractive', ...
    'The LCS support interpolation matrix is not contractive.');

Gamma = (eye(2)-MSS)\[1;0];
Xi = Gamma-[1;0];
inverseResidual = norm((eye(2)-MSS)*Gamma-[1;0],inf);
coefficientSignsOK = all(Gamma >= -par.auditTolerance) && ...
    all(Xi >= -par.auditTolerance);
assert(coefficientSignsOK,'QPAC:CoefficientEnvelopeSign', ...
    'Unexpected sign in Gamma or Xi.');
Gamma = max(Gamma,0);
Xi = max(Xi,0);

% Diagnostic for the configured derivative threshold. It does not determine
% whether other analytical branches are available.
dNplus = zeros(nOrient,1);
for io = 1:nOrient
    increments = phaseCoordinates(io,1) + ...
        (2*geom.n+1)*phaseCoordinates(io,2);
    dNplus(io) = min(distance_to_2pi(increments));
end
derivativeThresholdPassed = min(dNplus) >= par.derivativeThreshold;

%% 3. Self checks and curvature margin
selfIdentityError = zeros(par.L,4);
for ell = 1:par.L
    Cself = evaluate_source(support.range(ell),support.theta(ell), ...
        support.range(ell),support.theta(ell),geom,par);
    selfIdentityError(ell,:) = abs(Cself(1:4).'-[1,0,0,1]);
end
selfIdentitiesOK = max(selfIdentityError(:)) <= par.selfIdentityTolerance;
assert(selfIdentitiesOK,'QPAC:SelfIdentityFailure', ...
    'The normalized self-channel identities failed.');

[tauLo,tauHi] = tau_interval( ...
    par.rangeRows,par.thetaDomain,par.d);
qMin = quadratic_minimum(geom.G,tauLo,tauHi);
qMax = quadratic_maximum(geom.G,tauLo,tauHi);
sMin = min(sin(par.thetaDomain));
assert(sMin > 0,'QPAC:EndfireDomain', ...
    'The angular domain must stay inside (0,pi).');
sigmaLower2 = geom.kappa^2*sMin^2*qMin;

tauAbs = max(abs([tauLo,tauHi]));
hUniform = (abs(geom.x)+tauAbs*abs(geom.y))/sqrt(qMin);
UK2Self = 0;
UH2Self = 0;
profileCaps = cell(par.L,1);
for ir = 1:par.L
    profileCaps{ir} = global_profile_caps( ...
        geom,par.d/par.rangeRows(ir));
    UK2Self = max(UK2Self, ...
        sum(geom.b.*profileCaps{ir}.Q2));
    UH2Self = max(UH2Self, ...
        sum(geom.b.*hUniform.*profileCaps{ir}.Q2));
end

% The two support cells are singletons, so the complex-coefficient LCS
% calculation above also gives point-cell bounds for d2K and d2H.
UK2SS = max(lcsBound(:,5));
UH2SS = max(lcsBound(:,7));

Ecurv = 2*Xi(1)*UK2Self + 2*Xi(2)*UH2Self + ...
    2*(par.L-1)*(Gamma(1)*UK2SS+Gamma(2)*UH2SS);
mNear = 2*sigmaLower2-Ecurv;
mNearOK = isfinite(mNear) && mNear > 0;
assert(mNearOK,'QPAC:NonpositiveCurvatureMargin', ...
    'The curvature margin is not positive.');

%% 4. Near and far row-sum bounds
% Each range row contains one support at theta0. The Taylor interval on that
% row is [theta0-delta,theta0+delta].
nearGrid = linspace(par.theta0-par.delta, ...
                    par.theta0+par.delta,par.nNear);
farLeftGrid = linspace(par.thetaDomain(1), ...
                       par.theta0-par.delta,par.nFar);
farRightGrid = linspace(par.theta0+par.delta, ...
                        par.thetaDomain(2),par.nFar);

hNear = 2*par.delta/(par.nNear-1);
hFarLeft = (par.theta0-par.delta-par.thetaDomain(1))/(par.nFar-1);
hFarRight = (par.thetaDomain(2)-par.theta0-par.delta)/(par.nFar-1);

row = cell(par.L,1);
D2DirSumByRow = zeros(par.L,1);
D2AnalyticByRow = zeros(par.L,1);
D2SelectedByRow = zeros(par.L,1);
D2SelectedBranchByRow = strings(par.L,1);
D2AnalyticAuditByRow = false(par.L,1);
D3KernelByRow = zeros(par.L,1);
D3KernelAuditByRow = false(par.L,1);
D3DirectByRow = zeros(par.L,1);
D3SelectedByRow = zeros(par.L,1);
D3SelectedBranchByRow = strings(par.L,1);
D3SampleMaximumByRow = zeros(par.L,1);
farDirSumByRow = zeros(par.L,1);
farAnalyticByRow = zeros(par.L,1);
farSelectedByRow = zeros(par.L,1);
farSelectedBranchByRow = strings(par.L,1);
farAnalyticAuditByRow = false(par.L,1);
D2CoarseByRow = zeros(par.L,1);
D3CoarseByRow = zeros(par.L,1);
farCoarseByRow = zeros(par.L,1);

localSegment = [par.theta0-par.delta,par.theta0+par.delta];
localCellEdges = linspace(localSegment(1),localSegment(2), ...
    par.kernelBound.localCells+1);
coarseLocalCellEdges = localCellEdges(1:2:end);
if coarseLocalCellEdges(end) < localCellEdges(end)
    coarseLocalCellEdges(end+1) = localCellEdges(end); %#ok<AGROW>
end

farLeftCellEdges = linspace(par.thetaDomain(1), ...
    par.theta0-par.delta,par.kernelBound.farCellsPerSide+1);
farRightCellEdges = linspace(par.theta0+par.delta, ...
    par.thetaDomain(2),par.kernelBound.farCellsPerSide+1);
coarseFarLeftCellEdges = farLeftCellEdges(1:2:end);
coarseFarRightCellEdges = farRightCellEdges(1:2:end);
if coarseFarLeftCellEdges(end) < farLeftCellEdges(end)
    coarseFarLeftCellEdges(end+1) = farLeftCellEdges(end); %#ok<AGROW>
end
if coarseFarRightCellEdges(end) < farRightCellEdges(end)
    coarseFarRightCellEdges(end+1) = farRightCellEdges(end); %#ok<AGROW>
end

% Precompute normalized K/H coefficients and the nine-atom derivative
% dictionary. The Q2 construction uses its first five atoms.
analyticKernelCache = build_analytic_kernel_cache( ...
    support,geom,par);

for ir = 1:par.L
    rEval = par.rangeRows(ir);
    [~,F2Near,F3Near] = evaluate_common_rows( ...
        rEval,nearGrid,support,Gamma,geom,par);
    [F0Left,~,~] = evaluate_common_rows( ...
        rEval,farLeftGrid,support,Gamma,geom,par);
    [F0Right,~,~] = evaluate_common_rows( ...
        rEval,farRightGrid,support,Gamma,geom,par);

    [LF0,LF2] = row_lipschitz_constants( ...
        support,Gamma,geom,profileCaps{ir},par);

    % Each cell bound holds for every evaluation angle in the cell. Sum the
    % source contributions before maximizing over cells.
    D2KernelInfo = higher_derivative_kernel_row_bound( ...
        2,rEval,localCellEdges,support,Gamma, ...
        analyticKernelCache,geom,par);
    D3KernelInfo = higher_derivative_kernel_row_bound( ...
        3,rEval,localCellEdges,support,Gamma, ...
        analyticKernelCache,geom,par);

    farKernelInfo = normalized_kernel_far_row_bound( ...
        rEval,{farLeftCellEdges,farRightCellEdges}, ...
        support,Gamma,analyticKernelCache,geom,par);

    % Direct third-derivative bound on the full local interval.
    D3Direct = direct_third_derivative_row_bound( ...
        rEval,localSegment,support,Gamma,geom,par);
    D3Direct = D3Direct + ...
        floating_pad(D3Direct,par.roundoffFactor);

    D2Grid = max(F2Near);
    D3SampleMaximum = max(F3Near);
    farGridMaximum = max([F0Left,F0Right]);

    D2Pad = 0.5*hNear*LF2 + ...
        floating_pad(D2Grid,par.roundoffFactor);
    farPad = max(0.5*hFarLeft*LF0,0.5*hFarRight*LF0) + ...
        floating_pad(farGridMaximum,par.roundoffFactor);
    D2DirSum = D2Grid+D2Pad;
    farDirSum = farGridMaximum+farPad;

    [D2Selected,D2SelectedBranch] = select_enabled_bound( ...
        [D2KernelInfo.bound,D2DirSum], ...
        {'analytic-kernel','DirSum'}, ...
        [par.methods.enableAnalyticKernel,par.methods.enableDirSum]);

    [D3Selected,D3SelectedBranch] = select_enabled_bound( ...
        [D3KernelInfo.bound,D3Direct], ...
        {'analytic-kernel','direct-Q3'}, ...
        [par.methods.enableAnalyticKernel,par.methods.enableDirectQ3]);

    [farSelected,farSelectedBranch] = select_enabled_bound( ...
        [farKernelInfo.bound,farDirSum], ...
        {'analytic-kernel','DirSum'}, ...
        [par.methods.enableAnalyticKernel,par.methods.enableDirSum]);

    D2CoarseKernelInfo = higher_derivative_kernel_row_bound( ...
        2,rEval,coarseLocalCellEdges,support,Gamma, ...
        analyticKernelCache,geom,par);
    D3CoarseKernelInfo = higher_derivative_kernel_row_bound( ...
        3,rEval,coarseLocalCellEdges,support,Gamma, ...
        analyticKernelCache,geom,par);
    farCoarseKernelInfo = normalized_kernel_far_row_bound( ...
        rEval,{coarseFarLeftCellEdges,coarseFarRightCellEdges}, ...
        support,Gamma,analyticKernelCache,geom,par);

    D2DirSumByRow(ir) = D2DirSum;
    D2AnalyticByRow(ir) = D2KernelInfo.bound;
    D2SelectedByRow(ir) = D2Selected;
    D2SelectedBranchByRow(ir) = D2SelectedBranch;
    D2AnalyticAuditByRow(ir) = D2KernelInfo.pointAuditOK;
    D3KernelByRow(ir) = D3KernelInfo.bound;
    D3KernelAuditByRow(ir) = D3KernelInfo.pointAuditOK;
    D3DirectByRow(ir) = D3Direct;
    D3SelectedByRow(ir) = D3Selected;
    D3SelectedBranchByRow(ir) = D3SelectedBranch;
    D3SampleMaximumByRow(ir) = D3SampleMaximum;
    farDirSumByRow(ir) = farDirSum;
    farAnalyticByRow(ir) = farKernelInfo.bound;
    farSelectedByRow(ir) = farSelected;
    farSelectedBranchByRow(ir) = farSelectedBranch;
    farAnalyticAuditByRow(ir) = farKernelInfo.pointAuditOK;

    % Repeat the calculation on a nested coarser cover as a stability check.
    idxNear = 1:2:par.nNear;
    idxFar = 1:2:par.nFar;
    D2CoarseDirSum = max(F2Near(idxNear)) + hNear*LF2 + ...
        floating_pad(D2Grid,par.roundoffFactor);
    farCoarseDirSum = max([F0Left(idxFar),F0Right(idxFar)]) + ...
        max(hFarLeft*LF0,hFarRight*LF0) + ...
        floating_pad(farGridMaximum,par.roundoffFactor);
    D2CoarseByRow(ir) = select_enabled_bound( ...
        [D2CoarseKernelInfo.bound,D2CoarseDirSum], ...
        {'analytic-kernel','DirSum'}, ...
        [par.methods.enableAnalyticKernel,par.methods.enableDirSum]);
    D3CoarseByRow(ir) = select_enabled_bound( ...
        [D3CoarseKernelInfo.bound,D3Direct], ...
        {'analytic-kernel','direct-Q3'}, ...
        [par.methods.enableAnalyticKernel,par.methods.enableDirectQ3]);
    farCoarseByRow(ir) = select_enabled_bound( ...
        [farCoarseKernelInfo.bound,farCoarseDirSum], ...
        {'analytic-kernel','DirSum'}, ...
        [par.methods.enableAnalyticKernel,par.methods.enableDirSum]);

    row{ir} = struct( ...
        'range',rEval, ...
        'F2Near',F2Near,'F3Near',F3Near, ...
        'F0Left',F0Left,'F0Right',F0Right, ...
        'LF0',LF0,'LF2',LF2, ...
        'D2Pad',D2Pad, ...
        'D2DirSum',D2DirSum, ...
        'D2Analytic',D2KernelInfo.bound, ...
        'D2Selected',D2Selected, ...
        'D2SelectedBranch',D2SelectedBranch, ...
        'D2AnalyticCellBounds',D2KernelInfo.cellBounds, ...
        'D3Kernel',D3KernelInfo.bound, ...
        'D3KernelCellBounds',D3KernelInfo.cellBounds, ...
        'D3KernelPairK',D3KernelInfo.pairK, ...
        'D3KernelPairH',D3KernelInfo.pairH, ...
        'D3KernelPointAuditOK',D3KernelInfo.pointAuditOK, ...
        'D3KernelMaxPointRatioK',D3KernelInfo.maxPointRatioK, ...
        'D3KernelMaxPointRatioH',D3KernelInfo.maxPointRatioH, ...
        'D3Direct',D3Direct, ...
        'D3Selected',D3Selected, ...
        'D3SelectedBranch',D3SelectedBranch, ...
        'D3SampleMaximum',D3SampleMaximum, ...
        'farPad',farPad, ...
        'farDirSum',farDirSum, ...
        'farAnalytic',farKernelInfo.bound, ...
        'farSelected',farSelected, ...
        'farSelectedBranch',farSelectedBranch, ...
        'farAnalyticCellBounds',farKernelInfo.cellBounds);
end

selectedD2DominatesSamplesByRow = D2SelectedByRow + ...
    floating_pad(D2SelectedByRow,par.roundoffFactor) >= ...
    cellfun(@(z) max(z.F2Near),row);
selectedD2DominatesSamples = all(selectedD2DominatesSamplesByRow);
assert(selectedD2DominatesSamples,'QPAC:SelectedD2Underestimate', ...
    'The selected D2 envelope is below a sampled common-point row.');

selectedD3DominatesSamplesByRow = D3SelectedByRow + ...
    floating_pad(D3SelectedByRow,par.roundoffFactor) >= ...
    D3SampleMaximumByRow;
selectedD3DominatesSamples = all(selectedD3DominatesSamplesByRow);
assert(selectedD3DominatesSamples,'QPAC:SelectedD3Underestimate', ...
    ['The selected uniform third-derivative envelope is below a sampled ', ...
     'common-point third-derivative row.']);
assert(all(D3KernelAuditByRow),'QPAC:KernelD3PointAudit', ...
    ['A support-uniform d3K/d3H cell bound is below an endpoint or ', ...
     'midpoint audit value.']);
assert(all(D2AnalyticAuditByRow),'QPAC:KernelD2PointAudit', ...
    ['A support-uniform d2K/d2H cell bound is below an endpoint or ', ...
     'midpoint audit value.']);

selectedFarDominatesSamplesByRow = farSelectedByRow + ...
    floating_pad(farSelectedByRow,par.roundoffFactor) >= ...
    cellfun(@(z) max([z.F0Left,z.F0Right]),row);
selectedFarDominatesSamples = all(selectedFarDominatesSamplesByRow);
assert(selectedFarDominatesSamples,'QPAC:SelectedFarUnderestimate', ...
    'The selected far envelope is below a sampled common-point row.');
assert(all(farAnalyticAuditByRow),'QPAC:KernelFarPointAudit', ...
    ['A support-uniform K/H far-cell bound is below an endpoint or ', ...
     'midpoint audit value.']);

% Form each maximum after summing all source contributions at the same q.
D2 = max(D2SelectedByRow);
D3 = max(D3SelectedByRow);
etaFar = max(farSelectedByRow);

etaNear = 2*par.delta*D3/(3*mNear) + ...
    par.delta^2*D2^2/(2*mNear);
CQPAC = max([etaSS,etaNear,etaFar]);

D2Coarse = max(D2CoarseByRow);
D3Coarse = max(D3CoarseByRow);
etaFarCoarse = max(farCoarseByRow);
etaNearCoarse = 2*par.delta*D3Coarse/(3*mNear) + ...
    par.delta^2*D2Coarse^2/(2*mNear);

%% 5. Consistency checks
endpointScale = max(1,max(abs([par.thetaDomain,par.theta0,par.delta])));
endpointTolerance = 64*eps(endpointScale);
supportCellCoverageOK = numel(support.range) == par.L && ...
    numel(support.theta) == par.L; % singleton support cells
localSegmentCoverageOK = ...
    abs(nearGrid(1)-(par.theta0-par.delta)) <= endpointTolerance && ...
    abs(nearGrid(end)-(par.theta0+par.delta)) <= endpointTolerance;
kernelLocalCellCoverageOK = ...
    abs(localCellEdges(1)-localSegment(1)) <= endpointTolerance && ...
    abs(localCellEdges(end)-localSegment(2)) <= endpointTolerance && ...
    all(diff(localCellEdges) > 0);
kernelFarCellCoverageOK = ...
    abs(farLeftCellEdges(1)-par.thetaDomain(1)) <= endpointTolerance && ...
    abs(farLeftCellEdges(end)-(par.theta0-par.delta)) <= endpointTolerance && ...
    abs(farRightCellEdges(1)-(par.theta0+par.delta)) <= endpointTolerance && ...
    abs(farRightCellEdges(end)-par.thetaDomain(2)) <= endpointTolerance && ...
    all(diff(farLeftCellEdges) > 0) && ...
    all(diff(farRightCellEdges) > 0);
farCoverageOK = ...
    abs(farLeftGrid(1)-par.thetaDomain(1)) <= endpointTolerance && ...
    abs(farLeftGrid(end)-(par.theta0-par.delta)) <= endpointTolerance && ...
    abs(farRightGrid(1)-(par.theta0+par.delta)) <= endpointTolerance && ...
    abs(farRightGrid(end)-par.thetaDomain(2)) <= endpointTolerance;
% evaluate_common_rows uses the same q vector for every source term. The
% near and far intervals use the same value of par.delta.
commonPointNearByConstruction = true;
commonPointFarByConstruction = true;
sameRadiusByConstruction = true;
farBoundaryIncluded = farCoverageOK && kernelFarCellCoverageOK;

allFinite = all(isfinite([etaSS;Gamma;Xi;sigmaLower2;UK2Self; ...
    UH2Self;UK2SS;UH2SS;Ecurv;mNear;D2;D3;etaNear;etaFar;CQPAC]));
internalAuditOK = allFinite && flatEndOK && par.fresnelRegionOK && ...
    selfIdentitiesOK && lcsDominatesExact && ...
    supportUsesNontrivialLCS && lagExpansionAuditOK && ...
    selectedD2DominatesSamples && all(D2AnalyticAuditByRow) && ...
    selectedD3DominatesSamples && ...
    all(D3KernelAuditByRow) && ...
    selectedFarDominatesSamples && all(farAnalyticAuditByRow) && ...
    supportCellCoverageOK && ...
    localSegmentCoverageOK && kernelLocalCellCoverageOK && ...
    farCoverageOK && kernelFarCellCoverageOK && ...
    commonPointNearByConstruction && commonPointFarByConstruction && ...
    sameRadiusByConstruction;
qpacInequalitiesOK = etaSS < 1 && mNear > 0 && ...
    etaNear < 1 && etaFar < 1;
numericalQPACOK = internalAuditOK && qpacInequalitiesOK;
qpacSlack = 1-CQPAC;

usedAnalyticKernel = any(D2SelectedBranchByRow == "analytic-kernel") || ...
    any(D3SelectedBranchByRow == "analytic-kernel") || ...
    any(farSelectedBranchByRow == "analytic-kernel");
usedDirSum = any(D2SelectedBranchByRow == "DirSum") || ...
    any(farSelectedBranchByRow == "DirSum");
usedDirectQ3 = any(D3SelectedBranchByRow == "direct-Q3");
methodSelectionOK = ...
    (~usedAnalyticKernel || par.methods.enableAnalyticKernel) && ...
    (~usedDirSum || par.methods.enableDirSum) && ...
    (~usedDirectQ3 || par.methods.enableDirectQ3);
assert(methodSelectionOK,'QPAC:DisabledMethodSelected', ...
    'A disabled bound method was selected.');
if ~numericalQPACOK
    warning('QPAC:NumericalConditionNotMet', ...
        ['The QPAC checks did not pass. Inspect the saved ', ...
         'budgets and audit flags; no recovery claim is made by this run.']);
end

%% 6. Assemble and save results
result = struct();
result.parameters = par;
result.boundSelection = struct( ...
    'policy','minimum finite bound among enabled valid methods', ...
    'enabledMethods',par.methods);
result.support = support;
result.taper = struct('rho',rho,'b',geom.b,'n',geom.n);
result.geometry = geom;
result.phaseCoordinates = phaseCoordinates;
result.derivativeDiagnostic = struct( ...
    'dNplus',dNplus,'threshold',par.derivativeThreshold, ...
    'thresholdPassed',derivativeThresholdPassed);
result.LCS = struct( ...
    'channelNames',{channelNames}, ...
    'bound',lcsBound,'bestQ',lcsBestQ, ...
    'exact',exactSupport, ...
    'maxRelativeLagExpansionResidual',lcsRelativeLagResidual, ...
    'supportUsesNontrivialLCS',supportUsesNontrivialLCS, ...
    'lagExpansionAuditOK',lagExpansionAuditOK, ...
    'dominatesExact',lcsDominatesExact);
result.GSS = GSS;
result.MSS = MSS;
result.etaSS = etaSS;
result.Gamma = Gamma;
result.Xi = Xi;
result.inverseResidual = inverseResidual;
result.selfIdentityError = selfIdentityError;
result.sigmaLower2 = sigmaLower2;
result.qMin = qMin;
result.qMax = qMax;
result.UK2Self = UK2Self;
result.UH2Self = UH2Self;
result.UK2SS = UK2SS;
result.UH2SS = UH2SS;
result.Ecurv = Ecurv;
result.mNear = mNear;
result.D2 = D2;
result.D3 = D3;
result.etaNear = etaNear;
result.etaFar = etaFar;
result.CQPAC = CQPAC;
result.row = row;
result.padding = struct( ...
    'D2ByRow',cellfun(@(z) z.D2Pad,row), ...
    'farByRow',cellfun(@(z) z.farPad,row));
result.directD3 = struct( ...
    'segment',localSegment, ...
    'boundByRow',D3DirectByRow, ...
    'sampleMaximumByRow',D3SampleMaximumByRow, ...
    'dominatesSamplesByRow',D3DirectByRow + ...
        floating_pad(D3DirectByRow,par.roundoffFactor) >= ...
        D3SampleMaximumByRow, ...
    'usesFourthDerivative',false);
result.kernelD2 = struct( ...
    'method',['support-uniform five-atom analytical envelope ', ...
        '(Der/LCS/ResLin plus fast-cap alternative)'], ...
    'segment',localSegment, ...
    'cellEdges',localCellEdges, ...
    'analyticBoundByRow',D2AnalyticByRow, ...
    'dirSumBoundByRow',D2DirSumByRow, ...
    'selectedBoundByRow',D2SelectedByRow, ...
    'selectedBranchByRow',D2SelectedBranchByRow, ...
    'dominatesSamplesByRow',selectedD2DominatesSamplesByRow, ...
    'pointAuditOKByRow',D2AnalyticAuditByRow);
result.kernelD3 = struct( ...
    'method',['support-uniform nine-atom analytical envelope ', ...
        '(Der/LCS/ResLin plus fast-cap alternative)'], ...
    'segment',localSegment, ...
    'cellEdges',localCellEdges, ...
    'boundByRow',D3KernelByRow, ...
    'selectedBoundByRow',D3SelectedByRow, ...
    'selectedBranchByRow',D3SelectedBranchByRow, ...
    'sampleMaximumByRow',D3SampleMaximumByRow, ...
    'dominatesSamplesByRow',selectedD3DominatesSamplesByRow, ...
    'pointAuditOKByRow',D3KernelAuditByRow, ...
    'usesFourthDerivative',false, ...
    'QList',par.kernelBound.QList, ...
    'pList',par.kernelBound.pList);
result.farBounds = struct( ...
    'method',['support-uniform normalized K/H analytical envelope ', ...
        '(Der/LCS/ResLin) versus common-point DirSum'], ...
    'leftCellEdges',farLeftCellEdges, ...
    'rightCellEdges',farRightCellEdges, ...
    'analyticBoundByRow',farAnalyticByRow, ...
    'dirSumBoundByRow',farDirSumByRow, ...
    'selectedBoundByRow',farSelectedByRow, ...
    'selectedBranchByRow',farSelectedBranchByRow, ...
    'dominatesSamplesByRow',selectedFarDominatesSamplesByRow, ...
    'pointAuditOKByRow',farAnalyticAuditByRow);
result.coarseAudit = struct( ...
    'D2',D2Coarse,'D3',D3Coarse, ...
    'etaNear',etaNearCoarse,'etaFar',etaFarCoarse);
result.audit = struct( ...
    'flatEndOK',flatEndOK, ...
    'fresnelRegionOK',par.fresnelRegionOK, ...
    'supportCellCoverageOK',supportCellCoverageOK, ...
    'localSegmentCoverageOK',localSegmentCoverageOK, ...
    'kernelLocalCellCoverageOK',kernelLocalCellCoverageOK, ...
    'kernelFarCellCoverageOK',kernelFarCellCoverageOK, ...
    'farCoverageOK',farCoverageOK, ...
    'farBoundaryIncluded',farBoundaryIncluded, ...
    'commonPointNearByConstruction',commonPointNearByConstruction, ...
    'commonPointFarByConstruction',commonPointFarByConstruction, ...
    'sameRadiusByConstruction',sameRadiusByConstruction, ...
    'selfIdentitiesOK',selfIdentitiesOK, ...
    'supportUsesNontrivialLCS',supportUsesNontrivialLCS, ...
    'lagExpansionAuditOK',lagExpansionAuditOK, ...
    'selectedD2DominatesSamples',selectedD2DominatesSamples, ...
    'kernelD2PointAuditOK',all(D2AnalyticAuditByRow), ...
    'selectedD3DominatesSamples',selectedD3DominatesSamples, ...
    'kernelD3PointAuditOK',all(D3KernelAuditByRow), ...
    'selectedFarDominatesSamples',selectedFarDominatesSamples, ...
    'kernelFarPointAuditOK',all(farAnalyticAuditByRow), ...
    'coefficientSignsOK',coefficientSignsOK, ...
    'usedAnalyticKernel',usedAnalyticKernel, ...
    'usedDirSum',usedDirSum, ...
    'usedDirectQ3',usedDirectQ3, ...
    'methodSelectionOK',methodSelectionOK, ...
    'internalAuditOK',internalAuditOK, ...
    'qpacInequalitiesOK',qpacInequalitiesOK, ...
    'allFinite',allFinite);
result.numericalQPACOK = numericalQPACOK;
result.qpacSlack = qpacSlack;
result.arithmetic = 'MATLAB double precision';
result.continuumCoverage = 'analytical cell enclosures';
result.envelopeSystem = [ ...
    'fixed-point LCS for support; support-uniform analytical kernels for ', ...
    'D2/D3/far; common-row DirSum; direct uniform Q3'];
result.software = struct('function',mfilename,'matlabVersion',version);

thisDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(thisDir,'results');
figuresDir = fullfile(thisDir,'figures');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end
if ~exist(figuresDir,'dir')
    mkdir(figuresDir);
end

summary = make_summary_table(result);
csvPath = fullfile(resultsDir,'correlation_qpac_summary.csv');
matPath = fullfile(resultsDir,'correlation_qpac_results.mat');
texPath = fullfile(resultsDir,'correlation_qpac_table.tex');
writetable(summary,csvPath);
tableWritten = write_latex_table(texPath,result);
figureWritten = write_diagnostic_figure( ...
    figuresDir,nearGrid,farLeftGrid,farRightGrid,result);
result.outputStatus = struct( ...
    'csvWritten',exist(csvPath,'file') == 2, ...
    'tableWritten',tableWritten, ...
    'figureWritten',figureWritten, ...
    'csvPath',csvPath,'matPath',matPath,'texPath',texPath, ...
    'figuresDir',figuresDir);
save(matPath,'result');
print_summary(result,summary);
end

%% Configuration
function par = default_config()
par = struct();
par.c0 = 3e8;
par.fc = 100e9;
par.lambda = par.c0/par.fc;
par.k = 2*pi/par.lambda;
par.d = par.lambda/2;
par.Nr = 256;
par.L = 2;

par.theta0 = pi/2;
par.thetaHalfWidth = 0.15;
par.thetaDomain = par.theta0+[-par.thetaHalfWidth,par.thetaHalfWidth];
% Fixed localization radius for this reproducibility instance.
par.delta = 0.0066;

par.QList = 2:16;
par.nNear = 4001;
par.nFar = 8001;
par.chunkSize = 1000;
par.roundoffFactor = 4096;
par.auditTolerance = 1e-11;
par.selfIdentityTolerance = 2e-11;
par.lagResidualTolerance = 1e-8;
par.derivativeThreshold = 0.1;

% Support-uniform analytic kernel bounds for K/H and their second and
% third evaluation-angle derivatives.
par.kernelBound = struct();
par.kernelBound.QList = [2 8 16 32];
par.kernelBound.pList = 1:8;
par.kernelBound.localCells = 16;
% More cells give tighter far-region envelopes at a higher runtime.
par.kernelBound.farCellsPerSide = 64;
par.kernelBound.lcsSafety = 1e-12;
par.kernelBound.autoDeltaSafety = 0.98;
par.kernelBound.autoDeltaMin = 1e-12;
par.kernelBound.useFastCauchyCap = true;
par.kernelBound.useDerivativeBranch = true;
par.kernelBound.derivativeD0List = [ ...
    0.15 0.20 0.30 0.40 0.55 0.70 0.85 1.00 ...
    1.20 1.40 1.70 2.00 2.40 2.80 3.05];
par.kernelBound.derivativeRhoList = [ ...
    1e-8 3e-8 1e-7 3e-7 1e-6 3e-6 1e-5 3e-5 ...
    1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 0.07 0.15 ...
    0.3 0.6 1.0];
par.kernelBound.useResidueLinearBranch = true;
par.kernelBound.linearQList = [1 2 3 4 5 6 8 10 12 16 20 24 32];
par.kernelBound.linearAWindow = 0.12;
par.kernelBound.autoGammaSafety = 0.98;
par.kernelBound.autoGammaMin = 1e-12;

% Select the smallest finite enabled bound independently on each range row.
par.methods = struct();
par.methods.enableAnalyticKernel = true;
par.methods.enableDirSum = true;
par.methods.enableDirectQ3 = true;
end

function validate_config(par)
validateattributes(par.fc,{'double'},{'scalar','real','finite','positive'});
validateattributes(par.Nr,{'double'},{'scalar','integer','>=',10});
validateattributes(par.d,{'double'},{'scalar','real','finite','positive'});
validateattributes(par.delta,{'double'},{'scalar','real','finite','positive'});
validateattributes(par.nNear,{'double'},{'scalar','integer','>=',3});
validateattributes(par.nFar,{'double'},{'scalar','integer','>=',3});
validateattributes(par.chunkSize,{'double'}, ...
    {'scalar','integer','finite','positive'});
validateattributes(par.roundoffFactor,{'double'}, ...
    {'scalar','real','finite','positive'});
validateattributes(par.theta0,{'double'}, ...
    {'scalar','real','finite','>',0,'<',pi});
validateattributes(par.thetaDomain,{'double'}, ...
    {'vector','numel',2,'real','finite'});
assert(par.L == 2,'QPAC:SourceCount', ...
    'This fixed correlation example is implemented for exactly two sources.');
assert(mod(par.nNear,2) == 1 && mod(par.nFar,2) == 1, ...
    'QPAC:EvenGridSize','nNear and nFar must be odd.');
assert(par.thetaDomain(1) > 0 && par.thetaDomain(2) < pi && ...
    par.thetaDomain(1) < par.theta0-par.delta && ...
    par.theta0+par.delta < par.thetaDomain(2), ...
    'QPAC:InvalidAngularDomain', ...
    'The connected domain must contain nonempty near and far intervals.');
assert(all(par.QList >= 2 & par.QList <= par.Nr & ...
    par.QList == round(par.QList)), ...
    'QPAC:InvalidQList','Every residue denominator must be in [2,Nr].');
assert(all(par.kernelBound.QList >= 2 & ...
    par.kernelBound.QList <= par.Nr & ...
    par.kernelBound.QList == round(par.kernelBound.QList)), ...
    'QPAC:InvalidKernelQList', ...
    'Every kernel-bound residue denominator must be in [2,Nr].');
assert(all(par.kernelBound.pList >= 1 & ...
    par.kernelBound.pList == round(par.kernelBound.pList)), ...
    'QPAC:InvalidKernelPList', ...
    'Every Abel difference order must be a positive integer.');
validateattributes(par.kernelBound.localCells,{'double'}, ...
    {'scalar','integer','>=',1});
validateattributes(par.kernelBound.farCellsPerSide,{'double'}, ...
    {'scalar','integer','>=',1});
validateattributes(par.kernelBound.lcsSafety,{'double'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(par.kernelBound.autoDeltaSafety,{'double'}, ...
    {'scalar','real','finite','>',0,'<=',1});
validateattributes(par.kernelBound.autoDeltaMin,{'double'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(par.kernelBound.useFastCauchyCap,{'logical'}, ...
    {'scalar'});
validateattributes(par.kernelBound.useDerivativeBranch,{'logical'}, ...
    {'scalar'});
validateattributes(par.kernelBound.useResidueLinearBranch,{'logical'}, ...
    {'scalar'});
assert(all(isfinite(par.kernelBound.derivativeD0List)) && ...
    all(par.kernelBound.derivativeD0List > 0) && ...
    all(par.kernelBound.derivativeD0List < pi), ...
    'QPAC:InvalidDerivativeD0List', ...
    'Every derivative d0 threshold must lie strictly between 0 and pi.');
assert(all(isfinite(par.kernelBound.derivativeRhoList)) && ...
    all(par.kernelBound.derivativeRhoList >= 0) && ...
    all(par.kernelBound.derivativeRhoList <= 1), ...
    'QPAC:InvalidDerivativeRhoList', ...
    'Every derivative rho threshold must lie in [0,1].');
validateattributes(par.kernelBound.linearAWindow,{'double'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(par.kernelBound.autoGammaSafety,{'double'}, ...
    {'scalar','real','finite','>',0,'<=',1});
validateattributes(par.kernelBound.autoGammaMin,{'double'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(par.methods.enableAnalyticKernel,{'logical'}, ...
    {'scalar'});
validateattributes(par.methods.enableDirSum,{'logical'}, ...
    {'scalar'});
validateattributes(par.methods.enableDirectQ3,{'logical'}, ...
    {'scalar'});
assert(par.methods.enableAnalyticKernel || par.methods.enableDirSum, ...
    'QPAC:NoEnabledD2FarMethod', ...
    'Enable at least one D2/far bound method.');
assert(par.methods.enableAnalyticKernel || par.methods.enableDirectQ3, ...
    'QPAC:NoEnabledD3Method', ...
    'Enable at least one D3 bound method.');
assert(all(par.kernelBound.linearQList >= 1 & ...
    par.kernelBound.linearQList <= par.Nr & ...
    par.kernelBound.linearQList == round(par.kernelBound.linearQList)), ...
    'QPAC:InvalidLinearQList', ...
    'Every residue-linear denominator must be in [1,Nr].');
end

%% Geometry and channel evaluation
function [geom,rho] = build_geometry(par)
n = (0:par.Nr-1).';
rho = choose4(n).*choose4(par.Nr-1-n);
assert(sum(rho) > 0,'QPAC:ZeroTaper','The taper has zero mass.');
b = rho/sum(rho);
nbar = sum(b.*n);
n2 = n.^2;
n2bar = sum(b.*n2);
x = n-nbar;
y = n2-n2bar;
G = [sum(b.*x.^2),sum(b.*x.*y); ...
     sum(b.*x.*y),sum(b.*y.^2)];
assert(min(eig(G)) > 0,'QPAC:SingularTangentGram', ...
    'The tangent Gram matrix is singular.');

geom = struct();
geom.n = n;
geom.n2 = n2;
geom.b = b;
geom.x = x;
geom.y = y;
geom.G = G;
geom.kappa = par.k*par.d;
end

function v = choose4(m)
v = zeros(size(m));
mask = m >= 4;
z = m(mask);
v(mask) = z.*(z-1).*(z-2).*(z-3)/24;
end

function [omega1,omega2] = phase_coordinates_pair( ...
    rSource,thetaSource,rEval,thetaEval,par)
omega1 = par.k*par.d*(cos(thetaSource)-cos(thetaEval));
omega2 = 0.5*par.k*par.d^2*( ...
    sin(thetaEval)^2/rEval-sin(thetaSource)^2/rSource);
end

function h = tangent_vector(r,theta,geom,par)
z = geom.x+(par.d/r)*cos(theta)*geom.y;
q = sum(geom.b.*z.^2);
assert(q > 0,'QPAC:NonpositiveTangentNorm', ...
    'The tangent normalization is nonpositive.');
h = z/sqrt(q);
end

function coeff = fixed_channel_coefficients( ...
    rSource,thetaSource,rEval,thetaEval,geom,par)
hs = tangent_vector(rSource,thetaSource,geom,par);
he = tangent_vector(rEval,thetaEval,geom,par);

alpha = par.d/rEval;
c = cos(thetaEval);
s = sin(thetaEval);
z = geom.x+alpha*c*geom.y;
u = -geom.kappa*s*z;
ut = -geom.kappa*(c*geom.x+alpha*cos(2*thetaEval)*geom.y);
utt = geom.kappa*s*(geom.x+4*alpha*c*geom.y);
Q2 = -1i*ut-u.^2;
Q3 = -1i*utt-3*u.*ut+1i*u.^3;

b = geom.b;
coeff = { ...
    b, ...
    1i*b.*hs, ...
    -1i*b.*he, ...
    b.*he.*hs, ...
    b.*Q2, ...
    b.*Q3, ...
    1i*b.*hs.*Q2, ...
    1i*b.*hs.*Q3};
end

function C = evaluate_source( ...
    rSource,thetaSource,rEval,thetaEval,geom,par)
thetaEval = thetaEval(:).';
c = cos(thetaEval);
s = sin(thetaEval);

omega1 = geom.kappa*(cos(thetaSource)-c);
omega2 = 0.5*par.k*par.d^2*( ...
    s.^2/rEval-sin(thetaSource)^2/rSource);
phase = exp(1i*(geom.n*omega1+geom.n2*omega2));

hs = tangent_vector(rSource,thetaSource,geom,par);
zEval = geom.x+(par.d/rEval)*geom.y.*c;
qEval = sum(geom.b.*zEval.^2,1);
assert(all(qEval > 0),'QPAC:NonpositiveEvaluationTangentNorm', ...
    'An evaluation-side tangent normalization is nonpositive.');
hEval = zEval./sqrt(qEval);

u = -geom.kappa*s.*zEval;
ut = -geom.kappa*(c.*geom.x+ ...
    (par.d/rEval)*cos(2*thetaEval).*geom.y);
utt = geom.kappa*s.*(geom.x+ ...
    4*(par.d/rEval)*c.*geom.y);
Q2 = -1i*ut-u.^2;
Q3 = -1i*utt-3*u.*ut+1i*u.^3;

b = geom.b;
C = zeros(8,numel(thetaEval));
C(1,:) = abs(sum(b.*phase,1));
C(2,:) = abs(sum(1i*b.*hs.*phase,1));
C(3,:) = abs(sum(-1i*b.*hEval.*phase,1));
C(4,:) = abs(sum(b.*hEval.*hs.*phase,1));
C(5,:) = abs(sum(b.*Q2.*phase,1));
C(6,:) = abs(sum(b.*Q3.*phase,1));
C(7,:) = abs(sum(1i*b.*hs.*Q2.*phase,1));
C(8,:) = abs(sum(1i*b.*hs.*Q3.*phase,1));
end

%% Fixed-point lag-correlation bound
function out = lcs_fixed_point_bound( ...
    a,omega1,omega2,QList,roundoffFactor)
% For a singleton physical cell, the lag phase is a point.  Expanding each
% residue energy gives
%
% |S_s|^2 = sum_m |a_m|^2
%   + 2 sum_{h>=1} sum_m Re{a_{m+h} conj(a_m) exp(i DeltaPhi)}.
%
% The LCS upper bound is sum_s sqrt(upper enclosure of |S_s|^2).
a = a(:);
N = numel(a);
n = (0:N-1).';
l1Cap = sum(abs(a));
best = l1Cap+floating_pad(l1Cap,roundoffFactor);
bestQ = NaN; % NaN denotes the trivial triangle cap, not an LCS split.
maxRelativeResidualAtBest = NaN;

for Q = QList
    residueSum = 0;
    maxRelativeResidual = 0;
    for s = 0:Q-1
        idx = (s+1):Q:N;
        ni = n(idx);
        ai = a(idx);
        M = numel(idx);
        if M == 0
            continue;
        end

        lagEnergy = sum(abs(ai).^2);
        absoluteAccumulation = sum(abs(ai).^2);
        for h = 1:M-1
            for m = 1:M-h
                n0 = ni(m);
                n1 = ni(m+h);
                deltaPhase = (n1-n0)*omega1 + ...
                    (n1^2-n0^2)*omega2;
                lagProduct = ...
                    ai(m+h)*conj(ai(m))*exp(1i*deltaPhase);
                lagEnergy = lagEnergy+2*real(lagProduct);
                absoluteAccumulation = ...
                    absoluteAccumulation+2*abs(lagProduct);
            end
        end

        z = ai.*exp(1i*(ni*omega1+ni.^2*omega2));
        directEnergy = abs(sum(z))^2;
        relativeResidual = abs(lagEnergy-directEnergy)/ ...
            max([1,abs(lagEnergy),directEnergy]);
        maxRelativeResidual = max(maxRelativeResidual,relativeResidual);
        % Scale padding by the absolute accumulated lag mass and operation
        % count. This is a numerical safeguard, not directed rounding.
        energyPad = floating_pad(absoluteAccumulation,roundoffFactor)* ...
            max(1,M^2);
        if real(lagEnergy) < -energyPad
            error('QPAC:NegativeResidueEnergy', ...
                ['The correlation-expanded residue energy is negative ', ...
                 'beyond its arithmetic safeguard.']);
        end
        residueSum = residueSum+sqrt(max(real(lagEnergy)+energyPad,0));
    end

    candidate = min(l1Cap,residueSum);
    candidate = candidate+floating_pad(candidate,roundoffFactor);
    if candidate < best
        best = candidate;
        bestQ = Q;
        maxRelativeResidualAtBest = maxRelativeResidual;
    end
end

out = struct('bound',max(real(best),0), ...
    'bestQ',bestQ, ...
    'maxRelativeLagExpansionResidual',maxRelativeResidualAtBest);
end

%% Support-uniform analytic kernel bounds
function cache = build_analytic_kernel_cache( ...
    support,geom,par)
% The Q3 dictionary is exactly the one used in
% run_bounds_true_kernel_repro.m.  Coefficients are preprocessed once for
% every requested residue denominator and Abel difference order.
b = geom.b(:);
x = geom.x(:);
y = geom.y(:);
atoms = { ...
    b.*x; ...
    b.*y; ...
    b.*(x.^2); ...
    b.*(x.*y); ...
    b.*(y.^2); ...
    b.*(x.^3); ...
    b.*(x.^2.*y); ...
    b.*(x.*y.^2); ...
    b.*(y.^3)};

cache = struct();
cache.atoms = atoms;
cache.normalizedK = precompute_kernel_lcs_coefficient(b,par);
cache.normalizedH = cell(numel(support.range),1);
cache.K = cell(numel(atoms),1);
cache.H = cell(numel(support.range),numel(atoms));
for ia = 1:numel(atoms)
    cache.K{ia} = precompute_kernel_lcs_coefficient(atoms{ia},par);
end
for ell = 1:numel(support.range)
    hs = tangent_vector(support.range(ell),support.theta(ell),geom,par);
    cache.normalizedH{ell} = precompute_kernel_lcs_coefficient( ...
        1i*b.*hs,par);
    for ia = 1:numel(atoms)
        cache.H{ell,ia} = precompute_kernel_lcs_coefficient( ...
            1i*atoms{ia}.*hs,par);
    end
end
end

function info = higher_derivative_kernel_row_bound( ...
    order,rEval,cellEdges,support,Gamma,cache,geom,par)
% Bound
%   sum_l Gamma_K |d_theta^3 K(q,p_l)|
%       + Gamma_H |d_theta^3 H(q,p_l)|
% on every closed angular cell.  The source sum is taken before the
% maximum over cells, as required by the QPAC far-region condition.
assert(numel(cellEdges) >= 2 && all(diff(cellEdges) > 0), ...
    'QPAC:InvalidKernelCells', ...
    'The kernel-bound cell edges must be strictly increasing.');

nCells = numel(cellEdges)-1;
nSource = numel(support.range);
if order == 2
    nAtoms = 5;
    exactKRow = 5;
    exactHRow = 7;
elseif order == 3
    nAtoms = 9;
    exactKRow = 6;
    exactHRow = 8;
else
    error('QPAC:InvalidHigherOrder','Higher-kernel order must be 2 or 3.');
end
pairK = zeros(nSource,nCells);
pairH = zeros(nSource,nCells);
cellBounds = zeros(1,nCells);

for ic = 1:nCells
    thetaCell = [cellEdges(ic),cellEdges(ic+1)];
    gamma = higher_derivative_dictionary_gamma( ...
        order,rEval,thetaCell,geom,par);

    for ell = 1:nSource
        phaseBox = singleton_phase_box( ...
            support.range(ell),support.theta(ell), ...
            rEval,thetaCell,geom,par);

        BK = 0;
        BH = 0;
        for ia = 1:nAtoms
            BK = BK+gamma(ia)*kernel_scalar_analytic_cell_bound( ...
                cache.K{ia},phaseBox,par);
            BH = BH+gamma(ia)*kernel_scalar_analytic_cell_bound( ...
                cache.H{ell,ia},phaseBox,par);
        end

        if par.kernelBound.useFastCauchyCap
            [fastK,fastH] = higher_derivative_fast_caps( ...
                order,support.range(ell),support.theta(ell), ...
                rEval,thetaCell,geom,par);
            BK = min(BK,fastK);
            BH = min(BH,fastH);
        end

        pairK(ell,ic) = max(real(BK),0);
        pairH(ell,ic) = max(real(BH),0);
    end

    cellBounds(ic) = Gamma(1)*sum(pairK(:,ic)) + ...
        Gamma(2)*sum(pairH(:,ic));
end

% Check each cell envelope at its endpoints and midpoint.
pointAuditOK = true;
maxPointRatioK = 0;
maxPointRatioH = 0;
for ic = 1:nCells
    thetaAudit = [cellEdges(ic), ...
        0.5*(cellEdges(ic)+cellEdges(ic+1)),cellEdges(ic+1)];
    for ell = 1:nSource
        C = evaluate_source( ...
            support.range(ell),support.theta(ell), ...
            rEval,thetaAudit,geom,par);
        exactK = max(C(exactKRow,:));
        exactH = max(C(exactHRow,:));
        toleranceK = floating_pad(pairK(ell,ic),par.roundoffFactor);
        toleranceH = floating_pad(pairH(ell,ic),par.roundoffFactor);
        pointAuditOK = pointAuditOK && ...
            pairK(ell,ic)+toleranceK >= exactK && ...
            pairH(ell,ic)+toleranceH >= exactH;
        maxPointRatioK = max(maxPointRatioK, ...
            exactK/max(pairK(ell,ic),realmin));
        maxPointRatioH = max(maxPointRatioH, ...
            exactH/max(pairH(ell,ic),realmin));
    end
end

bound = max(cellBounds);
bound = bound+floating_pad(bound,par.roundoffFactor);
info = struct('bound',bound,'cellBounds',cellBounds, ...
    'pairK',pairK,'pairH',pairH,'cellEdges',cellEdges, ...
    'pointAuditOK',pointAuditOK, ...
    'maxPointRatioK',maxPointRatioK, ...
    'maxPointRatioH',maxPointRatioH);
end

function info = normalized_kernel_far_row_bound( ...
    rEval,edgeSets,support,Gamma,cache,geom,par)
% Support-uniform K/H envelope on the disconnected far set.  Each source
% is evaluated on the same physical cell; only after summation is the
% maximum over cells taken.
nSource = numel(support.range);
nCells = sum(cellfun(@(e) numel(e)-1,edgeSets));
pairK = zeros(nSource,nCells);
pairH = zeros(nSource,nCells);
cellBounds = zeros(1,nCells);
cellIntervals = zeros(nCells,2);
branchK = strings(nSource,nCells);
branchH = strings(nSource,nCells);

ic = 0;
for iset = 1:numel(edgeSets)
    edges = edgeSets{iset};
    assert(numel(edges) >= 2 && all(diff(edges) > 0), ...
        'QPAC:InvalidFarKernelCells', ...
        'Far analytical cell edges must be strictly increasing.');
    for jc = 1:(numel(edges)-1)
        ic = ic+1;
        thetaCell = [edges(jc),edges(jc+1)];
        cellIntervals(ic,:) = thetaCell;
        for ell = 1:nSource
            phaseBox = singleton_phase_box( ...
                support.range(ell),support.theta(ell), ...
                rEval,thetaCell,geom,par);
            [BK,labelK] = kernel_scalar_analytic_cell_bound( ...
                cache.normalizedK,phaseBox,par);
            [BH,labelH] = kernel_scalar_analytic_cell_bound( ...
                cache.normalizedH{ell},phaseBox,par);
            if BK > 1
                BK = 1;
                labelK = "unit-cap";
            end
            if BH > 1
                BH = 1;
                labelH = "unit-cap";
            end
            pairK(ell,ic) = max(real(BK),0);
            pairH(ell,ic) = max(real(BH),0);
            branchK(ell,ic) = labelK;
            branchH(ell,ic) = labelH;
        end
        cellBounds(ic) = Gamma(1)*sum(pairK(:,ic)) + ...
            Gamma(2)*sum(pairH(:,ic));
    end
end

% Check each cell envelope at its endpoints and midpoint.
pointAuditOK = true;
maxPointRatioK = 0;
maxPointRatioH = 0;
for ic = 1:nCells
    thetaAudit = [cellIntervals(ic,1),mean(cellIntervals(ic,:)), ...
        cellIntervals(ic,2)];
    for ell = 1:nSource
        C = evaluate_source( ...
            support.range(ell),support.theta(ell), ...
            rEval,thetaAudit,geom,par);
        exactK = max(C(1,:));
        exactH = max(C(2,:));
        toleranceK = floating_pad(pairK(ell,ic),par.roundoffFactor);
        toleranceH = floating_pad(pairH(ell,ic),par.roundoffFactor);
        pointAuditOK = pointAuditOK && ...
            pairK(ell,ic)+toleranceK >= exactK && ...
            pairH(ell,ic)+toleranceH >= exactH;
        maxPointRatioK = max(maxPointRatioK, ...
            exactK/max(pairK(ell,ic),realmin));
        maxPointRatioH = max(maxPointRatioH, ...
            exactH/max(pairH(ell,ic),realmin));
    end
end

bound = max(cellBounds);
bound = bound+floating_pad(bound,par.roundoffFactor);
info = struct('bound',bound,'cellBounds',cellBounds, ...
    'cellIntervals',cellIntervals,'pairK',pairK,'pairH',pairH, ...
    'branchK',branchK,'branchH',branchH, ...
    'pointAuditOK',pointAuditOK, ...
    'maxPointRatioK',maxPointRatioK, ...
    'maxPointRatioH',maxPointRatioH);
end

function gamma = higher_derivative_dictionary_gamma( ...
    order,rEval,thetaCell,geom,par)
% Coefficient envelopes for the nine Q3 dictionary atoms.  This is the
% order-three branch of the higher-derivative dictionary construction.
[cLo,cHi] = cos_interval(thetaCell);
[cAbsMin,cMax] = absolute_interval_range(cLo,cHi);
[~,s2Max] = sin2_interval(thetaCell);
sMax = sqrt(max(0,s2Max));
alpha = par.d/rEval;
kappa = geom.kappa;

p2m1Max = max_abs_affine_in_c2(cLo,cHi,-1,2);

if order == 2
    cm1c2Max = max_t_times_one_minus_t2(cAbsMin,cMax);
    c2m1c2Max = max_t2_times_one_minus_t2(cAbsMin,cMax);
    gamma = zeros(5,1);
    gamma(1) = kappa*cMax;
    gamma(2) = kappa*alpha*p2m1Max;
    gamma(3) = kappa^2*s2Max;
    gamma(4) = 2*kappa^2*alpha*cm1c2Max;
    gamma(5) = kappa^2*alpha^2*c2m1c2Max;
elseif order == 3
    p3m1Max = max_abs_affine_in_c2(cLo,cHi,-1,3);
    c2Max = cMax^2;
    c3Max = cMax^3;
    gamma = zeros(9,1);
    gamma(1) = kappa*sMax;
    gamma(2) = 4*kappa*alpha*cMax*sMax;
    gamma(3) = 3*kappa^2*cMax*sMax;
    gamma(4) = 3*kappa^2*alpha*p3m1Max*sMax;
    gamma(5) = 3*kappa^2*alpha^2*cMax*p2m1Max*sMax;
    gamma(6) = kappa^3*sMax^3;
    gamma(7) = 3*kappa^3*alpha*cMax*sMax^3;
    gamma(8) = 3*kappa^3*alpha^2*c2Max*sMax^3;
    gamma(9) = kappa^3*alpha^3*c3Max*sMax^3;
else
    error('QPAC:InvalidHigherOrder','Higher-kernel order must be 2 or 3.');
end

assert(cAbsMin >= 0 && all(gamma >= 0), ...
    'QPAC:NegativeDictionaryGamma', ...
    'A third-derivative dictionary envelope became negative.');
end

function box = singleton_phase_box( ...
    rSource,thetaSource,rEval,thetaCell,geom,par)
% Exact rectangular hull of (omega_1,omega_2) for a singleton source and
% fixed evaluation range over one angular cell.
[cLo,cHi] = cos_interval(thetaCell);
[s2Lo,s2Hi] = sin2_interval(thetaCell);
muSource = geom.kappa*cos(thetaSource);
muEval = geom.kappa*[cLo,cHi];
etaFactor = 0.5*par.k*par.d^2;
etaSource = etaFactor*sin(thetaSource)^2/rSource;
etaEval = etaFactor*[s2Lo,s2Hi]/rEval;
box = struct();
box.omega1 = [muSource-muEval(2),muSource-muEval(1)];
box.omega2 = [etaEval(1)-etaSource,etaEval(2)-etaSource];
end

function [BK,BH] = higher_derivative_fast_caps( ...
    order,rSource,thetaSource,rEval,thetaCell,geom,par)
% Fast l1/Cauchy caps provide alternatives to the dictionary bound.
[cLo,cHi] = cos_interval(thetaCell);
[c2Lo,c2Hi] = cos_interval(2*thetaCell);
[sLo,sHi] = cos_interval(thetaCell-pi/2);
sMax = max(abs([sLo,sHi]));
alpha = par.d/rEval;
tauEval = alpha*[cLo,cHi];
zEval = alpha*[c2Lo,c2Hi];

x = geom.x(:);
y = geom.y(:);
b = geom.b(:);
Aub = zeros(size(x));
Bub = zeros(size(x));
A4ub = zeros(size(x));
for nn = 1:numel(x)
    Aub(nn) = max_abs_affine_interval(x(nn),y(nn),tauEval);
    Bub(nn) = max_abs_bilinear_box( ...
        x(nn),y(nn),[cLo,cHi],zEval);
    A4ub(nn) = max_abs_affine_interval(x(nn),4*y(nn),tauEval);
end

if order == 2
    profileUpper = geom.kappa*Bub + ...
        geom.kappa^2*sMax^2.*Aub.^2;
elseif order == 3
    profileUpper = geom.kappa*sMax*A4ub + ...
        3*geom.kappa^2*sMax.*Aub.*Bub + ...
        geom.kappa^3*sMax^3.*Aub.^3;
else
    error('QPAC:InvalidHigherOrder','Higher-kernel order must be 2 or 3.');
end
hs = tangent_vector(rSource,thetaSource,geom,par);

K_l1 = sum(b.*profileUpper);
H_l1 = sum(b.*abs(hs).*profileUpper);
K_cauchy = sqrt(sum(b.*profileUpper.^2));
% sum_n b_n |h_s(n)|^2 = 1 by tangent normalization.
H_cauchy = K_cauchy;
BK = min(K_l1,K_cauchy);
BH = min(H_l1,H_cauchy);
end

function coeff = precompute_kernel_lcs_coefficient(a,par)
a = a(:);
qList = par.kernelBound.QList(:).';
pMax = max(par.kernelBound.pList);
coeff = struct();
coeff.a = a;
coeff.N = numel(a);
coeff.l1 = sum(abs(a));
coeff.diffNorms0to4 = kernel_diff_norms_0_to_4(a);
coeff.flatEndsR4 = coeff.N >= 9 && ...
    all(abs(a(1:4)) <= 1e-13) && ...
    all(abs(a(end-3:end)) <= 1e-13);
coeff.QList = qList;
coeff.Qstats = cell(numel(qList),1);

for iq = 1:numel(qList)
    Q = qList(iq);
    statsCell = cell(1,Q);
    for s = 1:Q
        statsCell{s} = kernel_coefficient_stats(a(s:Q:end),pMax);
    end
    coeff.Qstats{iq} = statsCell;
end
end

function stats = kernel_coefficient_stats(a,pMax)
a = a(:);
N = numel(a);
stats = struct();
stats.N = N;
stats.E = sum(abs(a).^2);
stats.C = zeros(max(N-1,0),1);
stats.W = cell(pMax,1);
for p = 1:pMax
    stats.W{p} = inf(max(N-1,0),1);
end

for h = 1:(N-1)
    corrSequence = a(1+h:N).*conj(a(1:N-h));
    [Csum,Wrow] = kernel_variation_profile(corrSequence,pMax);
    stats.C(h) = Csum;
    for p = 1:pMax
        stats.W{p}(h) = Wrow(p);
    end
end
end

function [Csum,Wrow] = kernel_variation_profile(c,pMax)
% Abel finite-difference variation constants W_{h,p}.
c = c(:);
Csum = sum(abs(c));
M = numel(c)-1;
Wrow = inf(1,pMax);

if isempty(c)
    Wrow(:) = 0;
    return;
end

differences = cell(pMax+1,1);
differences{1} = c;
for p = 1:pMax
    if numel(differences{p}) >= 2
        differences{p+1} = diff(differences{p});
    else
        differences{p+1} = [];
    end
end

for p = 1:pMax
    if p > M+1
        continue;
    end
    value = 0;
    valid = true;
    for r = 0:(p-1)
        dr = differences{r+1};
        if isempty(dr)
            valid = false;
            break;
        end
        value = value+2^(p-1-r)*(abs(dr(1))+abs(dr(end)));
    end
    dp = differences{p+1};
    if isempty(dp) && M-p+1 > 0
        valid = false;
    elseif ~isempty(dp)
        value = value+sum(abs(dp));
    end
    if valid
        Wrow(p) = value;
    end
end
end

function [Bbest,bestBranch] = kernel_scalar_analytic_cell_bound( ...
    coeff,phaseBox,par)
% Take the smallest valid bound among the l1, Der, LCS, and ResLin branches.
candidates = coeff.l1;
labels = "l1";

% Derivative branch.
if par.kernelBound.useDerivativeBranch && coeff.flatEndsR4
    dNLower = kernel_dN_plus_lower_box( ...
        coeff.N,phaseBox.omega1,phaseBox.omega2);
    rhoUpper = kernel_interval_sine_abs_max( ...
        phaseBox.omega2(1),phaseBox.omega2(2));
    for d0 = par.kernelBound.derivativeD0List(:).'
        if dNLower < d0-1e-14
            continue;
        end
        for rhoBar = par.kernelBound.derivativeRhoList(:).'
            if rhoUpper > rhoBar+1e-14
                continue;
            end
            beta = kernel_derivative_class_beta(d0,rhoBar);
            value = sum(beta.*coeff.diffNorms0to4);
            candidates(end+1) = min(coeff.l1,value); %#ok<AGROW>
            labels(end+1) = "Der"; %#ok<AGROW>
        end
    end
end

% LCS branch.
valueLCS = kernel_scalar_lcs_cell_bound(coeff,phaseBox,par);
candidates(end+1) = valueLCS;
labels(end+1) = "LCS";

% Residue-linear branch.
if par.kernelBound.useResidueLinearBranch
    w2 = sort(phaseBox.omega2(:).');
    for q = par.kernelBound.linearQList(:).'
        if q > coeff.N
            continue;
        end
        AList = kernel_linear_A_list(q,phaseBox,par);
        for Ares = AList
            epsBar = kernel_interval_distance_to_lattice_max( ...
                w2(1)-pi*Ares/q,w2(2)-pi*Ares/q,2*pi);
            if ~isfinite(epsBar) || epsBar > pi
                continue;
            end
            nuBar = kernel_interval_distance_to_lattice_max( ...
                q^2*(w2(1)-pi*Ares/q), ...
                q^2*(w2(2)-pi*Ares/q),2*pi);
            if ~isfinite(nuBar) || nuBar > pi
                continue;
            end
            gammaProfile = par.kernelBound.autoGammaSafety* ...
                kernel_gamma_profile_lower( ...
                q,Ares,epsBar,phaseBox);
            if any(gammaProfile <= par.kernelBound.autoGammaMin) || ...
                    any(~isfinite(gammaProfile))
                continue;
            end
            value = kernel_fixed_q_residue_linear_bound( ...
                coeff.a,q,nuBar,gammaProfile);
            candidates(end+1) = value; %#ok<AGROW>
            labels(end+1) = "ResLin"; %#ok<AGROW>
        end
    end
end

[Bbest,index] = min(candidates);
bestBranch = labels(index);
Bbest = max(real(Bbest+floating_pad(Bbest,par.roundoffFactor)),0);
end

function Bbest = kernel_scalar_lcs_cell_bound(coeff,phaseBox,par)
% Minimum of the l1 cap and every requested support-uniform LCS split.
Bbest = coeff.l1;
for iq = 1:numel(coeff.QList)
    BQ = kernel_fixed_q_lcs_cell(coeff,iq,phaseBox,par);
    Bbest = min(Bbest,BQ);
end
safety = par.kernelBound.lcsSafety + ...
    floating_pad(Bbest,par.roundoffFactor);
Bbest = max(real(min(coeff.l1,Bbest)+safety),0);
end

function B = kernel_fixed_q_lcs_cell(coeff,iq,phaseBox,par)
% Hybrid LCS bound: on each lag take the smaller of an interval-cosine
% real-part majorant and the Abel denominator magnitude majorant.
Q = coeff.QList(iq);
N = coeff.N;
if Q < 1 || Q > N
    B = coeff.l1;
    return;
end

w1 = sort(phaseBox.omega1(:).');
w2 = sort(phaseBox.omega2(:).');
Hq = ceil(N/Q)-1;
deltaProfile = zeros(1,max(Hq,0));
for h = 1:Hq
    deltaProfile(h) = par.kernelBound.autoDeltaSafety* ...
        kernel_interval_distance_to_lattice( ...
        h*Q^2*w2(1),h*Q^2*w2(2),pi);
end

statsCell = coeff.Qstats{iq};
Btotal = 0;
for s0 = 0:(Q-1)
    sub = coeff.a((s0+1):Q:N);
    M = numel(sub);
    if M == 0
        continue;
    end
    Fupper = sum(abs(sub).^2);
    stats = statsCell{s0+1};

    for h = 1:(M-1)
        Ucos = 0;
        for m0 = 0:(M-1-h)
            dlag = sub(m0+h+1)*conj(sub(m0+1));
            amplitude = abs(dlag);
            if amplitude == 0
                continue;
            end
            beta1 = h*Q;
            beta2 = 2*Q*s0*h+Q^2*(h^2+2*h*m0);
            J = kernel_affine_phase_interval( ...
                angle(dlag),beta1,w1,beta2,w2);
            Ucos = Ucos+amplitude* ...
                kernel_interval_cosine_max(J(1),J(2));
        end

        Umag = stats.C(h);
        dh = deltaProfile(h);
        if dh > par.kernelBound.autoDeltaMin && ...
                dh <= pi/2 && isfinite(dh)
            denominator = 2*sin(dh);
            for p = par.kernelBound.pList(:).'
                if p <= numel(stats.W)
                    Whp = stats.W{p}(h);
                    if isfinite(Whp)
                        Umag = min(Umag,Whp/denominator^p);
                    end
                end
            end
        end

        Fupper = Fupper+2*min(Ucos,Umag);
    end

    Btotal = Btotal+sqrt(max(real(Fupper),0));
    if Btotal >= coeff.l1
        Btotal = coeff.l1;
        break;
    end
end
B = max(real(min(coeff.l1,Btotal)),0);
end

function v = kernel_diff_norms_0_to_4(a)
a = a(:);
v = zeros(1,5);
v(1) = sum(abs(a));
for order = 1:4
    v(order+1) = sum(abs(diff(a,order,1)));
end
end

function beta = kernel_derivative_class_beta(d0,rhoBar)
s0 = sin(d0/2);
assert(s0 > 0 && isfinite(s0), ...
    'QPAC:InvalidDerivativeD0','Invalid derivative threshold d0.');
assert(rhoBar >= 0 && rhoBar <= 1 && isfinite(rhoBar), ...
    'QPAC:InvalidDerivativeRho','Invalid derivative rho bound.');
beta = zeros(1,5);
beta(1) = 105*rhoBar^4/(16*s0^8);
beta(2) = 105*rhoBar^3/(16*s0^7);
beta(3) = 45*rhoBar^2/(16*s0^6);
beta(4) = 5*rhoBar/(8*s0^5);
beta(5) = 1/(16*s0^4);
end

function dLower = kernel_dN_plus_lower_box(N,w1,w2)
w1 = sort(w1(:).');
w2 = sort(w2(:).');
dLower = inf;
for n0 = 0:(N-1)
    coefficient = 2*n0+1;
    values = [ ...
        w1(1)+coefficient*w2(1), ...
        w1(1)+coefficient*w2(2), ...
        w1(2)+coefficient*w2(1), ...
        w1(2)+coefficient*w2(2)];
    dLower = min(dLower,kernel_interval_distance_to_lattice( ...
        min(values),max(values),2*pi));
    if dLower <= 0
        dLower = 0;
        return;
    end
end
end

function value = kernel_interval_sine_abs_max(lo,hi)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if hi-lo >= pi
    value = 1;
    return;
end
values = [abs(sin(lo)),abs(sin(hi))];
for k = (floor((lo-pi/2)/pi)-1):(ceil((hi+pi/2)/pi)+1)
    t = pi/2+k*pi;
    if t >= lo && t <= hi
        values(end+1) = 1; %#ok<AGROW>
    end
end
value = max(values);
end

function AList = kernel_linear_A_list(q,phaseBox,par)
lo = phaseBox.omega2(1);
hi = phaseBox.omega2(2);
AList = [];
for Ares = 0:(2*q-1)
    distance = kernel_interval_distance_to_lattice( ...
        lo-pi*Ares/q,hi-pi*Ares/q,2*pi);
    if distance <= par.kernelBound.linearAWindow+1e-14
        AList(end+1) = Ares; %#ok<AGROW>
    end
end
if isempty(AList)
    center = 0.5*(lo+hi);
    distances = zeros(1,2*q);
    for Ares = 0:(2*q-1)
        distances(Ares+1) = kernel_distance_to_lattice( ...
            center-pi*Ares/q,2*pi);
    end
    [~,index] = min(distances);
    AList = index-1;
end
AList = unique(AList);
end

function gammaProfile = kernel_gamma_profile_lower( ...
    q,Ares,epsBar,phaseBox)
gammaProfile = zeros(1,q);
w1 = sort(phaseBox.omega1(:).');
for s = 0:(q-1)
    lo = q*w1(1)+pi*Ares*q-2*q*s*epsBar;
    hi = q*w1(2)+pi*Ares*q+2*q*s*epsBar;
    alphaLower = kernel_interval_distance_to_lattice(lo,hi,2*pi);
    gammaProfile(s+1) = max(0,min(1,sin(alphaLower/2)));
end
end

function B = kernel_fixed_q_residue_linear_bound( ...
    a,q,nuBar,gammaProfile)
a = a(:);
N = numel(a);
l1 = sum(abs(a));
nuBar = min(pi,abs(nuBar));
Bq = 0;
for s = 0:(q-1)
    sub = a((s+1):q:N);
    if isempty(sub)
        continue;
    end
    Ls = sum(abs(sub));
    M = numel(sub)-1;
    denominator = 2*gammaProfile(s+1);
    if denominator <= 0 || ~isfinite(denominator) || M == 0
        Bs = Ls;
    else
        variation = abs(sub(1))+abs(sub(end));
        for m = 0:(M-1)
            phaseJump = min(pi,nuBar*(2*m+1));
            variation = variation+kernel_scalar_phase_jump( ...
                sub(m+1),sub(m+2),phaseJump);
        end
        Bs = min(Ls,variation/denominator);
    end
    Bq = Bq+Bs;
    if Bq >= l1
        Bq = l1;
        break;
    end
end
B = max(real(min(l1,Bq)),0);
end

function value = kernel_scalar_phase_jump(a,b,T)
ra = abs(a);
rb = abs(b);
if ra == 0
    value = rb;
    return;
elseif rb == 0
    value = ra;
    return;
end
center = angle(b)-angle(a);
d0 = kernel_distance_to_lattice(center,2*pi);
dmax = min(pi,d0+abs(T));
value = sqrt(max(ra^2+rb^2-2*ra*rb*cos(dmax),0));
end

function J = kernel_affine_phase_interval(offset,beta1,I1,beta2,I2)
values = [ ...
    offset+beta1*I1(1)+beta2*I2(1), ...
    offset+beta1*I1(1)+beta2*I2(2), ...
    offset+beta1*I1(2)+beta2*I2(1), ...
    offset+beta1*I1(2)+beta2*I2(2)];
J = [min(values),max(values)];
end

function cmax = kernel_interval_cosine_max(lo,hi)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if hi-lo >= 2*pi
    cmax = 1;
    return;
end
values = [cos(lo),cos(hi)];
if ceil(lo/(2*pi)) <= floor(hi/(2*pi))
    values(end+1) = 1; %#ok<AGROW>
end
cmax = max(-1,min(1,max(values)));
end

function dmin = kernel_interval_distance_to_lattice(lo,hi,period)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if ceil(lo/period) <= floor(hi/period)
    dmin = 0;
    return;
end
k1 = round(lo/period);
k2 = round(hi/period);
dmin = min(abs(lo-k1*period),abs(hi-k2*period));
dmin = min(dmin,period/2);
end

function dmax = kernel_interval_distance_to_lattice_max(lo,hi,period)
if lo > hi
    tmp = lo; lo = hi; hi = tmp;
end
if hi-lo >= period
    dmax = period/2;
    return;
end
values = [kernel_distance_to_lattice(lo,period), ...
    kernel_distance_to_lattice(hi,period)];
for k = (floor((lo-period/2)/period)-1): ...
        (ceil((hi+period/2)/period)+1)
    t = (k+0.5)*period;
    if t >= lo && t <= hi
        values(end+1) = period/2; %#ok<AGROW>
    end
end
dmax = max(values);
end

function d = kernel_distance_to_lattice(x,period)
d = abs(mod(x+period/2,period)-period/2);
end

function [lo,hi] = sin2_interval(I)
a = min(I); b = max(I);
values = [sin(a)^2,sin(b)^2];
for k = (floor((a-pi/2)/pi)-1):(ceil((b-pi/2)/pi)+1)
    t = pi/2+k*pi;
    if t >= a && t <= b
        values(end+1) = 1; %#ok<AGROW>
    end
end
for k = (floor(a/pi)-1):(ceil(b/pi)+1)
    t = k*pi;
    if t >= a && t <= b
        values(end+1) = 0; %#ok<AGROW>
    end
end
lo = min(values);
hi = max(values);
end

function [lo,hi] = absolute_interval_range(a,b)
hi = max(abs(a),abs(b));
if a <= 0 && b >= 0
    lo = 0;
else
    lo = min(abs(a),abs(b));
end
end

function value = max_abs_affine_in_c2(cLo,cHi,a0,a1)
if cLo <= 0 && cHi >= 0
    tLo = 0;
else
    tLo = min(cLo^2,cHi^2);
end
tHi = max(cLo^2,cHi^2);
value = max(abs(a0+a1*[tLo,tHi]));
end

function value = max_t_times_one_minus_t2(tLo,tHi)
tLo = min(max(tLo,0),1);
tHi = min(max(tHi,tLo),1);
candidates = [tLo,tHi];
stationary = 1/sqrt(3);
if stationary >= tLo && stationary <= tHi
    candidates(end+1) = stationary; %#ok<AGROW>
end
value = max(candidates.*(1-candidates.^2));
end

function value = max_t2_times_one_minus_t2(tLo,tHi)
tLo = min(max(tLo,0),1);
tHi = min(max(tHi,tLo),1);
candidates = [tLo,tHi];
stationary = 1/sqrt(2);
if stationary >= tLo && stationary <= tHi
    candidates(end+1) = stationary; %#ok<AGROW>
end
value = max(candidates.^2.*(1-candidates.^2));
end

function value = max_abs_affine_interval(a,b,I)
I = sort(I(:).');
value = max(abs([a+b*I(1),a+b*I(2)]));
end

function value = max_abs_bilinear_box(x,y,cInterval,zInterval)
cInterval = sort(cInterval(:).');
zInterval = sort(zInterval(:).');
values = [ ...
    cInterval(1)*x+zInterval(1)*y, ...
    cInterval(1)*x+zInterval(2)*y, ...
    cInterval(2)*x+zInterval(1)*y, ...
    cInterval(2)*x+zInterval(2)*y];
value = max(abs(values));
end

%% Curvature and continuum padding
function [tauLo,tauHi] = tau_interval(ranges,thetaDomain,d)
[cLo,cHi] = cos_interval(thetaDomain);
invR = [1/max(ranges),1/min(ranges)];
values = d*[cLo*invR(1),cLo*invR(2), ...
            cHi*invR(1),cHi*invR(2)];
tauLo = min(values);
tauHi = max(values);
end

function [lo,hi] = cos_interval(I)
a = min(I); b = max(I);
values = [cos(a),cos(b)];
for k = (floor(a/pi)-1):(ceil(b/pi)+1)
    t = k*pi;
    if t >= a && t <= b
        values(end+1) = cos(t); %#ok<AGROW>
    end
end
lo = min(values);
hi = max(values);
end

function qMin = quadratic_minimum(G,tauLo,tauHi)
candidates = [tauLo,tauHi];
if G(2,2) > 0
    stationary = -G(1,2)/G(2,2);
    if stationary >= tauLo && stationary <= tauHi
        candidates(end+1) = stationary; %#ok<AGROW>
    end
end
q = G(1,1)+2*G(1,2)*candidates+G(2,2)*candidates.^2;
qMin = min(q);
assert(qMin > 0,'QPAC:NonpositiveQMinimum', ...
    'The tangent quadratic loses positivity on the hull.');
end

function qMax = quadratic_maximum(G,tauLo,tauHi)
q = G(1,1)+2*G(1,2)*[tauLo,tauHi]+ ...
    G(2,2)*[tauLo,tauHi].^2;
qMax = max(q);
end

function caps = global_profile_caps(geom,alpha)
absX = abs(geom.x);
absY = abs(geom.y);
kappa = geom.kappa;

caps.U = kappa*(absX+alpha*absY);
caps.Ut = kappa*(absX+alpha*absY);
caps.Utt = kappa*(absX+4*alpha*absY);
caps.Q2 = caps.Ut+caps.U.^2;
caps.Q3 = caps.Utt+3*caps.U.*caps.Ut+caps.U.^3;
end

function [LF0,LF2] = row_lipschitz_constants( ...
    support,Gamma,geom,caps,par)
LF0 = 0;
LF2 = 0;
for ell = 1:numel(support.range)
    hs = abs(tangent_vector( ...
        support.range(ell),support.theta(ell),geom,par));
    K1 = sum(geom.b.*caps.U);
    H1 = sum(geom.b.*hs.*caps.U);
    K3 = sum(geom.b.*caps.Q3);
    H3 = sum(geom.b.*hs.*caps.Q3);
    LF0 = LF0+Gamma(1)*K1+Gamma(2)*H1;
    LF2 = LF2+Gamma(1)*K3+Gamma(2)*H3;
end
end

function D3Row = direct_third_derivative_row_bound( ...
    rEval,thetaSegment,support,Gamma,geom,par)
%DIRECT_THIRD_DERIVATIVE_ROW_BOUND Uniform bound using only order three.
%
% For the evaluation phase, write
%
%   u   = -kappa*sin(theta)*(x+alpha*cos(theta)*y),
%   u_t = -kappa*(cos(theta)*x+alpha*cos(2*theta)*y),
%   u_tt=  kappa*sin(theta)*(x+4*alpha*cos(theta)*y).
%
% The third-derivative profile is bounded by
%
%   |Q3| <= |u_tt|+3|u||u_t|+|u|^3.
%
% Trigonometric interval extrema below make this bound uniform on the full
% Taylor segment.  No fourth derivative or grid-to-continuum padding enters.
[cLo,cHi] = cos_interval(thetaSegment);
[sLo,sHi] = cos_interval(thetaSegment-pi/2); % sin(theta)=cos(theta-pi/2)
[c2Lo,c2Hi] = cos_interval(2*thetaSegment);

cAbs = max(abs([cLo,cHi]));
sAbs = max(abs([sLo,sHi]));
c2Abs = max(abs([c2Lo,c2Hi]));

alpha = par.d/rEval;
absX = abs(geom.x);
absY = abs(geom.y);
kappa = geom.kappa;

U = kappa*sAbs*(absX+alpha*cAbs*absY);
Ut = kappa*(cAbs*absX+alpha*c2Abs*absY);
Utt = kappa*sAbs*(absX+4*alpha*cAbs*absY);
Q3 = Utt+3*U.*Ut+U.^3;

K3 = sum(geom.b.*Q3);
D3Row = 0;
for ell = 1:numel(support.range)
    hs = abs(tangent_vector( ...
        support.range(ell),support.theta(ell),geom,par));
    H3 = sum(geom.b.*hs.*Q3);
    D3Row = D3Row+Gamma(1)*K3+Gamma(2)*H3;
end
end

function [F0,F2,F3] = evaluate_common_rows( ...
    rEval,thetaEval,support,Gamma,geom,par)
thetaEval = thetaEval(:).';
M = numel(thetaEval);
F0 = zeros(1,M);
F2 = zeros(1,M);
F3 = zeros(1,M);

for first = 1:par.chunkSize:M
    last = min(M,first+par.chunkSize-1);
    q = thetaEval(first:last);
    f0 = zeros(size(q));
    f2 = zeros(size(q));
    f3 = zeros(size(q));
    for ell = 1:numel(support.range)
        C = evaluate_source(support.range(ell),support.theta(ell), ...
            rEval,q,geom,par);
        f0 = f0+Gamma(1)*C(1,:)+Gamma(2)*C(2,:);
        f2 = f2+Gamma(1)*C(5,:)+Gamma(2)*C(7,:);
        f3 = f3+Gamma(1)*C(6,:)+Gamma(2)*C(8,:);
    end
    F0(first:last) = f0;
    F2(first:last) = f2;
    F3(first:last) = f3;
end
end

function d = distance_to_2pi(x)
d = abs(mod(x+pi,2*pi)-pi);
end

function pad = floating_pad(value,factor)
pad = factor*eps(max(1,abs(value)));
end

function [value,label] = select_enabled_bound(values,labels,enabled)
%SELECT_ENABLED_BOUND Return the tightest finite enabled enclosure.
values = values(:).';
enabled = logical(enabled(:).');
assert(numel(values) == numel(labels) && ...
    numel(values) == numel(enabled), ...
    'QPAC:BoundSelectionSize', ...
    'Bound values, labels, and enable flags must have equal length.');
admissible = enabled & isfinite(values) & values >= 0;
assert(any(admissible),'QPAC:NoEnabledBound', ...
    'At least one finite nonnegative enabled bound is required.');
work = values;
work(~admissible) = inf;
[value,index] = min(work);
label = string(labels{index});
end

%% Output helpers
function summary = make_summary_table(result)
quantity = { ...
    'eta_SS';'sigma_lower_squared';'E_curv';'m_near'; ...
    'D2_analytic_upper';'D2_DirSum_upper';'D2_selected_upper'; ...
    'D3_analytic_upper';'D3_direct_Q3_upper'; ...
    'D3_selected_upper';'eta_near_upper'; ...
    'eta_far_analytic_upper';'eta_far_DirSum_upper'; ...
    'eta_far_selected_upper'; ...
    'C_QPAC_upper';'QPAC_slack';'internal_audit_OK'; ...
    'numerical_QPAC_OK'};
value = [result.etaSS;result.sigmaLower2;result.Ecurv;result.mNear; ...
    max(result.kernelD2.analyticBoundByRow); ...
    max(result.kernelD2.dirSumBoundByRow);result.D2; ...
    max(result.kernelD3.boundByRow); ...
    max(result.directD3.boundByRow);result.D3; ...
    result.etaNear;max(result.farBounds.analyticBoundByRow); ...
    max(result.farBounds.dirSumBoundByRow);result.etaFar; ...
    result.CQPAC;result.qpacSlack;double(result.audit.internalAuditOK); ...
    double(result.numericalQPACOK)];
summary = table(quantity,value);
end

function written = write_latex_table(filename,result)
written = false;
fid = fopen(filename,'w');
assert(fid >= 0,'QPAC:CannotWriteTable','Could not create the LaTeX table.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
p = result.parameters;

fprintf(fid,'%% Auto-generated by run_lcs_example_qpac.m\n');
fprintf(fid,'\\begin{tabular}{c|c}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'QPAC quantity & Numerical value or bound\\\\\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'$f_c$ & $%.12g\\,\\mathrm{GHz}$\\\\\n',p.fc/1e9);
fprintf(fid,'$N_r$ & $%d$\\\\\n',p.Nr);
fprintf(fid,'$r_1$ & $%.12g\\,\\mathrm{m}$\\\\\n',p.rangeRows(1));
fprintf(fid,'$r_2$ & $%.12g\\,\\mathrm{m}$\\\\\n',p.rangeRows(2));
fprintf(fid,'$\\delta$ & $%.12g$\\\\\n',p.delta);
fprintf(fid,'$\\eta_{\\rm SS}$ & $%.12g$\\\\\n',result.etaSS);
fprintf(fid,'$m_{\\rm near}$ & $%.12g$\\\\\n',result.mNear);
fprintf(fid,'$D_2^{\\rm analytic}(\\delta)$ & $%.12g$\\\\\n', ...
    max(result.kernelD2.analyticBoundByRow));
fprintf(fid,'$D_2^{\\rm DirSum}(\\delta)$ & $%.12g$\\\\\n', ...
    max(result.kernelD2.dirSumBoundByRow));
fprintf(fid,'$D_2(\\delta)$ (selected) & $%.12g$\\\\\n',result.D2);
fprintf(fid,'$D_3^{\\rm analytic}(\\delta)$ & $%.12g$\\\\\n', ...
    max(result.kernelD3.boundByRow));
fprintf(fid,'$D_3^{\\rm direct}(\\delta)$ & $%.12g$\\\\\n', ...
    max(result.directD3.boundByRow));
fprintf(fid,'$D_3(\\delta)$ (selected) & $%.12g$\\\\\n',result.D3);
fprintf(fid,'$\\eta_{\\rm near}(\\delta)$ & $%.12g$\\\\\n',result.etaNear);
fprintf(fid,'$\\eta_{\\rm far}^{\\rm analytic}(\\delta)$ & $%.12g$\\\\\n', ...
    max(result.farBounds.analyticBoundByRow));
fprintf(fid,'$\\eta_{\\rm far}^{\\rm DirSum}(\\delta)$ & $%.12g$\\\\\n', ...
    max(result.farBounds.dirSumBoundByRow));
fprintf(fid,'$\\eta_{\\rm far}(\\delta)$ (selected) & $%.12g$\\\\\n', ...
    result.etaFar);
fprintf(fid,'$\\mathfrak C_{\\rm QPAC}$ & $%.12g$\\\\\n',result.CQPAC);
fprintf(fid,'QPAC slack $1-\\mathfrak C_{\\rm QPAC}$ & $%.12g$\\\\\n', ...
    result.qpacSlack);
fprintf(fid,'Numerical QPAC check & %d\\\\\n',result.numericalQPACOK);
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}\n');
written = true;
end

function written = write_diagnostic_figure( ...
    figuresDir,nearGrid,farLeftGrid,farRightGrid,result)
written = false;
fig = [];
try
    fig = figure('Visible','off','Color','w', ...
        'Position',[100 100 1200 360], ...
        'PaperPositionMode','auto');
    cleanup = onCleanup(@() close_if_valid(fig)); %#ok<NASGU>

    subplot(1,3,1);
    exactMax = max(result.LCS.exact(:,1:4),[],1);
    boundMax = max(result.LCS.bound(:,1:4),[],1);
    bar([exactMax(:),boundMax(:)]);
    set(gca,'XTick',1:4,'XTickLabel',{'K','H','dK','dH'});
    ylabel('channel magnitude');
    legend({'exact','LCS bound'},'Location','northwest');
    title('Support--support correlation');
    grid on;

    subplot(1,3,2); hold on;
    colors = lines(numel(result.row));
    assert(result.D2 > 0 && result.D3 > 0, ...
        'QPAC:ZeroLocalBound','D2 and D3 must be positive for plotting.');
    localHandles = gobjects(2*numel(result.row),1);
    localLabels = cell(2*numel(result.row),1);
    ih = 0;
    for ir = 1:numel(result.row)
        ih = ih+1;
        localHandles(ih) = plot(nearGrid, ...
            result.row{ir}.F2Near/result.D2, ...
            '-','Color',colors(ir,:),'LineWidth',1.1);
        localLabels{ih} = sprintf('D_2 row %d',ir);
        ih = ih+1;
        localHandles(ih) = plot(nearGrid, ...
            result.row{ir}.F3Near/result.D3, ...
            '--','Color',colors(ir,:),'LineWidth',1.1);
        localLabels{ih} = sprintf('D_3 row %d',ir);
    end
    xlabel('$\theta$ (rad)','Interpreter','latex');
    ylabel('row sum / reported bound');
    title('Common-point local rows');
    legend(localHandles,localLabels,'Location','best');
    grid on;

    subplot(1,3,3); hold on;
    for ir = 1:numel(result.row)
        plot(farLeftGrid,result.row{ir}.F0Left, ...
            '-','Color',colors(ir,:),'LineWidth',1.1);
        plot(farRightGrid,result.row{ir}.F0Right, ...
            '-','Color',colors(ir,:),'LineWidth',1.1);
    end
    axisLimits = xlim;
    plot(axisLimits,[1 1],'k--','LineWidth',1.0);
    xlabel('$\theta$ (rad)','Interpreter','latex');
    ylabel('$\sum_\ell(\Gamma_K|K|+\Gamma_H|H|)$', ...
        'Interpreter','latex');
    title(sprintf('Common-point far row; bound %.4f',result.etaFar));
    grid on;

    pdfName = fullfile(figuresDir,'correlation_qpac_diagnostics.pdf');
    pngName = fullfile(figuresDir,'correlation_qpac_diagnostics.png');
    if exist('exportgraphics','file') == 2 || ...
            exist('exportgraphics','builtin') == 5
        exportgraphics(fig,pdfName,'ContentType','vector');
        exportgraphics(fig,pngName,'Resolution',300);
    else
        print(fig,pdfName,'-dpdf','-painters','-bestfit');
        print(fig,pngName,'-dpng','-r300');
    end
    written = exist(pdfName,'file') == 2 && exist(pngName,'file') == 2;
catch ME
    warning('QPAC:FigureNotWritten', ...
        'Diagnostic figure was not written: %s',ME.message);
end
end

function close_if_valid(fig)
if ~isempty(fig) && isgraphics(fig)
    close(fig);
end
end

function print_summary(result,summary)
p = result.parameters;
fprintf('\nCorrelation-route QPAC example\n');
fprintf('------------------------------\n');
fprintf('N_r                            = %d\n',p.Nr);
fprintf('f_c                            = %.12g GHz\n',p.fc/1e9);
fprintf('lambda                         = %.12g m\n',p.lambda);
fprintf('aperture                       = %.12g m\n',p.aperture);
fprintf('radiating Fresnel interval     = [%.12g, %.12g] m\n', ...
    p.reactiveBoundary,p.fraunhoferBoundary);
fprintf('range rows                     = [%.12g, %.12g] m\n', ...
    p.rangeRows(1),p.rangeRows(2));
fprintf('common support angle           = %.12g rad\n',p.theta0);
fprintf('connected angular domain       = [%.12g, %.12g] rad\n', ...
    p.thetaDomain(1),p.thetaDomain(2));
fprintf('common localization radius     = %.12g rad\n',p.delta);
fprintf('derivative threshold           = %.12g\n',p.derivativeThreshold);
fprintf('min d_N^+                      = %.12g\n', ...
    min(result.derivativeDiagnostic.dNplus));
fprintf('derivative threshold passed    = %d\n', ...
    result.derivativeDiagnostic.thresholdPassed);
fprintf('G_SS =\n');
disp(result.GSS);
fprintf('Gamma = [%.12g, %.12g]^T\n', ...
    result.Gamma(1),result.Gamma(2));
fprintf(['budget selection               = rowwise minimum over enabled ', ...
    'valid bounds\n']);
fprintf('analytic kernel method enabled = %d\n', ...
    p.methods.enableAnalyticKernel);
fprintf('DirSum method enabled          = %d\n', ...
    p.methods.enableDirSum);
fprintf('direct-Q3 method enabled       = %d\n', ...
    p.methods.enableDirectQ3);
fprintf('analytic D2 row bounds         = [%.12g, %.12g]\n', ...
    result.kernelD2.analyticBoundByRow(1), ...
    result.kernelD2.analyticBoundByRow(2));
fprintf('DirSum D2 row bounds           = [%.12g, %.12g]\n', ...
    result.kernelD2.dirSumBoundByRow(1), ...
    result.kernelD2.dirSumBoundByRow(2));
fprintf('selected D2 branches           = [%s, %s]\n', ...
    char(result.kernelD2.selectedBranchByRow(1)), ...
    char(result.kernelD2.selectedBranchByRow(2)));
fprintf('analytic D2 point audits       = [%d, %d]\n', ...
    result.kernelD2.pointAuditOKByRow(1), ...
    result.kernelD2.pointAuditOKByRow(2));
fprintf('analytic D3 row bounds         = [%.12g, %.12g]\n', ...
    result.kernelD3.boundByRow(1),result.kernelD3.boundByRow(2));
fprintf('direct-Q3 D3 row bounds        = [%.12g, %.12g]\n', ...
    result.directD3.boundByRow(1),result.directD3.boundByRow(2));
fprintf('selected D3 row bounds         = [%.12g, %.12g]\n', ...
    result.kernelD3.selectedBoundByRow(1), ...
    result.kernelD3.selectedBoundByRow(2));
fprintf('selected D3 branches           = [%s, %s]\n', ...
    char(result.kernelD3.selectedBranchByRow(1)), ...
    char(result.kernelD3.selectedBranchByRow(2)));
fprintf('max sampled D3 rows            = [%.12g, %.12g]\n', ...
    result.directD3.sampleMaximumByRow(1), ...
    result.directD3.sampleMaximumByRow(2));
fprintf('kernel cell point audits       = [%d, %d]\n', ...
    result.kernelD3.pointAuditOKByRow(1), ...
    result.kernelD3.pointAuditOKByRow(2));
fprintf('analytic far row bounds        = [%.12g, %.12g]\n', ...
    result.farBounds.analyticBoundByRow(1), ...
    result.farBounds.analyticBoundByRow(2));
fprintf('DirSum far row bounds          = [%.12g, %.12g]\n', ...
    result.farBounds.dirSumBoundByRow(1), ...
    result.farBounds.dirSumBoundByRow(2));
fprintf('selected far branches          = [%s, %s]\n', ...
    char(result.farBounds.selectedBranchByRow(1)), ...
    char(result.farBounds.selectedBranchByRow(2)));
fprintf('analytic far point audits      = [%d, %d]\n', ...
    result.farBounds.pointAuditOKByRow(1), ...
    result.farBounds.pointAuditOKByRow(2));
disp(summary);
fprintf('fine/coarse eta_near           = %.12g / %.12g\n', ...
    result.etaNear,result.coarseAudit.etaNear);
fprintf('fine/coarse eta_far            = %.12g / %.12g\n', ...
    result.etaFar,result.coarseAudit.etaFar);
fprintf('QPAC slack (1-C_QPAC)          = %.12g\n',result.qpacSlack);
fprintf('internalAuditOK                = %d\n',result.audit.internalAuditOK);
fprintf('numericalQPACOK                = %d\n',result.numericalQPACOK);
fprintf('results MAT-file                = %s\n', ...
    result.outputStatus.matPath);
fprintf('summary CSV                     = %s\n', ...
    result.outputStatus.csvPath);
fprintf('LaTeX table                     = %s\n', ...
    result.outputStatus.texPath);
fprintf('diagnostic figure written       = %d\n', ...
    result.outputStatus.figureWritten);
if result.numericalQPACOK
    fprintf(['Interpretation: all implemented QPAC inequalities and ', ...
        'internal consistency checks pass.\n']);
else
    fprintf(['Interpretation: at least one QPAC inequality or internal ', ...
        'audit does not pass; inspect the reported quantities.\n']);
end
fprintf(['Numerics: MATLAB double precision; continuum coverage uses ', ...
    'analytical cell bounds.\n\n']);
end
    
