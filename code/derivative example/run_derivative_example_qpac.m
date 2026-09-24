% Self-contained reproducibility code for the derivative QPAC example.
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
% This single MATLAB file contains the complete LCS example implementation:
% configuration, support-uniform LCS bound routines, deterministic sampling
% diagnostics, QPAC constant assembly, and LaTeX table generation.
% No external MATLAB files or input data files are required.
% See the repository README for the complete terms of use.

function out = run_derivative_example_qpac(opt)
%RUN_DERIVATIVE_TWO_RANGE_QPAC Two-range, continuous-window QPAC example.
%
%   out = run_derivative_two_range_qpac();
%   out = run_derivative_two_range_qpac(struct('farPoints',8001));
%
% The support has one source on each of two *distinct* range rows.  Its
% angles vary independently over two closed intervals.  The script bounds:
%   (i)  support interactions using the fourth-order derivative branch;
%   (ii) second derivatives at distinct support points;
%   (iii) complete near sets by common-point row sums plus Lipschitz padding;
%   (iv) complete far sets, including the opposite range row, by common-point
%        row sums plus Lipschitz padding.
%
% The padding is analytical. Floating-point values have an additional small
% numerical reserve; for a machine-verified proof, repeat the calculations
% with directed-rounding interval arithmetic. MATLAB base functions suffice.

if nargin < 1, opt = struct(); end
assert(isstruct(opt) && isscalar(opt),'Options must be a scalar struct.');
def = struct('supportPoints',101,'nearPoints',601, ...
    'farPoints',4001,'batchSize',64,'roundingReserve',1e-9);
names = fieldnames(opt);
for j=1:numel(names)
    assert(isfield(def,names{j}),'Unknown option: %s',names{j});
    def.(names{j}) = opt.(names{j});
end
opt = def;
assert(all([opt.supportPoints,opt.nearPoints,opt.farPoints]>=3) && ...
    all(fix([opt.supportPoints,opt.nearPoints,opt.farPoints]) == ...
        [opt.supportPoints,opt.nearPoints,opt.farPoints]));
assert(opt.batchSize>=1 && opt.batchSize==fix(opt.batchSize));
assert(opt.roundingReserve>=0 && isfinite(opt.roundingReserve));

%% Physical model and non-singleton support class
p = struct();
p.Nr = 128;
p.c0 = 3e8;
p.fc = 10e9;
p.lambda = p.c0/p.fc;
p.d = p.lambda/2;
p.kappa = 2*pi*p.d/p.lambda;
p.ranges = [10,100];
p.alpha = p.d./p.ranges;
p.windows = [-.401,-.399;.399,.401]; % theta = pi/2 + t
p.domain = [-.401,.401];
p.delta = .012;
p.derivativeThreshold = 1.9;
p.curvatureSinThreshold = .002;
p.aperture = (p.Nr-1)*p.d;
p.reactiveBoundary = .62*sqrt(p.aperture^3/p.lambda);
p.fraunhoferBoundary = 2*p.aperture^2/p.lambda;
assert(all(p.ranges>p.reactiveBoundary & ...
    p.ranges<p.fraunhoferBoundary));
assert(p.windows(1,2)<p.windows(2,1));
assert(all(p.windows(:,1)>=p.domain(1) & ...
    p.windows(:,2)<=p.domain(2)));

g = geometry(p);
M = derivative_size_caps(p,g);
source = cell(2,1);
tan = cell(2,1);
grid = cell(2,1);
for j=1:2
    tan{j} = tangent_cell(p.windows(j,:),p.alpha(j),g);
    grid{j} = linspace(p.windows(j,1),p.windows(j,2), ...
        opt.supportPoints).';
    source{j} = source_grid(grid{j},p.alpha(j),g);
end
supportFill = max(diff(p.windows,1,2))/(opt.supportPoints-1)/2;

%% Support-support derivative branch, in both ordered directions
u = zeros(1,4); % K,H,dK,dH
dLower = inf;
rhoUpper = 0;
for evalRow=1:2
    srcRow = 3-evalRow;
    phase = phase_box(p.windows(srcRow,:),p.alpha(srcRow), ...
        p.windows(evalRow,:),p.alpha(evalRow),g);
    dLower = min(dLower,phase.dLower);
    rhoUpper = max(rhoUpper,phase.sinCurvatureUpper);
    assert(phase.dLower>0,'Derivative branch fails on a support cell.');

    hS=tan{srcRow}.h; hE=tan{evalRow}.h;
    C_H=bsxfun(@times,g.b,hS);
    C_dK=bsxfun(@times,g.b,hE);
    C_dH=zeros(g.N,4);
    col=0;
    for ie=1:2
        for is=1:2
            col=col+1;
            C_dH(:,col)=g.b.*hE(:,ie).*hS(:,is);
        end
    end
    here=[derivative_bound(g.b,phase), ...
        tan{srcRow}.arc*max(derivative_bound(C_H,phase)), ...
        tan{evalRow}.arc*max(derivative_bound(C_dK,phase)), ...
        tan{srcRow}.arc*tan{evalRow}.arc* ...
        max(derivative_bound(C_dH,phase))];
    u=max(u,here);
end
u=up(u,opt); dLower=down(dLower,opt);
rhoUpper=up(rhoUpper,opt);
assert(dLower>p.derivativeThreshold && ...
    rhoUpper<p.curvatureSinThreshold, ...
    'Prescribed derivative support-class inequalities fail.');
G=[u(1),u(2);u(3),u(4)];
etaSS=up(max(abs(eig(G))),opt); % L-1=1
assert(etaSS<1,'Support interpolation is not certified.');
Gamma=(eye(2)-G)\[1;0];
Gamma=up(Gamma,opt);
Xi=max(0,Gamma-[1;0]);

%% Support curvature: exact sampled derivatives plus analytic padding
u2SS=zeros(1,2);
for evalRow=1:2
    srcRow=3-evalRow;
    q=grid{evalRow};
    best=zeros(1,2);
    for j=1:opt.batchSize:numel(q)
        idx=j:min(j+opt.batchSize-1,numel(q));
        [K2,H2]=pair_matrices(q(idx),p.alpha(evalRow), ...
            source{srcRow},2,g);
        best=max(best,[max(abs(K2(:))),max(abs(H2(:)))]);
    end
    % One evaluation grid and one independent source grid cover each cell.
    best=best+[M.M3*supportFill+M.M2*M.M1*supportFill, ...
               M.M3*supportFill+M.M2*M.hPrime*supportFill];
    u2SS=max(u2SS,best);
end
u2SS=up(u2SS,opt);
UK2self=up(g.kappa^2*M.qMax,opt);
UH2self=up(M.M2,opt);
Ecurv=up(2*(Xi(1)*UK2self+Xi(2)*UH2self)+ ...
    2*(Gamma(1)*u2SS(1)+Gamma(2)*u2SS(2)),opt);
mNear=down(2*M.sigmaMin2-Ecurv,opt);
assert(mNear>0,'The local curvature margin is not positive.');

%% Near: cover every Taylor segment for both continuous support windows
near2=aggregate_grid(p,g,source,Gamma,M,opt,2);
near3=aggregate_grid(p,g,source,Gamma,M,opt,3);
D2=near2.upper;
D3=near3.upper;
etaNear=up(2*p.delta*D3/(3*mNear)+ ...
    p.delta^2*D2^2/(2*mNear),opt);

%% Far: q can lie on either row; only same-row source is excluded locally
far=aggregate_grid(p,g,source,Gamma,M,opt,0);
etaFar=far.upper;
CQ=up(max([etaSS,etaNear,etaFar]),opt);

out=struct('parameters',p,'options',opt,'GSS',G,'uSS',u, ...
    'derivativeLower',dLower,'sinCurvatureUpper',rhoUpper, ...
    'etaSS',etaSS,'Gamma',Gamma,'Xi',Xi,'u2SS',u2SS, ...
    'UK2self',UK2self,'UH2self',UH2self,'Ecurv',Ecurv, ...
    'mNear',mNear,'D2',D2,'D3',D3,'near2',near2, ...
    'near3',near3,'etaNear',etaNear,'far',far, ...
    'etaFar',etaFar,'CQPAC',CQ,'QPACpass',CQ<1, ...
    'derivativeCaps',M);

fprintf('\nTwo-range derivative-route QPAC example\n');
fprintf('N_r = %d, f_c = %.5g GHz, d = %.5g m\n', ...
    p.Nr,p.fc/1e9,p.d);
fprintf('Aperture = %.9g m; radiating Fresnel interval = [%.9g, %.9g] m\n', ...
    p.aperture,p.reactiveBoundary,p.fraunhoferBoundary);
fprintf('Range rows = [%.9g, %.9g] m\n',p.ranges);
fprintf('Angular domain = pi/2 + [%.9g, %.9g] rad\n',p.domain);
fprintf('Support windows = pi/2 + [%.9g, %.9g], [%.9g, %.9g] rad\n', ...
    p.windows(1,:),p.windows(2,:));
fprintf('Localization radius = %.9g rad\n',p.delta);
fprintf('d_N^+ lower       = %.12g (threshold %.6g)\n', ...
    dLower,p.derivativeThreshold);
fprintf('|sin omega_2| max = %.12g (threshold %.6g)\n', ...
    rhoUpper,p.curvatureSinThreshold);
fprintf('uSS [K H dK dH]  = [%.12g %.12g %.12g %.12g]\n',u);
fprintf('eta_SS            = %.12g\n',etaSS);
fprintf('m_near lower      = %.12g\n',mNear);
fprintf('D2 upper          = %.12g (grid %.12g, padding %.12g)\n', ...
    D2,near2.gridMax,near2.padding);
fprintf('D3 upper          = %.12g (grid %.12g, padding %.12g)\n', ...
    D3,near3.gridMax,near3.padding);
fprintf('eta_near          = %.12g\n',etaNear);
fprintf('eta_far           = %.12g (grid %.12g, padding %.12g)\n', ...
    etaFar,far.gridMax,far.padding);
fprintf('C_QPAC            = %.12g; strict margin = %.12g\n', ...
    CQ,1-CQ);
fprintf('QPAC pass         = %d\n',out.QPACpass);
assert(out.QPACpass,'The QPAC recovery inequalities did not all pass.');
end

function g=geometry(p)
n=(0:p.Nr-1).';
rho=choose4(n).*choose4(p.Nr-1-n);
b=rho/sum(rho);
g=struct('N',p.Nr,'n',n,'b',b,'kappa',p.kappa, ...
    'nBar',sum(b.*n),'n2Bar',sum(b.*n.^2));
g.x=n-g.nBar;
g.y=n.^2-g.n2Bar;
g.A=sum(b.*g.x.^2);
g.B=sum(b.*g.x.*g.y);
g.C=sum(b.*g.y.^2);
assert(g.C>0 && all(rho(1:4)==0) && all(rho(end-3:end)==0));
end

function v=choose4(z)
v=zeros(size(z)); mask=z>=4; z=z(mask);
v(mask)=z.*(z-1).*(z-2).*(z-3)/24;
end

function q=qpoly(t,g)
q=g.A+2*g.B*t+g.C*t.^2;
end

function M=derivative_size_caps(p,g)
% Uniform Cauchy bounds for the first through fourth angular derivatives.
% u, u', u'', u''' are the four gauged phase derivatives.
C=max(abs(sin(p.domain)));
alpha=max(p.alpha);
tau=alpha*C;
vertex=min(tau,max(-tau,-g.B/g.C));
qMin=min(qpoly([-tau,tau,vertex],g));
qMax=max(qpoly([-tau,tau],g));
sMin=cos(max(abs(p.domain)));
sigmaMin2=g.kappa^2*sMin^2*qMin;
assert(qMin>0 && sigmaMin2>0);
M1=g.kappa*sqrt(qMax);
U=g.kappa*(abs(g.x)+tau*abs(g.y));
U1=g.kappa*(C*abs(g.x)+alpha*abs(g.y));
U2=g.kappa*(abs(g.x)+4*tau*abs(g.y));
U3=g.kappa*(C*abs(g.x)+4*alpha*abs(g.y));
M2=sqrt(sum(g.b.*(U1.^2+U.^4)));
M3=sqrt(sum(g.b.*((U2+U.^3).^2+9*U.^2.*U1.^2)));
M4=sqrt(sum(g.b.*((U.^4+4*U.*U2+3*U1.^2).^2+ ...
    (U3+6*U.^2.*U1).^2)));
M=struct('qMin',qMin,'qMax',qMax,'sigmaMin2',sigmaMin2, ...
    'M1',M1,'M2',M2,'M3',M3,'M4',M4, ...
    'hPrime',2*M2/sqrt(sigmaMin2));
end

function T=tangent_cell(w,alpha,g)
tau=-alpha*sin(fliplr(w));
H=zeros(g.N,2);
for j=1:2
    H(:,j)=(g.x+tau(j)*g.y)/sqrt(qpoly(tau(j),g));
end
c=sum(g.b.*H(:,1).*H(:,2));
assert(c>-1,'Tangent endpoints are antipodal.');
T=struct('h',H,'arc',sqrt(2/(1+min(1,c))));
end

function P=phase_box(ws,as,we,ae,g)
% Ordered pair: source (ws,as), then evaluation (we,ae).
cs=-sin(fliplr(ws)); ce=-sin(fliplr(we));
s2s=sin2_range(ws); s2e=sin2_range(we);
w1=g.kappa*[cs(1)-ce(2),cs(2)-ce(1)];
w2=g.kappa/2*[ae*s2e(1)-as*s2s(2), ...
               ae*s2e(2)-as*s2s(1)];
dn=inf;
for n=0:g.N-1
    dn=min(dn,lattice_lower(w1+(2*n+1)*w2));
end
rho=min(1,max(abs(w2))); % |sin(omega_2)| <= min(1,|omega_2|)
s=sin(dn/2);
if s<=0, beta=inf(1,5);
else
    beta=[105*rho^4/(16*s^8),105*rho^3/(16*s^7), ...
          45*rho^2/(16*s^6),5*rho/(8*s^5),1/(16*s^4)];
end
P=struct('dLower',dn,'sinCurvatureUpper',rho,'beta',beta, ...
    'omega1Box',w1,'omega2Box',w2);
end

function r=sin2_range(w)
v=cos(w).^2;
r=[min(v),max(v)];
if w(1)<=0 && w(2)>=0, r(2)=1; end
end

function d=lattice_lower(a)
if ceil(a(1)/(2*pi))<=floor(a(2)/(2*pi)), d=0;
else
    d=min(abs(mod(a+pi,2*pi)-pi));
end
end

function B=derivative_bound(a,P)
assert(all(all(a(1:4,:)==0)) && all(all(a(end-3:end,:)==0)));
if isvector(a), a=a(:); end
B=zeros(1,size(a,2));
for j=1:5
    B=B+P.beta(j)*sum(abs(a),1);
    if j<5, a=diff(a,1,1); end
end
end

function S=source_grid(t,alpha,g)
[psi,h]=gauged_atom(t,alpha,g);
S=struct();
S.t=t;
S.Wpsi=bsxfun(@times,psi,g.b.').';
S.Wh=bsxfun(@times,h,g.b.').';
end

function [psi,h,u,ut,utt]=gauged_atom(t,alpha,g)
t=t(:); theta=pi/2+t;
co=cos(theta); si=sin(theta);
n=g.n.';x=g.x.';y=g.y.';
phase=g.kappa*(co*n-alpha/2*(si.^2)*(n.^2));
chi=g.kappa*(g.nBar*co+alpha/4*g.n2Bar*cos(2*theta));
psi=exp(1i*bsxfun(@minus,phase,chi));
v=bsxfun(@plus,x,alpha*co*y);
tau=alpha*co;
h=-1i*bsxfun(@rdivide,v,sqrt(qpoly(tau,g))).*psi;
u=-g.kappa*bsxfun(@times,si,v);
ut=-g.kappa*(co*x+alpha*cos(2*theta)*y);
utt=g.kappa*bsxfun(@times,si,bsxfun(@plus,x,4*alpha*co*y));
end

function [K,H]=pair_matrices(q,alphaEval,S,order,g)
[psi,~,u,ut,utt]=gauged_atom(q,alphaEval,g);
E=conj(psi);
switch order
    case 0
        E=conj(psi);
    case 2
        E=E.*(-1i*ut-u.^2);
    case 3
        E=E.*(-1i*utt-3*u.*ut+1i*u.^3);
    otherwise
        error('Only orders 0, 2, and 3 are used.');
end
K=E*S.Wpsi;
H=E*S.Wh;
end

function out=aggregate_grid(p,g,source,Gamma,M,opt,order)
% Continuum majorant from a common-evaluation tensor grid. The near tube
% deliberately includes every admissible source-dependent Taylor segment.
if order==0, count=opt.farPoints;
else, count=opt.nearPoints;
end
srcFill=max(diff(p.windows,1,2))/(opt.supportPoints-1)/2;
gridMax=-inf;
qFill=0;
for row=1:2
    if order==0, interval=p.domain;
    else
        interval=[max(p.domain(1),p.windows(row,1)-p.delta), ...
                  min(p.domain(2),p.windows(row,2)+p.delta)];
    end
    q=linspace(interval(1),interval(2),count).';
    thisFill=max(diff(q))/2;
    qFill=max(qFill,thisFill);
    for j=1:opt.batchSize:count
        ix=j:min(j+opt.batchSize-1,count);
        qChunk=q(ix);
        total=zeros(numel(qChunk),1);
        for srcRow=1:2
            [K,H]=pair_matrices(qChunk,p.alpha(row), ...
                source{srcRow},order,g);
            V=Gamma(1)*abs(K)+Gamma(2)*abs(H);
            if order==0 && srcRow==row
                % A far point on this row is outside the delta-neighborhood
                % of its own source. The nearest-grid triple remains in the
                % relaxed constraint below.
                relaxed=p.delta-thisFill-srcFill;
                assert(relaxed>0,'Increase the far/source grid size.');
                mask=abs(bsxfun(@minus,qChunk,source{srcRow}.t.')) ...
                    >= relaxed;
                V(~mask)=-inf;
            end
            total=total+max(V,[],2);
        end
        gridMax=max(gridMax,max(total));
    end
end
assert(isfinite(gridMax),'No feasible far grid triple was found.');
switch order
    case 0
        evalLip=M.M1;
        srcLip=Gamma(1)*M.M1+Gamma(2)*M.hPrime;
    case 2
        evalLip=M.M3;
        srcLip=M.M2*(Gamma(1)*M.M1+Gamma(2)*M.hPrime);
    case 3
        evalLip=M.M4;
        srcLip=M.M3*(Gamma(1)*M.M1+Gamma(2)*M.hPrime);
end
padding=2*sum(Gamma)*evalLip*qFill+2*srcLip*srcFill;
out=struct('gridMax',gridMax,'padding',padding, ...
    'upper',up(gridMax+padding,opt),'evaluationFill',qFill, ...
    'sourceFill',srcFill,'order',order);
end

function z=up(z,opt)
z=z+opt.roundingReserve*max(1,abs(z));
end

function z=down(z,opt)
z=z-opt.roundingReserve*max(1,abs(z));
end
