function [seq varBound] = PLDSLaplaceInference_sp(params,seq)
%
% function seq = PLDSlpinf(params,seq)
%
%
% simplified version of PLDSLaplaceInference, for a specific case of use.
%
% call simplified functions
% use one line code instead of call PoissonBaseMeasure
% remove unuseful branch useS
% Sile Hu 2016-11-14
% changes for call PLDSLaplaceInferenceCore_sp instead of PLDSLaplaceInferenceCore
% remove Tmax
% Sile Hu 2016-11-14
computeVarBound = true; %!!! make this nice

Trials      = numel(seq);
[yDim xDim] = size(params.model.C);
varBound    = 0;

mps = params.model;
mps.E = [];
mps.D = [];
mps.initV = mps.Q0;
mps.initx = mps.x0;

runinfo.nStateDim = xDim;
runinfo.nObsDim   = yDim;

infparams.initparams = mps;
infparams.runinfo    = runinfo;
infparams.notes      = params.model.notes;

%-- indat is the data sequence for PLDSLaplaceInferenceCore_sp input --

T = size(seq.y,2);
indat = seq;
%should set other data structure independent flag for this
if isfield(seq,'posterior') && isfield(seq.posterior,'xsm')
    indat.xxInit = seq.posterior.xsm;
else
    indat.xxInit = getPriorMeanLDS_sp(params,T,seq);%seq.T=201;seq.y = 12*201 double
end

seqRet = PLDSLaplaceInferenceCore_sp(indat, infparams);%indat.T=201;indat.y = 12*201 double;indat.xxInit=1*201 double

seq.posterior.xsm      = seqRet.x; % PLDS_filter_smoother smooth_state
seq.posterior.Vsm      = seqRet.V; % PLDS_filter_smoother stateCov_smooth
seq.posterior.VVsm     = seqRet.VV; % PLDS_filter_smoother sqrtm( stateCov_smooth(:,:,t) * stateCov_smooth(:,:,t+1) );
seq.posterior.lamOpt   = exp(vec(seqRet.Ypred)); % compute predicted lambda; PLDS_filter_smoother didn't record this

%---------------------------var inference---------------------------------
computeVarBound = 0;
if computeVarBound %PLDS_filter_smoother varBound
    
    VarInfparams    = params.model;
    VarInfparams.CC = zeros(xDim,xDim,yDim);
    for yy=1:yDim%yDim=12
        VarInfparams.CC(:,:,yy) = params.model.C(yy,:)'*params.model.C(yy,:);%CC =1*1*12 double;params.model.C=12*1 double
    end
    VarInfparams.CC = reshape(VarInfparams.CC,xDim^2,yDim);
    
    Cl = {}; for t=1:T; Cl = {Cl{:} params.model.C}; end
    Wmax = sparse(blkdiag(Cl{:}));
    
    % iterate over trials
    
    %-- VarInfparams is the data sequence for VariationalInferenceDualCost_sp input --
    T = size(seq.y,2);
    VarInfparams.d = repmat(params.model.d,T,1);
    
    VarInfparams.mu         = vec(getPriorMeanLDS_sp(params,T,seq));
    VarInfparams.W          = Wmax(1:yDim*T,1:xDim*T);
    VarInfparams.y          = seq.y;
    VarInfparams.Lambda     = buildPriorPrecisionMatrixFromLDS(params,T);%--- need to be understand
    VarInfparams.WlamW      = sparse(zeros(xDim*T));
    VarInfparams.dualParams = [];
    
    if isfield(params.model,'baseMeasureHandle')
        %PoissonBaseMeasure, log y!
        VarInfparams.DataBaseMeasure = -sum(log(gamma(vec(seq.y)+1)));
        seq.posterior.DataBaseMeasure = VarInfparams.DataBaseMeasure;
    end
    
    lamOpt = seq.posterior.lamOpt;
    
    [DualCost, ~, varBound] = VariationalInferenceDualCost_sp(lamOpt,VarInfparams);  
    seq.posterior.varBound = varBound;
    
%     varBound = 0;
%     varBound = varBound + seq.posterior.varBound;
   
end

binsize = 0.05;
y = seq.y; 
A = params.model.A; % m-by-m
C = params.model.C; % C-by-m
d = params.model.d; % C-by-1
Q = params.model.Q; % m-by-m
smooth_state = seqRet.x;
stateCov_smooth = seqRet.V;
% correlation estimation (Smith & Brown, 2003)
QQ1 = zeros(1,T);
QQ2 = zeros(1,T);

for t=1:T
    QQ1(t) = stateCov_smooth(t) + smooth_state(t) * smooth_state(t)';% Wk
end
for t=1:T-1
    tem = A * stateCov_smooth(t+1);
    QQ2(t) = tem + smooth_state(t) * smooth_state(t+1)';% Wk,k+1
end

% compute the Q function
Qfun = 0;
% compute E[log p(N0,t|x,theta*)]
for t=1:T
    rate = exp(C * smooth_state(:,t) + d);
    tem1 = sum(y(:,t) .* (C * smooth_state(:,t) + d + log(binsize)) - rate * binsize);   % sum wrt channel
    Qfun = Qfun + tem1; % sum wrt time
end    


%NOTE:  a'*b*a = trace(a*a'*b) 

% noise covariance Q has to be full rank
% compute E[log p(x|Rho,alpha,sigma^2)] neglecting last two terms
for t=1:T-1
    tem2 = -0.5 * trace( (QQ1(t+1) - 2*A*QQ2(t) + A*QQ1(t)*A') * inv(Q));  
    Qfun = Qfun + tem2 - 0.5*log(det(Q)) -0.5*log(2*pi);
end
seq.posterior.varBound = Qfun;
varBound = Qfun;
