function [f, df] = PLDSMStepObservationCost_sp(vecCd,seq,params)
%
% function [f, df] = PLDSMStepObservationCost(vecCd,seq,params)
%
% Mstep for observation parameters C,d for standard PLDS with exp-Poisson observations
%
% Input:
%	- convention Cd = [C d]  and vecCd = vec(Cd)
%
% to do: 
%
%       0) analyze run time
%
%
% (c) L Buesing 2014
%
% simplified version of PLDSMStepObservationCost, for a specific case of use.
% 
% remove useCmask/useS branch,xDim,Trials
% change matrix operation to scalar
% Sile Hu 2016-11-21

yDim    = size(seq(1).y,1);

CdMat   = reshape(vecCd,yDim,2);
C       = CdMat(:,1);
d       = CdMat(:,2);
CC      = zeros(yDim,1);% C^2

for yd=1:yDim
  CC(yd) = C(yd)'*C(yd);
end

f   = 0;				% current value of the cost function = marginal likelihood
df  = zeros(size(C));			% derivative wrt C
dfd = zeros(yDim,1);			% derivative wrt d
 
y    = seq.y;%observation
m    = seq.posterior.xsm; %latent variable
Vsm  = seq.posterior.Vsm';%latent variable covariance

h    = bsxfun(@plus,C*m,d);% each columb:C*m(i)+d
rho  = CC*Vsm;% each columb:CC*Vsm(i)

yhat = exp(h+rho/2);
f    = f+sum(vec(y.*h-yhat));% sum(y*(C*m+d) - exp((C*m+d)+C^2*cov/2))

TT   = yhat*Vsm';%exp((C*m+d)+C^2*cov/2)*cov
TT   = bsxfun(@times,TT,vec(C));% exp((C*m+d)+C^2*cov/2)*cov*C

df   = df  + (y-yhat)*m'-TT;% (y - exp((C*m+d)+C^2*cov/2))*m - exp((C*m+d)+C^2*cov/2)*cov*C
dfd  = dfd + sum((y-yhat),2);% sum over time (y - exp((C*m+d)+C^2*cov/2))

f  = -f;
df = -vec([df dfd]);
