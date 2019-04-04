function Mu = getPriorMeanLDS_sp(params,T,seq)
% 
% simplified version of getPriorMeanLDS, for a specific case of use.
%
% remove unuseful branches of useB and ~isempty(seq)
% remove useless variable xAdd
% replace varargin input para with simple seq
% Sile Hu 2016-11-14
%
A    = params.model.A;
x0   = params.model.x0;

xDim = size(params.model.A,1);

Mu = zeros(xDim,T);
Mu(:,1) = x0;
% forward filtering to get initxx
for t=2:T
  Mu(:,t) = A*Mu(:,t-1);
end