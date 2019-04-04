function params = PLDSMStepObservation_sp_fixd(params,seq)
%
% function params = PLDSMStepObservation(params,seq)
%
%
% simplified version of PLDSMStepObservation, for a specific case of use.
% 
% remove useCmask branch
% Sile Hu 2016-11-21
%

minFuncOptions = params.opts.algorithmic.MStepObservation.minFuncOptions;

[yDim xDim] = size(params.model.C);

%CdInit = vec([params.model.C params.model.d]); % warm start at current parameter values -- what is warm start mean?
CdInit = vec([params.model.C]); % warm start at current parameter values -- what is warm start mean?

%MStepCostHandle = @PLDSMStepObservationCost_sp;
MStepCostHandle = @PLDSMStepObservationCost_sp_fixd;

%%% optimization %%%
% for test
%CdInit = [-0.471283;0.1619;0.0902355;0.362439;0.22494;0.746859;0.0614833;-0.651503;-1.40789;-0.766934;-0.748916;-0.909273;-0.921149;-0.595436];

CdOpt = minFunc_sp(MStepCostHandle,CdInit,minFuncOptions,seq,params);
%CdOpt = minFunc(MStepCostHandle,CdInit,minFuncOptions,seq,params);
CdOpt = reshape(CdOpt,yDim,xDim);

params.model.C = CdOpt(:,1:xDim);
%params.model.d = CdOpt(:,end);


