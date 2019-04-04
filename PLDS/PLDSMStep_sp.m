function [params seq] = PLDSMStep_sp(params,seq)
%
% params = PLDSMStep(params,seq) 
%
% 
% simplified version of PLDSMStep, for a specific case of use.
% 
% Sile Hu 2016-11-21
% remove LDSTransformParams
% Sile Hu 2016-11-22

params = LDSMStepLDS_sp(params,seq);
%params = PLDSMStepObservation_sp(params,seq);
params = PLDSMStepObservation_sp_fixd(params,seq);