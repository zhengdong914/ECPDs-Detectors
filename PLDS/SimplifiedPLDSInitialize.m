function params = SimplifiedPLDSInitialize(seq,xDim,params)
%
% simplified version of PLDSInitialize() function
% Initializing PLDS parameters using Nuclear Norm Minimization
% Sile Hu 2016

% set initial params
yDim       = size(seq(1).y,1);
Trials     = numel(seq);
params     = PLDSsetDefaultParameters(params,xDim,yDim);					% set standard parameter values
% **this function will be concerned later, when I can figure out which part
% of the parameters are necessary or not needed for simplified version


% -----Nuclear Norm Minimization initialization--------
disp('Initializing PLDS parameters using Nuclear Norm Minimization')
        
dt = params.opts.algorithmic.NucNormMin.dt; 
options = params.opts.algorithmic.NucNormMin.options;
%  get param from default param, options are for MODnucnrmminWithd

seqRebin.y = [seq.y]; seqRebin = rebinRaster(seqRebin,dt);
%the result of rebin is sum every 10 samples to 1 bin, so the bins number
%become T/dt

Y  = [seqRebin.y];
options.lambda = options.lambda*sqrt(size(Y,1)*size(Y,2));
%  scale lamda by size of Y

Yext = [];
[Y,Xu,Xs,Xv,d] = MODnucnrmminWithd( Y, options , 'Yext', Yext );
%** to be considered later, need to understand math equations

params.model.d = d-log(dt);%** why?

params.model.C = Xu(:,1:xDim)*Xs(1:xDim,1:xDim);
u = [];
params.model = LDSObservedEstimation(Xv(:,1:xDim)',params.model,dt,u);
params.model.Xpca = Xv(:,1)';
params.model.Xs   = Xs(1,1);
params.modelInit = params.model;
