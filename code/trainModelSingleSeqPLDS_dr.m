% train single seq model
% author: Sile Hu, yuehusile@gmail.com
% date: 2017-3-11
% 2017-7-31 Sile Hu, ignore all zero unit for training, then set the C and
% d of that unit with arbitrary value(zero_unit_C=min(C)/100,zero_unit_d=average(d))
function model = trainModelSingleSeqPLDS_dr(seqACC,seqS1,options)

seq = seqACC;
for i=1:length(seqACC)
    seq(i).y = [seq(i).y;seqS1(i).y];
    seq(i).nACC=size(seqACC(i).y,1);
    seq(i).nS1=size(seqS1(i).y,1);
end

for trial = 1:length(seq)  %  length(seqACC) must equal to length(seqS1)
    
    %ignore all zero units
    noSpike = [];
    for i=1:size(seq(trial).y,1)
        if sum(seq(trial).y(i,:))<1
            noSpike = [noSpike i];
        end
    end
    
    noSpike_1 = length(find(noSpike<=seq(trial).nACC));
    noSpike_2 = length(noSpike)-noSpike_1;
    if (noSpike_1==seq(trial).nACC || noSpike_2==seq(trial).nS1)
        model(trial).params = [];
        model(trial).newseq.posterior.xsm = zeros(options.xdim,ceil((options.TPost+options.TPre)/options.binsize)+1);
        model(trial).newseq.posterior.Vsm = zeros(options.xdim * (ceil((options.TPost+options.TPre)/options.binsize)+1),options.xdim);
        model(trial).varBound = [];
        model(trial).infTime = [];
        continue;
    end
        
    seq(trial).y(noSpike,:) = [];
    % initialize
    params = [];

    paramsInit = PLDSInitialize(seq(trial),options.xdim,options.initMethod,options.doubleRegion,params); % xDim = 1
    
    % EM iteration - model training
    paramsInit.opts.algorithmic.EMIterations.maxIter     = 50;
    paramsInit.opts.algorithmic.EMIterations.maxCPUTime  = 30; %unit: s    
    [params newseq varBound infTime] = PLDS_EM(paramsInit,seq(trial),options.binsize);
    
    % set the C and d of that unit with arbitrary value
    noSpikeC = min(abs(params.model.C))/100;
    noSpiked = mean(params.model.d);
    C = params.model.C;
    d = params.model.d;
    for i=1:length(noSpike)
        C = [C(1:noSpike(i)-1,:);noSpikeC;C(noSpike(i):end,:)];
        d = [d(1:noSpike(i)-1);noSpiked;d(noSpike(i):end)];
    end
    params.model.C = C;
    params.model.d = d;
    
    model(trial).params = params;
    model(trial).newseq = newseq;
    model(trial).varBound = varBound;
    model(trial).infTime = infTime;
end