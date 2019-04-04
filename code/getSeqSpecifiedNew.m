% get Seq from data with specified neurons(idx), stimulus type, alignment, binsize,
% and length before and after alignment(TPre,TPost)
% input:
%        data      - n * 1 structure, n is number of neurons
%        idx       - idx of neurons that will be included in the seq
%        stimulus  - stimulus type that will be label/aligned in the seq
%        options   - seq options;
% output:
%        seq       - sequence structure
% hsl 2017-3-13
%function seq = getSeqSpecifiedNew(data, idx, sField, alignment, binsize, TT1, TT2)
function [seq,testseq] = getSeqSpecifiedNew(data, idx, stimulus, teststimulus, options)
stimulusField = ['laserOn', stimulus];
withdrawField = ['withdraw', stimulus];
teststimulusField = ['laserOn', teststimulus];
alignment = options.alignment;
binsize = options.binsize;
TPre = options.TPre;
TPost = options.TPost;

if isfield(data,stimulusField) && eval(['~isempty(data(1).' stimulusField ')'])
    edges = [-TPre:binsize:TPost];
    NoTrials = eval(['length(data(1).' stimulusField ')']);
    seq = struct('T', length(edges), 'y',zeros(length(idx),length(edges)));
    for trial= 1:NoTrials
        Spikecount = zeros(length(idx),length(edges));
        for cc=1:length(idx)
            if (strcmp('stimulus',alignment) || eval(['isempty(data(idx(cc)).' withdrawField ')']) || strcmp(stimulus,'50') || strcmp(stimulus,'VF'))
                trigger = eval(['data(idx(cc)).' stimulusField]); % stimulus
            elseif(strcmp('withdraw',alignment))
                trigger = eval(['data(idx(cc)).' withdrawField]); % withdraw
            end
            ind = find(data(idx(cc)).spikes>=trigger(trial)-TPre & data(idx(cc)).spikes<=trigger(trial)+TPost);
%             ind = find(data(idx(cc)).spikes>=trigger(trial)-TPre-TPost-1 & data(idx(cc)).spikes<=trigger(trial)-1);
            count = histc(data(idx(cc)).spikes(ind)-trigger(trial), edges);
            Spikecount(cc,:) = count;
        end
        seq = setfield(seq,{trial},'T', size(Spikecount,2));
        seq = setfield(seq,{trial},'y', Spikecount);
    end
    if(NoTrials>1)
        testseq = struct();
        for trial= 1:NoTrials         
            if (strcmp('stimulus',alignment) || eval(['isempty(data(idx(1)).' withdrawField ')']) || strcmp(stimulus,'50') || strcmp(stimulus,'VF'))
                trigger = eval(['data(idx(1)).' stimulusField]); % stimulus
                testtrigger = eval(['data(idx(1)).' teststimulusField]);
            elseif(strcmp('withdraw',alignment))
                trigger = eval(['data(idx(1)).' withdrawField]); % withdraw
                testtrigger = eval(['data(idx(1)).' teststimulusField]);
            end
            if trial~=NoTrials
                edges = [binsize:binsize:trigger(trial+1)-trigger(trial)];
                testind= testtrigger>trigger(trial) & testtrigger<trigger(trial+1);
            else
                edges = [binsize:binsize:max(trigger(NoTrials),testtrigger(end))+TPost-trigger(NoTrials)];
                testind= testtrigger>trigger(NoTrials);
            end
            Spikecount = zeros(length(idx),length(edges));
            mark = floor((testtrigger(testind)-trigger(trial)-TPost)/binsize)+1;
            for cc=1:length(idx)
                if trial~=NoTrials
                    ind = find(data(idx(cc)).spikes>trigger(trial)+TPost & data(idx(cc)).spikes<=trigger(trial+1)+TPost);
                else
                    ind = find(data(idx(cc)).spikes>trigger(NoTrials)+TPost & data(idx(cc)).spikes<=max(trigger(NoTrials),testtrigger(end))+TPost);
                end
                count = histc(data(idx(cc)).spikes(ind)-trigger(trial)-TPost, edges);
                Spikecount(cc,:) = count;
            end
            testseq = setfield(testseq,{trial},'T', size(Spikecount,2));
            testseq = setfield(testseq,{trial},'y', Spikecount);
            testseq = setfield(testseq,{trial},'Test', mark);
        end
        
    else
        testseq=[];
    end
else
    disp(['no such field, or empty field' stimulusField]);
    seq = [];
    testseq = [];
end