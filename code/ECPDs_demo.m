% dual recording decoding demo for new Plexon data
% trial-by-trial structure
% Zhengdong Xiao 2017-3-11
% update : Zhengdong Xiao 2017-7-25

% NOTE: before running, please make sure the paths are correct 
% when set loadModel=1, please make sure the model of current data exist
% in the model directory
%% --- add path & clear space ---
% PLDS code path
plds_core_path = '../PLDS';
addpath(genpath(plds_core_path));
% LDS code path
lds_core_path = '../lds';
addpath(genpath(lds_core_path));
% minFunc path
minFunc_path = '../minFunc'; %
addpath(genpath(minFunc_path));
data_path = './';
addpath(genpath(data_path));
draw_path = '../draw_figure/'; %path of draw figure functions
addpath(draw_path);
clear;clc;
close all;
%% -- run options and load data --
% 0 - disable
% 1 - enable
saveModel = 1;
% directory for saving trained models
model_dir = '../model';
loadModel = 0;% used trained model when set it as 1; training new model when set it as 0
train_lds = 0;% not using lds model 
dataOpt = 'data';
load([dataOpt '.mat']);
%% --- set options ---
options.binsize = 0.05; %unit: s
options.TPre = 5;       %unit: s, before stimulus
options.TPost = 5;     %unit: s, after stimulus
options.alignment = 'stimulus';
options.initMethod = 'ExpFamPCA';%'NucNormMin';%'ExpFamPCA'; % initial all parameters of PLDS model(C,A,Q...); don't change this
options.xdim = 1;% input dim is 1, first order model , don't change this 
options.doubleRegion = 0; % for AR-2 model, don't change 
options.sim = 0;  % simulation for AR-2 model,dont' change
options.PP = 0;% if stimulus is pin prick, change it to 1
%% --- get trial sequence --- 
% get all the unit raster in every trial during 10s, seqACC250
if isfield(data,'allUnits')
    if isfield(data.allUnits, 'spike')
        [seq,testseq] = getRecSeq(data.allUnits.spike,options);
    else
        [seq,testseq] = getRecSeq(data.allUnits,options);
    end
else
    [seq,testseq] = getRecSeq(data.allUnitsAdd,options);
end
%% --- training phase ---
idx_minus = find(dataOpt=='-');
dataOpt(idx_minus) = '_';
if loadModel
    filename = [model_dir '/' dataOpt '.mat'];  
    load(filename);
    model = eval([dataOpt '.model']);
    options = eval([dataOpt '.options']);
    if train_lds
        model_lds = eval([dataOpt '.model_lds']);
    end
else
    model = trainModelPLDS_dr(seq,options);
    if train_lds
        model_lds = trainModelLDS(seq,options);
    end
    if saveModel
        filename = [model_dir '/' dataOpt '.mat'];  
        eval([dataOpt '.model=model;']);
        eval([dataOpt '.options=options;']);
        if train_lds
            eval([dataOpt '.model_lds=model_lds;']);
        end
        save(filename, dataOpt);
    end
end
if train_lds
    model.lds = model_lds;
end
options.region = seq.region;
options.stimulus = seq.stimulus;
for i=1:length(options.stimulus)
    eval(['options.ntrial' seq.stimulus{i} '=seq.ntrial' seq.stimulus{i}]);
end
%% testing phase and evaluation
options.n_premodels = 3;% using previous three models to predict current trial; % odd value for majority vote;
options.n_ccmodels = 3; % using previous three ccf models to predict current trial; 
onlineDecodeResult = cell(2,options.n_premodels);
for n_decode=1:2 %  loop 1 is for detecting positive stimulus trials ;loop 2 is for detecting negative stimulus trials
    %
    % --- online decode ---
    % ***********************
    options.trial250 = -1;
    options.trial50 = -1;
    seq.stimulus = options.stimulus;
    % ***********************
    options.baseline = [1 4];% [start end] unit:s, 0 is the start of seq; model free need baseline to compute
    %options.currentTrial = 1; % decode training seq, should ignore first data point
    options.currentTrial = 0;% disable current model predict current trial.
    options.preTrial=1;% using previous model to predict current trial
    if(n_decode==2)% using positive model(250mW) to predict negative trials (50mw)
        options.currentTrial = 0;
        options.preTrial= 0;
    end
    % 'PLDS_FB' - forward backward PLDS, the result get while training
    % 'PLDS_forward' - forward filter PLDS, for online decoding
    % 'LDS'          - transfered LDS filter, for online decoding
    % 'ModelFree'    - model free decoding method
    % 'PLDS_PF'      - particle filter decoding method for PLDS online decoding
    % 'GT'           - ground truth of latent variable when using simulated data
    % 'PLDS_PF_SM2'   - particle filter PLDS smoothing, based on PLDS_PF v2
    options.decoders = {'PLDS_forward'};
    % options.decoders = {'GT', 'PLDS_forward','PLDS_PF'};

    % PF params, for silas paper
    % best for sim2_4 trial=9, np=5000, eps=0.005, qp=0.9
    options.pf.n_particles = 1000;
    options.pf.eps = 0.05;
    options.pf.qp = 0.9; % proportion of q0
    options.pf.version = 3;
    plotTitle = ['pf' num2str(options.pf.version) ': eps=' num2str(options.pf.eps) ' qp=' num2str(options.pf.qp) ' np=' num2str(options.pf.n_particles)];

    % --- compute z-score ---
    options.confidenceInterval = 2;% confidence interval = options.confidenceInterval * standard deviation
    options.ccconfidenceInterval = 2; % for ccf
    % --- compute detection ---
    options.threshold = 1.65; % 95%
    options.rou = [0.2,0.4,0.6,0.8]; % ccf facting factor 

    z_scores = cell(1,options.n_premodels); % zscore of the z vlaue(latent variable)
    for n_models=1:options.n_premodels
        options.n_models = n_models;
        onlineDecodeResult{n_decode,n_models} = onlineDecodeTrial(seq,testseq,model,options);    
        if(options.n_ccmodels==options.n_premodels)
            if n_decode==1
                z_scores{n_models} = computePvalueTrial(onlineDecodeResult{n_decode,n_models},seq,options);
            else
                z_scores{n_models} = computePvalueTrial(onlineDecodeResult{2,n_models},seq,options,onlineDecodeResult{2,n_models});
            end
        end
    end
    
    if(options.n_ccmodels~=options.n_premodels) %ccf value
        temp=onlineDecodeResult{n_decode,1};
        for n_ACCmodels=1:options.n_premodels
            for n_S1models=1:options.n_premodels
                eval(['temp.' options.decoders{1} '.ACC' seq.stimulus{1} '=onlineDecodeResult{' num2str(n_ACCmodels) '}.' options.decoders{1} '.ACC' seq.stimulus{1}]);
                eval(['temp.' options.decoders{1} '.S1' seq.stimulus{1} '=onlineDecodeResult{' num2str(n_S1models) '}.' options.decoders{1} '.S1' seq.stimulus{1}]);
                z_scores{(n_ACCmodels-1)*options.n_premodels+n_S1models} = computePvalueTrial(temp,seq,options);
            end
        end
%         options.n_premodels=options.n_premodels^2;
    end
        % --- draw plots ---
        % settings
        %plotOpt.raster = 1;
    plotOpt.region ={'ACC','S1'};
    plotOpt.decoders = {'all'};%{'PLDS_FB','PLDS_forward'};
    %plotOpt.cross = [1 1]; TODO
    plotOpt.trials = {'all'};%[1 2];
    plotOpt.baseline = 1;
    plotOpt.legend = 1;
    plotOpt.stimulus = {'all'};%{'50' 'VF'};
    plotOpt.withdraw = 1;
    plotOpt.zscoreColor = ['r','b','m','c','y','k'];
    plotOpt.baseColor = 'b';
    plotOpt.stimulusColor = 'g';
    plotOpt.withdrawalColor = 'k';
    plotOpt.threshColor = 'k';
    plotOpt.plotTitle = plotTitle;
    plotOpt.shadedCI = 1; % whether the confidence interval is drew as shaded area or not
    plotOpt.ylim = [-8 8];
    plotOpt.dirSaveFigure = dataOpt;
    plotOpt.regionsInOne = 1;
    options.plotOpt = plotOpt;
    options.trialPP = -1;
    options.trial250 = -1;
    options.trial50 = -1;
    options.trial150 = -1;
    if (n_decode==2)&&size(seq.stimulus,2)>1
        seq.stimulus = {options.stimulus{1}};
        options.teststimulus = {options.stimulus{2}};
    else
        seq.stimulus = {options.stimulus{1}};
        options.teststimulus = {options.stimulus{1}};
    end
    
    if 0 % plot the figure
       drawZscorePlots_Correlation(z_scores,seq,options); 
    end
    % --- statistics analysis --- evaluate the result, TP or FP
    options.areaThr=1.6; % ccf area threshold 
    options.cc=3; % rou=0.6, forget factor 
    if(n_decode==1)
        detectionRange_s = [0 3]; % unit=second, 0 is alignment event
    else
        detectionRange_s = [0 3];
    end
    detectionRange = (detectionRange_s + options.TPre)/options.binsize;
    detectionRange = detectionRange(1):detectionRange(end);
    % detectInfo = getDetectInfo(z_score, detectionRange, Thr, binsize, tpre);
    if ~options.PP
%         th_range=[10:-0.05:0];
        th_range = 1.65;
    else
        th_range=[10:-0.05:0.1,0.05:-0.01:0]; % for AUROC plot
    end
%     th_range=[1.65];
%     th_range=[3.08 2.57 2.33 1.65 1.28];
if 1
if 0
    w_range=[0:1:60];% decition making window for majority voting (1:60 bins)
    for ii=1:size(w_range,2)
        options.ACC_w=w_range(ii);
        options.S1_w=w_range(ii);
        ACC_r = zeros(eval(['length(z_scores{1}.' options.decoders{1} '.' plotOpt.region{1} options.stimulus{1} ')']),options.n_premodels);
        S1_r = ACC_r;
        for n_models=1:options.n_premodels
            detectStat = getDetectStat(z_scores{n_models},detectionRange,options);
            for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
                ACC_r(i,n_models) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).isDetected']);
                S1_r(i,n_models) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{2} '.allTrials(i).isDetected']);
            end
            ACCrate(n_decode,ii,n_models) = sum(ACC_r((1+options.n_premodels):end,n_models))/(size(ACC_r(:,n_models),1)-options.n_premodels);
            S1rate(n_decode,ii,n_models) = sum(S1_r((1+options.n_premodels):end,n_models))/(size(S1_r(:,n_models),1)-options.n_premodels);
        end
        detectStat = getDetectStatCorr(z_scores,detectionRange,options); 
        ACCensemble_r=zeros(size(ACC_r,1),1);
        S1ensemble_r=ACCensemble_r;
        for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
            ACCensemble_r(i)= eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).isDetected']);
            S1ensemble_r(i) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{2} '.allTrials(i).isDetected']);
        end
        ACCmajorityrate(n_decode,ii)=sum(ACCensemble_r((1+options.n_premodels):end))/(size(ACCensemble_r,1)-options.n_premodels);
        S1majorityrate(n_decode,ii)=sum(S1ensemble_r((1+options.n_premodels):end))/(size(S1ensemble_r,1)-options.n_premodels);
    end
else
    options.ACC_w=60;  %60 bins 
    options.S1_w=3; % 3 bins 
    for ii=1:size(th_range,2)
        options.threshold = th_range(ii);
        ACC_r = zeros(eval(['length(z_scores{1}.' options.decoders{1} '.' plotOpt.region{1} options.stimulus{1} ')']),options.n_premodels);
        S1_r = ACC_r;
        for n_models=1:options.n_premodels
            detectStat = getDetectStat(z_scores{n_models},detectionRange,options);
            for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
                ACC_r(i,n_models) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).isDetected']);
                S1_r(i,n_models) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{2} '.allTrials(i).isDetected']);
            end

            ACCrate(n_decode,ii,n_models) = sum(ACC_r((1+options.n_premodels):end,n_models))/(size(ACC_r(:,n_models),1)-options.n_premodels);            
            S1rate(n_decode,ii,n_models) = sum(S1_r((1+options.n_premodels):end,n_models))/(size(S1_r(:,n_models),1)-options.n_premodels);
%             rate(n_decode,ii,n_models) = sum(ACC_r(:,n_models) | S1_r(:,n_models))/size(ACC_r(:,n_models),1);
        end
        detectStat = getDetectStatCorr(z_scores,detectionRange,options); 
        ACCensemble_r=zeros(size(ACC_r,1),1);
        S1ensemble_r=ACCensemble_r;
        for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
            ACCensemble_r(i)= eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).isDetected']);
            S1ensemble_r(i) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{2} '.allTrials(i).isDetected']);
        end

        ACCmajorityrate(n_decode,ii)=sum(ACCensemble_r((1+options.n_premodels):end))/(size(ACCensemble_r,1)-options.n_premodels);
        S1majorityrate(n_decode,ii)=sum(S1ensemble_r((1+options.n_premodels):end))/(size(S1ensemble_r,1)-options.n_premodels);
    end  
end
end

if 1 % for ccf
if 1
    if ~options.PP
%         area_range=[120:-20:20,10:-0.2:0.2];
        area_range=[120:-20:20,10:-0.2:0.2,0.01:-0.001:0];
    else
        area_range=[120:-20:20,10:-0.2:0.2,0.01:-0.001:0];
    end
%     area_range=[0.8];
    options.cc=20; 
    for jj=1:size(area_range,2)
        options.areaThr=area_range(jj);
        CCACC_r = zeros(size(options.rou,2),eval(['options.ntrial' options.stimulus{n_decode}]),options.n_premodels);
        CC_result = CCACC_r;
        for n_models=1:options.n_premodels
            detectStat = getDetectStat(z_scores{n_models},detectionRange,options);
            for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
                CCACC_r(:,i,n_models) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).cc_area']);
                CC_result(:,i,n_models) = CCACC_r(:,i,n_models)>area_range(jj);
            end
        end
        CCrate(n_decode,jj,:,:) = permute(sum(CC_result(:,(1+options.n_premodels):end,:),2)/(size(CC_result,2)-options.n_premodels),[1,3,2]);
        %majority vote
        detectStat = getDetectStatCorr(z_scores,detectionRange,options);
        CCensemble_r=zeros(size(CCACC_r,2));
        for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
            CCensemble_r(i) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).cc_isDetected']);
        end
        CCmajorityrate(n_decode,jj,:) = sum(CCensemble_r((1+options.n_premodels):end))/(size(CCensemble_r,1)-options.n_premodels);
    end
else    
    options.areaThr=1.6;
    cc_range=[0:20];
    for jj=1:size(cc_range,2)
        options.cc=cc_range(jj);
        CCACC_r = zeros(size(options.rou,2),eval(['options.ntrial' options.stimulus{n_decode}]),options.n_premodels);
        CC_result = CCACC_r;
        for n_models=1:options.n_premodels
            detectStat = getDetectStat(z_scores{n_models},detectionRange,options);
            for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
                CCACC_r(:,i,n_models) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).cc_area']);
                CC_result(:,i,n_models) = CCACC_r(:,i,n_models)>options.areaThr;
            end
        end
        CCrate(n_decode,jj,:,:) = permute(sum(CC_result(:,(1+options.n_premodels):end,:),2)/(size(CC_result,2)-options.n_premodels),[1,3,2]);
        %majority vote
        detectStat = getDetectStatCorr(z_scores,detectionRange,options);
        CCensemble_r=zeros(size(CCACC_r,2));
        for i=1:size(eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials']),2)
            CCensemble_r(i) = eval(['detectStat.' options.decoders{1} '.S' options.stimulus{1} '.' plotOpt.region{1} '.allTrials(i).cc_isDetected']);
        end
        CCmajorityrate(n_decode,jj,:) = sum(CCensemble_r((1+options.n_premodels):end))/(size(CCensemble_r,1)-options.n_premodels);
    end
end
end
end
%% plot AUROC curve
exp1.name = 'CC';
exp1.rates = CCrate(:,:,3,1);
exp2.name = 'CC(ensemble)';
exp2.rates = CCmajorityrate;
exp3.name = 'ACC';
exp3.rates = ACCrate(:,:,1);
exp4.name = 'ACC(ensemble)';
exp4.rates = ACCmajorityrate;
exp5.name = 'S1';
exp5.rates = S1rate(:,:,1);
exp6.name = 'S1(ensemble)';
exp6.rates = S1majorityrate;
area = auroc('1.jpg',exp3,exp5,exp1,exp4,exp6,exp2);