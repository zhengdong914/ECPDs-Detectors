clear;clc;close all
path(path,'..\')
path(path,'minFunc')
path(path,'minFunc\minFunc')
path(path,'minFunc\minFunc\mex')
path(path,'minFunc\minFunc\compiled')
path(path,'minFunc\logisticExample')
path(path,'minFunc\autoDif')
path(path,'..\cfa')
%--------data options------------
data_opt = 1;%2;%1;
if data_opt == 1
    load pain_S1_population.mat;
    C = 12;       % channel numbers
    xDim = 1;     % hidden variable demension 
    NoTrials = 40;% trial numbers
    Laser_Withdrawal_delay = [6.4644    4.7855   10.8595    4.4506    2.3268    7.6650    3.7381    5.5232    4.7824    5.8435    8.1605    6.9696    6.6968    5.2428 ...
                              10.7525    6.9056   30.3126    4.5093    5.2475    7.5384    5.7005    5.0687    4.6128    5.2260    4.7582    4.9953    6.5654    3.7573 ...    
                              8.0604    8.7205    6.2504    5.9855    5.9615    5.4247    3.5108    2.1469    3.5210    5.4032    5.9143    2.5045];
elseif data_opt == 2
    load ../newACCdata.mat;
    C = 15; 
    NoTrials = 31; 
    data = newACCdata; 
    clear newACCdata
    laser = 1000 * [0.0664    0.1220    0.2578    0.3125    0.4056    0.5184    0.5733    0.6280    0.6822    0.7360    0.7928    0.8472 ...
        0.9016    0.9548    1.0090    1.0796    1.1339    1.2409    1.2968    1.3511    1.4260    1.4792    1.5376    1.5980 ...
        1.6563    1.7097    1.7639    1.8175    1.8732    1.9272    1.9799];
    
    withdrawal = 1000*[0.0672    0.1228    0.2597    0.3138    0.4068    0.5202    0.5740    0.6286    0.6828    0.7366    0.7947    0.8476 ...
        0.9023    0.9555    1.0096    1.0805    1.1348    1.2419    1.2979    1.3519    1.4270    1.4803    1.5386    1.5990 ...
        1.6572    1.7110    1.7647    1.8184    1.8741    1.9280    1.9814];
end

%--------allocate data structure------
binsize = 0.05; %50ms
TT1 = 5;
TT2 = 5;
edges = [-TT1:binsize:TT2];
%T: number of sample points
%y: C dimension(channel) data of T samples.
if data_opt == 1
    seq = struct('T', length(edges), 'y',zeros(C,length(edges)));
    for j=1:NoTrials
        seq = setfield(seq,{j},'T', AllTrialSeq(j).T);
        seq = setfield(seq,{j},'y', AllTrialSeq(j).y);
    end
    yDim = C;
elseif data_opt == 2
    ind_C = [1:C];
    xDim = 1;
    seq = struct('T', length(edges), 'y',zeros(length(ind_C),length(edges)));
    for trial= 1:NoTrials
        yDim = length(ind_C); 
        Spikecount = zeros(length(ind_C),length(edges));
        for c=1:length(ind_C);
            cc = ind_C(c);
            %trigger = data(c).laserOn; % laser
            trigger = data(cc).recFlick; % withdrawal
            ind = find(data(cc).spikes>=trigger(trial)-TT1 & data(cc).spikes<=trigger(trial)+TT2);
            count = histc(data(cc).spikes(ind)-trigger(trial), edges);
            Spikecount(c,:) = count;
        end
        seq = setfield(seq,{trial},'T', size(Spikecount,2));
        seq = setfield(seq,{trial},'y', Spikecount);
    end 
end

% for test
C_rec = [];
d_rec = [];
x_rec = [];
a_rec = [];
pi_rec = [];
q_rec = [];
q0_rec = [];
x0_rec = [];

% linear filter
z_score_rec = [];
z_var_rec = [];
% nonlinear
z_score_rec1 = [];
z_var_rec1 = [];
base2_rec = [];

%--------Iterate training and decoding over all trials--------
maxIter = 50;
for trial= 1: NoTrials  
    % for test
%     if trial==18 || trial==29 || trial==38
%         C_rec = [C_rec zeros(size(params.model.C)) zeros(size(params.model.C))];
%         continue;
%     end
    display(trial);
    h = figure(trial);
    % -------- Model Parameter Initialization
    params = [];
    
    params = PLDSInitialize(seq(trial),xDim,'ExpFamPCA',params);
    params1 = PLDSInitialize(seq(trial),xDim,'NucNormMin',params);
    %params = PLDSInitialize(seq(trial),xDim,'params',params);
    params.model.inferenceHandle = @PLDSLaplaceInference;                           % comment out for using variational infernce
    params.opts.algorithmic.EMIterations.maxIter     = maxIter;						% setting max no of EM iterations
    params.opts.algorithmic.EMIterations.maxCPUTime  = 30;			   		        % setting max CPU time for EM to 600s

     C_rec = [C_rec [params.model.C params1.model.C]];
%     d_rec = [d_rec [params.model.d params1.model.d]];
%     x_rec = [x_rec [params.model.Xpca' params1.model.Xpca']];
%     a_rec = [a_rec; [params.model.A params1.model.A]];
%     pi_rec = [pi_rec; [params.model.Pi]];%params1.model.Pi = 0.05, pi=sum(x^2)/T
%     q_rec = [q_rec; [params.model.Q params1.model.Q]];
%     q0_rec = [q0_rec; [params.model.Q0]];% params1.model.Q0 = 0.05
    %x0_rec = [x0_rec [params.model.x0 params1.model.x0]];% all 0
    
    [params newseq varBound] = PLDS_EM(params,seq(trial),binsize);

    %---------online filtering----------------
    if trial>1
        [pred_state, var] = PLDS_online_filtering(params,seq(trial).y,last_z0,binsize);
        %pred_Zscore = (pred_state-base1)/base2;  % previous trial stat (base1, base2) 
    end
    last_z0 = newseq.posterior.xsm(end);
    
    % add for verify modify correctness
%     if (trial == 2)
%         load('trial2_z.mat', 'zz');
%         if sum(abs(zz-newseq.posterior.xsm))>1e-4
%             display('something wrong, the result is not exactly same');
%         else
%             display('all right, the result is exactly the same');
%         end
%     end
    
    %--------assement/draw figures------------ ** modify at last
    subplot(311);
    imagesc(edges, [1:C],newseq.y);axis xy
    colormap(flipud(bone));
    ylabel('Cell','fontsize',16);
    title('Single Trial','fontsize',20);
    set(gca,'fontsize',16)

    range = [1: 80];
    base1 = mean(newseq.posterior.xsm(range));
    base2 = std(newseq.posterior.xsm(range));
    
    temporal_Z = (newseq.posterior.xsm - base1) / base2 ;
    subplot(313);
    hold on;
    tem1 = temporal_Z * sign(median(params.model.C));
    tem2 = newseq.posterior.Vsm(1:xDim:end);
    plot(edges,tem1 ,'r-','linewidth',2);% red curve middle: temporal_Z
    hold on;plot(edges,tem1'+2*sqrt(tem2)/base2 ,'r--','linewidth',1);% red curve up: temporal_Z + 2*std
    hold on;plot(edges,tem1'-2*sqrt(tem2)/base2 ,'r--','linewidth',1);% red curve down: temporal_Z - 2*std
    
    hold on;plot([edges(1) edges(end)],[1.65 1.65],'k--','linewidth',1);% black line up: 1.65
    hold on;plot([edges(1) edges(end)],[-1.65 -1.65],'k--','linewidth',1);% black line down: -1.65
    
    hold on;plot([edges(1) edges(end)],[3.09 3.09],'b--','linewidth',1);% blue line up: 3.09
    hold on;plot([edges(1) edges(end)],[-3.09 -3.09],'b--','linewidth',1);% blue line down: -3.09
    
    % detect cross thresh
    % mean in the middle
    
    % plot detection result
%     [y1,index1] = max(tem1(TT1/binsize: TT1/binsize + 5/binsize));  % 6 second window
%     [y2,index2] = min(tem1(TT1/binsize: TT1/binsize + 5/binsize));
    [y1,index1] = max(tem1((TT1-1)/binsize: TT1/binsize));  % 6 second window
    [y2,index2] = min(tem1((TT1-1)/binsize: TT1/binsize));
    
    %confident interval level
    ci = 1;
    
    if y1 > -y2 % peak > trough
        temp0 = tem1'-ci*sqrt(tem2)/base2;  % positive - lower CI
        %Stat1(trial) = max(temp0(TT1/binsize: TT1/binsize + 5/binsize));
        Stat1(trial) = max(temp0((TT1-1)/binsize: TT1/binsize));%-1s-0s
    elseif -y2 > y1 % trough > peak
        temp0 = tem1'+ ci*sqrt(tem2)/base2; % negative + higher CI
        %Stat1(trial) = abs(min(temp0(TT1/binsize: TT1/binsize + 5/binsize)));
        Stat1(trial) = abs(min(temp0((TT1-1)/binsize: TT1/binsize)));%-1s-0s
    end
    
%PLDS online filter:
    if trial > 1
       hold on; 
       %plot(edges, pred_Zscore,'k-','linewidth',2); % black curve: linear filter result
       z_score = (pred_state - mean(pred_state(range)))/ std(pred_state(range));
       plot(edges, z_score, 'm-','linewidth',2);% mageta curve: linear filter result
       z_score_rec = [z_score_rec; z_score];
       z_var_rec = [z_var_rec; var'];
    end
    z_score_rec1 = [z_score_rec1; tem1];
    z_var_rec1 = [z_var_rec1; tem2'];
    base2_rec = [base2_rec; base2];
    
    ylabel('Z-score','fontsize',16);
    xlabel('Time (s)','fontsize',16);
    title('PLDS','fontsize',20);
    set(gca,'fontsize',16)
    %plotname = ['directwithQfun' num2str(trial)];
    %plotname = ['fb' num2str(trial)];
    %saveas(h,plotname,'fig');    
end

% save('z_linear_expPca.mat','z_score_rec');
% save('z_var_linear_expPca.mat','z_var_rec');
% save('z_nonlinear_expPca.mat','z_score_rec1');
% save('z_var_nonlinear_expPca.mat','z_var_rec1');
% save('z_base2_expPca.mat','base2_rec');

% save('z_linear_Nuc.mat','z_score_rec');
% save('z_var_linear_Nuc.mat','z_var_rec');
% save('z_nonlinear_Nuc.mat','z_score_rec1');
% save('z_var_nonlinear_Nuc.mat','z_var_rec1');%for test
%save('z_base2_Nuc.mat','base2_rec');
%[length(find(Stat1>1.65)), length(find(Stat1>2.33)), length(find(Stat1>2.58)), length(find(Stat1>3.09)), length(find(Stat1>3.32)), length(find(Stat1>3.52))]
detection165 = find(Stat1>1.65)';
detection233 = find(Stat1>2.33)';
detection258 = find(Stat1>2.58)';
detection309 = find(Stat1>3.09)';
detection332 = find(Stat1>3.32)';
detection352 = find(Stat1>3.52)';

a=0;
