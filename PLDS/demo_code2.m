
clear;clc;close all

data_opt = 2

if data_opt == 1
    load /Users/zhechen/Documents/MATLAB/TobiData/continuous/ACC_250mW.mat;
    C = 11; NoTrials = 32; data = ACC_250mW; clear ACC_250mW
    laser  = 1000*[0.0490    0.0884    0.1250    0.1621    0.1991    0.2353    0.2724    0.3082    0.3843    0.4213    0.4597    0.5730 ...
        0.6093    0.6462    0.6883    0.7243    0.7951    0.8306    0.8681    0.9398    1.0145    1.0530    1.1243    1.1611 ...
        1.1978    1.2355    1.3464    1.4193    1.4569    1.4954    1.5342    1.5715];
    withdrawal = 1000*[0.0495    0.0888    0.1256    0.1627    0.1999    0.2361    0.2728    0.3090    0.3855    0.4215    0.4604    0.5734 ...
        0.6098    0.6468    0.6887    0.7247    0.7952    0.8310    0.8685    0.9401    1.0154    1.0535    1.1252    1.1615 ...
        1.1982    1.2360    1.3469    1.4198    1.4574    1.4957    1.5346    1.5718];
    
elseif data_opt == 2
    %load /Users/zhechen/Documents/MATLAB/TobiData/continuous/newACCdata.mat;
    
    load ../newACCdata.mat;
    C = 15; NoTrials = 31; data = newACCdata; clear newACCdata
    laser = 1000 * [0.0664    0.1220    0.2578    0.3125    0.4056    0.5184    0.5733    0.6280    0.6822    0.7360    0.7928    0.8472 ...
        0.9016    0.9548    1.0090    1.0796    1.1339    1.2409    1.2968    1.3511    1.4260    1.4792    1.5376    1.5980 ...
        1.6563    1.7097    1.7639    1.8175    1.8732    1.9272    1.9799];
    
    withdrawal = 1000*[0.0672    0.1228    0.2597    0.3138    0.4068    0.5202    0.5740    0.6286    0.6828    0.7366    0.7947    0.8476 ...
        0.9023    0.9555    1.0096    1.0805    1.1348    1.2419    1.2979    1.3519    1.4270    1.4803    1.5386    1.5990 ...
        1.6572    1.7110    1.7647    1.8184    1.8741    1.9280    1.9814];
    
    
elseif data_opt == 3
    load /Users/zhechen/Documents/MATLAB/TobiData/continuous/newdata/ratA250mw_ACC_SUA.mat
    xDim = 1;
    C = 8;    NoTrials = 34; data = ratA250mw_ACC_SUA; clear ratA250mw_ACC_SUA
    laser = 1000 * [0.0669    0.1153    0.1631    0.2127    0.2609    0.3090    0.3567    0.4071    0.4543    0.5030    0.5521 ...
         0.5984    0.6460    0.6958    0.7441    0.7955    0.8456    0.8941    0.9442    0.9931    1.0401    1.0879 ...
         1.1377    1.1873    1.2353    1.2833    1.3320    1.3820    1.4309    1.4785    1.5265    1.5745    1.7377 1.8630];

    withdrawal = 1000*[0.0681    0.1161    0.1637    0.2136    0.2639    0.3105    0.3581    0.4085    0.4559    0.5038    0.5529 ...
         0.5991    0.6470    0.6967    0.7450    0.7972    0.8468    0.8952    0.9455    0.9942    1.0409    1.0886 ...
         1.1394    1.1887    1.2359    1.2843    1.3332    1.3831    1.4317    1.4791    1.5270    1.5768    1.7387 1.8642];
end

Laser_Withdrawal_delay  = withdrawal - laser;
       
  
    
binsize = 0.05; % 50 ms
TT1 = 5; TT1 = 5;


TT2 = TT1;
edges = [-TT1:binsize:TT2];  
  
  
xDim = 1; ind_C = [1:C];

path(path,'..\')
path(path,'minFunc')
path(path,'minFunc\minFunc')
path(path,'minFunc\minFunc\mex')
path(path,'minFunc\minFunc\compiled')
path(path,'minFunc\logisticExample')
path(path,'minFunc\autoDif')
path(path,'..\cfa')
    
maxIter = 50;
for trial= 1:5 % NoTrials
    
    seq = struct('T', length(edges), 'y',zeros(length(ind_C),length(edges)));
    yDim = length(ind_C);
    
    
    Spikecount = zeros(length(ind_C),length(edges));
    for c=1:length(ind_C);
        cc = ind_C(c);
        trigger = data(c).laserOn; % laser
        % trigger = data(cc).recFlick; % withdrawal
        
        ind = find(data(cc).spikes>=trigger(trial)-TT1 & data(cc).spikes<=trigger(trial)+TT2);
        
        count = histc(data(cc).spikes(ind)-trigger(trial), edges);
        Spikecount(c,:) = count;
    end
  
     
    seq = setfield(seq,{trial},'T', size(Spikecount,2));   
    seq = setfield(seq,{trial},'y', Spikecount); 
     
    
    
    figure(trial);
    
    
    params = [];
    
    % -------- LDS 
    path(path,'/Users/zhechen/Documents/MATLAB/NichoData/lds');
    y1 = log(Spikecount+1); % log-transformation
    y2 = sqrt(Spikecount); % square root
    ssm1 = lds(y1',1,length(edges),300,1e-6);
    %ssm2 = lds(y2',1,length(edges),200,1e-6);
    
    % -------- PLDS Parameter Initialization
    
    params = PLDSInitialize(seq(trial),xDim,'NucNormMin',params);
    %params = PLDSInitialize(seq(trial),xDim,'ExpFamPCA',params);
    
    params.model.inferenceHandle = @PLDSLaplaceInference;                           % comment out for using variational infernce
    
    params.opts.algorithmic.EMIterations.maxIter     = maxIter;						% setting max no of EM iterations
    params.opts.algorithmic.EMIterations.maxCPUTime  = 30;			   		        % setting max CPU time for EM to 600s
    
    %[params newseq varBound] = PopSpikeEM(params,seq(trial));
    %[params newseq Qfun] = PLDS_EM(params,newseq);
    
    [params newseq Qfun] = PLDS_EM(params,seq(trial));
    
 
    
    if trial>1
        [pred_state,pred_stateCov] = PLDS_online_filtering(params,Spikecount,last_z0,last_Q0,binsize);
        pred_Zscore = (pred_state-base1)/base2;  % previous trial stat (base1, base2)
        
        % Kalman filtering using LDS 
        [pred_stateKF,pred_stateCovKF] = KF_filter(ssm1,Spikecount,last_z0KF,last_Q0KF); 
    end
    
    
    %% assessment
    
    
    
    
    h1=subplot(311);
    imagesc(edges, [1:length(ind_C)],newseq.y);axis xy
    colormap(flipud(bone));
    ylabel('Cell','fontsize',16);
    title(num2str(trial),'fontsize',24);
    set(gca,'fontsize',16)
    
    set(gca,'Xticklabel',[]);
    
    
    
    range = [1: 4/binsize];
    base1 = mean(newseq.posterior.xsm(range));
    base2 = std(newseq.posterior.xsm(range));
    
    
    temporal_Z = (newseq.posterior.xsm - base1) / base2 ;
    
%     p1 = get(h1, 'pos');
%     % 4-element vector: [left, bottom, width, height]
%     h2=subplot(212);
%     p2 = p1;
%     p2(2) = p1(2) - 0.38;
%     set(h2, 'pos', p2);
    
    subplot(312);
    hold on;
    tem1 = temporal_Z; % * sign(max(params.model.C)); % * sign(median(params.model.C));
    tem2 =  newseq.posterior.Vsm(1:xDim:end);
    
    
    ratio1 =   (newseq.posterior.xsm ) / (base2);
    
    if trial > 1
        base3 = pred_state-mean(pred_state(range));
        ratio2 =   base3   / std(pred_state(range));
    end
    
    CUSUM = zeros(1,length(ratio1));
    CUSUM2 = zeros(1,length(ratio1));
    CUSUM31 = zeros(C,length(ratio1));
    CUSUM32 = zeros(C,length(ratio1));
    CUSUM3 = zeros(1,length(ratio1));
    CUSUM41 = zeros(C,length(ratio1));
    CUSUM4 = zeros(1,length(ratio1));
    
    for t=2:length(ratio1)
        % Poisson 
        omega = exp(newseq.posterior.xsm(t))*binsize - exp(base2)*binsize;
        temp = exp(newseq.posterior.xsm(t)) *binsize * (newseq.posterior.xsm(t)-base2) - omega;
        CUSUM(t) = max(0, CUSUM(t-1)+ temp);
        
        
        
        
        % Poisson on raw data
        lam0 = mean(seq(trial).y(:,range),2)+ 0.1;
        lam11 = lam0 + 2*sqrt(lam0);             
        temp11 = seq(trial).y(:,t) .* log(lam11./lam0) - ( lam11-lam0);
        CUSUM31(ind_C,t) = max(0, CUSUM31(ind_C,t-1)+ temp11);
        CUSUM3(t) = max(CUSUM31(:,t));
        %lam12 = max(0.1, lam0- 2*sqrt(lam0));
        % temp12 = seq(trial).y(:,t) .* log(lam12./lam0) - ( lam12-lam0);
        %CUSUM32(t) = min(0, CUSUM3(t-1)+ min(temp12));
        
        % Gauss on raw data
        mean0 = mean(seq(trial).y(:,range),2);
        var0 = var(seq(trial).y(:,range),0,2);
        var1 = var0;
        
        temp11 = ( var1*log(2*pi) + var1.*log(var1) + (seq(trial).y(:,t)-mean0-sqrt(var0)).^2 ) ./ ...
                 ( var0*log(2*pi) + var0.*log(var0) + (seq(trial).y(:,t)-mean0).^2 );
        temp11 =  temp11 .* (var0./var1); 
             
        CUSUM41(ind_C,t) = max(0, CUSUM41(ind_C,t-1)+ temp11);     
        CUSUM4(t) = max(CUSUM41(:,t));
        
%         if CUSUM3(t) > 10
%             CUSUM31(:,t) = 0;
%         end
    end
    
    hold on;
    fill([edges'; flipdim(edges',1)], [tem1'+sqrt(tem2)/base2; flipdim(tem1'-sqrt(tem2)/base2,1)], [6 6 6]/8,'EdgeColor',[6 6 6]/8);  
    plot(edges,tem1 ,'r-','linewidth',2);
    %hold on;plot(edges,tem1'+1*sqrt(tem2)/base2 ,'r--','linewidth',1);
    %hold on;plot(edges,tem1'-1*sqrt(tem2)/base2 ,'r--','linewidth',1);
    
 
    hold on;plot([edges(1) edges(end)],[1.65 1.65],'k--','linewidth',1);
    hold on;plot([edges(1) edges(end)],[-1.65 -1.65],'k--','linewidth',1);
    
    
    
    tem = squeeze(ssm1.x);
    tem_Z1 = (tem - mean(tem(range))) / std(tem(range));
    hold on;
    fill([edges'; flipdim(edges',1)], [tem_Z1+sqrt(squeeze(ssm1.P)); flipdim(tem_Z1-sqrt(squeeze(ssm1.P)),1)], [7 7 7]/8,'EdgeColor',[7 7 7]/8);

    plot(edges, tem_Z1,'b-','linewidth',2);
    %plot(edges, tem_Z1 + sqrt(squeeze(ssm1.P)),'k--','linewidth',1);
    %plot(edges, tem_Z1 - sqrt(squeeze(ssm1.P)),'k--','linewidth',1);
    
    %tem = squeeze(ssm2.x);
    %tem_Z2 = (tem - mean(tem(range))) / std(tem(range));
    %plot(edges, tem_Z2,'k-','linewidth',2);
    
    
    if trial > 1
        hold on;
        tem1 = (pred_state - mean(pred_state(range)))/ std(pred_state(range));
        tem2 =  pred_stateCov / std(pred_state(range));
        
        plot(edges, tem1, 'm-','linewidth',2);      
        %hold on;plot(edges,  tem1'+sqrt(tem2)/base2 ,'m--','linewidth',1);
        %hold on;plot(edges,  tem1'-sqrt(tem2)/base2 ,'m--','linewidth',1);
        
        tem3 = (pred_stateKF - mean(pred_stateKF(range)))/ std(pred_stateKF(range));
        tem4 =  pred_stateCovKF / std(pred_stateKF(range));
        
        plot(edges, tem3, 'g-','linewidth',2); 
        
    end
    
    ylabel('Z-score','fontsize',16);
    set(gca,'fontsize',16)
    
    %hold on;plot([0 0],[-5 10],'k-','linewidth',2)
    ylim([-5 5]);
    
    
    
    subplot(313);
    hold on;plot(edges, CUSUM3,'k','linewidth',2);
    %hold on;plot(edges, CUSUM4,'y','linewidth',2);
    
    hold on;plot([edges(1) edges(end)],[3.4 3.4],'k--','linewidth',1);
    ylabel('Cumulative Stat','fontsize',16);
    xlabel('Time (s)','fontsize',16);
    set(gca,'fontsize',16);
    ylim([-1 10]);
    
    % for next trial 
    last_z0 = newseq.posterior.xsm(end); 
    tem2 =  newseq.posterior.Vsm(1:xDim:end);last_Q0 = tem2(end-xDim+1:end);
    tem1 = squeeze(ssm1.x); tem2 = squeeze(ssm1.P);
    last_z0KF = tem1(end); last_Q0KF = tem2(end);
    
    
    if trial == 1 % cell pruning 
        %[ssm1.C params.model.C]
        [temp2,ind2]=sort(abs(params.model.C),'descend');
        ind_C = ind2(1:8);
        
    end
end