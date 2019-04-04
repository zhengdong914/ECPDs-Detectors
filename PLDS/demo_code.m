
clear;clc;close all


path(path,'..\')
path(path,'minFunc')
path(path,'minFunc\minFunc')
path(path,'minFunc\minFunc\mex')
path(path,'minFunc\minFunc\compiled')
path(path,'minFunc\logisticExample')
path(path,'minFunc\autoDif')
path(path,'..\cfa')

load pain_S1_population.mat; % 150 mW
C = 12;
xDim = 1;
NoTrials = 40;
Laser_Withdrawal_delay = [6.4644    4.7855   10.8595    4.4506    2.3268    7.6650    3.7381    5.5232    4.7824    5.8435    8.1605    6.9696    6.6968    5.2428 ...
    10.7525    6.9056   30.3126    4.5093    5.2475    7.5384    5.7005    5.0687    4.6128    5.2260    4.7582    4.9953    6.5654    3.7573 ...
    8.0604    8.7205    6.2504    5.9855    5.9615    5.4247    3.5108    2.1469    3.5210    5.4032    5.9143    2.5045];

% time 0 is laser / withdrawal 

binsize = 0.05;
TT1 = 5; TT1 = 5;


TT2 = TT1;
edges = [-TT1:binsize:TT2];
seq = struct('T', length(edges), 'y',zeros(C,length(edges)));

%%

  

maxIter = 50;
for trial=  1:3 % [1: NoTrials]
    
    binsize = 0.05;
    
    for j=1:NoTrials
        seq = setfield(seq,{j},'T', AllTrialSeq(j).T);
        seq = setfield(seq,{j},'y', AllTrialSeq(j).y); yDim = C;
    end
    
    
    figure(trial);
    
    
    params = [];
    
    % -------- LDS 
    %path(path,'/Users/zhechen/Documents/MATLAB/NichoData/lds');
    path(path,'../lds');
    y1 = log(seq(trial).y+1); % log-transformation
    y2 = sqrt(seq(trial).y); % square root
    ssm1 = lds(y1',1,length(edges),200,1e-6);
    %ssm2 = lds(y2',1,length(edges),200,1e-6);
    
    % -------- PLDS Parameter Initialization
    
    params = PLDSInitialize(seq(trial),xDim,'NucNormMin',params);
    
    params.model.inferenceHandle = @PLDSLaplaceInference;                           % comment out for using variational infernce
    
    params.opts.algorithmic.EMIterations.maxIter     = maxIter;						% setting max no of EM iterations
    params.opts.algorithmic.EMIterations.maxCPUTime  = 30;			   		        % setting max CPU time for EM to 600s
    
    %[params newseq varBound] = PopSpikeEM(params,seq(trial));
    %[params newseq Qfun] = PLDS_EM(params,newseq);
    
    [params newseq Qfun] = PLDS_EM(params,seq(trial));
    
    if trial>1
        [pred_state,pred_stateCov] = PLDS_online_filtering(params,seq(trial).y,last_z0,last_Q0,binsize);
        pred_Zscore = (pred_state-base1)/base2;  % previous trial stat (base1, base2)
        
        % Kalman filtering using LDS 
        [pred_stateKF,pred_stateCovKF] = KF_filter(ssm1,seq(trial).y,last_z0KF,last_Q0KF); 
    end
    
    
    %% assessment
    
    
    
    
    h1=subplot(311);
    imagesc(edges, [1:C],newseq.y);axis xy
    colormap(flipud(bone));
    ylabel('Cell','fontsize',16);
    title(num2str(trial),'fontsize',24);
    set(gca,'fontsize',16)
    
    set(gca,'Xticklabel',[]);
    
    
    
    range = [1: 100];
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
    tem1 = temporal_Z * sign(max(params.model.C)); % * sign(median(params.model.C));
    tem2 =  newseq.posterior.Vsm(1:xDim:end);
    
    
    ratio1 =   (newseq.posterior.xsm ) / (base2);
    
    if trial > 1
        base3 = pred_state-mean(pred_state(range));
        ratio2 =   base3   / std(pred_state(range));
    end
    
    CUSUM = zeros(length(ratio1),1);
    CUSUM2 = zeros(length(ratio1),1);
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
        CUSUM31(:,t) = max(0, CUSUM31(:,t-1)+ temp11);
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
             
        CUSUM41(:,t) = max(0, CUSUM41(:,t-1)+ temp11);     
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
end