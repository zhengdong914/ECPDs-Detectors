% draw z-score plots based on different options
% one figure for each trial
% author: Sile Hu
% date: 2017-3-13
function drawZscorePlots_Correlation(z_scores,seq,options)
if length(z_scores)>1
    zscores=z_scores{1};
else
    zscores=z_scores;
end
if strcmp(options.plotOpt.stimulus,'all')
    stimulus = seq.stimulus;
else
    stimulus = options.plotOpt.stimulus;
end

if strcmp(options.plotOpt.decoders,'all')
    decoders = options.decoders;
else
    decoders = options.plotOpt.decoders;
end

if strcmp(options.plotOpt.region,'all')
    region = seq.region;
else
    region = options.plotOpt.region;
end

legends = {};
edges = [-options.TPre:options.binsize:options.TPost];
for s=1:length(stimulus)
    if strcmp(options.plotOpt.trials,'all')
        if options.currentTrial || options.preTrial
            trials = 1:eval(['seq.ntrial' stimulus{s}]);
        else
            trials = 1:eval(['length(zscores.' options.decoders{1} '.' seq.region{1} seq.stimulus{1} ')']);
        end
    else
        trials = options.plotOpt.trials;
    end
    st = 1;
    et = length(trials);
    if isfield(options,['trial' seq.stimulus{s}])
        single_trial = eval(['options.trial' seq.stimulus{s}]);
        if (single_trial>=1 && single_trial<=et)
            trials = single_trial;
            et = 1;
        end
    end
    if strcmp(stimulus{s},'50') || strcmp(stimulus{s},'150') || strcmp(stimulus{s},'250')
        mw = 'mW';
    else
        mw = [];
    end
    for t=st:et
        figure('name',[options.teststimulus{s} ' trial' num2str(trials(t))]);

        subplotIdx = '411';
        shift = eval(['seq.ntrial' options.teststimulus{s} '- length(zscores.' options.decoders{1} '.' seq.region{1} seq.stimulus{1} ')']);
        seqall = eval(['seq.' seq.region{1} options.teststimulus{s} '(trials(t)+shift)']);
        titleStr = ['test stimulus=' options.teststimulus{s} mw ',trial' num2str(trials(t)) ' ' region{1} '=1:' num2str(size(seqall.y,1))];        
        if(length(seq.region)==1)
            seqall.y = [seqall.y;eval(['seq.' seq.region{1} options.teststimulus{s} '(trials(t)+shift).y'])];
        else
            seqall.y = [seqall.y;eval(['seq.' seq.region{2} options.teststimulus{s} '(trials(t)+shift).y'])];
        end
        
        drawRasterTrialNew(seqall,titleStr,edges,subplotIdx);          
        
        % lower plot
        subplotIdx = '412';
        % draw baseline
        legends = drawBaseline(options,edges,legends,subplotIdx);
        legends = drawThresh(options,edges,legends,subplotIdx);

        for r=1:length(region)
            if (r==2)
                for d=1:length(decoders)
                    if ~(strcmp(decoders{d},'GT'))
                        if (length(seq.region)==1)
                            z_score = eval(['zscores.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t));']);
                        else
                            z_score = eval(['zscores.' decoders{d} '.' seq.region{2} stimulus{s} '(trials(t));']);
                        end
                        l_z = length(z_score.zscore) - 1;
                        decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                        if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')) && options.currentTrial
                            %decoder_edges = decoder_edges(2:end);
                            decoder_edges = decoder_edges + options.binsize;
                        end
                        color = options.plotOpt.zscoreColor(d+r-1);
                        legends = drawZscoreTrial(z_score,decoder_edges,legends,decoders{d},region{r},options,color,subplotIdx);
                    end
                end
            else
                for d=1:length(decoders)
                    z_score = eval(['zscores.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t));']);
                    l_z = length(z_score.zscore) - 1;
                    decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                    if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')) && options.currentTrial
                        %decoder_edges = decoder_edges(2:end);
                        decoder_edges = decoder_edges + options.binsize;
                    end
                    color = options.plotOpt.zscoreColor(d+r-1);
                    legends = drawZscoreTrial(z_score,decoder_edges,legends,decoders{d},region{r},options,color,subplotIdx);
                end
            end
        end

        % draw event
        if isfield(seq, ['withdrawal' options.teststimulus{s} ])
            withdrawEvent = eval(['seq.withdrawal' stimulus{s} ';']);
        else
            if isfield(seq, ['withdraw' options.teststimulus{s} ])
                withdrawEvent = eval(['seq.withdraw' options.teststimulus{s} ';']);
            else
                withdrawEvent = [];
            end
        end
        legends = drawEvent(withdrawEvent,options,stimulus{s},legends,trials(t),subplotIdx);

        %legend = drawCross();

        % font, lengend, limit
        ylabel('Z-score','fontsize',24);
        set(gca,'fontsize',20);
        ylim(options.plotOpt.ylim);
%         xlabel('time(s)','fontsize',24);
        xlim([-options.TPre,options.TPost])
        subplot(subplotIdx)
%         lgd = legend(legends);
%         lgd.FontSize = 10;
        legends = [];
%         title(options.plotOpt.plotTitle);

        subplotIdx = '413';         
        drawBaseline(options,edges,legends,subplotIdx);       
        for d=1:length(decoders)
            if ~(strcmp(decoders{d},'GT'))
                for rou=1:length(options.rou)-1
                    cct = eval(['zscores.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t)).ccvalue(:,:,rou);']);                
                    options.C_thresh=eval(['zscores.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t)).ccthresh(:,:,rou);']);
                    inx1 = find(cct(options.TPre/options.binsize:(options.TPre+2)/options.binsize)>options.C_thresh(2) | cct(options.TPre/options.binsize:(options.TPre+2)/options.binsize)<options.C_thresh(1));
                    inx2 = find(cct((options.TPre+2)/options.binsize:(options.TPre+options.TPost)/options.binsize)>options.C_thresh(1) & cct((options.TPre+2)/options.binsize:(options.TPre+options.TPost)/options.binsize)<options.C_thresh(2));
                    if ~isempty(inx1)&&~isempty(inx2)
                        options.crossvalue=[(inx1(1)-1)*options.binsize (inx2(1)-1)*options.binsize+2];
                    else
                        options.crossvalue=[];
                    end
                    color = options.plotOpt.zscoreColor(4-rou);
                    legends = drawCorrelationPlots(cct,options,edges,legends,subplotIdx,color,options.rou(rou));
                    if(options.rou(rou)==0.6)
                        area=zeros(size(cct));
                        for tt=2:length(cct)
                            if(cct(tt)>options.C_thresh(2))
                                area(tt)=area(tt-1)+cct(tt)-options.C_thresh(2);
                            elseif(cct(tt)<options.C_thresh(1))
                                area(tt)=area(tt-1)+options.C_thresh(1)-cct(tt);
                            end
                        end
                    end
                end
            end
        end
        ylabel('CCF','fontsize',24);
        set(gca,'fontsize',20);
        xlim([-options.TPre,options.TPost])
        ylim([-2,2]);
%         xlabel('time (s)','fontsize',24);
        subplot(subplotIdx)
        lgd = legend(legends);
        lgd.FontSize = 10;
        legends = [];
        if ~isempty(options.crossvalue)
%             title(['Cross Thresh  On:' num2str(options.crossvalue(1)) 's   Off:' num2str(options.crossvalue(2)) 's     Confidential Interval:' num2str(options.confidenceInterval)]);
        end
        
        subplotIdx = '414';
        subplot(subplotIdx);
        hold on;
        p_area=plot(edges,area,[color '-'],'linewidth',2,'DisplayName',['area above threshold']);
        legends=[p_area];
        p_thr1 = plot([edges(1) edges(end)],[1.6 1.6],[options.plotOpt.threshColor '--'],'linewidth',1,'DisplayName',['Area threshold:' num2str(1)]);
        ylabel(['Area above';'threshold '],'FontSize',24);
        set(gca,'fontsize',20);
        xlim([-options.TPre,options.TPost])
        xlabel('Time (s)','fontsize',24);
        
        % model free
%         ytest = 11 - seqall.y(9,:);
%         lam0 = mean(ytest(:,[1/options.binsize : 4/options.binsize]),2)+ 0.1;
%         lam11 = lam0 - 2*sqrt(lam0);   
% %         lam11 = lam0*0.5;
%         CUSUM31 = zeros(size(ytest,1),size(ytest,2));
%         CUSUM3 = zeros(1,size(ytest,2));
%         for t=2:length(ytest)
%             temp11 = ytest(:,t) .* log(lam11./lam0) - ( lam11-lam0);
%             temp33 = CUSUM31(:,t-1)+temp11;
%             [M,I] = max(temp33);
%             
%             CUSUM31(:,t) = max(0, M);
%             temp22=CUSUM31(:,t);
%             CUSUM3(t) = max(temp22);
% %             if CUSUM3(t)>200
% %                 CUSUM31(:, t) = 0;
% %             end
%                 
%         end
%         ModelFreeZ = max(CUSUM3,0);
%         plot(edges,ModelFreeZ ,[color '-'],'linewidth',2,'DisplayName',[' Model-free']);
        
        
        annotation('textbox',[0.0635000000000001 0.950000000000000 0.0109047619047619 0.0265339966832503],...
        'String',{'A'},'LineStyle','none','FontSize',40,'FitBoxToText','off');
        annotation('textbox',[0.0635000000000001 0.750000000000000 0.0109047619047619 0.0265339966832503],...
        'String',{'B'},'LineStyle','none','FontSize',40,'FitBoxToText','off');
        annotation('textbox',[0.0635000000000001 0.550000000000000 0.0109047619047619 0.0265339966832503],...
        'String',{'C'},'LineStyle','none','FontSize',40,'FitBoxToText','off');
        annotation('textbox',[0.0635000000000001 0.350000000000000 0.0109047619047619 0.0265339966832503],...
        'String',{'D'},'LineStyle','none','FontSize',40,'FitBoxToText','off');
    
        figureDir = ['../../figures/' options.plotOpt.dirSaveFigure];
        if ~exist(figureDir,'dir')
            mkdir(figureDir);
        end
%        savefig(['R:\jwanglab\jwanglabspace\Qiaosheng Zhang\EP recording\Earth 10\5_7_2018\offline_analysis\' options.plotOpt.dirSaveFigure '/' stimulus{s} '_trial' num2str(trials(t)) 'it200.fig']);
        %savefig(['R:\jwanglab\jwanglabspace\Qiaosheng Zhang\EP recording\Earth 10\4_30_2018\offline_analysis\trial_' num2str(trials(t)) '.fig']);
    end
end
end