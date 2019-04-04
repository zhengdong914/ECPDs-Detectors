% draw z-score plots based on different options
% one figure for each trial
% author: Sile Hu
% date: 2017-3-13
function drawZscorePlots_SIM(z_scores1,z_scores2,seq,options)

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
        if options.currentTrial || strcmp(seq.stimulus{s},'SIM')
            trials = 1:eval(['seq.ntrial' stimulus{s}]);
        else
            trials = 1:eval(['seq.ntrial' options.teststimulus{s} '-1']);
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
        figure('name',[stimulus{s} ' trial' num2str(trials(t))]);

        subplotIdx = '211';
        seqall = eval(['seq.' seq.region{1} stimulus{s} '(trials(t))']);
        if options.currentTrial || strcmp(seq.stimulus{s},'SIM')
            titleStr = ['stimulus=' stimulus{s} mw ',trial' num2str(trials(t)) ' ' region{1} '=1:' size(seqall.y,1)];
        else
            titleStr = ['stimulus=' options.teststimulus{s} mw ',trial' num2str(trials(t)) ' ' region{1} '=1:' size(seqall.y,1)];
        end
        if options.currentTrial || strcmp(seq.stimulus{s},'SIM')
            if(length(seq.region)==1)
                seqall.y = [seqall.y;eval(['seq.' seq.region{1} stimulus{s} '(trials(t)).y'])];
            else
                seqall.y = [seqall.y;eval(['seq.' seq.region{2} stimulus{s} '(trials(t)).y'])];
            end
        else
            if(length(seq.region)==1)
                seqall.y = [seqall.y;eval(['seq.' seq.region{1} options.teststimulus{s} '(trials(t)).y'])];
            else
                seqall.y = [seqall.y;eval(['seq.' seq.region{2} options.teststimulus{s} '(trials(t)).y'])];
            end
        end
        drawRasterTrialNew(seqall,titleStr,edges,subplotIdx);            
%         subplotIdx = '312';         
% 
%         for d=1:length(decoders)
%             if ~(strcmp(decoders{d},'GT'))
%                 cct = eval(['z_scores2.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t)).ccvalue;']);                
%                 options.C_thresh=eval(['z_scores2.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t)).ccthresh;']);
%                 inx1 = find(cct(options.TPre/options.binsize:(options.TPre+2)/options.binsize)>options.C_thresh(2) | cct(options.TPre/options.binsize:(options.TPre+2)/options.binsize)<options.C_thresh(1));
%                 inx2 = find(cct((options.TPre+2)/options.binsize:(options.TPre+10)/options.binsize)>options.C_thresh(1) & cct((options.TPre+2)/options.binsize:(options.TPre+10)/options.binsize)<options.C_thresh(2));
%                 if ~isempty(inx1)&&~isempty(inx2)
%                     options.crossvalue=[(inx1(1)-1)*options.binsize (inx2(1)-1)*options.binsize+2];
%                 else
%                     options.crossvalue=[];
%                 end
%                 drawCorrelationPlots(cct,options,edges,legends,subplotIdx);
%             end
%         end
%         ylabel('C-Value','fontsize',24);
%         set(gca,'fontsize',20);
%         xlim([-options.TPre,options.TPost])
%         xlabel('time(s)','fontsize',24);
%         subplot(subplotIdx)
%         lgd = legend(legends);
%         lgd.FontSize = 10;
%         legends = [];
%         if ~isempty(options.crossvalue)
%             title(['Cross Thresh  On:' num2str(options.crossvalue(1)) 's   Off:' num2str(options.crossvalue(2)) 's     Confidential Interval:' num2str(options.confidenceInterval)]);
%         end
        % lower plot
        subplotIdx = '212';
        % draw baseline
        legends = drawBaseline(options,edges,legends,subplotIdx);
        legends = drawThresh(options,edges,legends,subplotIdx);

        for r=1:length(region)
            if (r==2)
                for d=1:length(decoders)
                    if ~(strcmp(decoders{d},'GT'))
                        if (length(seq.region)==1)
                            z_score = eval(['z_scores2.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t));']);
                        else
                            z_score = eval(['z_scores2.' decoders{d} '.' seq.region{2} stimulus{s} '(trials(t));']);
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
                    z_score = eval(['z_scores1.' decoders{d} '.' seq.region{1} stimulus{s} '(trials(t));']);
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
        if isfield(seq, ['withdrawal' stimulus{s} ])
            withdrawEvent = eval(['seq.withdrawal' stimulus{s} ';']);
        else
            if isfield(seq, ['withdraw' stimulus{s} ])
                withdrawEvent = eval(['seq.withdraw' stimulus{s} ';']);
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
        xlabel('time(s)','fontsize',24);
        xlim([-options.TPre,options.TPost])
        subplot(subplotIdx)
        lgd = legend(legends);
        lgd.FontSize = 10;
        legends = [];
        title(options.plotOpt.plotTitle);

        figureDir = ['../../figures/' options.plotOpt.dirSaveFigure];
        if ~exist(figureDir,'dir')
            mkdir(figureDir);
        end
%         savefig(['../../figures/' options.plotOpt.dirSaveFigure '/' stimulus{s} '_trial' num2str(trials(t)) 'it200.fig']);
    end
end
end