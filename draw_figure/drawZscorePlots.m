% draw z-score plots based on different options
% one figure for each trial
% author: Sile Hu
% date: 2017-3-13
function drawZscorePlots(z_scores,seq,options)

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
        trials = 1:eval(['seq.ntrial' stimulus{s}]);
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
%         order=[7 10 4 1 5 6 2 8 9 3 11 12];
%         seq.PAINSIM(t).y=seq.PAINSIM(t).y(order,:);
%         figure('name',[stimulus{s} ' trial' num2str(trials(t))]);
        figure;
        if options.plotOpt.regionsInOne
            %-----------------------------------------
            % draw different region in the same figure
            %-----------------------------------------
            % draw raster of all units
            subplotIdx = '211';
            seqall = eval(['seq.' region{1} stimulus{s} '(trials(t))']);
            titleStr = ['stimulus=' stimulus{s} mw ',trial' num2str(trials(t)) ' ' region{1} '=1:' size(seqall.y,1)];
            for r=2:length(region)
                seqall.y = [seqall.y;eval(['seq.' region{r} stimulus{s} '(trials(t)).y'])];
            end
            drawRasterTrialNew(seqall,titleStr,edges,subplotIdx);
            
            % lower plot
            subplotIdx = '212';
            % draw baseline
           
            legends = drawBaseline(options,edges,legends,subplotIdx);
            legends = drawThresh(options,edges,legends,subplotIdx);
           
            for r=1:length(region)
                for d=1:length(decoders)
                    z_score = eval(['z_scores.' decoders{d} '.' region{r} stimulus{s} '(trials(t));']);
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
        else
            %-------------------------------------------
            % draw different region in different figures
            %-------------------------------------------
            for r=1:length(region)
                % draw raster
                subplotIdx = ['2' num2str(length(region)) num2str(r)];
                titleStr = ['region=' region{r} ',stimulus=' stimulus{s} mw ',trial' num2str(trials(t))];
                drawRasterTrialNew(eval(['seq.' region{r} stimulus{s} '(trials(t))']),titleStr,edges,subplotIdx)
                
                subplotIdx = ['2' num2str(length(region)) num2str(r+length(region))];
                % draw baseline
                legends = drawBaseline(options,edges,legends,subplotIdx);
                legends = drawThresh(options,edges,legends,subplotIdx);
                for d=1:length(decoders)
                    z_score = eval(['z_scores.' decoders{d} '.' region{r} stimulus{s} '(trials(t));']);
                    l_z = length(z_score.zscore)-1;
                    decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                    if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')) && options.currentTrial
                        decoder_edges = decoder_edges(2:end);
                    end
                    color = options.plotOpt.zscoreColor(d);
                    legends = drawZscoreTrial(z_score,decoder_edges,legends,decoders{d},options,color,subplotIdx);
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
                %withdrawEvent = eval(['seq.withdrawal' stimulus{s} ';']);
                legends = drawEvent(withdrawEvent,options,stimulus{s},legends,trials(t),subplotIdx);
                %legend = drawCross();
                ylabel('Z-score','fontsize',24);
                set(gca,'fontsize',20);
                ylim(options.plotOpt.ylim);
                xlim([-options.TPre,options.TPost])
                xlabel('time(s)','fontsize',24);
                subplot(subplotIdx)
                lgd = legend(legends);
                lgd.FontSize = 10;
                legends = [];
                title(options.plotOpt.plotTitle);
            end
        end
%         figureDir = ['../../figures/' options.plotOpt.dirSaveFigure];
%         if ~exist(figureDir,'dir')
%             mkdir(figureDir);
%         end
        %savefig(['../../figures/' options.plotOpt.dirSaveFigure '/' stimulus{s} '_trial' num2str(trials(t)) '.fig']);
    end
end
end