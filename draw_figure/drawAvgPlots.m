function drawAvgPlots(lowerBoundAvg,options)

if strcmp(options.plotOpt.stimulus,'all')
    stimulus = options.stimulus;
else
    stimulus = options.plotOpt.stimulus;
end

if strcmp(options.plotOpt.decoders,'all')
    decoders = options.decoders;
else
    decoders = options.plotOpt.decoders;
end

if strcmp(options.plotOpt.region,'all')
    region = options.region;
else
    region = options.plotOpt.region;
end

legends = {};
edges = [-options.TPre:options.binsize:options.TPost];
for s=1:length(stimulus)
    if isfield(options,['trial' options.stimulus{s}])
        single_trial = eval(['options.trial' options.stimulus{s}]);
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
    figure('name',[stimulus{s} ' trial average']);
    
    %-----------------------------------------
    % draw different region in the same figure
    %-----------------------------------------
    % draw baseline
    legends = drawBaseline(options,edges,legends,'');
    legends = drawThresh(options,edges,legends,'');
    for r=1:length(region)
        for d=1:length(decoders)
            edges = [-options.TPre:options.binsize:options.TPost];
            lbavg_st = eval(['lowerBoundAvg.' decoders{d} '.S' stimulus{s} '.' region{r}  ';']);
            if isfield(lbavg_st,'lbavg')
                if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')) && options.currentTrial
                    edges = edges(2:end);
                end
                color = options.plotOpt.zscoreColor(d+r-1);
                lgd0 = plot(edges,lbavg_st.lbavg ,[color '-'],'linewidth',2,'DisplayName',[region{r} ' ' decoders{d} ' avg lb']);
                legends = [legends lgd0];
            end
        end
    end
    hold on;
    lgd0 = plot(edges,lowerBoundAvg.gt ,['r' '-'],'linewidth',2,'DisplayName',['z ground truth']);
    legends = [legends lgd0];
    
    % font, lengend, limit
    ylabel('Z-score','fontsize',24);
    set(gca,'fontsize',20);
    ylim(options.plotOpt.ylim);
    xlabel('time(s)','fontsize',24);
    %subplot(subplotIdx)
    lgd = legend(legends);
    lgd.FontSize = 10;
    legends = [];
    title(options.plotOpt.plotTitle);
    
    figureDir = ['../../figures/' options.plotOpt.dirSaveFigure];
    if ~exist(figureDir,'dir')
        mkdir(figureDir);
    end
    savefig(['../../figures/' options.plotOpt.dirSaveFigure '/' stimulus{s} '_lbavg_pf' num2str(options.pf.version) '.fig']);
end
end