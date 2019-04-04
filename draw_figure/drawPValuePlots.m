% draw z-score plots based on different options
% one figure for each trial
% author: Sile Hu
% date: 2017-3-13
function drawPValuePlots(z_scores,seq,options)

if(options.doubleRegion~=1)
    % load('AR1sim2.mat');
    load('AR1result.mat');
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
    for t=st:et
        figure('name',[stimulus{s} ' trial' num2str(trials(t))]);
        for r=1:length(region)
            edges = [-options.TPre:options.binsize:options.TPost];
            % draw raster
            subplotIdx = ['2' num2str(length(region)) num2str(r)];
            if strcmp(stimulus{s},'50') || strcmp(stimulus{s},'150') || strcmp(stimulus{s},'250')
                mw = 'mW';
            else
                mw = [];
            end
            titleStr = ['region=' region{r} ',stimulus=' stimulus{s} mw ',trial' num2str(trials(t))];
            drawRasterTrialNew(eval(['seq.' region{r} stimulus{s} '(trials(t))']),titleStr,edges,subplotIdx)
            
            subplotIdx = ['2' num2str(length(region)) num2str(r+length(region))];
            % TODO: draw baseline
            legends = drawBaseline(options,edges,legends,subplotIdx);
            legends = drawPValueThresh(options,edges,legends,subplotIdx);
            for d=2:2%1:length(decoders)
                z_score = eval(['z_scores.' decoders{d} '.' region{r} stimulus{s} '(trials(t));']);
                l_z = length(z_score.pValue) - 1;
                decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')||strcmp(decoders{d},'GT')) && options.currentTrial
                    decoder_edges = decoder_edges + options.binsize;
                end
                color = options.plotOpt.zscoreColor(d);
                legends = drawPValueTrialAR2(z_score,decoder_edges,legends,decoders{d},options,color,subplotIdx);
            end
            for d=2:2
                z_score1 = eval(['AR1result.' decoders{d} '.' region{r} stimulus{s} '(trials(t));']);
                l_z = length(z_score1.pValue) - 1;
                decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')) && options.currentTrial
                     decoder_edges = decoder_edges + options.binsize;
                end
                color = options.plotOpt.zscoreColor(d+3);
                legends = drawPValueTrial(z_score1,decoder_edges,legends,decoders{d},options,color,subplotIdx);
            end
            withdrawEvent = eval(['seq.withdraw' stimulus{s} ';']);
            legends = drawEvent(withdrawEvent,options,stimulus{s},legends,trials(t),subplotIdx);
            %legend = drawCross();
            ylabel('P-Value','fontsize',16);
            set(gca,'fontsize',16);
            ylim([1e-4 1]);
            xlim([-options.TPre,options.TPost])
            xlabel('time(s)','fontsize',16);
            subplot(subplotIdx)
            lgd = legend(legends);
            lgd.FontSize = 10;
            legends = [];
        end

%         savefig(['../../figures/sim_a3__2_u10_iu4_t10_l201_test2intensity/' stimulus{s} '_trial' num2str(t) '.fig']);
%         savefig(['../../figures/Neptune3_6_28_2017_15050_online/' stimulus{s} '_trial' num2str(t) '.fig']);

    end
end