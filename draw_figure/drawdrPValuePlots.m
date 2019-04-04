% draw z-score plots based on different options
% one figure for each trial
% author: Sile Hu
% date: 2017-3-13
function drawdrPValuePlots(z_scores,seq,options)
if(options.sim==1)
    load('AR1s0.mat');
    load('AR1s2.mat');
else
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
        edges = [-options.TPre:options.binsize:options.TPost];
        % draw raster
        subplotIdx = '211';
        if strcmp(stimulus{s},'50') || strcmp(stimulus{s},'150') || strcmp(stimulus{s},'250')
            mw = 'mW';
        else
            mw = [];
        end
        titleStr = ['region=ACC&S1' ',stimulus=' stimulus{s} mw ',trial' num2str(trials(t))];
        drawRasterTrialNew(eval(['seq.dr' stimulus{s} '(trials(t))']),titleStr,edges,subplotIdx)

        subplotIdx = '212';
        % TODO: draw baseline
        legends = drawBaseline(options,edges,legends,subplotIdx);
        legends = drawPValueThresh(options,edges,legends,subplotIdx);
        for d=1:length(decoders)
            z_score = eval(['z_scores.' decoders{d} '.dr' stimulus{s} '(trials(t));']);
            l_z = length(z_score.pValue) - 1;
            decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
            if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')||strcmp(decoders{d},'GT')) && options.currentTrial
                decoder_edges = decoder_edges + options.binsize;
            end
            color = options.plotOpt.zscoreColor(d);
            legends = drawPValueTrialAR2(z_score,decoder_edges,legends,decoders{d},options,color,subplotIdx);
        end
        if(options.sim==1)
            for d=1:2
                if(d==1)
                    z_score1 = eval(['AR1s0.' decoders{1} '.' region{1} stimulus{s} '(trials(t));']);
                else
                    z_score1 = eval(['AR1s2.' decoders{1} '.' region{1} stimulus{s} '(trials(t));']);
                end
                l_z = length(z_score1.pValue) - 1;
                decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                if ~(strcmp(decoders{1},'PLDS_FB')||strcmp(decoders{1},'ModelFree')) && options.currentTrial
                     decoder_edges = decoder_edges + options.binsize;
                end
                color = options.plotOpt.zscoreColor(d+1);
                legends = drawPValueTrial(z_score1,decoder_edges,legends,decoders{1},region{1},color,subplotIdx);
            end            
        else
            for d=1:1
                for j=1:length(region)
                    z_score1 = eval(['AR1result.' decoders{d} '.' region{j} stimulus{s} '(trials(t));']);
                    l_z = length(z_score1.pValue) - 1;
                    decoder_edges = [-options.TPre:options.binsize:l_z*options.binsize - options.TPre];
                    if ~(strcmp(decoders{d},'PLDS_FB')||strcmp(decoders{d},'ModelFree')) && options.currentTrial
                         decoder_edges = decoder_edges + options.binsize;
                    end
                    color = options.plotOpt.zscoreColor(j+1);
                    legends = drawPValueTrial(z_score1,decoder_edges,legends,decoders{d},region{j},color,subplotIdx);
                end
            end
        end
        
        withdrawEvent = eval(['seq.withdraw' stimulus{s} ';']);
        legends = drawEvent(withdrawEvent,options,stimulus{s},legends,trials(t),subplotIdx);
        %legend = drawCross();
        ylabel('P-Value','fontsize',16);
        set(gca,'fontsize',16);
        ylim([1e-10 1]);
        xlim([-options.TPre,options.TPost])
        xlabel('time(s)','fontsize',16);
        subplot(subplotIdx)
        lgd = legend(legends);
        lgd.FontSize = 10;
        legends = [];

%         savefig(['../../figures/allUnits Mars_50_250_02092017/' stimulus{s} '_trial' num2str(t) '.fig']);
%         savefig(['../../figures/sim_a3__2_u10_iu4_t10_l201_test2intensity/' stimulus{s} '_trial' num2str(t) '.fig']);
%         savefig(['../../figures/Neptune3_6_28_2017_15050_online/' stimulus{s} '_trial' num2str(t) '.fig']);

    end
end