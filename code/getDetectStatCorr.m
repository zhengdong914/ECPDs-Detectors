% 2017-7-27 Sile Hu, yuehusile@gmail.com
% description:
% get detection infomation and statistics of all the
% decoders/regions/stimulus/trials configured in 'options', with
% detection range specified, based on z_scores
%
% input:
%         z_scores       - z-score structure array
%         detectionRange - range of detection (index)
%         options        - data options
% output:
%         detectStat - detection information and statistics
function detectStat = getDetectStatCorr(z_scores,detectionRange,options)
for d=1:length(options.decoders)
    for r=1:length(options.region)
        for s=1:1%length(options.stimulus)
            if options.currentTrial || options.preTrial
                trials = eval(['options.ntrial' options.stimulus{s}]);
            else
                trials = eval(['length(z_scores{1}.' options.decoders{1} '.' options.region{r} options.stimulus{s} ')']);
            end
            st = 1;
            et = trials;
            if isfield(options,['trial' options.stimulus{s}])
                single_trial = eval(['options.trial' options.stimulus{s}]);
                if (single_trial>=1 && single_trial<=trials)
                    st = single_trial;
                    et = single_trial;
                end
            end
            eval(['detectStat.' options.decoders{d} '.S' options.stimulus{s} '.' options.region{r} '.ntrials=trials']);
            totalDetect = 0;
            for t = st:et
                z_score=cell(size(z_scores));
                for i=1:length(z_score)
                    z_score{i} = eval(['z_scores{i}.' options.decoders{d} '.' options.region{r} options.stimulus{s} '(' num2str(t) ')']);
                end
                % getDetectInfo tobeDone
                eval(['detectStat.' options.decoders{d} '.S' options.stimulus{s} '.' options.region{r} '.allTrials(' num2str(t) ')=getDetectInfoCorr(z_score,detectionRange,options.threshold, options, options.TPre,options.areaThr,r);']);
                totalDetect = totalDetect + eval(['detectStat.' options.decoders{d} '.S' options.stimulus{s} '.' options.region{r} '.allTrials(' num2str(t) ').isDetected']);
            end
            eval(['detectStat.' options.decoders{d} '.S' options.stimulus{s} '.' options.region{r} '.totalDetect=totalDetect;']);
        end
    end
end

