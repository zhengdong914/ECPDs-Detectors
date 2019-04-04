% get stimulus seq based on data option
% author: Sile Hu
% date: 2017-3-13
%--------------------------------------------------------------------------
% modification note:
% created: 2017-3-13 v0
% modified: 
% 2017-7-19 v0 Sile Hu, fix non-cell adaptation
% 20170817  v0 Sile Hu, using strsplit() to handle non-cell adaptation,
%                       adapation to old recording files
%-------------------------------------------------------------------------- 
function [seq,testseq] = getStimulusSeq(allUnits,idx,region,options)
if (options.PP==1)
    seq.stimulus = {'PP','VF'};
    seq.teststimulus = {'VF'};
else
    seq.stimulus = allUnits(1).intensity;
    if ~iscell(seq.stimulus)
        seq.stimulus = strrep(seq.stimulus,'50mw','50');
        seq.stimulus = strsplit(seq.stimulus);
    end
    if(size(seq.stimulus,2)==1)
        seq.teststimulus = seq.stimulus;
    else
        seq.teststimulus = {seq.stimulus{2}};
    end
end
seq.region = region;
if ~iscell(seq.region)
    seq.region = strsplit(seq.region);
end
for j=1:length(seq.stimulus)
    for i=1:length(seq.region)
        eval(['[seq.' seq.region{i} seq.stimulus{j} ',testseq.' seq.region{i} seq.stimulus{j} ']= getSeqSpecifiedNew(allUnits, idx{' num2str(i) '}, ''' seq.stimulus{j} '''' ',' '''' seq.teststimulus{1} ''' , options);']);
    end
    eval(['seq.ntrial' seq.stimulus{j} '= length(allUnits(1).laserOn' seq.stimulus{j} ');']);
    if isfield(allUnits(1),['withdraw' seq.stimulus{j} ])
        eval(['seq.withdraw' seq.stimulus{j} ' = getRelativeEventTime(allUnits(1).withdraw' seq.stimulus{j} ',allUnits(1).laserOn' seq.stimulus{j} ',options);']);
    else
        if ~strcmp(seq.stimulus{j}, '50') && ~strcmp(seq.stimulus{j}, 'VF') 
            if isfield(allUnits(1),'recFlick')
                eval(['seq.withdraw' seq.stimulus{j} ' = getRelativeEventTime(allUnits(1).recFlick, allUnits(1).laserOn' seq.stimulus{j} ',options);']);
            else
                eval(['seq.withdraw' seq.stimulus{j} ' = [];']);
            end
        else
            eval(['seq.withdraw' seq.stimulus{j} ' = [];']);
        end
    end
end