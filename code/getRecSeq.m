% get seq from recording data based on data option
% author: Sile Hu
% date: 2017-3-13
function [seq,testseq] = getRecSeq(allUnits,options)

region = {};
C = length(allUnits);
for i=1:C
    region{i} = allUnits(i).region;
end

region = unique(region);

idx = {};
for i=1:length(region)
    idxTemp = [];
    for cc=1:C
        if strcmp(allUnits(cc).region,region{i})
            idxTemp = [idxTemp cc];
        end
    end
    idx{i} = idxTemp;
end

[seq,testseq] = getStimulusSeq(allUnits,idx,region,options);

if(options.doubleRegion==1)
    for j=1:2
        eval(['seq.dr' seq.stimulus{j} '=seq.' region{1} seq.stimulus{j}]);
        for i=1:length(eval(['seq.' region{1} seq.stimulus{j}]))
            eval(['seq.dr' seq.stimulus{j} '(i).y=[' 'seq.dr' seq.stimulus{j} '(i).y;' 'seq.' region{2} seq.stimulus{j} '(i).y]']);
        end
    end
end