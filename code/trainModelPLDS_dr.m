% train PLDS model based on data options 
% author: Sile Hu
% date: 2017-3-11
function model = trainModelPLDS_dr(seq,options)
if (options.sim==1)
    eval(['model.dr' seq(1).stimulus{1} '=trainModelSingleSeqPLDS_dr(seq(1).' seq(1).region{1} seq(1).stimulus{1} ', seq(2).' seq(2).region{1} seq(2).stimulus{1} ', options )']);
else
    for j=1:length(seq.stimulus)
        if(options.doubleRegion==1)
            eval(['model.dr' seq.stimulus{j} '=trainModelSingleSeqPLDS_dr(seq.' seq.region{1} seq.stimulus{j} ', seq.' seq.region{2} seq.stimulus{j} ', options )']);
        else
            for i=1:length(seq.region)
                eval(['model.' seq.region{i} seq.stimulus{j} '=trainModelSingleSeqPLDS(seq.' seq.region{i} seq.stimulus{j} ',options )']);
            end
        end
    end
end