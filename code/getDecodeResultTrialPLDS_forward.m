% get forward backward PLDS decoding result(the result get while training)
% author: Sile Hu, yuehusile@gmail.com
% date: 2017-3-13
% 2017-08-07 add single trial option
function decodeResult = getDecodeResultTrialPLDS_forward(seq,testseq,model,options)
for j=1:length(seq.stimulus)
    trials = eval(['seq.ntrial' seq.stimulus{j}]);
    st = 1;
    et = trials;
    if isfield(options,['trial' seq.stimulus{j}])
        single_trial = eval(['options.trial' seq.stimulus{j}]);
        if (single_trial>=1 && single_trial<=trials)
            st = single_trial;
            et = single_trial;
        end
    end
    for i=1:length(seq.region)
        numtest = 1;  % index of test trials
        for k=st:et
            params = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').params;']);
            if (options.currentTrial)
                ytest = eval(['seq.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').y;']);
                last_z0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').newseq(end).posterior.xsm(:,1);']);
                if(size(last_z0,1)==2)
                    last_Q0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').newseq(end).posterior.Vsm(1:2,:);']);
                else
                    last_Q0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').newseq(end).posterior.Vsm(1,:);']);
                end
            else
                if  ~options.preTrial &&~isempty(testseq)
                    kk=k;
                    ytest = eval(['testseq.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').y;']);
                    if kk>(options.n_models-1)
                        kk=k-options.n_models+1;
                    else
                        kk=et+k-options.n_models+1;
                    end                     
                    params = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(kk) ').params;']);
                    last_z0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(kk) ').newseq(end).posterior.xsm(:,end);']);
                    last_Q0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(kk) ').newseq(end).posterior.Vsm(end,:);']);
                else
                    kk=k;
                    ytest = eval(['seq.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').y;']);
                    if kk>options.n_models
                        kk=k-options.n_models;
                    else
%                         kk=et+k-options.n_models;
                    end                    
                    params = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(kk) ').params;']);
                    last_z0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(kk) ').newseq(end).posterior.xsm(:,1);']);
                    last_Q0 = eval(['model.' seq.region{i} seq.stimulus{j} '('  num2str(kk) ').newseq(end).posterior.Vsm(1,:);']);
                end
            end
            [online_state, online_var] = PLDS_online_filtering(params,ytest,last_z0,last_Q0,options.binsize);
            if ~options.preTrial && ~isempty(testseq) && ~options.currentTrial && isfield(eval(['testseq.' seq.region{i} seq.stimulus{j} '('  num2str(k) '),']), 'Test')
                test = eval(['testseq.' seq.region{i} seq.stimulus{j} '('  num2str(k) ').Test;']);
                if 1
                    if ~isempty(test)
                        for ii=1:length(test)
                            edges=[test(ii)-options.TPre/options.binsize:test(ii)+options.TPost/options.binsize];
                            if(edges(1)<=0)
                                edges = 1;
                            end
                            eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(numtest) ').x =online_state(:,edges);']);
                            eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(numtest) ').V =online_var(:,edges);']);
                            numtest = numtest + 1;
                        end
                    else
                        eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(numtest) ').x =[];']);
                        eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(numtest) ').V =[];']);
                    end
                else          
                    for tt=1:options.timesneg
                        tmp=size(online_state,2)-((21-(tt-1)*4)/options.binsize);
                        edges=[tmp-options.TPre/options.binsize:tmp+options.TPost/options.binsize];
                        eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(numtest) ').x =online_state(:,edges(2:end));']);
                        eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(numtest) ').V =online_var(:,edges(2:end));']);
                        numtest = numtest + 1;                   
                    end 
                end
            else
                eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(k) ').x =online_state;']);
                eval(['decodeResult.' seq.region{i} seq.stimulus{j} '(' num2str(k) ').V =online_var;']);
            end
        end
    end
end