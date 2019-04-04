% online decoding based on trained model, trial structure
% author: Sile Hu
% date: 2017-3-13
function onlineDecodeResult = onlineDecodeTrial(seq,testseq,model,options)
for i=1:length(options.decoders)
    decoder = options.decoders{i};
    if (options.doubleRegion==1)
        decodeResult = eval(['getDecodeResultTrial' decoder '_dr(seq,testseq,model,options);']);
    else
        if strcmp(decoder,'GT')
            decodeResult = eval(['getDecodeResultTrial' decoder '(testseq,model,options);']);
        else
            decodeResult = eval(['getDecodeResultTrial' decoder '(seq,testseq,model,options);']);
        end
    end
    eval(['onlineDecodeResult.' decoder ' = decodeResult;']);
end