% compute event(withdraw) time - stimulus time
% author: Sile Hu
% date: 2017-3-11
% output: relaitveEvent - 2*length(event), [idx event time - stimulus time]
%--------------------------------------------------------------------------
% modification note:
% created: 2017-3-11 v0
% modified: 
% 2017-7-19 v0 Sile Hu, fix redundent event time bug
%-------------------------------------------------------------------------- 
function relaitveEvent = getRelativeEventTime(event,stimulus,options)
relaitveEvent = [];
for i=1:length(event)
    evDiff = abs(stimulus - event(i));
    idx = find(evDiff == min(evDiff));
    if (~isempty(relaitveEvent))
        while (~isempty(find(relaitveEvent(:,1)==idx)))
            idx = idx + 1;
        end
    end
    relaitveEvent = [relaitveEvent; [idx event(i)-stimulus(idx)]];
end