% 2017-7-28 Sile Hu, yuehusile@gmail.com
% description:
% a cross threshold is defined by cross threshold by z-score or confidence intervals,
% at least two consecutive points, to suppress false positives
%
% input:
%         z_score - z-score structure, including upper bound and lower
%         bound if using PLDS or LDS
%         range   - range of detection
%         Thr     - threshold
%         binsize - bin size
%         tpre    - time before alignment event (laser on etc.)
% output:
%         firstDetect - first cross threshold (detection) position
%         isDetected - 1=detected,0=not detected
function detectInfo = getDetectInfo(z_score, range, Thr, binsize, tpre)

z = z_score.zscore(range);
if isfield(z_score,'lowerBound')
    findRangeLower = z_score.lowerBound(range);
else
    findRangeLower = z_score.zscore(range);
end
if isfield(z_score,'upperBound')
    findRangeHigher = z_score.upperBound(range);
else
    findRangeHigher = z_score.zscore(range);
end
if isfield(z_score,'ccvalue')
    findRangeCC = z_score.ccvalue(1,range,:);
    findRangeTh = z_score.ccthresh;
end

L = [];
H = [];
for i=2:length(findRangeLower)-1
    if (findRangeLower(i+1)>Thr && findRangeLower(i)>Thr && findRangeLower(i-1)<=Thr)
        L = [L i];
    end
    if (findRangeHigher(i+1)<-Thr  && findRangeHigher(i)<-Thr && findRangeHigher(i-1)>=-Thr)
        H = [H i];
    end
end
if isempty(L) && findRangeLower(1)>Thr && findRangeLower(2)>Thr
    L = 1;
end
if isempty(H) && findRangeHigher(1)<-Thr && findRangeHigher(2)<-Thr
    H = 1;
end
isDetected = 0;
if isfield(z_score,'ccvalue')
    CCtemp=zeros(size(findRangeCC));
    area=CCtemp;
    for i=1:size(findRangeCC,3)
        for j=1:size(findRangeCC,2)
            if(findRangeCC(1,j,i)>findRangeTh(1,2,i) || findRangeCC(1,j,i)< findRangeTh(1,1,i))
                CCtemp(1,j,i)= max(findRangeCC(1,j,i)-findRangeTh(1,2,i),findRangeTh(1,1,i)-findRangeCC(1,j,i));
                if j==1
                    area(1,j,i)= CCtemp(1,j,i);
                else
                    area(1,j,i)= area(1,j-1,i)+ CCtemp(1,j,i);
                end
            else
                CCtemp(1,j,i)=0;
                area(1,j,i)=0;
            end
        end
    end
    cc_area=permute(max(area,[],2),[3,1,2]);
else
    cc_area=[];
end

firstDetect = length(z);
%display('offline');
if ~isempty(L)
    firstDetect = min(L(1),firstDetect);
    %display(L);
end
if ~isempty(H)
    firstDetect = min(H(1),firstDetect);
    %display(H);
end

if firstDetect ~= length(z)
    isDetected = 1;
end

% if(isDetected==1)
% %    detect_range = [max(firstDetect-10,1):min(firstDetect+10,length(z))];
% %    if(sum(findRangeCC(detect_range)>findRangeTh(2))||sum(findRangeCC(detect_range)<findRangeTh(1)))
%     if isfield(z_score,'ccvalue')
%        if(findRangeCC(firstDetect)>findRangeTh(2)||findRangeCC(firstDetect)<findRangeTh(1))
%            cc_isDetected = 1;   
%        else
%            cc_isDetected = 0;
%        end
%     end
% end

detectInfo.firstDetect = firstDetect + range(1) - 1;
firstDetectTime = detectInfo.firstDetect * binsize - tpre;
detectInfo.firstDetectTime = firstDetectTime;
detectInfo.isDetected = isDetected;
% detectInfo.cc_isDetected = cc_isDetected;
detectInfo.cc_area = cc_area;