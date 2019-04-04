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
function detectInfo = getDetectInfoCorr(z_scores, range, Thr, options, tpre, areaThr, region)
findRangeLower=cell(1,length(z_scores));
findRangeHigher=cell(1,length(z_scores));
cc_area=cell(1,length(z_scores));
L = cell(1,length(z_scores));
H = cell(1,length(z_scores));
CC = cell(1,length(z_scores));
if isfield(z_scores{1},'ccvalue')
    CCtemp=zeros(length(z_scores),size(range,2),size(z_scores{1}.ccvalue,3));
    area=CCtemp;
end
for m=1:length(z_scores)
    z_score = z_scores{m};
    z = z_score.zscore(range);
    if isfield(z_score,'lowerBound')
        findRangeLower{m} = z_score.lowerBound(range);
    else
        findRangeLower{m} = z_score.zscore(range);
    end
    if isfield(z_score,'upperBound')
        findRangeHigher{m} = z_score.upperBound(range);
    else
        findRangeHigher{m} = z_score.zscore(range);
    end
    if isfield(z_score,'ccvalue')
        findRangeCC = z_score.ccvalue(1,range,:);
        findRangeTh = z_score.ccthresh;
    end
    for i=2:length(range)-1
        if (findRangeLower{m}(i+1)>Thr && findRangeLower{m}(i)>Thr && findRangeLower{m}(i-1)<=Thr)
            L{m} = [L{m} i];
        end
        if (findRangeHigher{m}(i+1)<-Thr  && findRangeHigher{m}(i)<-Thr && findRangeHigher{m}(i-1)>=-Thr)
            H{m} = [H{m} i];
        end
    end
    if isempty(L{m}) && findRangeLower{m}(1)>Thr && findRangeLower{m}(2)>Thr
        L{m} = 1;
    end
    if isempty(H{m}) && findRangeHigher{m}(1)<-Thr && findRangeHigher{m}(2)<-Thr
        H{m} = 1;
    end

    if isfield(z_score,'ccvalue')
        for i=1:size(findRangeCC,3)
            for j=1:size(findRangeCC,2)
                if(findRangeCC(1,j,i)>findRangeTh(1,2,i) || findRangeCC(1,j,i)< findRangeTh(1,1,i))
                    CCtemp(m,j,i)= max(findRangeCC(1,j,i)-findRangeTh(1,2,i),findRangeTh(1,1,i)-findRangeCC(1,j,i));
                    if j==1
                        area(m,j,i)= CCtemp(m,j,i);
                    else
                        area(m,j,i)= area(m,j-1,i)+ CCtemp(m,j,i);
                    end
                else
                    CCtemp(m,j,i)=0;
                    area(m,j,i)=0;
                end
            end
        end
%         cc_area{m}=permute(max(area,[],2),[3,1,2]);
%     else
%         cc_area{m}=[];
        CC{m} = find(area(m,:,3)>=areaThr);
    end
    
end

isDetected = 0;
if(region==1)
    det_window=options.ACC_w;
else
    det_window=options.S1_w;
end
for m=1:length(z_scores)
    LH=[L{m},H{m}];
    if ~isempty(LH)
        tmp=1:length(z_scores);
        tmp(m)=[];
        for l=1:length(LH)
            count=1;
            for n=1:length(tmp)
                tmp1=[L{tmp(n)},H{tmp(n)}];
                if ~isempty(find(tmp1<=(LH(l)+det_window)&tmp1>=LH(l), 1))
                    count=count+1;
                end
            end
            if(count>length(z_scores)/2)
                isDetected=1;
            end
            if isDetected
                break;
            end
        end
        if isDetected
           break;
        end
    end
end
cc_isDetected=0;
cc_det_window=options.cc;
for m=1:length(z_scores)
    if ~isempty(CC{m})
        tmp=1:length(z_scores);
        tmp(m)=[];
        for cc=1:length(CC{m})
            count=1;
            for n=1:length(tmp)
                if ~isempty(find(CC{tmp(n)}<=(CC{m}(cc)+cc_det_window)&CC{tmp(n)}>=CC{m}(cc), 1))
                    count=count+1;
                end
            end
            if(count>length(z_scores)/2)
                cc_isDetected=1;
            end
            if cc_isDetected
                break;
            end
        end
        if cc_isDetected
           break;
        end
    end
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

detectInfo.isDetected = isDetected;
detectInfo.cc_isDetected = cc_isDetected;
% detectInfo.cc_area = cc_area;