function lgd = drawPValueTrialAR2(z_score,edges,lgd,decoder,options,color,subplotIdx)
subplot(subplotIdx);
decoder = strrep(decoder,'_','\_');
lgd0 = plot(edges,z_score.pValue ,[color '-'],'linewidth',2,'DisplayName',['AR(2)-' decoder ' p-value']);
set(gca,'YLim',[1e-4 1]);
set(gca,'YScale','log');
lgd = [lgd lgd0];
if isfield(z_score,'upperBound')
hold on;
% lgd1 = plot(edges,z_score.upperBound ,[color '--'],'linewidth',1,'DisplayName',[decoder ' confident interval']);
lgd1=[];
hold on;
plot(edges,z_score.lowerBound ,[color '--'],'linewidth',1);
lgd = [lgd lgd1];
end