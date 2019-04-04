% draw raster for a trial on the plot handled by h
% author: Sile Hu
% date: 2017-3-13
% input: seq - spike sequence
%        sub   - subplot index str
function drawRasterY(y,edges)
% -- raster --
figure,
imagesc(edges, [1:size(y,1)],y);axis xy
colormap(flipud(bone));
ylabel('Cell','fontsize',16);
xlabel('time(s)','fontsize',16);
set(gca,'fontsize',16)