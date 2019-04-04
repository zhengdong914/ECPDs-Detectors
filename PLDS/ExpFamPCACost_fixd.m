function [f df] = ExpFamPCACost_fixd(CXd,Y,xDim,lambda,s,CposConstrain,fix_d)
%
% [f df] = ExpFamPCACost(CXd,Y,xDim,lambda)
%
% (c) L Buesing 2014

[yDim T] = size(Y);


d  = fix_d;
C  = CXd(1:yDim,:);
if CposConstrain; C = exp(C); end;
X  = CXd(yDim+1:end,:)';

nu = bsxfun(@plus,C*X+s,d);
Yhat = exp(nu);

f = sum(vec(-Y.*nu+Yhat))+lambda/2*(norm(C,'fro')^2+norm(X,'fro'));

YhatmY = Yhat-Y;

gX = C'*YhatmY+lambda*X;
gC = YhatmY*X'+lambda*C;
if CposConstrain; gC = gC.*C; end
gd = sum(YhatmY,2);

%df = [vec([gC;gX']);gd];
df = [vec([gC;gX'])];

