% function [L,Ph,LL]=fa(X,K,cyc,tol);
% 
% Maximum Likelihood Factor Analysis using EM
%
% X - data matrix
% K - number of factors
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% L - factor loadings 
% Ph - diagonal uniquenesses matrix
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM 
%
% Use ffa.m, which is considerably faster!


function [L,Ph,LL]=fa(X,K,cyc,tol);

if nargin<4  tol=0.0001; end;
if nargin<3  cyc=100; end;

N=length(X(:,1));
D=length(X(1,:));
tiny=exp(-700);

X=X-ones(N,1)*mean(X);
XX=X'*X/N;

cX=cov(X);
scale=det(cX)^(1/D);
L=randn(D,K)*sqrt(scale/K);

%for test
%L = [0.1846   -0.4657   -0.0282    0.6299   -0.1941   -0.5298    0.0624]';
Ph=diag(cX);

I=eye(K);

lik=0; LL=[];

const=-D/2*log(2*pi);

for i=1:cyc;

  %%%% E Step %%%%
  Phd=diag(1./Ph);
  LP=Phd*L;
  MM=Phd-LP*inv(I+L'*LP)*LP';
  dM=sqrt(det(MM));
  beta=L'*MM;
  XM=X*MM;
  EZ=XM*L;
  EZZ=I-beta*L +beta*XX*beta';

  %%%% Compute log likelihood %%%%
  
  oldlik=lik;
  H=const+log(dM)-0.5*sum(XM.*X,2);
  lik=sum(H); 
  fprintf('cycle %i lik %g \n',i,lik);
  LL=[LL lik];
  
  %%%% M Step %%%%

  L=X'*EZ*inv(EZZ)/N;
  Ph=diag(XX-L*EZ'*X/N);

  if (i<=2)    
    likbase=lik;
  elseif (lik<oldlik)     
    disp('violation');
  elseif ((lik-likbase)<(1+tol)*(oldlik-likbase)|~isfinite(lik))  
    break;
  end;

end
