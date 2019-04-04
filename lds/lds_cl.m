%
% function lik=lds_cl(net,X,K,T);
% 
% Linear Dynamical System  - Calculate Likelihood
%
% X - N x p data matrix
% K - size of state space (default 2)
% T - length of each sequence (N must evenly divide by T, default T=N)
% net - a structure consisting of the network parameters
%
% lik - log likelihoods of X under the LDS model specified by net

function lik=lds_cl(net,X,K,T);

% net is  [A,C,Q,R,x0,P0,Mu,LL,LM]
% net=struct('type','lds','A',[],'C',[],'Q',[],'R',[],'x0',[],'P0',[],'Mu',[],'LL',[],'LM',[]);

p=length(X(1,:));
N=length(X(:,1));

if nargin<4   T=N; end;
if nargin<3   K=2; end;

A=net.A; C=net.C; Q=net.Q;
R=net.R; x0=net.x0; P0=net.P0; Mu=net.Mu; 

X=X-ones(N,1)*Mu;

if (rem(N,T)~=0)
  disp('Error: Data matrix length must be multiple of sequence length T');  return; 
end;
N=N/T;

lik=0;

Y=reshape(X,T,N,p);
Y=permute(Y,[2 3 1]); % Y is (N,p,T)

[lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3]=kalmansmooth(A,C,Q,R,x0,P0,Y);

  

