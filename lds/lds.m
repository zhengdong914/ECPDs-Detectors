%
% function net=lds(X,K,T,cyc,tol);
% 
% Adaptive Linear Dynamical System 
%
% X - N x p data matrix
% K - size of state space (default 2)
% T - length of each sequence (N must evenly divide by T, default T=N)
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% net is a structure consisting of:
%
% A - state transition matrix
% C - observation (output) matrix 
% Q - state noise covariance 
% R - observation noise covariance
% x0 - initial state mean
% P0 - initial state covariance
% Mu - output mean
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM

% add x,P output by hsl 2016-12-2

function net=lds(X,K,T,cyc,tol);

% net is  [A,C,Q,R,x0,P0,Mu,LL,LM]

net=struct('type','lds','A',[],'C',[],'Q',[],'R',[],'x0',[],'P0',[],'Mu',[],'LL',[],'LM',[],'x',[],'P',[]);

p=length(X(1,:));
N=length(X(:,1));

if nargin<5   tol=0.0001; end;
if nargin<4   cyc=100; end;
if nargin<3   T=N; end;
if nargin<2   K=2; end;

Mu=mean(X);
X=X-ones(N,1)*Mu;

if (rem(N,T)~=0)
  disp('Error: Data matrix length must be multiple of sequence length T');
  return;
end;
N=N/T;

  if (K<=p) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize with Factor Analysis
    
    fprintf('\nInitializing with Factor Analysis...\n');
    [L,Ph,LM]=fa(X,K,100,0.001);
    i=length(LM);
    C=L;
    R=Ph;
    Phi=diag(1./R);
    temp1=Phi*L;
    temp2=Phi-temp1*inv(eye(K)+L'*temp1)*temp1';
    temp1=X*temp2*L;
    x0=mean(temp1)';
    Q=cov(temp1);
    P0=Q;
    t1=temp1(1:N*T-1,:);
    t2=temp1(2:N*T,:);
    
    A=inv(t1'*t1+Q)*t1'*t2;
    clear temp1 temp2 Phi t1 t2;
    fprintf('FA log likelihood %f after %i steps\nInitialized.\n',LM(i),i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
    % AR1 model
    LM=[];
    Q=eye(K);
    P0=Q;
    cX=cov(X);
    R=0.5*diag(diag(cX));
    [u,s,v]=svd(cX-R);
    C=u*[sqrtm(s) zeros(p,K-p)];
    C=C+randn(p,K)*max(max(abs(C)))/20;
    beta=C'*inv(C*C'+R);
    t1=X*beta'+0.001*randn(N*T,K);
    x0=mean(t1)';

    t1b=t1(1:N*T-1,:);
    t1c=t1(2:N*T,:);
    A=inv(t1b'*t1b+Q)*t1b'*t1c;
    clear t1 t1b t1c cX;
    R=diag(R);
  end;

lik=0;
LL=[];

Y=reshape(X,T,N,p);
Y=permute(Y,[2 3 1]); % Y is (N,p,T), analogously to X

YY=sum(X.*X)'/(T*N); 

for cycle=1:cyc
  
  % E STEP
  
  oldlik=lik;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3]=kalmansmooth(A,C,Q,R,x0,P0,Y);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  LL=[LL lik];
  fprintf('cycle %g lik %g\n',cycle,lik);
  
  if (cycle<=2)
    likbase=lik;
  elseif (lik<oldlik) 
    fprintf(' violation');
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~isfinite(lik)) 
    fprintf('\n');
    break;
  end;
  fprintf('\n');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % M STEP
  
  % Re-estimate A,C,Q,R,x0,P0;
  
  x0=sum(Xfin(:,:,1),1)'/N; 
  T1=Xfin(:,:,1)-ones(N,1)*x0';
  P0=Pfin(:,:,1)+T1'*T1/N; 
  
  C=YX*inv(Ptsum)/N;
  
  R=YY-diag(C*YX')/(T*N);
    
  A=A1*inv(A2); 
    
  Q=(1/(T-1))*diag(diag((A3-A*(A1')))); 
  if (det(Q)<0) 
    fprintf('Q problem\n');
  end;
  fprintf('cycle %g x0 %g q0%g\n',cycle,x0,P0);
end;

net.A=A;
net.C=C;
net.Q=Q;
net.R=R;
net.x0=x0;
net.P0=P0;
net.Mu=Mu;
net.LL=LL;
net.LM=LM;
net.x = Xfin;
net.P = Pfin;
