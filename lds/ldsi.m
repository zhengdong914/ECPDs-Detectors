% function [A,B,C,Q,R,x0,P0,Mu,LL]=ldsi(X,D,K,T,cyc,tol,fast);
% 
% Adaptive Linear Dynamical System with Inputs
%
% X - N x p data matrix
% D - number of inputs (rest assumed to be outputs; default 0)
% K - size of state space (default 2)
% T - length of each sequence (N must evenly divide by T, default T=N)
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
% fast - flag for using faster algorithm that detects approx convergence of
%        Riccati equations (fast==1 use faster algorithm, default 0)
%
% A - state transition matrix
% B - input matrix
% C - observation (output) matrix 
% Q - state noise covariance 
% R - observation noise covariance
% x0 - initial state covariance
% P0 - initial state mean
% Mu - input/output mean
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM

function [A,B,C,Q,R,x0,P0,Mu,LL]=ldsi(X,D,K,T,cyc,tol,fast);

p=length(X(1,:));
N=length(X(:,1));
tiny=exp(-700);
problem=0;

if nargin<7   fast=0; end;
if nargin<6   tol=0.0001; end;
if nargin<5   cyc=100; end;
if nargin<4   T=N; end;
if nargin<3   K=2; end;
if nargin<2   D=0; end;

if (D==0)
  [A,C,Q,R,x0,P0,Mu,LL]=ldso(X,K,T,cyc,tol,fast); % call old LDS code
  B=[];
else
  Mu=mean(X);
  X=X-ones(N,1)*Mu;

  if (rem(N,T)~=0)
    disp('Error: Data matrix length must be multiple of sequence length T');
    return;
  end;
  N=N/T;
  
  % divide into inputs and outputs 
  
  U=X(:,1:D); % NT*D
  Y=X(:,D+1:p);
  p=p-D;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize with Linear Regression
  
  fprintf('\nInitializing with Linear Regression...\n');
  beta=U\Y; % D x p answer
  [t1,t2,t3]=svd(beta'); % p x D
  t3=t3';
  
  if K<=D,
    if K<=p,
      C=t1*t2(:,1:K);
    elseif p<=D & K>p,
      C=t1*[t2(:,1:p) 0.1*randn(p,K-p)];
    elseif p>D & K>p,
      C=t1*[t2(1:K,1:K); 0.1*randn(K-p,K)];
    end;	
    B=t3(1:K,:);
  else    
    if p<=D,
      C=t1*[t2(1:p,1:p) 0.1*randn(p,K-p)];
    elseif p>D,
      C=t1*[t2(1:D,1:D) 0.1*randn(D,K-D); 0.1*randn(p-D,K)];
    end;
    B=[t3; 0.1*randn(K-D,D)];
  end;

  Xhat=U*B';
  Yhat=Xhat*C';
  Ydiff=Y-Yhat;
  R=sum(Ydiff.*Ydiff,1)/(N*T);
  A=eye(K);
  x0=mean(Xhat)';
  t4=cov(Xhat,1);
  Q=t4+eye(K)*max(eig(t4))*0.01; % ill conditioned.
  P0=Q;
  
  clear t1 t2 t3 t4 Xhat Yhat Ydiff beta X;

  fprintf('\nInitialized.\n');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  I=eye(K);
  Y=Y';
  U=U';
  lik=0;
  LL=[];
  const=(2*pi)^(-p/2);

  Yalt=reorder(Y,T);
  Ualt=reorder(U,T);

  YY=sum(Yalt.*Yalt,2)/(T*N);

  Xpre=zeros(K,T*N);   % P(x_t | y_1 ... y_{t-1})
  Xcur=zeros(K,T*N);   % P(x_t | y_1 ... y_t)
  Xfin=zeros(K,T*N);   % P(x_t | y_1 ... y_T)    given all outputs

  Ppre=zeros(T*K,K);
  Pcur=zeros(T*K,K);
  Pfin=zeros(T*K,K); % Var(x_t|{y})
  Pt=zeros(T*K,K); % E(x_t x_t | {y})
  Pcov=zeros(T*K,K); % Var(x_t,x_t-1|{y}) used to be x_t,x_{t+1}
  Pc=zeros(T*K,K); % E(x_t x_t-1 | {y}) used to be x_t,x_{t+1}
  Kcur=zeros(T*K,p);
  invP=zeros(T*p,p);
  J=zeros(T*K,K);

  for cycle=1:cyc

    % E STEP
    
    % FORWARD PASS
    
    Xpre(:,1:N)=x0*ones(1,N);
    Ppre(1:K,:)=P0;
    invR=diag(1./(R+(R==0)*tiny));
    
    for t=1:T
      tk=(t-1)*K+1:t*K;    tN=(t-1)*N+1:t*N;  tp=(t-1)*p+1:t*p;
      if (K<p)
	T1=invR*C;
	invP(tp,:)=invR-T1*Ppre(tk,:)*inv(I+C'*T1*Ppre(tk,:))*T1';
      else
	invP(tp,:)=inv(diag(R)+C*Ppre(tk,:)*C');
      end;
      Kcur(tk,:)=Ppre(tk,:)*C'*invP(tp,:);
      Xcur(:,tN)=Xpre(:,tN)+Kcur(tk,:)*(Yalt(:,tN)-C*Xpre(:,tN)); 
      Pcur(tk,:)=Ppre(tk,:)-Kcur(tk,:)*C*Ppre(tk,:);
      if (t<T)
	Xpre(:,tN+N)=A*Xcur(:,tN)+B*Ualt(:,tN);
	Ppre(tk+K,:)=A*Pcur(tk,:)*A'+Q;
      end;
    end;  
    
    % BACKWARD PASS
    
    t=T; tk=(t-1)*K+1:t*K;  tN=(t-1)*N+1:t*N;

    Xfin(:,tN)=Xcur(:,tN);
    Pfin(tk,:)=Pcur(tk,:); 
    Pt(tk,:)=Pfin(tk,:) + Xfin(:,tN)*Xfin(:,tN)'/N; 
    
    for t=(T-1):-1:1
      tk=(t-1)*K+1:t*K;  tN=(t-1)*N+1:t*N;
      J(tk,:)=Pcur(tk,:)*A'*inv(Ppre(tk+K,:));
      Xfin(:,tN)=Xcur(:,tN)+J(tk,:)*(Xfin(:,tN+N)-Xpre(:,tN+N));
      Pfin(tk,:)=Pcur(tk,:)+J(tk,:)*(Pfin(tk+K,:)-Ppre(tk+K,:))*J(tk,:)';
      Pt(tk,:)=Pfin(tk,:) + Xfin(:,tN)*Xfin(:,tN)'/N; % E(x_t x_t)
    end;
    
    t=T;  
    tk=(t-1)*K+1:t*K; 
    tN=(t-1)*N+1:t*N;
    Pcov(tk,:)=(I-Kcur(tk,:)*C)*A*Pcur(tk-K,:);
    Pc(tk,:)=Pcov(tk,:)+Xfin(:,tN)*Xfin(:,tN-N)'/N;
    
    for t=(T-1):-1:2
      tk=(t-1)*K+1:t*K; 
      tN=(t-1)*N+1:t*N;
      Pcov(tk,:)=Pcur(tk,:)*J(tk-K,:)'+J(tk,:)*(Pcov(tk+K,:)-A*Pcur(tk,:))*J(tk-K,:)';
      Pc(tk,:)=Pcov(tk,:)+Xfin(:,tN)*Xfin(:,tN-N)'/N;
    end;    
    
    % Calculate Likelihood

    oldlik=lik;
    lik=0;
    for t=1:T % Using innovations form
      tN=(t-1)*N+1:t*N; 
      tp=(t-1)*p+1:t*p;
      MM=invP(tp,:);
      dM=sqrt(det(MM));
      if (isreal(dM) & dM>0)
	Ydiff=Yalt(:,tN)-C*Xpre(:,tN);
	lik=lik+N*log(dM)-0.5*sum(sum(Ydiff.*(MM*Ydiff)));
      else
	problem=1;
      end;
    end;
    if problem 
      fprintf(' problem '); problem=0;
    end;

    lik=lik+N*T*log(const);
    LL=[LL lik];

    fprintf('cycle %g lik %g',cycle,lik);

    if (cycle<=2)
      likbase=lik;
    elseif (lik<oldlik) 
      fprintf(' violation');
    elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~finite(lik)) 
      fprintf('\n');
      break;
    end;
    fprintf('\n');
    
    % M STEP
    
    % Re-estimate A,B,C,Q,R,x0,P0;
    
    x0old=x0;
    x0=sum(Xfin(:,1:N),2)/N; 
    T1=Xfin(:,1:N)'-ones(N,1)*mean(Xfin(:,1:N)');
    P0=Pfin(1:K,:)+T1'*T1/N; 
    
    Ptsum=zeros(K);
    for i=1:T
      tk=(i-1)*K+1:i*K;
      Ptsum=Ptsum+Pt(tk,:);
    end;

    YX=Yalt*Xfin';
    C=YX*pinv(Ptsum)/N;
    
    R=YY-diag(C*YX')/(T*N);
    
    PC1=zeros(K); 
    PT1=zeros(K);
    PT2=zeros(K);
    UU=zeros(D);
    UX1=zeros(D,K);
    UX2=zeros(D,K);
    for t=2:T 
      tk=(t-1)*K+1:t*K;      tN=(t-1)*N+1:t*N; 
      UX1=UX1+Ualt(:,tN-N)*Xfin(:,tN-N)'/N;
      UX2=UX2+Ualt(:,tN-N)*Xfin(:,tN)'/N;
      PC1=PC1+Pc(tk,:); %P_{t,t-1}
      PT1=PT1+Pt(tk,:);
      PT2=PT2+Pt(tk-K,:);
      UU=UU+Ualt(:,tN-N)*Ualt(:,tN-N)'/N;
    end;

    T1=[PT2 UX1'; UX1 UU];
    T2=[PC1 UX2'];
    T3=T2*inv(T1);
    A=T3(:,1:K);
    B=T3(:,K+1:K+D);
    
    Q1=PT1-A*PC1'-B*UX2;
    Q=diag(diag(Q1))/(T-1);

    if (det(Q)<0 | sum(sum(abs(Q-Q')))>0.001) 
      fprintf('Q problem\n');
    end;
    Q=(Q+Q')/2;

  end;
end;

