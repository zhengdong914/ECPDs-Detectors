% function [A,C,Q,R,x0,P0,LL,LM]=lds_f(X,K,T,cyc,tol);
% 
% Adaptive Linear Dynamical System Using Fast Riccati Approximation
%
% X - N x p data matrix
% T - length of each sequence (N must evenly divide by T, default T=N)
% K - size of state space (default 2)
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% A - state transition matrix
% C - output matrix 
% Q - state noise covariance 
% R - ouput noise covariance
% x0 - initial state covariance
% P0 - initial state mean
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM

function [A,C,Q,R,x0,P0,LL,LM]=lds_f(X,K,T,cyc,tol);

p=length(X(1,:));
N=length(X(:,1));
tiny=exp(-700);
thres=1e-8;
problem=0;

if nargin<5   tol=0.0001; end;
if nargin<4   cyc=100; end;
if nargin<3   T=N; end;
if nargin<2   K=2; end;

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
  Q=cov(temp1,1);
  P0=Q;
  t1=temp1(1:N*T-1,:);
  t2=temp1(2:N*T,:);
  
  A=inv(t1'*t1+Q)*t1'*t2
  clear temp1 temp2 Phi t1 t2;
  fprintf('FA log likelihood %f after %i steps\nInitialized.\n',LM(i),i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  % AR1 model
  Q=eye(K);
  P0=Q;
  cX=cov(X,1);
  R=0.5*diag(cX);
  [u,s,v]=svd(cX-R);
  C=u*[sqrtm(s) zeros(p,K-p)];
  C=C+randn(p,K)*max(max(abs(C)))/20;
  cX
  C*C'+R
  beta=C'*inv(C*C'+R);
  t1=X*beta'+0.001*randn(N*T,K);
  x0=mean(t1)';

  t1b=t1(1:N*T-1,:);
  t1c=t1(2:N*T,:);
  A=inv(t1b'*t1b+Q)*t1b'*t1c
  clear t1 t1b t1c cX;
end;

I=eye(K);
Y=X';
lik=0;
LL=[];
const=(2*pi)^(-p/2);

Yalt=reorder(Y,T);

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
  
  tset=T+1;
  settle=0;
  invR=diag(1./R);
    
  for t=1:T
    tk=(t-1)*K+1:t*K;    tN=(t-1)*N+1:t*N;  tp=(t-1)*p+1:t*p;
    if (~settle)
      temp1=invR*C;
      invP(tp,:)=invR-temp1*Ppre(tk,:)*inv(I+C'*temp1*Ppre(tk,:))*temp1';
      Kcur(tk,:)=Ppre(tk,:)*C'*invP(tp,:);
      Pcur(tk,:)=Ppre(tk,:)-Kcur(tk,:)*C*Ppre(tk,:);
    else
      Kcur(tk,:)=Kcurset;
      Pcur(tk,:)=Pcurset;
    end;
    Xcur(:,tN)=Xpre(:,tN)+Kcur(tk,:)*(Yalt(:,tN)-C*Xpre(:,tN)); 
    if (t<T)
      Xpre(:,tN+N)=A*Xcur(:,tN);
      if (~settle)
	Ppre(tk+K,:)=A*Pcur(tk,:)*A'+Q;
      end;
    end;
    if (t>5 & ~settle)
      n1=norm(Kcur(tk,:)-Kcur(tk-K,:))/norm(Kcur(tk,:));
      if (n1<thres) 
	settle=1;
	Ppreset=Ppre(tk,:);
	invPset=invP(tp,:);
	Kcurset=Kcur(tk,:);
	Pcurset=Pcur(tk,:);
	tset=t;
      end;
    end;
  end;  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % BACKWARD PASS
  
  t=T; tk=(t-1)*K+1:t*K;  tN=(t-1)*N+1:t*N;

  Xfin(:,tN)=Xcur(:,tN);
  if (t>tset)
    Pfin(tk,:)=Pcurset;
  else
    Pfin(tk,:)=Pcur(tk,:); 
  end;
  Pt(tk,:)=Pfin(tk,:) + Xfin(:,tN)*Xfin(:,tN)'/N;
  
  settle2=0; % settling of backward recurions
  if (tset<T)
    Jset=Pcurset*A'*inv(Ppreset);
  end;
  
  for t=(T-1):-1:1
    tk=(t-1)*K+1:t*K;  tN=(t-1)*N+1:t*N;

    if (t<tset | ~settle2)
      if (t>tset)
	J(tk,:)=Jset;
	Xfin(:,tN)=Xcur(:,tN)+J(tk,:)*(Xfin(:,tN+N)-A*Xcur(:,tN));
	Pfin(tk,:)=Pcur(tk,:)+J(tk,:)*(Pfin(tk+K,:)-Ppreset)*J(tk,:)';
      else
	J(tk,:)=Pcur(tk,:)*A'*inv(Ppre(tk+K,:));
	Xfin(:,tN)=Xcur(:,tN)+J(tk,:)*(Xfin(:,tN+N)-A*Xcur(:,tN));
	Pfin(tk,:)=Pcur(tk,:)+J(tk,:)*(Pfin(tk+K,:)-Ppre(tk+K,:))*J(tk,:)';
      end;
      Pt(tk,:)=Pfin(tk,:) + Xfin(:,tN)*Xfin(:,tN)'/N;
    else
      J(tk,:)=Jset;
      Xfin(:,tN)=Xcur(:,tN)+Jset*(Xfin(:,tN+N)-A*Xcur(:,tN));
      Pfin(tk,:)=Pfinset;
      Pt(tk,:)=Pfin(tk,:) + Xfin(:,tN)*Xfin(:,tN)'/N;
    end;      
    
    if (t<T-5 & ~settle2)
      n2=norm(Pfin(tk,:)-Pfin(tk+K,:))/norm(Pfin(tk,:));
      if (n2<thres)
	settle2=1;
	Pfinset=Pcurset+Jset*(Pfin(tk,:)-Ppreset)*Jset';
	tset2=t;
      end;
    end;
    
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  t=T;  tk=(t-1)*K+1:t*K; tN=(t-1)*N+1:t*N;
  if (t>tset)
    Pcov(tk,:)=(I-Kcurset*C)*A*Pcurset;
  else
    Pcov(tk,:)=(I-Kcur(tk,:)*C)*A*Pcur(tk-K,:);
  end;
  Pc(tk,:)=Pcov(tk,:)+Xfin(:,tN)*Xfin(:,tN-N)'/N;
  
  settle3=0;
  for t=(T-1):-1:2
    tk=(t-1)*K+1:t*K; tN=(t-1)*N+1:t*N;
    if (t>tset+1 & settle3)
      Pcov(tk,:)=Pcovset;
    elseif (t>tset+1)
      Pcov(tk,:)=Pcurset*Jset'+Jset*(Pcov(tk+K,:)-A*Pcurset)*Jset'; 
    else    
      Pcov(tk,:)=Pcur(tk,:)*J(tk-K,:)'+J(tk,:)*(Pcov(tk+K,:)-A*Pcur(tk,:))*J(tk-K,:)'; 
    end;
    if (t<tset2 & ~settle3)
      n3=norm(Pcov(tk,:)-Pcov(tk+K,:))/norm(Pcov(tk,:));
      if (n3<thres)
	settle3=1;
	Pcovset=Pcurset*Jset'+Jset*(Pcov(tk+K,:)-A*Pcurset)*Jset';
	tset3=t;
      end;
    end;
    Pc(tk,:)=Pcov(tk,:)+Xfin(:,tN)*Xfin(:,tN-N)'/N;
  end;    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Calculate Likelihood

  oldlik=lik;
  lik=0;
  for t=1:T % Using innovations form
    tk = (t-1)*K+1:t*K; tN=(t-1)*N+1:t*N; tp=(t-1)*p+1:t*p;
    if (t<tset)
      MM=invP(tp,:);
      dM=sqrt(det(MM));
    end;
    if (isreal(dM) & dM>0)
      Ydiff=Yalt(:,tN)-C*Xpre(:,tN);
      lik=lik+N*log(dM)-0.5*sum(sum(Ydiff.*(MM*Ydiff)));
    else 
      problem=1;
    end;
  end;
  if problem 
    fprintf(' problem'); problem=0;
  end;

  lik=lik+N*T*log(const);
  LL=[LL lik];

  fprintf('cycle %g lik %g',cycle,lik);
  
  if (cycle<=2)
    likbase=lik;
  elseif (lik<oldlik) 
    fprintf(' violation');
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~finite(lik)) 
    break;
  end;
  fprintf('\n');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % M STEP
  
  % Re-estimate A,C,Q,R,x0,P0;
  
  x0=sum(Xfin(:,1:N),2)/N; 
  P0=Pfin(1:K,:)+cov(Xfin(:,1:N)',1); 
  
  Ptsum=zeros(K);
  for i=1:T
    tk=(i-1)*K+1:i*K;
    Ptsum=Ptsum+Pt(tk,:);
  end;

  YX=Yalt*Xfin';
  C=YX*inv(Ptsum)/N;
  
  R=YY-diag(C*YX')/(T*N); 
  
  A1=zeros(K);
  A2=zeros(K);
  A3=zeros(K);
  for t=2:T 
    tk=(t-1)*K+1:t*K;
    A1=A1+Pc(tk,:); 
    A2=A2+Pt(tk-K,:);
    A3=A3+Pt(tk,:);
  end;
  A=A1*inv(A2); 
  if(abs(eig(A))>1.001) 
    disp('Warning - unstable system');
  end;
  
  Q=(1/(T-1))*diag(diag(A3-A*(A1'))); 
  if (det(Q)<0) 
    fprintf('Q problem\n');
  end;
  
  P0=(P0+P0')/2;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;


