% function [lik, likv]=ldsi_cl(X,D,K,T,A,B,C,Q,R,x0,P0,Mu);
%
% Calculate Likelihood for Linear Dynamical System with Inputs
%
% X - N x p data matrix
% K - size of state space 
% T - length of each sequence (N must evenly divide by T)
% A - state transition matrix
% B - input matrix
% C - observation matrix 
% Q - state noise covariance 
% R - observation noise covariance
% x0 - initial state covariance
% P0 - initial state mean
% Mu - input/output mean
% 
% lik - log likelihood of X 
% likv - vector of log likelihoods
%

function [lik, likv]=ldsi_cl(X,D,K,T,A,B,C,Q,R,x0,P0,Mu);

p=length(X(1,:));
N=length(X(:,1));
tiny=exp(-700);

X=X-ones(N,1)*Mu;

if (rem(N,T)~=0)
  disp('Error: Data matrix length must be multiple of sequence length T');
  return;
end;
N=N/T;
I=eye(K);

% divide into inputs and outputs 
  
U=X(:,1:D); % NT*D
X=X(:,D+1:p);
p=p-D;
  
Y=X';
U=U';
const=(2*pi)^(-p/2);

Xpre=zeros(K,T*N);   % P(x_t | y_1 ... y_{t-1})
Xcur=zeros(K,T*N);   % P(x_t | y_1 ... y_t)
Ppre=zeros(K,T*K);
Pcur=zeros(K,T*K);
Kcur=zeros(p,T*K);

Yalt=reorder(Y,T);
Ualt=reorder(U,T);

% forward pass

Xpre(:,1:N)=x0*ones(1,N);
Ppre(:,1:K)=P0;
invR=diag(1./R);
temp1=invR*C;

lik=0;
likv=zeros(N,T);
settle=0;
thres=1e-9;

for t=1:T
  tk=(t-1)*K+1:t*K;    tN=(t-1)*N+1:t*N;
  if (~settle)
    Ppretk=Ppre(:,tk)';
    MM=invR-temp1*Ppretk*inv(I+C'*temp1*Ppretk)*temp1';
    dM=sqrt(det(MM));
    Kcurtk=Ppretk*C'*MM; 
    Pcurtk=Ppretk-Kcurtk*C*Ppretk;
    
    XpretN=Xpre(:,tN);
    Ydiff=Yalt(:,tN)-C*XpretN;
    XcurtN=XpretN+Kcurtk*Ydiff;
    
    Xcur(:,tN)=XcurtN;
    Pcur(:,tk)=Pcurtk';
    Kcur(:,tk)=Kcurtk';
    if (t<T)
      Xpre(:,tN+N)=A*XcurtN+B*Ualt(:,tN);
      Ppre(:,tk+K)=(A*Pcurtk*A'+Q)';
    end;
    if (t>5)
      if (norm(Kcurtk-Kcur(:,tk-K)')/norm(Kcurtk)<thres) 
	settle=1;
      end;
    end;
  else
    XpretN=Xpre(:,tN);
    Ydiff=Yalt(:,tN)-C*XpretN;
    XcurtN=XpretN+Kcurtk*Ydiff;
    if (t<T)
      Xpre(:,tN+N)=A*XcurtN+B*Ualt(:,tN);
    end;
  end;
  
  if (isreal(dM)& dM>0)
    likv(:,t)=log(const)+log(dM)-0.5*sum(Ydiff.*(MM*Ydiff),1)';
    lik=lik+sum(likv(:,t));
  else
    disp('problem');
  end;
end;
likv=likv';
likv=likv(:);

