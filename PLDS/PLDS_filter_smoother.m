
function [smooth_state, stateCov_smooth, newseq, Qfun] = PLDS_filter_smoother(params,seq, z0,Q0, binsize) 


 
if nargin < 5; 
    binsize = 0.05; % 50 ms
end
if nargin < 4; 
    Q0 = params.model.Q0;    
end
if nargin < 3; 
    z0 = params.model.x0;
end



A = params.model.A; % m-by-m
C = params.model.C; % C-by-m
d = params.model.d; % C-by-1
Q = params.model.Q; % m-by-m

y = seq.y; 
[NoUnit,T] = size(y);

% when z0 is univariate, m=1
m = size(z0,1); 
state = zeros(m,T);
pred_state = zeros(m,T);
smooth_state = zeros(m,T);

stateCov_pred = zeros(m,m,T);
stateCov_filt = zeros(m,m,T);
stateCov_smooth = zeros(m,m,T);

% assuming single-trial structure 
for t=1:T
    % prediction
    z_pred = A * z0; % m-by-1
    Q_pred = A * Q0 * A' + Q; % m-by-m
    rate = exp(C * z_pred + d); % C-by-1
    y_pred = rate * binsize; % C-by-1 
    
    % filtering
    invQ = inv(Q_pred) + (C' * diag(y_pred) * C); %m-by-m
    Q_filt = inv(invQ);    % m-by-m
    z_filt = z_pred + Q_filt * (C' * (y(:,t) - y_pred));  % m-by-1 
    

%     z_filt = z_pred + Q_pred * (C' * (y(:,t) - y_pred));  % m-by-1
%     
%     invQ = inv(Q_pred) + (C' * diag(y_pred) * C); %m-by-m
%     Q_filt = inv(invQ);    % m-by-m
    
    
 
    % save result 
    pred_state(:,t) = z_pred; 
    state(:,t) = z_filt;
    stateCov_pred(:,:,t) = Q_pred;
    stateCov_filt(:,:,t) = Q_filt;
    % update
    z0 = z_filt; 
    Q0 = Q_filt;
 
  
end

smooth_state(:,T) = state(:,T);
stateCov_smooth(:,:,T) = stateCov_filt(:,:,T);

% smoothing 
for t=T-1:-1:1
    
    B_t = A * stateCov_filt(:,:,t) * inv(stateCov_pred(:,:,t+1)); % m-by-m
    smooth_state(:,t) = state(:,t) + B_t * (smooth_state(:,t+1) - pred_state(:,t+1)); 
    
    stateCov_smooth(:,:,t) = stateCov_filt(:,:,t) + B_t * (stateCov_smooth(:,:,t+1) - stateCov_pred(:,:,t+1)) * B_t';
    
end


 
newseq = seq;  

newseq.posterior.xsm = smooth_state; % posterior mean over latent
newseq.posterior.Vsm = zeros(T*m,m); % posterior covariance over latent
for t=1:T
    newseq.posterior.Vsm((t-1)*m+1: t*m, :) = stateCov_smooth(:,:,t);
end
s0 = squeeze(stateCov_smooth(1,1,:));
newseq.posterior.VVsm = zeros((T-1)*m,m); % posterior time-lag 1 covariance
for t=1:T-1
    newseq.posterior.VVsm((t-1)*m+1: t*m, :) = sqrtm( stateCov_smooth(:,:,t) * stateCov_smooth(:,:,t+1) );
end

% correlation estimation (Smith & Brown, 2003)
QQ1 = zeros(m,m,T);
QQ2 = zeros(m,m,T);

for t=1:T
    QQ1(:,:,t) = stateCov_smooth(:,:,t) + smooth_state(:,t) * smooth_state(:,t)';% Wk
end
for t=1:T-1
    B_t = A * stateCov_filt(:,:,t) * inv(stateCov_pred(:,:,t+1));
    tem = B_t * stateCov_smooth(:,:,t+1);
    QQ2(:,:,t) = tem + smooth_state(:,t) * smooth_state(:,t+1)';% Wk,k+1
end

% compute the Q function
Qfun = 0;
% compute E[log p(N0,t|x,theta*)]
for t=1:T
    rate = exp(C * smooth_state(:,t) + d);
    tem1 = sum(y(:,t) .* (C * smooth_state(:,t) + d + log(binsize)) - rate * binsize);   % sum wrt channel
    Qfun = Qfun + tem1; % sum wrt time
end    


%NOTE:  a'*b*a = trace(a*a'*b) 

% noise covariance Q has to be full rank
% compute E[log p(x|Rho,alpha,sigma^2)] neglecting last two terms
for t=1:T-1
    tem2 = -0.5 * trace( (QQ1(:,:,t+1) - 2*A*QQ2(:,:,t) + A*QQ1(:,:,t)*A') * inv(Q));  
    Qfun = Qfun + tem2 - 0.5*log(det(Q)) -0.5*m*log(2*pi);
end


newseq.posterior.varBound = Qfun; 