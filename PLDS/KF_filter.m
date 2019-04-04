

function [state, stateCov] = KF_filter(ssm,y, z0, Q0) 


if nargin < 4; 
    Q0 = ssm.P0;    
end
if nargin < 3; 
    z0 = ssm.x0;    
end

% z0 is the latent state estimate from the last time point prior to y(1:T) 
% Q0 is the latent state variance from the last time point prior to y(1:T)

A = ssm.A; % m-by-m
C = ssm.C; % C-by-m
d = ssm.Mu'; % C-by-1

Q = ssm.Q; % m-by-m
R = diag(ssm.R); % C-by-C 


 
[NoUnit,T] = size(y);
% when z0 is univariate, m=1
m = size(z0,1); 
state = zeros(m,T);
stateCov = zeros(m,m,T);

% assuming single-trial structure 
for t=1:T
    % prediction
    z_pred = A * z0; % m-by-1
    Q_pred = A * Q0 * A' + Q; % m-by-m
    
    y_pred = C * z_pred + d; 
    % filtering
    K = Q_pred * C' * inv(C * Q_pred * C'+ R); %m-by-C  % Kalman gain 
     
    z_filt = z_pred + K * (y(:,t) - y_pred);  % m-by-1 
    Q_filt = (eye(m) - K*C) * Q_pred; % m-by-m
    
  
   
    % save result 
    state(:,t) = z_filt;
    stateCov(:,:,t) = Q_filt;
    
    % update
    z0 = z_filt; 
    Q0 = Q_filt; 
    

end

if m==1
    stateCov = squeeze(stateCov);
end