function smooth = PLDSLaplaceInferenceCore_sp(data, params)
%
% simplified version of PLDSLaplaceInferenceCore, for a specific case of use.
% Compute laplace-approximated posterior means and covariances for PLDS
%
% output:
% x - posterior latent variable Xk|K
% V - posterior latent covariance Wk
% VV - posterior latent variable Wk,k+1
% Ypred - Eta
% LL - likelihood
%
% remove unuseful branches of useB
% remove useless H, ds, because it keep to be all zero
% change expressions from matrix basis to scalar basis
% Sile Hu 2016-11-14
% change logdet() to log() for scalar case
% Sile Hu 2016-11-18
% changes for call sym_blk_tridiag_inv_v1_sp instead of sym_blk_tridiag_inv_v1
% change variable name to make it clear
% Sile Hu 2016-11-19
% remove chron and blkdiag
% Sile Hu 2016-11-22

runinfo = params.runinfo;
mps = params.initparams;

% initialize latent states to 0
% xInit = randn(runinfo.nStateDim, size(data.y,2)); % for subsequent EM iterations in the poisson case we'll want to initialize with previous x's.

xx = data.xxInit;%initial value of latent state
%compute inverse for latter use
Qinv = 1/mps.Q;
Q0inv = 1/mps.initV;
AQiA = mps.A'*Qinv*mps.A;% A^2/Q

LL = -inf;
LLprev = -inf;
while 1
    
    XX = xx(:,2:end) - mps.A*xx(:,1:end-1);% XX = x(k)-A*x(k-1): erro between state and its prediction via previous state
    % temporal variable
    QiXX = zeros(size(data.xxInit));
    QiXX(:,2:end) = - Qinv*XX;% -(x(k)-A*x(k-1))/Q erro normalized by covariance
    QiXX(:,1) = - Q0inv*( xx(:,1) - mps.initx );
    
    T = size(data.y,2);
    
    % sum over s and u (C and T)
    % compute Ypred from latent variable
    Ypred = bsxfun(@plus, mps.C*xx, mps.d);% C*x + d, Ypred is Eta in paper
    
    %% Latent-state grad & hessian (this should be the same for both poisson and gaussian likelihoods)
    
    lat_grad =  [ QiXX(:,1:end-1) - mps.A'*QiXX(:,2:end),   QiXX(:,end) ];%latent gradient

    lat_hess_diag = -diag([Q0inv+AQiA, ones(1,T-2) * (Qinv + AQiA), Qinv]);
    
    %this equation result in a sparse diagnal matrix(a vector) for scalar x    
    II_c = circshift(speye(T), [0 1]); II_c(end,1) = 0;% 1 right shifted identity matrix
    lat_hess_off_diag = mps.A'*Qinv * II_c; 
    lat_hess_off_diag = lat_hess_off_diag + lat_hess_off_diag';
    
    lat_hess = lat_hess_diag  + lat_hess_off_diag;
    
    %% Poisson Observation gradient and hessian
    
    Lambda = exp(Ypred);% Ypred is Eta=C*x + d in paper
    YL = data.y-Lambda;% y-lambda    
    YC = zeros(runinfo.nObsDim, size(data.y,2), runinfo.nStateDim);
    YC(:,:) = bsxfun(@times, YL, mps.C);%(y-lambda).*C
    poiss_hess_blk = -mps.C'*bsxfun(@times, Lambda, mps.C);%sum over channels (-lambda*C^2)
    
    poiss_grad = sum(YC,1);
    poiss_hess = diag(poiss_hess_blk);%diagnal is poisson hess
    
    %% Add the latent and observation hessian & gradient
    
    Hess = poiss_hess + lat_hess; % diagnal is poisson hess, left and right of diagnal is latent hess
    
    Grad = poiss_grad + lat_grad;
    
    %% Compute newton step, perform line search
    
    %break when converge
    if LL - LLprev < 1e-10 % TODO: put tolerance in setup file
        break
    end
    
    xold = xx;
    %--------there's a matrix inverse should be handled
    UPDATE  = (Hess \ Grad(:))';%compute Newton's method update
    
    LLprev = LL;
    ii = 0;
    
    LOGDETS = log(mps.initV) + (T-1)*log(mps.Q);%=log for scalar log(Q0)+(T-1)*log(Q);
    
    %this is line scan loop of Newton's method
    while 1
        dx = 2^(-ii); ii = ii + .1;
        xx = xold - dx*UPDATE;% dx is step size for update, it is getting smaller by a ratio of 2^-0.1(back tracking)
        
        %% Recompute just the likelihood @ dx step
        Ypred = bsxfun(@plus, mps.C*xx, mps.d);
        Lambda = exp(Ypred);
        
        if dx < .001 % step size too small
            break
        end
        
        % Let's see if we can't compute the likelihood.
        % fprintf('\ncomputing likelihood:\n')
        
        XX = xx(:,2:end) - mps.A*xx(:,1:end-1);
        X0 = (xx(:,1) - mps.initx);
        
        %----------why LL is like this ? ----
        PLL = sum(sum(data.y.*Ypred-Lambda));%likelihood = sum(y*Eta-Lambda)?
        GnLL = LOGDETS + sum(sum(XX.*(Qinv*XX))) + sum(sum(X0.*(Q0inv*X0)));
        
        LL = PLL - GnLL/2;
        %% Finish recomputing the likelihood
        
        if LL > LLprev % found a larger likelihood
            break
        end
    end
    
    %fprintf('%0.2f --> %0.3f', dx, LL);
    
end

%%

smooth = [];
smooth.VV = [];

Hess_diag = -[-(Q0inv + AQiA) + poiss_hess_blk(1), -(Qinv + AQiA) + poiss_hess_blk(2:end-1), -Qinv + poiss_hess_blk(end)];

%compute smoothed covariance
[smooth.V,smooth.VV]=sym_blk_tridiag_inv_v1_sp(Hess_diag, -(mps.A'*Qinv), T);

smooth.x = xx;
smooth.loglik = LL;
smooth.Ypred = Ypred;

%% With minfunc
%
% opt = struct('Method','csd', 'maxFunEvals', 500, 'Display', 'on'); % TODO: Make this an option in opt struc
% xInitmin = minFunc( @(x) poiss_obj_fun(x, data, params), xInit(:), opt);