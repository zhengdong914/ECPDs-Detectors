echo on;
clf;

% Linear Dynamical System EM Demo
% Copyright (c) by Zoubin Ghahramani, University of Toronto
% 
% A very simple demo of the Matlab 5 code for training a Linear Dynamical 
% System using the EM algorithm. Based on the following paper:
% 
%  Ghahramani, Z. and Hinton, G.E. (1996) 
%     Parameter estimation for linear dynamical systems 
%     University of Toronto Technical Report CRG-TR-96-2, 6 pages. 
%
%    ftp://ftp.cs.toronto.edu/pub/zoubin/tr-96-2.ps.gz
%    http://www.cs.toronto.edu/~zoubin/
% 
% Look at Contents.m for a description of the Matlab files
% 
% I recommend that you try to compile the CMEX code, i.e. using 'make kalman',
% so that you have the faster executable version of the Kalman smoother for
% your machine.


% generate some random data 

T=100;
X=randn(T,4);  % a sequence of 100 4-d vectors

plot(X);

% Hit any key to continue 

pause;

fprintf('\n\n\n\n');

% train LDS with 2-d state vector for 30 cycles of EM or until convergence

net=lds(X,2,T,100,0.00001);

% plot log likelihood (log L) per sample in bits

subplot(121)
plot(net.LL/(100*log(2)));             
ylabel('Log likelihood per sample (bits)'); 
xlabel('Iterations of EM'); 

subplot(122)
semilogy(diff(net.LL));
grid;
ylabel('Convergence'); 
xlabel('Iterations of EM'); 

% convergence can be quite slow towards the end...

% show some of the parameters

net.A
net.C



