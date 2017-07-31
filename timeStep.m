%% Cleaning 
clear all, close all, clc;

%% Load problem files and solution
addpath('../ProblemFiles');
addpath('../SolutionDatabase');

%% Spatial coordinates: x direction
nx = 1000; Lx = 100; hx = 2*Lx/nx; x = -Lx +[0:nx-1]'*hx;
iR = [1:nx]'; iS = nx+iR; idx = [iR iS];

%% Specify initial conditions and parameters
% p0(1)  = 1500;  % Lambda  
% p0(2)  = 3;      % k
% p0(3)  = 4;   % Betas   
% p0(4)  = 4;   % BetaR1  
% p0(5)  = 15;  % BetaR2  
% p0(6)  = 20;      % D       
% p0(7)  = 0.04;    % F       
% p0(8)  = 0.8;      % G       
% p0(9)  = 100;  % C       
% p0(10) = 1;   % B       
% 
% sigma    = @(I) 1-exp(-(I/p0(1)).^p0(2));
% r0 = 1600*ones(size(x)) + 100*cos(0.04*x); s0 = sigma(r0);
% r0=580*rand(size(r0));s0=0.1*rand(size(s0));
% u0 = [r0; s0];

 % Specify initial conditions and parameters
  sol=load('/Users/pmxtw2/Dropbox/Cities/Tim/Continuation/p0Continuation_2bumps_4GMRES/solution_0000100.mat');
  %sol = load('../Continuation/SOLUTION.mat');

  u0  = sol.u; u0(iR)=u0(iR)+200*rand(size(iR));u0(iS)=u0(iS)+0.1*rand(size(iS));
  p0  = sol.p;


%% Synaptic kernel
betas   = p0(3);
betaR1  = p0(4);
betaR2  = p0(5);
kern = @(x,beta) sqrt(1/(2*pi*beta^2)) * exp( -(1/(2*beta^2))*x.^2);
ws  = kern(x,betas);
wr1 = kern(x,betaR1);
wr2 = kern(x,betaR2);

%figure; plot(x,ws,x,wr1,x,wr2);

PlotSolution(x,u0,p0,[],idx,false);

%% Define handle to right-hand side, jacobian and time output function
prob    = @(t,u) CitiesModel(u,p0,ws,wr1,wr2,x,Lx,idx);
timeOut = @(t,u,flag) TimeOutput(t,u,flag);
  
%% Integrate ODE
options = odeset('OutputFcn',timeOut);
tSpan = [0:100:4000]; 
[t,U] = ode113(prob,tSpan,u0,options);
    %ode45 gives some errors?
 
%% Plot final solution
PlotSolution(x,U(end,:)',p0,[],idx,false);
% 
% % Plot history
PlotHistory(x,t,U,p0,[],idx);
% 
%% Save
  p = p0;
  u = U(end,:)';
  save('../Continuation/SOLUTION.mat','u','p');

%% Plot conserved quantity
figure;
plot(t,2*Lx*mean(U(:,iR),2),'.-');
xlabel('t'); ylabel('total population');


