% Cleaning 
clear all, close all, clc;

%% Load folders
addpath('../SecantContinuation/0.2/src');
addpath('../ProblemFiles');
addpath('../SolutionDatabase');


%% Spatial coordinates: x direction
nx = 1000; Lx = 100; hx = 2*Lx/nx; x = -Lx +[0:nx-1]'*hx;
iR = [1:nx]'; iS = nx+iR; idx = [iR iS];

%% Specify initial conditions and parameters
% p0(1)  = 3000;  % Lambda  
% p0(2)  = 3;      % k
% p0(3)  = 4;   % Betas   
% p0(4)  = 4;   % BetaR1  
% p0(5)  = 10;  % BetaR2  
% p0(6)  = 5;      % D       
% p0(7)  = 0.2;    % F       
% p0(8)  = 1;      % G       
% p0(9)  = 3000;  % C       
% p0(10) = 1e-4;   % B       
% 
% sigma    = @(I) 1-exp(-(I/p0(1)).^p0(2));
% r0 = 1500*ones(size(x)) + 100*cos(x); s0 = sigma(r0);
% u0 = [r0; s0];

% % Load initial conditions and parameters
 %sol = load('SOLUTION.mat');
 %sol = load('data/solution_0000001.mat');
 %sol=load('/Users/pmxtw2/Dropbox/Cities/Tim/Continuation/p0Continuation_3bumpsGMRES/solution_0000000.mat');
 sol=load('/Users/pmxtw2/Dropbox/Cities/Tim/Continuation/p0Continuation_2bumps_4GMRES/solution_0000000.mat');
 u0  = sol.u;
 p0  = sol.p;

%% Synaptic kernel
betas   = p0(3);
betaR1  = p0(4);
betaR2  = p0(5);
kern = @(x,beta) sqrt(1/(2*pi*beta^2)) * exp( -(1/(2*beta^2))*x.^2);
ws  = kern(x,betas);
wr1 = kern(x,betaR1);
wr2 = kern(x,betaR2);

% Plot initial solution
% PlotSolution(x,u0,p0,[],false);

% Define handle to right-hand side, jacobian and time output function
prob     = @(u,p) CitiesModel(u,p,ws,wr1,wr2,x,Lx,idx,[]);

jac      = @(u,p,v) CitiesJacobian(u,p,ws,wr1,wr2,x,Lx,idx,v);
plotSol  = @(u,p,parent) PlotSolution(x,u,p,parent,idx,false);
solMeas  = @(step,u,p) SolutionMeasures(step,u,p,x,idx,'periodic');
compSpec = @(u,p) ComputeSpectrum(u,p,ws,wr1,wr2,x,Lx,idx);
plotSpec = @(lambda,p,parent) PlotSpectrum(lambda,p,parent);

%% Assign problem 
stepPars.iContPar                         = 3;
stepPars.pMin                             = 0;
stepPars.pMax                             = 1e5;

stepPars.s0                               = 2;%-0.1;
stepPars.sMin                             = 0.1;%0.001;
stepPars.sMax                             = 500;%0.2;
stepPars.maxSteps                         = 20000;

stepPars.nPrint                           = 1;
stepPars.nSaveSol                         = 1;
stepPars.finDiffEps                       = 1e-7;
stepPars.optNonlinIter                    = 8;

 stepPars.NewtonGMRESOptions.nonlinTol     = 1e-3;
 stepPars.NewtonGMRESOptions.nonlinMaxIter = 20;
 stepPars.NewtonGMRESOptions.linTol        = 1e-2;%1e-3;
 stepPars.NewtonGMRESOptions.linrestart    = 200;
 stepPars.NewtonGMRESOptions.linmaxit      = 30;
 stepPars.NewtonGMRESOptions.damping       = 1.0;
 stepPars.NewtonGMRESOptions.display       = 0;

% stepPars.fsolveOptions        = optimset('Display','iter',...
%                                          'TolFun',1e-3,...
%                                          'MaxIter',20,...
%                                          'Jacobian','off');

stepPars.dataFolder                       = 'Data';
stepPars.PlotSolution                     = plotSol;
stepPars.BranchVariables                  = solMeas;
stepPars.PlotBranchVariableId             = 6;
stepPars.ComputeEigenvalues               = compSpec;
stepPars.PlotSpectrum                     = plotSpec;

%% Run
 branch = SecantContinuationNewtonGMRES(prob,jac,u0,p0,stepPars);
 %branch = SecantContinuation(prob,u0,p0,stepPars);
 
%Rename data file:
%
%!mv Data/ p0Continuation_2bumps_4GMRES
