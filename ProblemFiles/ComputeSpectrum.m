% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function [V,D] = ComputeSpectrum(u,p,ws,wr1,wr2,x,Lx,idx);

  %% Linear function handle
  jhandle = @(v) CitiesJacobian(u,p,ws,wr1,wr2,x,Lx,idx,v);

  %% Call eigenvalue solver
  eigsOpts = []; 
  %eigOpts.maxit=500;
  %eigOpts.tol=1e-6;

  [V,D,flag] = eigs(jhandle,size(u,1),10,'lr',eigOpts);

  if flag ~= 0
    warning('Not all eigenvalues converged');
  end

end
