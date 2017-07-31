function F = SolutionMeasures(step,u,p,x,idx,bConds)

  %% Rename parameters
  Lx = max(x);
  hx = abs(x(2)-x(1));
  nx = size(x,1);

  %% Compute quadrature weights and l2norm
  w = ones(nx,1); 
  if ~strcmp(bConds,'periodic')
    w(1) = 0.5; w(end) = 0.5; 
  end
  
  %% Compute L2-norm
  v = zeros(nx,1);
  for k = 1:size(idx,2)
    v = v + u(idx(:,k)).^2;
  end
  l2Norm = sqrt(sum( hx * w .* v)/(2*Lx));

  %% Compute L2-norm of first component
  l2NormU1 = sqrt(sum( hx * w .* u(idx(:,1)).^2)/(2*Lx));

  %% Width of U1
  fr = @(u) 1./( 1+ exp( -20*( u - 0.5 ) ) );
  WidthU1 = sum( hx * w.* fr( u(idx(:,1)) ) );

  %% Max of U1
  maxU1 = max(u(idx(:,1)));

  %% Min of U1
  minU1 = min(u(idx(:,1)));
  
  %% Max of U2
  maxU2 = max(u(idx(:,2)));
  
  %% Min of U2
  minU2 = min(u(idx(:,2)));
  
  %% L2-norm of U2
  l2NormU2 = sqrt(sum( hx * w .* u(idx(:,2)).^2)/(2*Lx));

  %% Allocate
  F = zeros(1,8);

  %% Assign branch variables
  F = [l2Norm maxU1 minU1 l2NormU1 WidthU1 maxU2 minU2 l2NormU2];
  
end
