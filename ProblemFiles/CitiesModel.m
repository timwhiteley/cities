function [F,Jv] = CitiesModel(u,p,ws,wr1,wr2,x,Lx,idx,v);
 
% Rename parameters
  lambda  = p(1);
  k       = p(2);
  % betas   = p(3);
  % betaR1  = p(4);
  % betaR2  = p(5);
  D       = p(6);
  f       = p(7);
  g       = p(8);
  c       = p(9);
  b       = p(10);
  nx      = size(u,1)/2;

  % Split variable
  iR = idx(:,1); iS = idx(:,2);
  r = u(iR);  s = u(iS);

  % Auxiliary functions
  compConv = @(v1,v2) (2*Lx/nx)*ifftshift(real(ifft(fft(v1) .* fft(v2))));
  AFun     = @(I,s) (1-s).*I;
  sigma    = @(I) 1-exp(-(I/lambda).^k);
  sigmaprime = @(I) (k./lambda)*(I./lambda).^(k-1) .* exp(-((I./lambda).^k));

  % Auxiliary vector
  Is  = compConv(wr1,s);
  Ir  = compConv(ws,r);
  Ir2 = compConv(wr2,r);
  A   = AFun(Is,s);
  IA  = compConv(wr2,A);

  % Right-hand side
  F = zeros(size(u));
  F(iR) = b*r .* (1-r/c) + D *( A .* Ir2 - r .* IA );
  F(iS) = (f + g*s) .* ( sigma( Ir) - s );
  

 %Compute Jacobian action
 if nargout > 1 && ~isempty(v)
 

    
 Jv=zeros(size(u));
 Jv(iS)=(f + g.*s).*(compConv(ws,v(iR)).*(sigmaprime(Ir))-v(iS));%...
         %+(f+g.*(s+v(iS))).*(sigma(Ir)-s);%This term is 0 if at ss?? 
 
 A_tilde=compConv(wr1,v(iS)).*(1-s) - v(iS).*(Is);
 Jv(iR)=D*(A.*compConv(wr2,v(iR))+A_tilde.*Ir2 - v(iR).*IA - r.*(compConv(wr2,A_tilde))) + b.*(v(iR).*(1-2*r./c));

 end
 
end
