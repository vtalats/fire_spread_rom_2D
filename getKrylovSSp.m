function [ Vr, Wr ] = getKrylovSSp( s, A, b, c, E, lf )
%GETKRYLOVSSP Generates Rational Krylov subspace given sigma.
%   The following Rational Krylov subspaces are created
%   
%   Vr = [(s_1E-A)\b (s_2E-A)\b ... (s_rE - A)\b]
%
%   Wr' = [c/(s_1E-A); c/(s_2E-A); ... ;c/(s_rE-A)]


k = 1;

n = size(A,1);
r = length(s);

Vr = zeros(n,r);
Wr = zeros(r,n);

while k <= r
    Ak = ((s(k).*E) - A);
  
  if (abs(imag(s(k)))/abs(s(k))) < 1e-08 
    lf.pmsg(lf.PED, 'sigma_%d = %f is real.', k, s(k));
    Vr(:,k) = real(Ak\b);
    Wr(k,:) = real(c/Ak);
    k = k+1;
  else
    lf.pmsg(lf.PED, 'sigma_%d = %f%+fi is complex.', k, real(s(k)), imag(s(k)));
    vtmp = Ak\b;
    wtmp = c/Ak;
    Vr(:,k:k+1) = [real(vtmp) imag(vtmp)];
    Wr(k:k+1,:) = [real(wtmp);imag(wtmp)];
    k = k+2;
  end
  
end
  
[Vr,~]=qr(Vr,0);

[Wr,~] = qr(Wr',0); 


end

