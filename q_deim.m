function [ S, M ] = q_deim( U )
  %Q_DEIM - Calculate the DEIM selection matrix.
  % This is based on the method developed in 
  %     Zlatko Drmac and Serkan Gugercin. A New Selection Operator for the 
  %     Discrete Empirical Interpolation Method: improved a priori error bound 
  %     and extensions. arXiv preprint arXiv:1505.00370, 2015.
  %
  % Input : U n x m with orthonormal columns
  % Output : S selection of m row indices with guaranteed upper bound
  %             norm(inv(U(S,:))) <= sqrt(nxm+1) * O(2^m).
  %        : M the matrix U*inv(U(S,:));
  %
  % The Q?DEIM projection of an nx1 vector f is M*f(S). 
  % Coded by Zlatko Drmac, April 2015.
  % 
  [n,m] = size(U) ;

  if nargout == 1
    [~,~,P] = qr(U','vector') ; 
    S = P(1:m);
  else
    [~,R,P] = qr(U','vector') ; 
    S = P(1:m) ; 
    M = [eye(m) ; (R(:,1:m)\R(:,m+1:n))'] ; 
    Pinverse(P) = 1 : n ; 
    M = M(Pinverse,:) ;
  end
end
