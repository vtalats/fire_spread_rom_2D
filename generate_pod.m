  function [POD,POD_eig,POD_energy] = generate_pod(snapshots,M,r)
%-------------------------------------------------------------------------------
%  generate_pod - Compute POD basis vectors from simulation snapshots.  
%
%  Given snapshot data, a weighting matrix M (usually the mass matrix in FEM),
%  an (optional) maximum number of basis vectors, and an (optional) set of 
%  times where snapshots are provided (the default assumes they are uniformly
%  spaced), this function returns the solution to the discretized POD problem:
%
%        POD(:,i) = argmin_{POD(:,i)'*M*POD(:,j)=0, j<i}  E( POD(:,i) )
%
%  where E( POD(:,i) ) == \sum_{k=1}^n \| snapshots(k)-POD(:,i)*a(k) \|_M^2 
%        a(k) = < snapshots(k), POD(:,i) >_M,
%        < POD(:,i), POD(:,i) >_M = 1, and < a, b >_M == a'*M*b.
%                      
%  Copyright (c) 2012, Jeff Borggaard, Virginia Tech
%  Version: 2.0
%
%  Usage:   [POD,pod_energy,total_energy] = generate_pod(snapshots,M,r,time)
%
%  Variables:  POD
%                        the POD basis  (output)
%              snapshots
%                        solution snapshots   size (N,q)
%              M
%                        a positive definite matrix    size (N*m,N*m)
%                        (used to define a weighted inner product, e.g. 
%                        the finite element mass matrix)
%
%                        if m>1, and M is N*N, then the same M is used 
%                        to weigh each component.
%              r
%                        number of POD basis vectors computed
%                        (optional: default is q, the # of snapshots)
%
%  Example usage:
%              snapshots = [ u; v; w ];  % u is N x q (q timesteps)
%
%  Note that either
%              MM                            = spalloc(3*N,3*N,3*nnz(M));
%              MM(    1:  N,    1:  N)       = M;
%              MM(  N+1:2*N,  N+1:2*N)       = M;
%              MM(2*N+1:3*N,2*N+1:3*N)       = M;
%              [POD,pod_energy,total_energy] = generate_pod(snapshots,MM,q,time)
%  or 
%              [POD,pod_energy,total_energy] = generate_pod(snapshots,M,q,time)
%
%  produce the same output.
%-------------------------------------------------------------------------------
  small_problem = 5000;

  if ( nargin<2 & nargin>3 )
    error('call as generate_pod(snapshots,MM) or add r option');
  end

  [nm, ~] = size(M);
  [n , q] = size(snapshots);

  if ( mod(n,nm)~=0 )
    error('snapshots and M have incompatible dimensions')
  end
  
  if ( nargin<3 )
    r = q;
  end
 
  
  if ( n <= small_problem )                    % Solve "Small Problems" Directly
    if ( nargin==3 )
      %-------------------------------------------------------------------------
      %  Weight the velocity data
      %-------------------------------------------------------------------------
      % get the dimension if same M is used for each direction...
      dim = n/nm;
    
      [R,~,S] = chol(M);
  %    [R] = chol(M);
  %    S = speye(size(M));
      clear M
    
      rows = 1:nm:n+1;
      ss = snapshots;
      for i=1:dim
        r_range = rows(i):rows(i+1)-1;
        snapshots(r_range,:) = R*S'*snapshots(r_range,:);
      end
    end

    [U,Sigma,V] = svd(snapshots,0);
%     [U,Sigma,V] = svds(snapshots,p);
    clear snapshots
  
    [n1,n2]  = size(U);
  
    if ( nargin==2 )
      POD = U(:,1:r);
    else
      for i=1:dim
        POD(rows(i):rows(i+1)-1,:) = S* (R\U(rows(i):rows(i+1)-1, 1:r));
      end
    end
  
    POD_eig   = diag(Sigma(1:min(n1,n2),1:min(n1,n2))).^2;
    tmp       = cumsum(POD_eig);
    POD_energy = tmp/tmp(end);

  else                                          % We have a "Large Size Problem"
    dim  = n/nm;
    rows = 1:nm:n+1;
  
    % Compute the temporal autocorrelation matrix
    Rt = zeros(q,q);
    for i=1:dim
      Rt = Rt + snapshots(rows(i):rows(i+1)-1,:)'*M*...
                snapshots(rows(i):rows(i+1)-1,:);
    end
  
    Rt = 0.5*(Rt + Rt');
    [Psi,D2] = eigs(Rt,r);
    [~,index] = sort(diag(D2),'descend');
  
    %sqrt(diag(D2(index,index)))
  
    POD = snapshots*Psi(:,index);
  
    for i=1:r
      normP2 = 0;
      for j=1:dim
        normP2 = normP2 + POD(rows(j):rows(j+1)-1,i)'*M*...
                          POD(rows(j):rows(j+1)-1,i);
      end
      POD(:,i) = POD(:,i) /sqrt( normP2 );
    end
    POD_eig   = diag( D2(index,index) );
    tmp       = cumsum(POD_eig);
    POD_energy = tmp/trace(Rt);
  end
  
  % Depending on the way eigs returns we could get an e-vector in
  % the opposite direction than a previous run with the same data.  To
  % prevent this, we force the first term in the vector to always be
  % positive.
  for k=1:size(POD,2)
    if POD(1,k) < 0
      POD(:,k) = -POD(:,k);
    end
  end
end % function generate_pod
