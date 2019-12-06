function [ Vr, Wr, siter ] = irka_siso( A, b, c, d, E, r, sysinit, loglevel, Stol,lf )
%IRKA_SISO Reduces an nxn SISO dynamical system to an rxr SISO reduced system using IRKA.
%
%   Usage:
%     
%   Inputs:
%     A, b, c, d, E - state space matrices 
%             A and E - n x n
%                   b - n x 1
%                   c - 1 x n
%                   d - 1 x 1
%     r - size of the reduced order model.  Typically, this should be
%         r << n, and most certainly r < n.  If r >= n, then no reduction 
%         is performed and the original system is returned with no error.
%     loglevel - controls how much gets logged
%     tol - tolerence for sigma difference
%
%   Outputs:
%     Ar, Br, Cr, Dr - reduced system state space matrices 
%            Ar and Er --> r x r
%                   br --> r x 1
%                   cr --> 1 x r
%                   dr --> 1 x 1
%     err - The H_2 error between the system and reduced models
% 
%     
%

if nargin < 10
  logName = [datestr(now,'mmddyyyy') '.irkas'];
  lf = Msgcl(loglevel,logName);
  if nargin < 9
    Stol = 1e-5;
    if nargin < 7
      loglevel = 2;
    end
  end
  old_ll = lf.logLevel;
else
  old_ll = lf.logLevel;
  lf.logLevel = loglevel;
end

n = size(A,1);

if size(A,2) ~= n
  error('The ''A'' matrix must be square.');
elseif (size(b,2) ~= 1) || (size(c,1) ~= 1)
  error('Too many inputs or outputs.  This is only for SISO systems');
elseif (size(d,1) ~= 1) || (size(d,2) ~= 1)
  error('The ''D'' matrix is the wrong size.  It should be (num outouts) x (num inputs).');
end

calc_Hinf_err = 0;
Htol = 1e-10;

lf.pmsg(lf.ALL,'==========================================================');
lf.pmsg(lf.ALL,'* Beginning IRKA SISO model reduction');
lf.pmsg(lf.ALL,'*      n = %d ==> r = %d',n,r);
lf.pmsg(lf.ALL,'*   Log level is set to %d.',loglevel);


if n < r
  lf.pmsg(lf.ALL, 'Since r > n, the original system is being returned.');
  Vr = eye(n);
  Wr = eye(n);
else
  if sysinit
    lf.pmsg(lf.WARN,'Setting initial sigma based on eigs of A.');
    sigma = -eigs(E\A,r);
  elseif exist('initialrom.mat','file')
    lf.pmsg(lf.WARN,'Found initialization file.  Using to set sigma');
    load InitialROM;
    sigma = -eig(Er\Ar);
  else
    lf.pmsg(lf.PED, 'Creating %d sigma values using logspace.',r);
    sigma = logspace(-3,-1,r)';
  end
  sigma = sort(sigma);
  
  [Vr,Wr] = getKrylovSSp(sigma,A,b,c,E,lf);

  lf.pmsg(lf.PED,'  Creating initial reduced matrices.');
  Ar = Wr'*A*Vr;
  Er = Wr'*E*Vr;
  Br = Wr'*b;
  Cr = c*Vr;
  
  if calc_Hinf_err
    sysr = dss(Ar,Br,Cr,full(d),Er);
  end
  
  k = 1;
  error_in_s = inf;
  error_in_H = inf;
  min_err    = inf;
  min_Herr   = inf;
  above_tol  = 1;
  
 
  while (k <= 100) && (above_tol)
    lf.pmsg(lf.WARN,'Starting irka iteration %d',k);
    siter(:,k) = sigma;
    % sigma_new = sort(-eig(Er\Ar));
    [X,S]=eig(Er\Ar);

    lf.pmsg(lf.PED,'  Build new shifts and directions.');
    sigma_new = -diag(S);   
    neg_ent = find(sigma_new<0);
    sigma_new(neg_ent) = -sigma_new(neg_ent);
        
    lf.pmsg(lf.PED, '  Calculating the relative error of the sigma shift vector.');
    error_in_s = norm(sigma_new-sigma)/norm(sigma);
    sigma = sigma_new;
    
    [Vr,Wr] = getKrylovSSp(sigma,A,b,c,E,lf);
    
    lf.pmsg(lf.PED,'  Creating reduced matrices for this iteration');
    Ar = Wr'*A*Vr;
    Er = Wr'*E*Vr;
    Br = Wr'*b;
    Cr = c*Vr;

    if calc_Hinf_err
      lf.pmsg(lf.PED, 'Calculating the relative error of the transfer function.');
      sysr_new = dss(Ar,Br,Cr,full(d),Er);
      error_in_H = norm(sysr-sysr_new,inf)/norm(sysr,inf);
      sysr = sysr_new;
    end


    if error_in_s < min_err
      min_err = error_in_s;
    end
    
    if error_in_H < min_Herr
      min_Herr = error_in_H;
    end
    
    lf.pmsg(lf.WARN, '  Sigma error = %e',error_in_s);
    if calc_Hinf_err
      lf.pmsg(lf.WARN, '      H error = %e',error_in_H);
    end;
    
    if error_in_s <= Stol && k>2
      above_tol = 0;
    elseif error_in_H <= Htol && k>2
      above_tol = 0;
    end
    
    k = k+1;
  end

%   br = Wr'*b;
%   cr = c*Vr;
%   dr = d;
  if k < 100
    lf.pmsg(lf.ALL,'*  Converged.');
    lf.pmsg(lf.ALL,'*  Total iterations = %d',k-1);
    lf.pmsg(lf.ALL,'*  Final sigma error = %e',error_in_s);
    if calc_Hinf_err    
      lf.pmsg(lf.ALL,'*  Final H_inf error = %e',error_in_H);
    end
  else
    lf.pmsg(lf.ALL,'*  WARNING: Failed to converge!!');
    lf.pmsg(lf.ALL,'*  Minimum sigma error = %e',min_err);
    if calc_Hinf_err    
      lf.pmsg(lf.ALL,'*  Minimum H_inf error = %e',min_Herr);
    end
  end
  
end

lf.pmsg(lf.ALL,'*  End IRKA model reduction');
lf.pmsg(lf.ALL,'==========================================================');
lf.logLevel = old_ll;
end

