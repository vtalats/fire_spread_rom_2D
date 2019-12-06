function [nl_out] = f(TS,n,B)
	%F Calculate the combustion nonlinearity
	%   The basic fire-spread model has an Arrhenius kinetics term
	%   that couples the fuel mass fraction to the energy equation.
	%			F(T,S) := S*exp{-B/T} for all T > 0
	%
	%		Inputs
	%			TS - Temperature and Mass Fraction [T;S] where size(T)=size(S)
	%     n  - lenght(T)
	%			B  - Beta constant
	%		Output
	%			nl_out - output of the nonlinear function F(S,T)
	%
	%--------------------------------------------------------------------------
	% Author: Alan Lattimer, JENSEN HUGHES
	% Date: June 2016
	%--------------------------------------------------------------------------

  T = TS(1:n,:);
  S = TS(n+1:end,:);
  aurK = zeros(size(T));
  pos_idx = find(T>0);
  aurK(pos_idx) = exp(-B./(T(pos_idx)));
  nl_out     = S.*aurK;
end


