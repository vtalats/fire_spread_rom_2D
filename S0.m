function [ S ] = S0( x, init_cond )
  %S0 Sets various initial conditions for fuel mass fraction
  %   Used to test different input conditions
  % Inputs
  %   x         : domain location
  %   init_cond : initial condition number
  % Output
  %   S         : fuel mass fraction
  %
  %------------------------------------------------------------
  % Alan Lattimer, VT, April 2016
  %------------------------------------------------------------
  
  b = 1000;
  switch init_cond
    case 0
      S = ones(size(x));
    case 1
      S = 0.85 + 0.15.*cos(8.*pi.*x./b);
    case 2
      S = ones(size(x));
      S(x<(0.2*b)&x>(0.19*b)) = 0.1;
      S(x<(0.81*b)&x>(0.8*b)) = 0.1;
    case 3
      S = 1 - 0.2*rand(size(x));
    case 4
      S = (1/(0.25*b)).*x;
      S(S>1) = 1;
      S(x>(0.75*b)) = 0.8;
    case 5
      S = (0.4/b).*x + 0.6;
    case 6
      S = 1 - ((0.4/b).*x);
    otherwise 
      S = ones(size(x));
  end
  
  
end

