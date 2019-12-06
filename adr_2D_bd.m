function [ T, S, t, xyi, A, tfin ] = adr_2D_bd( bounds, n, T0_fcn, S0_fcn, tspan, param, bc, lf )
%ADR_2D_BD Solution to the 1D advection-diffusion-reaction model using finite backward differences.
%   Basic Equations
%       T_t = k(T_{xx}) - v.T_x + A[S*exp{-B/(T-Ta)} - C(T-Ta)]
%       S_t = -Cs*S*exp{-B/(T-Ta)}, T>Ta
%     with conditions:
%       T(x,0) = T0_fcn(x), T(0,t) = T(L,t) = Ta
%
%   In order to simplify the boundary conditions the following shifted
%   problem is actually solved:
%       T_t = k(T_{xx}) - v.T_x + A[S.exp{-B/T} - CT]
%       S_t = -Cs*S*exp{-B/T}, T>0
%     with conditions:
%       T(x,0) = T0_fcn(x), T(0,t) = T(L,t) = 0
%
%
%--------------------------------------------------------------------------
%   Usage
%     [ T, S, t, x, A, tfin ] = adr_1D_bd( bounds, n, T0_fcn, tspan, param, lf )
%
%   Inputs
%     bounds    : end points of the boundary (Typically this should be [0,L;0,L])
%     n         : number of intervals to discretize the domain [nx,ny]
%			T0_fcn    : function handle defining the initial condition for T;
%			            should have form T0(x,y)
%			S0_fcn    : function handle defining the initial condition for S;
%			            should have form S0(x,y)
%			tspan     : Time span to solve the system
%     param     : System parameters
%                 k  -> thermal diffusivity
%                 v  -> wind speed
%                 A  -> temp rise/sec at max burning
%                 B  -> Arrhenius reaction rate
%                 C  -> coefficient for heat transfer to environment
%                 Cs -> fuel disappearance rate
%                 Ta -> ambient temp (K)
%     bc        : boundary condition
%                 'dirichlet' or 'neumann'
%     lf        : Logfile for output messages - default is to log to
%                 MMDDYYY.adrbd at log level 2
%
%   Outputs
%
%   Requires
%     MsgCls - class for logging error messages
%
%   Alan Lattimer, JENSEN HUGHES, 2016
%--------------------------------------------------------------------------

if nargin < 7
  error('Not enough inputs.')
end
if nargin < 8
  loglevel = 2;
  logName = [datestr(now,'mmddyyyy') '.adr2d'];
  lf = Msgcl(loglevel,logName);
end


lf.pmsg(lf.ALL,'********************************************************');
lf.pmsg(lf.ALL,'* Advection-Diffusion-Reaction');
lf.pmsg(lf.ALL,'* Finite Backward Difference Solver');

% Generate mesh and A matrix
[ A, xy, inodes, ~] = build_adr_mats2D( bounds, n, param, bc, lf );

% N is number of internal nodes
N = length(inodes);
% x and y are the internal node values
x = xy(inodes,1);
y = xy(inodes,2);

TS0          = zeros(2*N,1);
TS0(1:N)     = T0_fcn(x,y);         % Set the initial temperature conditions
TS0(N+1:end) = S0_fcn(x,y);         % Set initial mass fractions for fuel

% figure('Name','Test Mass Fraction');
% plot(x,TS0(N+1:end));


diffeq = @(t,TS) A*TS + [param.A.*f(TS,N,param.B);-param.Cs.*f(TS,N,param.B)];
% diffeq = @(t,T) A*T;
lf.pmsg(lf.ERR,'Solving the full ODE...');
tic;
[t,TS] = ode45(diffeq,tspan,TS0);
tfin = toc;
lf.pmsg(lf.ERR,'  => completed in %f seconds.',tfin);

% Nfull = length(xy);
%
% T = zeros(size(t,1),Nfull);
% S = zeros(size(t,1),Nfull);
% T(:,inodes) = TS(:,1:N);
% T = T + param.Ta;
% S(:,inodes) = TS(:,N+1:end);
%

xyi = xy(inodes,:);
T = TS(:,1:N) + param.Ta;
S = TS(:,N+1:end);

end



