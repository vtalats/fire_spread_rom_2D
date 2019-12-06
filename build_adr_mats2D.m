function [ A, xy, inodes, e_conn] = build_adr_mats2D( bounds, n, param, bc, lf )
  %BUILD_ADR_MATS Matrices for the 2D advection-diffusion-reaction model using finite backward differences.
  %   Basic Equations
  %       T_t = k(T_{xx}) - v.T_x + A[S*exp{B/(T-Ta)} - C(T-Ta)]
  %       S_t = -Cs*S*exp{Bp/(T-Ta)}, T>Ta
  %     with conditions:
  %       T(x,0) = T0_fcn(x), T(0,t) = T(L,t) = Ta
  %
  %   In order to simplify the boundary conditions the following shifted
  %   problem is actually solved:
  %       T_t = k(T_{xx}) - v.T_x + A[S.exp{Bp/T} - CpT]
  %       S_t = -Cs*S*exp{Bp/T}, T>0
  %     with conditions:
  %       T(x,0) = T0_fcn(x), T(0,t) = T(L,t) = 0
  %
	%		This fuction returns the mesh and the A matrix for
	%     \dot{u} = Au + F(u) where u = [T;S]
  %
  %--------------------------------------------------------------------------
  %   Usage
  %     [ x, A ] = build_adr_mats2D( bounds, n, param, lf )
  %
  %   Inputs
  %     bounds    : domain of the 2D boundary (Typically this should be [0,L;0,L])
  %     n         : number of intervals to discretize the domain n = [nx,ny]
  %     param     : System parameters
  %                 k  -> thermal diffusivity
  %                 v  -> wind speed
  %                 A  -> temp rise/sec at max burning
  %                 B  -> Arrhenius reaction rate
  %                 C  -> coefficient for heat transfer to environment
  %                 Cs -> fuel disappearance rate
  %                 Ta -> ambient temp (K)
  %     lf        : Logfile for output messages - default is to log to
  %                 MMDDYYY.adrbd at log level 2
  %
  %   Outputs
	%     A         : A matrix from above
	%			xy        : mesh points
	%			e_conn    : connectivity matrix
	%			inodes    : internal nodes
  %
  %   Requires
  %     MsgCls - class for logging error messages
  %
  %   Alan Lattimer, Virginia Tech, 2015
  %--------------------------------------------------------------------------

  if nargin < 4
    error('Not enough inputs.')
  end
  if nargin < 5
    loglevel = 2;
    logName = [datestr(now,'mmddyyyy') '.badr'];
    lf = Msgcl(loglevel,logName);
  end

  switch bc
    case 'dirichlet'
      % Solve on internal nodes only
      N = prod(n-2);
      nx = n(1)-2;
    case 'neumann'
      N = prod(n);
      nx = n(1);
    otherwise
      % default to Dirichlet bc
      N = prod(n-2);
      nx = n(1)-2;
  end


	[xy, e_conn, bnodes] = twod_mesh_Q(1, n(1), n(2), bounds);

	inodes = find(~bnodes);
	h = xy(2,1) - xy(1,1);


  lf.pmsg(lf.ALL,'*-------------------------------------------------------')
  lf.pmsg(lf.ALL,'* Generate ADR Finite Backward Difference Matrices');
  lf.pmsg(lf.ALL,'*-------------------------------------------------------')
  lf.pmsg(lf.ALL,'* Parameters:');
  lf.pmsg(lf.ALL,'*   k  = %e',param.k);
  lf.pmsg(lf.ALL,'*   v  = (%4.2f,%4.2f)',param.v)
  lf.pmsg(lf.ALL,'*   A  = %e',param.A);
  lf.pmsg(lf.ALL,'*   B  = %e',param.B);
  lf.pmsg(lf.ALL,'*   C  = %e',param.C)
  lf.pmsg(lf.ALL,'*   Cs = %e',param.Cs)
  lf.pmsg(lf.ALL,'*   Ta = %6.2f',param.Ta);
  lf.pmsg(lf.ALL,'* Discretization:');
  lf.pmsg(lf.ALL,'*   domain: [%d,%d] x [%d,%d]',bounds');
  lf.pmsg(lf.ALL,'*   nx = %d, ny = %d',n);
  lf.pmsg(lf.ALL,'*   h  = %f',h);
  lf.pmsg(lf.ALL,'*-------------------------------------------------------')
  lf.pmsg(lf.ALL,'* Alan Lattimer, JENSEN HUGHES, 2016');
	lf.pmsg(lf.ALL,'* NOTE: Currently the convection term is turned OFF.');
  lf.pmsg(lf.ALL,'*-------------------------------------------------------')
  lf.pmsg(lf.ERR,'Building system matrices...');
  % Create helper vector to generate system matrices
  d  = ones(N,1);

  % Set nx to the number of internal nodes in the x direction

  % Set indices used to zero out portions of the matrix related to boundary
  % nodes.
	z1 = nx:nx:(N - nx);
	z2 = z1 + 1;

  AT1 = (param.k/h^2) .* spdiags([d d -4*d d d],[-nx -1 0 1 nx], N, N);
  u = param.v(1).*d;
  v = param.v(2).*d;
  AT2 = AT1 + ((-1/(2*h)).* spdiags([-v -u u v],    [-nx -1 1 nx],   N, N));
  AT2(z1,z2) = 0;
  AT2(z2,z1) = 0;
  AT3 = -param.C.*param.A.*speye(N);
  AT  = AT2 + AT3;
  AS = sparse(N,N);
  A  = blkdiag(AT,AS);


  lf.pmsg(lf.PED,'  => completed.');

end

