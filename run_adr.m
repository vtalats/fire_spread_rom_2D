%RUN_ADR This is a helper script to run the coupled advection-diffusion-reaction model
%   Basic Equations
%       T_t = k(T_{xx}) - v.T_x + A[S*exp{B/(T-Ta)} - C(T-Ta)]
%       S_t = -Cs*S*exp{Bp/(T-Ta)}, T>Ta
%     with conditions:
%       T(x,0) = T0_fcn(x), T(0,t) = T(L,t) = Ta
%
%   This script was developed to support a proposal with B. Lattimer, S.
%   Gugercin, and J. Borggaard.
%   Alan Lattimer, Virginia Tech, 2015
%--------------------------------------------------------------------------

%% Full Model
run_full_model = 1;
close all;

% Full Model
if run_full_model
  
  clc;
  clear;
  close all;

  logName = [datestr(now,'mmddyyyy') '.romtst'];
  loglevel = 99;
  lf = Msgcl(loglevel,logName);
  
  fig_dir = 'latex/Figures/';
  if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
  end

  a = 0; b = 1000;
  c = 0; d = 1000;
  bounds = [a,b;c,d];
  n = [501,501];
  tf = 3000;
  nt = 800;
  tspan = linspace(0,tf,nt);
  Tc = 1200;
  Ta = 300;
  sigma = 10*sqrt(2);
  xc = (a+b)/2;
  yc = (c+d)/2;
  bc = 'dirichlet';
  movie_time = 10; % seconds
  frate  = floor(nt/movie_time);
  make_movies = 0;

  param.k   = 2.1360e-01;
  param.v   = [0,0];
  param.A   = 1.8793e02;
  param.B   = 5.5849e02;
  param.C   = 4.8372e-05;
  param.Cs  = 1.6250e-01;
  param.Ta  = 300; % degrees K

  T0_fcn = @(x,y) Tc.*(exp(-(((x-xc)./sigma).^2)).*exp(-(((y-yc)./sigma).^2)));
  S0_fcn = @(x,y) 0.8 + 0.2.*cos(12.*pi.*x./b).*sin(12.*pi./d);
%   S0_fcn = @(x,y) ones(size(x));

  tic;
  [ T, S, t, xy, A, tf ] = adr_2D_bd( bounds, n, T0_fcn, S0_fcn, tspan, param, bc, lf );
  texec = toc;
  lf.pmsg(lf.ERR,'Completed execution in %f seconds.',texec);

  x = xy(:,1);
  y = xy(:,2);
  
  N = length(xy);
  I = speye(N);
  
  if make_movies
    gen_movie_2D(T,xy,n-2,frate,'temperature');
  end


end
bw_plot = 0;

if bw_plot
  colormap(mycolormap(1));
else
  colormap(mycolormap(2));
end



%% POD Model
lf.pmsg(lf.ERR,'*-------------------------------------------------------');
lf.pmsg(lf.ERR,'POD Model:');

rT = 250;
rS = 150;
lf.pmsg(lf.ERR,'  ROM Sizes:');
lf.pmsg(lf.ERR,'    => r_T = %d',rT);
lf.pmsg(lf.ERR,'    => r_S = %d',rS);

lf.pmsg(lf.ERR,'Generating POD basis vectors.');
[PODT,SigmaT,~] = generate_pod(T'-param.Ta,I,rT);
[PODS,SigmaS,~] = generate_pod(S',I,rS);

sing_fig = figure('Name', 'Singular Values');
if bw_plot
  semilogy(1:rT,SigmaT./SigmaT(1),'k.',1:rS,SigmaS./SigmaS(1),'k*');
else
  semilogy(1:rT,SigmaT./SigmaT(1),'k.',1:rS,SigmaS./SigmaS(1),'b.');
end
title('Singular Values');
xlabel('number');
ylabel('Normalized \sigma');
legend('POD_T','POD_S');


V = blkdiag(PODT, PODS);
Ar = V'*A*V; 

TS0          = zeros(2*N,1);
TS0(1:N)     = T0_fcn(x,y);
TS0(N+1:end) = S0_fcn(x,y);

TS0r = V'*TS0;
lf.pmsg(lf.ERR,'** Skiping POD solve.');

%   ode = @(t,TSr) Ar*TSr + V'*f(V*TSr,n,param);
% ode = @(t,TSr) Ar*TSr + V'*[param.A.*f(V*TSr,N,param.B);-param.Cs.*f(V*TSr,N,param.B)];
% 
% lf.pmsg(lf.ERR,'  Solving reduced model...');
% tic;
% [tr,TSr] = ode45(ode,tspan,TS0r);
% tred = toc;
% lf.pmsg(lf.ERR,'    => completed in %f seconds.',tred);
% lf.pmsg(lf.ERR,'       => %4.2f times faster than the full model.',tf/tred);
% 
% TSfull = TSr*V';
% Tred = TSfull(:,1:N)+param.Ta;
% Sred = TSfull(:,N+1:end);

% figure('Name','Temperature (POD)');
% surf(x,t,TSfull(:,1:n)+param.Ta,'Edgecolor','none')
% if bw_plot
%   colormap(map)
%   fname = sprintf('bw_temp_pod_%d_%d',rT,rS);
% else
%   fname = sprintf('temp_pod_%d_%d',rT,rS);
% end
% xlabel('x (m)');
% ylabel('time (s)');
% zlabel('Temp (K)');
% title(sprintf('Temperature (POD : r_T=%d, r_S=%d)',rT,rS));   
% saveas(gcf,[fig_dir fname '.png'])
% % saveas(gcf,[fig_dir fname '.eps'],'epsc')
% % savefig([fig_dir fname '.fig']);
% 
% 
% errT = norm(T-Tred)/norm(T);
% errS = norm(S-Sred)/norm(S);
% 
% lf.pmsg(lf.ERR,'  Errors:');
% lf.pmsg(lf.ERR,'    => The temperature relative error is %e.',errT);
% lf.pmsg(lf.ERR,'    => The mass fraction relative error is %e.',errS);
% lf.pmsg(lf.ERR,'*-------------------------------------------------------');
% 
% %% DEIM Model
lf.pmsg(lf.ERR,'POD/DEIM Model:');
rD = 300;
lf.pmsg(lf.ERR,'  ROM Sizes:');
lf.pmsg(lf.ERR,'    => r_T    = %d',rT);
lf.pmsg(lf.ERR,'    => r_S    = %d',rS);
lf.pmsg(lf.ERR,'    => r_DEIM = %d',rD);

lf.pmsg(lf.ERR,'  => Calculate the nonlinear term.');
fsnaps = f([T';S'],N,param.B);

% [Uf,Sigmaf,~] = svd(fsnaps,0);
% U = Uf(:,1:rD);
lf.pmsg(lf.ERR,'  => Generating POD basis vectors.');
[U,Sigmaf,~] = generate_pod(fsnaps,I,rD);

figure(sing_fig);
hold on
semilogy(1:rD,Sigmaf./Sigmaf(1),'r.');
hold off
title('Singular Values');
xlabel('number');
ylabel('Normalized \sigma');
legend('POD_T','POD_S','DEIM');
fname = sprintf('singular');
saveas(gcf,[fig_dir fname '.eps'],'epsc')
savefig([fig_dir fname '.fig']);

% p = deim(U,rD);

lf.pmsg(lf.ERR,'  => Building constant matrices.');

p2 = q_deim(U);

P = I(:,p2);

PTV = blkdiag(P',P')*V;

ET =  param.A.*(PODT'*U)/(P'*U);
ES = -param.Cs.*(PODS'*U)/(P'*U);
E  = blkdiag(ET,ES);

odeDEIM = @(t,TSr) Ar*TSr + E*repmat(f(PTV*TSr,rD,param.B),2,1);

lf.pmsg(lf.ERR,'  => Solving DEIM reduced model...');
tic;
[trd,TSrd] = ode45(odeDEIM,tspan,TS0r);
tdeim = toc;
lf.pmsg(lf.ERR,'     => completed in %f seconds.',tdeim);
lf.pmsg(lf.ERR,'        => %4.2f times faster than the full model.',tf/tdeim);

lf.pmsg(lf.ERR,'     => Projecting ROM solution to approximate full-order solution.');
TSdeim = TSrd*V';
TredD = TSdeim(:,1:N)+param.Ta;
SredD = TSdeim(:,N+1:end);


if make_movies
  lf.pmsg(lf.ERR,'  => Generating a movie.');
  gen_movie_2D(TredD,xy,n-2,frate,'temperature_deim');
end

% 
% figure('Name','Temperature (POD/DEIM)');
% surf(x,t,TSdeim(:,1:n)+param.Ta,'Edgecolor','none')
% if bw_plot
%   colormap(map)
%   fname = sprintf('bw_temp_pod_deim_%d_%d_%d',rT,rS,rD);
% else
%   fname = sprintf('temp_pod_deim_%d_%d_%d',rT,rS,rD);
% end
% xlabel('x (m)');
% ylabel('time (s)');
% zlabel('Temp (K)');
% title(sprintf('Temperature Profile over Time - POD/DEIM')); 
% saveas(gcf,[fig_dir fname '.png'])
% % saveas(gcf,[fig_dir fname '.eps'],'epsc')
% % savefig([fig_dir fname '.fig']);
% 

errTd = norm(T-TredD)/norm(T);
errSd = norm(S-SredD)/norm(S);

lf.pmsg(lf.ERR,'  => The temperature relative error is %e.',errTd);
lf.pmsg(lf.ERR,'  => The mass fraction relative error is %e.',errSd);
lf.pmsg(lf.ERR,'*-------------------------------------------------------');


% 
% err_PDT = 100*(errTd-errT)/errT;
% err_PDS = 100*(errSd-errS)/errS;
% 
% lf.pmsg(lf.ERR,'  => %% Change in temperature error is %4.1f%%.',err_PDT);
% lf.pmsg(lf.ERR,'  => %% Change in mass fraction error is %4.1f%%.',err_PDS);
% lf.pmsg(lf.ERR,'*-------------------------------------------------------');
% 
% 
% %% Graph Temp Profile
% tidx = 667;
% lf.pmsg(lf.ERR,'Generating temperature profile at t = %d.',floor(t(tidx)));
% figure('Name','Temperature Profile')
% if bw_plot
%   plot(x,T(1,:),'-.','Color',[0.5,0.5,0.5]);
%   hold on
%   plot(x,T(tidx,:),'Color',[0.5,0.5,0.5]);
%   plot(x,TredD(tidx,:),'k--');
% %  plot(x(30:60:end),TredD(tidx,30:60:end),'k*');
%   hold off
%   fname = sprintf('bw_temp_profile_t_%d_%d_%d_%d',floor(t(tidx)),rT, rS, rD);
% else
%   plot(x,T(1,:),x,T(tidx,:),x,Tred(tidx,:),x,TredD(tidx,:));
%   fname = sprintf('temp_profile_t_%d_%d_%d_%d',floor(t(tidx)),rT, rS, rD);
% end
% xlabel('x (m)');
% ylabel('Temperature (K)');
% %title(sprintf('Temp profile at t=%d (POD - r_T=%d, r_S=%d; DEIM m=%d)',round(t(tidx)),rT,rS,rD));
% title(sprintf('Temperature Profile at t = %d s',floor(t(tidx))));
% legend('Initial Condition','Full Model', 'POD/DEIM', 'Location', 'best');
% saveas(gcf,[fig_dir fname '.png'])
% % saveas(gcf,[fig_dir fname '.eps'],'epsc')
% % savefig([fig_dir fname '.fig']);
% 
% %% Graph Wavefront
% lf.pmsg(lf.ERR,'Generating wavefront plot over time.');
% midn = floor(n/2);
% [~,maxj]     = max(T(:,1:midn),[],2);
% [~,maxjred]  = max(Tred(:,1:midn),[],2);
% [~,maxjredD] = max(TredD(:,1:midn),[],2);
% maxTloc      = x(maxj);
% maxTloc_red  = x(maxjred);
% maxTloc_redD = x(maxjredD);
% [~,maxj2]     = max(T(:,midn+1:end),[],2);
% [~,maxjred2]  = max(Tred(:,midn+1:end),[],2);
% [~,maxjredD2] = max(TredD(:,midn+1:end),[],2);
% maxTloc2      = x(maxj2+midn);
% maxTloc_red2  = x(maxjred2+midn);
% maxTloc_redD2 = x(maxjredD2+midn);
% figure('Name','Wavefront');
% if bw_plot
%   plot(t,maxTloc,'Color',[0.4, 0.4, 0.4])
%   hold on;
%   plot(t,maxTloc_red,'k--')
%   plot(t(1:30:nt),maxTloc_redD(1:30:nt),'ko');
%   plot(t,maxTloc2,'Color',[0.4, 0.4, 0.4])
%   plot(t,maxTloc_red2,'k--')
%   plot(t(1:30:nt),maxTloc_redD2(1:30:nt),'ko');
%   hold off;
%   fname = sprintf('bw_wavefront_%d_%d_%d',rT,rS,rD);
% else
%   plot(t,maxTloc,'k',t,maxTloc_red,'b',t(1:30:nt),maxTloc_redD(1:30:nt),'ro');
%   hold on;
%   plot(t,maxTloc2,'k',t,maxTloc_red2,'b',t(1:30:nt),maxTloc_redD2(1:30:nt),'ro');
%   hold off;
%   fname = sprintf('wavefront_%d_%d_%d',rT,rS,rD);
% end
% xlabel('time (s)');
% ylabel('x position (m)');
% title('Location of Maximum Temperature over Time');
% legend('Full Model','POD', 'POD/DEIM', 'Location', 'best');
% saveas(gcf,[fig_dir fname '.png'])
% % saveas(gcf,[fig_dir fname '.eps'],'epsc')
% % savefig([fig_dir fname '.fig']);
% 
lf.pmsg(lf.ALL,'* Completed.');
lf.pmsg(lf.ALL,'*-------------------------------------------------------');
