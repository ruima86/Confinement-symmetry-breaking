%Initial state in the absence of activity
clear;
close all;
%%%%%%%%% geometric parameters %%%%%%%%%%%%%%%%%%%%%%%%%
ri = 1e-3;                             % minimum radius to avoid singularity

% geometry of the cell
R0 = 3;                                % R_0 in Table S3. Spherical radius
vol0 = 2/3*R0^3;                       % V_0. Conserved volume of the cell/(2*pi)  

%%%%%%%%%% viscosity and friction %%%%%%%%%%%%%%%%%%%%%%
mu = 1;                                % related with etas and etab in Table S3 via  mu = etas
mub = 0;                               % related with etas and etab in Table S3 via  mub = (etab-etas)/2
Ga2 = 0.01;                            % Gamma in Table S3. Tangential friction coefficient
%Rh = sqrt(etab/Ga2);

%%%%%%%%%% motor dynamics and active stress %%%%%%%%%%%                 
ksig = 5;                              % k in Table S3. Turnover rate
sigeq0 = 1;                            % c_0 in Table S3. Steady state concentration
sigeqM = 5;                            % related with Delta_c in Eq. (14) via sigeqM = (Delta_c+1)*c_0
Sp = 50;                               % related with w in Table S3 via Sp = 1/w. Profile width
xi = 150;                              % (xi*Delta mu)_0 in Table S3. Active contractility
Sig0 = 10;                             % c_s in Table S3. Saturation concentration 
Di = 1;                                % D in Table S3. Diffusion constant

%%%%%%%%% elastic properties of the surface %%%%%%%%%%%
ka = 1;                                % kappa in Table S3. Bending rigidity of the surfacce  
la = 0;                                % gamma in Table S3. Passive tension

%%%%%%%%% bvp4c parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTol = 1e-5;                           % relative tolerence
Nmax = 50000;                          % maximum iteration steps
dt = 0.01;                             % simulation time step
T = 5;                                 % total simulation time
%%%%%%%%% initial guess of the lagrangian multipliers %%%%
h0 = R0*pi;                            % total arclength
v0 = 0;                                % translational velocity
fa = xi*sigeq0^2/(sigeq0^2+ Sig0^2);     % active stress
pre0 = 2*(la+fa)/R0;                   % pressure difference (volume conservation)
H = 2*R0;

%%%%%%%%  Confinement   %%%%%%%%%%%%%%%%%%%%%%%
Vext = 0.1;                            % P_0 in Table S3
Vr = 10;                               
sp = 50;                               % w_p in Table S3
b0 = 3.5;                              % minor axis of confinement
a0 = 1.3.*b0;                          % major axis of confinement

% Load the steady state solution as the beginning the dynamics
filename = ('./data/Shape_ConstraintSteadystate_Isotropic_varShell__R0_3_ka_1_sigeq0_1_sigeqM_1_mu_1_mub_0_Ga_0.01_ksig_5_xi_0_Sig0_10_Di_1.mat');
load(filename);
[~,ind] = min(abs(cell2mat(sol_all(:,1))-a0)+abs(cell2mat(sol_all(:,2))-b0));
a = sol_all{ind,1};
b = sol_all{ind,2};
sol = sol_all{ind,4};
sol.parameters(3) = [];

% Save the data in the following file
DIR = './data/';
filename = [DIR,'Shape_ConstraintIso_dynamics',...
    '_ka_',num2str(ka),...
    '_R0_',num2str(R0),...
    '_mu_',num2str(mu),...
    '_mub_',num2str(mub),...
    '_Ga_',num2str(Ga2),...
    '_xi_',num2str(xi),...
    '_Di_',num2str(Di),...
    '_Sig0_',num2str(Sig0),...
    '_sigeq0_',num2str(sigeq0),...
    '_sigeqM_',num2str(sigeqM),...
    '_ksig_',num2str(ksig),...
    '_a0_',num2str(a),...
    '_b0_',num2str(b),'.mat'];
sol_all = cell(2,3);
j = 0;
Tr = 0:0.01:T;

%%%%%% Dynamic simulation with active contractility %%%%%%%
for t = 0:dt:T
    t
    para = sol.parameters;
    uc = sol.x;
    r0c = sol.y(3,:);
    z0c = sol.y(4,:);
    p0c = sol.y(1,:);
    h0c = sol.y(12,:);
    sig0c = sol.y(9,:);
    sol0 = sol;
    
    yeq = @(u,y,para) shapefv(u,y,para,la,ksig,mu,mub,Ga2,Sig0,xi,Di,ka,sigeq0,sigeqM,Sp,sp,Vext,Vr,a,b,uc,r0c,z0c,p0c,h0c,sig0c,dt);
    ybc = @(ya,yb,para) twobcfv(ya,yb,para,ka,ri,vol0);
    jach = @(u,y,para) jac(u,y,para,la,ksig,mu,mub,Ga2,Sig0,xi,Di,ka,sigeq0,sigeqM,Sp,sp,Vext,a,b,uc,r0c,z0c,p0c,h0c,sig0c,dt);
%    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax,'FJacobian',jach,'Vectorized','on');
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
    sol = bvp5c(yeq,ybc,sol,opts);
    uc = sol.x;
    yc = deval(sol,uc);
    rc = yc(3,:);
    zc = yc(4,:);
    pc = yc(1,:);
    vu = yc(7,:);
    sigma = yc(9,:);
    y0 = deval(sol0,uc);
    r0c = y0(3,:);
    z0c = y0(4,:);
    dr = (rc - r0c)/dt;
    dz = (zc - z0c)/dt;
    psi = pc;
    vn = cos(psi).*dz + sin(psi).*dr;
       
    if min(abs(Tr-t))<0.5*dt && sol.stats.maxerr < RTol
        j = j+1;
        sol_all{j,1} = t;
        sol_all{j,2} = sol;
        sol_all{j,3} = vn;
        save(filename,'sol_all');
    end

    if sol.stats.maxerr>RTol
        break;
    end
  
end





