%Initial state in the absence of activity
clear;
close all;
%%%%%%%%% geometric parameters %%%%%%%%%%%%%%%%%%%%%%%%%
ri = 1e-3;                     % minimum radius to avoid singularity

% geometry of the cell
R0 = 3;                        % R_0 in Table S3.    Spherical radius
vol0 = 2/3*R0^3;               % V_0                Conserved volume of the cell/(2*pi)  

%%%%%%%%%% viscosity and friction %%%%%%%%%%%%%%%%%%%%%%
mu = 1;                        % related with etas and etab in Table S3 via  mu = etas
mub = 0;                       % related with etas and etab in Table S3 via  mub = (etab-etas)/2
Ga2 = 0.01;                    % Gamma in Table S3        Tangential friction coefficient
%Rh = sqrt(etab/Ga2);

%%%%%%%%%% motor dynamics and active stress %%%%%%%%%%%                 
ksig = 5;                      % k in Table S3                 Turnover rate
sigeq0 = 1;                    % c_0 in Table S3               Steady state concentration
sigeqM = 1;                    % related with Delta_c in Eq. (14) via     sigeqM = (Delta_c+1)*c_0
Sp = 50;                       % related with w in Table S3 via Sp = 1/w;   profile width
xi = 0;                        % (xi*Delta mu)_0 in Table S3   Active contractility
Sig0 = 10;                     % c_s in Table S3               Saturation concentration 
Di = 1;                        % D in Table S3                 Diffusion constant

%%%%%%%%% elastic properties of the surface %%%%%%%%%%%
ka = 1;                        % kappa in Table S3.  Bending rigidity of the surfacce  
la = 0;                        % gamma in Table S3. Passive tension

%%%%%%%%% bvp4c parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTol = 1e-5;                   % relative tolerence
Nmax = 50000;                  % maximum iteration steps

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

%%%%%%%% Initial guess of the surface shape and motor concentration %%%%%%
%%%%%%%% The shape is spherical. The motor concentration is uniform %%%%%%
yeq = @(u,y,para)shapef(u,y,para,ka,la,ksig,mu,mub,Ga2,Sig0,xi,Di,sigeq0,sigeqM,Sp,sp,Vext,a0,b0,Vr);
ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,ri,vol0);
opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
yinit = @(u) guessf(u,ka,R0,ri,sigeq0);
solinit = bvpinit(linspace(0,1,101),yinit,[h0,pre0,v0]);

% Find the first solution
sol = bvp5c(yeq,ybc,solinit,opts);

DIR = './data/';
filename = [DIR,'Shape_ConstraintSteadystate_Isotropic_varShell_',...
    '_R0_',num2str(R0),...
    '_ka_',num2str(ka),...
    '_sigeq0_',num2str(sigeq0),...
    '_sigeqM_',num2str(sigeqM),...    
    '_mu_',num2str(mu),...
    '_mub_',num2str(mub),...
    '_Ga_',num2str(Ga2),...
    '_ksig_',num2str(ksig),...
    '_xi_',num2str(xi),...
    '_Sig0_',num2str(Sig0),...
    '_Di_',num2str(Di),'.mat'];

j = 0;
sol_all = cell(1,4);
for ratio = 1:-0.0002:0.79
    a = a0*ratio;
    b = b0*ratio;
    [ratio,a,b]
    yeq = @(u,y,para)shapef(u,y,para,ka,la,ksig,mu,mub,Ga2,Sig0,xi,Di,sigeq0,sigeqM,Sp,sp,Vext,a,b,Vr);
    ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,ri,vol0);
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
    sol0 = sol;
    sol = bvp5c(yeq,ybc,sol,opts);
    if sol.stats.maxerr < RTol
        j = j+1;
        sol_all{j,1} = a;
        sol_all{j,2} = b;
        H = sol.y(4,1)-sol.y(4,end);
        sol_all{j,3} = H;
        sol_all{j,4} = sol;
        save(filename,'sol_all');
    else
        sol = sol0;
        break;
    end
end