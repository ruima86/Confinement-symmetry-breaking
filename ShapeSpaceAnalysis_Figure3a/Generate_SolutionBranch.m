% Exemplary computation of a shape solution branch
% Red curve with V_cell/V_shell = 0.638 in Figure 3b
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
sigeqM = 5;                    % related with Delta_c in Eq. (14) via     sigeqM = (Delta_c+1)*c_0
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
b0 = 3.1936;                              % minor axis of confinement
a0 = 1.3.*b0;                          % major axis of confinement

%%%%%%%% The first solution of branch %%%%%%%%%%
filename = ('./data/Shape_ConstraintSteadystate_Isotropic_varShell__R0_3_ka_1_sigeq0_1_sigeqM_5_mu_1_mub_0_Ga_0.01_ksig_5_xi_0_Sig0_10_Di_1.mat');
load(filename);
[~,ind] = min(abs(cell2mat(sol_all(:,1))-a0)+abs(cell2mat(sol_all(:,2))-b0));
a = sol_all{ind,1};
b = sol_all{ind,2};
sol = sol_all{ind,4};

%%%%%%%% Save the data in the following file
DIR = pwd;
filename = [DIR,'/data/Shape_ConstraintSteadystate_Isotropic_varActivity',...
                     '_R0_',num2str(R0),...
                     '_ka_',num2str(ka),...
                     '_sigeq0_',num2str(sigeq0),...
                     '_sigeqM_',num2str(sigeqM),...
                     '_mu_',num2str(mu),...
                     '_mub_',num2str(mub),...
                     '_Ga_',num2str(Ga2),...
                     '_ksig_',num2str(ksig),...
                     '_Sig0_',num2str(Sig0),...
                     '_Di_',num2str(Di),...
                     '_a_',num2str(a),...
                     '_b_',num2str(b),'.mat'];

j = 0;
sol_all = cell(1,5);
% varying xi
Mc = 0:0.1:30.32;
[sol,sol_all,j] = iterate_xi(filename,Mc,sol,sol_all,j,la,ksig,mu,mub,Ga2,Sig0,sigeqM,Di,ka,sigeq0,Sp,RTol,Nmax,ri,vol0,sp,Vext,a,b,Vr);

%varying H
xi = sol_all{j,1};
H = sol_all{j,3};
H0 = sol_all{j-1,3};
d = 0.0002;
if H>H0
    Hc = H+d:d:8.15756;
else
    Hc = H-d:-d:H-3;
end
[sol,sol_all,j] = iterateH(filename,Hc,sol,sol_all,j,sigeqM,la,ksig,mu,mub,Ga2,Sig0,xi,Di,ka,sigeq0,Sp,RTol,Nmax,ri,vol0,sp,Vext,a,b,Vr);

% varying xi
xi = sol_all{j,1};
xi0 = sol_all{j-1,1};
d = 0.05;
if xi > xi0
    Mc = xi+d:d:150;
else
    Mc = xi-d:-d:0;
end
[sol,sol_all,j] = iterate_xi(filename,Mc,sol,sol_all,j,la,ksig,mu,mub,Ga2,Sig0,sigeqM,Di,ka,sigeq0,Sp,RTol,Nmax,ri,vol0,sp,Vext,a,b,Vr);

%%%%%% Draw the solution branch in parameter space (Pe-H) %%%%%
figure(1)
fc = sigeq0.^2/(sigeq0.^2+ Sig0.^2); 
xi_scale = Di*mu/(R0.^2*fc);
xi = cell2mat(sol_all(:,1));
H = cell2mat(sol_all(:,3));
plot(xi./xi_scale,H./R0,'color','r','LineWidth',3);
box on;
set(gca,'FontName','Times New Roman','FontSize',36,'LineWidth',3.5,'Layer','top');
xlabel('$Pe$','FontSize',48,'Interpreter','latex');
ylabel('$\widetilde{H}$','FontSize',48,'Interpreter','latex');






function [y1,y2,y3] = iterate_xi(filename,Mc,sol,sol_all,j,la,ksig,mu,mub,Ga2,Sig0,sigeqM,Di,ka,sigeq0,Sp,RTol,Nmax,ri,vol0,sp,Vext,a,b,Vr)
% varying xi
if length(sol.parameters)>3
    sol.parameters(4) = [];
end
for xi = Mc
    yeq = @(u,y,para)shapef(u,y,para,ka,la,ksig,mu,mub,Ga2,Sig0,xi,Di,sigeq0,sigeqM,Sp,sp,Vext,a,b,Vr);
    ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,ri,vol0);
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
    sol0 = sol;
    sol = bvp5c(yeq,ybc,sol,opts);
    if sol.stats.maxerr < RTol
        j = j+1;
        sol_all{j,1} = xi;
        sol_all{j,2} = sigeqM;
        H = sol.y(4,1)-sol.y(4,end);
        sol_all{j,3} = H;
        sol_all{j,4} = max(sol.y(3,:));
        sol_all{j,5} = sol;
        [xi,H]
        save(filename,'sol_all');
    else
        sol = sol0;
        break;
    end
end
y1 = sol;
y2 = sol_all;
y3 = j;
end



function [y1,y2,y3] = iterateH(filename,Hc,sol,sol_all,j,sigeqM,la,ksig,mu,mub,Ga2,Sig0,xi,Di,ka,sigeq0,Sp,RTol,Nmax,ri,vol0,sp,Vext,a,b,Vr)
sol.parameters(4) = xi;
for H = Hc
    yeq = @(u,y,para)shapeHxi(u,y,para,la,ksig,mu,mub,Ga2,Sig0,sigeqM,Di,ka,sigeq0,Sp,sp,Vext,a,b,Vr);
    ybc = @(ya,yb,para)twobcH(ya,yb,para,ka,ri,H,vol0);
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
    sol0 = sol;
    sol = bvp5c(yeq,ybc,sol,opts);
    if sol.stats.maxerr < RTol
        j = j+1;
        xi = sol.parameters(4);
        sol_all{j,1} = xi;
        sol_all{j,2} = sigeqM; 
        sol_all{j,3} = H;
        sol_all{j,4} = max(sol.y(3,:));
        sol_all{j,5} = sol;
        [xi,H]
        save(filename,'sol_all');
    else
        sol = sol0;
        break;
    end
end
y1 = sol;
y2 = sol_all;
y3 = j;
end
