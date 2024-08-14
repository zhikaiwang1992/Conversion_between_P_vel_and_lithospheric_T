%*************************************************************************%
% Script: Main_test_influence_of_anelasticity.m
%
%
% Z.K. WANG @ IPGP. Email: zwang@ipgp.fr
%*************************************************************************%
%% Global configuration
clear; clc; close all;
set(0,'defaultfigurecolor','w');
set(0,'DefaultAxesFontSize',16);

%% Path
path('/home/zkw/Software/Raytracing_tomo/matab_scripts/',path);
path('/home/zkw/Software/Raytracing_tomo/create_vm/',path);
path('/home/zkw/IPGPwork/ILAB/Convert_vel_to_temperature/',path);


%% Parameters
Dwater = 1000; % water density. Unit: kg/m^3
Hwater = 4000; % depth of seafoor, Unit m
Dsedim = 1800; % sediment density. Unit: kg/m^3; 
Hsedim = 1000; % thickness of sediment. Unit: m
Dcrust = 2980; % crustal density. Unit: kg/m^3
Hcrust = 5500; % thickness of crust. Unit: m
Htarg  = 30000; % 

g      = 9.8;  % gravity unit. Unit: m/s^2
Trange = 200:10:1250;

out  = ['temperature_velocity_olivine65_Orthopyroxene31_A148.txt'];


%%
Pref   = 101325;    % Unit: Pa,  1Pa = 1kg/m/s^2
Tref   = 300;       % Unit: K
ncomp  = 5;         % Number of total types of minerals.
Perc   = [ 0.75, 0.21 ,0.035, 0, 0.005];

Dref0  = [3222 , 3198 , 3280 , 3578 , 3565 ];   % Density of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
Kref0  = [129E9, 111E9, 105E9, 198E9, 173E9];   % Compressional module of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
Uref0  = [82E9 , 81E9 , 67E9 , 108E9, 92E9 ];   % Shear module of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
dK_dT0 = [-16E6, -12E6, -13E6, -28E6, -21E6];   % Temperature derivatives of Compressional module of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
dU_dT0 = [-14E6, -11E6, -10E6, -12E6, -10E6];   % Temperature derivatives of shear module of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
dK_dP0 = [4.2  , 6.0  , 6.2  , 5.7  , 4.9  ];   % Pressure derivatives of Compressional module of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
dU_dP0 = [1.4  , 2.0  , 1.7  , 0.8  , 1.4  ];   % Pressure derivatives of shear module of Oivine, Orthopyroxene, Clinopyroxene, Spinel, Garnet.
alpha0 = [0.2010E-4, 0.3817E-4, 0.3206E-4, 0.6969E-4, 0.0991E-4];

dMPTdenp = 0;
Anelas   = 1;          % Whether consider anelasticity. 
                       % 1: Yes; 0: No.

A   = 0.148;  % Unit: dimensionless, from (0.148, Porter et al., 2019; 0.049, Shapiro e al., 2004)
a   = 0.15;    % Unit: dimensionless, from Porter et al., 2019
H   = 500E3;   % Unit: J/mol, from Porter et al., 2019
V   = 20E-6;    % Unit: m^3/mol, from Porter et al., 2019
d   = 0.001;   % Unit: m, from Porter et al., 2019
R   = 8.314;   % Unit: J/(K*mol), From Wiki
fd  = 1.0;     % Unit: Hz. Dominant frequency for approximation

vps  = 0;       % Unit: dimensionless, vp/vs ratio
                % if you want to update the vp/vs ratio in calculation, set 'vps' to 0 or negative.
                % If you want to use a constant vp/vs, set 'vps' to the true vp/vs ratio.

%%
omega = 2*pi*fd;
Trange_K = Trange+273;
Hmantle = Htarg - Hwater - Hsedim -Hcrust;
Ptarg =  g * (Hwater * Dwater + Hsedim * Dsedim + Hcrust * Dcrust + Hmantle * sum(Perc.*Dref0));


%% VRHaveraging


%% 
Vsyn = zeros(size(Trange));
NT   = length(Trange);

for ii = 1:1:NT
    Ttarg = Trange_K(ii);
    
    [Kref, Uref, Dref, dK_dP, dK_dT, dU_dP, dU_dT, alpha]=VRHaveraging(ncomp, Pref, Ptarg, Tref, Ttarg, Perc, Kref0, Uref0, Dref0, dK_dP0, dK_dT0, dU_dP0, dU_dT0, alpha0);
    
    [Vsyn(ii)] = Goes_convert_temperature_to_vel(Pref, Tref, Kref, Uref, Dref, Ptarg, Ttarg, alpha, ...
                                                                  dK_dP, dK_dT, dU_dP, dU_dT, A, a, V, H, R, omega, vps, dMPTdenp, Anelas);
end

figure()
    plot(Trange,Vsyn,'LineWidth',2)
    axis([min(Trange),max(Trange),min(Vsyn)+100,max(Vsyn)+100])
  
fid = fopen(out,'w');
for i=1:1:NT
    fprintf(fid,'%f\t%f\n',Trange(i),Vsyn(i));
end
fclose(fid);
