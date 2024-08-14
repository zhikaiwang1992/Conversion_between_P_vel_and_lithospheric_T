function [Vsyn] = Porter_convert_temperature_to_vel(Pref, Tref, Kref, Uref, Dref, Ptarg, Tinit, alpha0, ...
                                                                  dK_dP, dK_dT, dU_dP, dU_dT, A, a, d, V, H, R, omega, vps, dMPTdenp, Anelas)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Goes_convert_vel_to_temperature.m
%%
%%  Function: Convert the temperature to velocity at a specific depth given a reference velocity at the
%%            reference pressure and temperature using the model proposed by Goes and Govers (2000,JGR).
%%
%%  Zhikai WANG, zwang@ipgp.fr 
%%
%% last change at 30/03/2019
%%
%%  Reference: S. Goes and R. Govers, 2000, Shallow mantle temperature under Europe from P and S wave tomography,
%%             Journal of Geophysical Research, 105(B5), 11153-11169.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qk = 1000; % Q for bulk wave, from Cammarano et al., 2003
% Qk = 479;  % Q for bulk wave, from Anderson and Given, 1982
if(dMPTdenp == 1)
    dK_dP = dK_dP * exp(alpha0 * (Tinit-Tref));
    dU_dP = dU_dP * exp(alpha0 * (Tinit-Tref));    
end
%% Step 1:  Update compressional and shear moduli and density at (Ptarg, Tinit)
K = Kref + (Tinit-Tref) * dK_dT + (Ptarg-Pref) * dK_dP;
U = Uref + (Tinit-Tref) * dU_dT + (Ptarg-Pref) * dU_dP;
D = Dref * (1 - alpha0*(Tinit-Tref) + (Ptarg-Pref)/K);

%% Step 2: Calculate shear and compressional anelasticity Qu and Qp at temperature (Ptarg, Tinit)
if(Anelas == 1)
    E = H + Ptarg * V;
    B = A * d^a;
    Qu = B * omega^a * exp((a*E)/(R*Tinit));
    
    if(vps > 0)   % vps is constant during iteration
        L = 4/3 / vps^2;
    else          % vp/vs varies during iteration
        Vp = sqrt((K + 4/3 * U)/D);
        Vs = sqrt(U/D);
        vps = Vp/Vs;
        L  = 4/3 / vps^2;
    end
    Qp_1 = (1-L)/Qk + L/Qu;
    Qp   = 1 / Qp_1;
    %% Step 3: calculate the anharmonic and anelastic velocity at (Ptarg, Tinit)
    Vanh  = sqrt((K + 4/3 * U)/D);
    Vsyn  = Vanh * (1 - 2/Qp/tan(pi*a/2));   %% calculated velocity in this iteration
    
    %% Step 4: calculate the anharmonic and anelastic velocity derivatives with respect to temperature
    dV_dT_anh  = 0.5/D/Vanh * (dK_dT + 4/3*dU_dT + Vanh * Vanh * (Dref * alpha0 + Dref * (Ptarg-Pref)/K^2 * dK_dT));
    dV_dT_anel = 1/Qp * a * H / (2 * R * Tinit * Tinit * tan(pi*a/2));
    dV_dT = dV_dT_anh + dV_dT_anel;
else
    %% This part neglect the anelasticity effect.
    Vanh  = sqrt((K + 4/3 * U)/D);    
    dV_dT_anh  = 0.5/D/Vanh * (dK_dT + 4/3*dU_dT + Vanh * Vanh * (Dref * alpha0 + Dref * (Ptarg-Pref)/K^2 * dK_dT));
    dV_dT = dV_dT_anh;
    Vsyn = Vanh;
end

end
