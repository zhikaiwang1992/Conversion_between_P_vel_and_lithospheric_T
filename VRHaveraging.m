function [Kref, Uref, Dref, dK_dP, dK_dT, dU_dP, dU_dT, alpha]=VRHaveraging(N, Pref, P, Tref, T, Perc, Kref0, Uref0, Dref0, dK_dP0, dK_dT0, dU_dP0, dU_dT0, alpha0)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VRHaveraging.m
%%
%%  Function: Calculate elastic parameters for a combination of mineral using 
%%            laboratory parameters for a pure singe mineral.
%%
%%  Zhikai WANG, zwang@ipgp.fr 
%%
%% last change at 30/03/2019
%%
%%  Reference: S. Goes and R. Govers, 2000, Shallow mantle temperature under Europe from P and S wave tomography,
%%             Journal of Geophysical Research, 105(B5), 11153-11169.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( nargin ~=14 )
    error('The input parameters are not enough');
end

if( (length(Kref0)~=N) || (length(Uref0)~=N) || (length(Dref0)~=N) )
    error('Elements of Kref or Uref or Dref are not N')
end

if( (length(dK_dP0)~=N) || (length(dK_dT0)~=N) || (length(dU_dP0)~=N) || (length(dU_dT0)~=N) )
    error('Elements of derivatives are not N')
end

if( (length(alpha0)~=N) )
    error('Elements of alpha0 are not N');
end


Nperc = length(Perc);
if(Nperc > N)
    Perc = Perc(1:N);
elseif(Nperc < N)
    Perc(Nperc+1:N) = 0;
end
Nsum = sum(Perc);
if(Nsum ~=1)
    warning('The summation of is not 1');
end

%%
alpha = sum(alpha0 .* Perc);

Dref_PT = Dref0 .* (1-alpha0*(T-Tref)+(P-Pref)./Kref0);
Dref = sum(Dref_PT .* Perc);

%%
Kvoigt = sum(Kref0 .* Perc);
Kreuss = 1/sum(Perc./Kref0);
Uvoigt = sum(Uref0 .* Perc);
Ureuss = 1/sum(Perc./Uref0);
Kref = (Kvoigt + Kreuss)/2;
Uref = (Uvoigt + Ureuss)/2;

%%
dK_dP = sum(dK_dP0.*Perc) + 1/Kreuss/Kreuss * sum(Perc.*dK_dP0./Kref0./Kref0);
dK_dT = sum(dK_dT0.*Perc) + 1/Kreuss/Kreuss * sum(Perc.*dK_dT0./Kref0./Kref0);
dU_dP = sum(dU_dP0.*Perc) + 1/Ureuss/Ureuss * sum(Perc.*dU_dP0./Uref0./Uref0);
dU_dT = sum(dU_dT0.*Perc) + 1/Ureuss/Ureuss * sum(Perc.*dU_dT0./Uref0./Uref0);

end