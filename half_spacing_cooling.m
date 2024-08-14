function [T] = half_spacing_cooling(h, t, T0, Tm, K0)

    %% Convert K from unit m^2/s to km^2/Ma
    K = K0 * 365 * 24 * 3600;       %% Thermal Diffusivity, unit in km^2/Ma
    
    %% Calculate temperature using half-spacing cooling model
    T = T0 + erf(0.5*h/sqrt(K*t))*(Tm-T0);
end