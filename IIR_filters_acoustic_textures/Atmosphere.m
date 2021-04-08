function [AtmosphereCoeffa, AtmosphereCoeffb] = Atmosphere(NrFreqAirAbs,fnCalc,NormFreqAir,fo,hAir,Ttemp,distance1)
% Initialize T0
T0 = 293.15;
e = exp(1); % Eulerian number

% Oxygen relation frequency (Equation 2.3)
FrO = 24 + 4.04 * 10^4 * hAir * ((0.02+hAir)/(0.391+hAir));

% Nitrogen relation frequency (Equation 2.4)
FrN = ((Ttemp)/(T0))^-1/2 * (9 + 280 * hAir * e^(-4.17 * ((Ttemp)/(T0)^-1/3)-1));

AirAbsorb = zeros(1,3);

for ii = 1:NrFreqAirAbs
    AirAbsorb(ii) = (869 * fnCalc(ii)^2 * 1.84 * 10^-11 * ((Ttemp)/(T0))+((Ttemp)/(T0))^-(2/5) * (0.01275 * ((e^-2239.1/Ttemp)/FrO + 1000^2/FrO)) + 0.1068 * ((e^-3352/Ttemp)/(FrN + 1000^2/FrN))* distance1) / 100;
end

% Calculate air absorption coefficient for frequency freqAir (Equation 2.2)
% AirAbsorb1000 = (869 * 1000^2 * 1.84 * 10^-11 * ((Ttemp)/(T0))+((Ttemp)/(T0))^-(2/5) * (0.01275 * ((e^-2239.1/Ttemp)/FrO + 1000^2/FrO)) + 0.1068 * ((e^-3352/Ttemp)/(FrN + 1000^2/FrN))* distance1) / 100;
% AirAbsorb2000 = (869 * 2000^2 * 1.84 * 10^-11 * ((Ttemp)/(T0))+((Ttemp)/(T0))^-(2/5) * (0.01275 * ((e^-2239.1/Ttemp)/FrO + 2000^2/FrO)) + 0.1068 * ((e^-3352/Ttemp)/(FrN + 2000^2/FrN))* distance1) / 100;
% AirAbsorb4000 = (869 * 4000^2 * 1.84 * 10^-11 * ((Ttemp)/(T0))+((Ttemp)/(T0))^-(2/5) * (0.01275 * ((e^-2239.1/Ttemp)/FrO + 4000^2/FrO)) + 0.1068 * ((e^-3352/Ttemp)/(FrN + 4000^2/FrN))* distance1) / 100;

% Linear amplitudes for each frequency
amps = [0 AirAbsorb(1,1) AirAbsorb(1,2) AirAbsorb(1,3) 1];

% Syntax for function
[AtmosphereCoeffa, AtmosphereCoeffb] = yulewalk(fo,NormFreqAir,amps);