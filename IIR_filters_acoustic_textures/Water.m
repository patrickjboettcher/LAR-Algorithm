function [WaterCoeffa, WaterCoeffb] = Water(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.008 0.008 0.013 0.015 0.02 0.025 1];

% Syntax for function
[WaterCoeffa, WaterCoeffb] = yulewalk(fo,NormFreqVal,amps);