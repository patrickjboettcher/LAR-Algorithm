function [WoodCoeffa, WoodCoeffb] = Wood(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.15 0.11 0.1 0.07 0.06 0.07 1];

% Syntax for function
[WoodCoeffa, WoodCoeffb] = yulewalk(fo,NormFreqVal,amps);
