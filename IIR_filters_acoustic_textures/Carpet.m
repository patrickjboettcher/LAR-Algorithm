function [CarpetCoeffa, CarpetCoeffb] = Carpet(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.01 0.02 0.06 0.015 0.25 0.45 1];

% Syntax for function
[CarpetCoeffa, CarpetCoeffb] = yulewalk(fo,NormFreqVal,amps);
