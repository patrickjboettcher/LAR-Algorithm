function [GlassCoeffa, GlassCoeffb] = Glass(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.18 0.06 0.04 0.03 0.02 0.02 1];

% Syntax for function
[GlassCoeffa, GlassCoeffb] = yulewalk(fo,NormFreqVal,amps);
