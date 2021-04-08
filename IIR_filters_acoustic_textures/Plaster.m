function [PlasterCoeffa, PlasterCoeffb] = Plaster(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.01 0.02 0.02 0.03 0.04 0.05 1];

% Syntax for function
[PlasterCoeffa, PlasterCoeffb] = yulewalk(fo,NormFreqVal,amps);