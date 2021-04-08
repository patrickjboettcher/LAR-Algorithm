function [BrickCoeffa, BrickCoeffb] = Brick(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.03 0.03 0.03 0.04 0.05 0.07 1];

% Syntax for function
[BrickCoeffa, BrickCoeffb] = yulewalk(fo,NormFreqVal,amps);
