function [MarbleCoeffa, MarbleCoeffb] = Marble(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [00.01 0.01 0.01 0.01 0.02 0.02 1];

% Syntax for function
[MarbleCoeffa, MarbleCoeffb] = yulewalk(fo,NormFreqVal,amps);