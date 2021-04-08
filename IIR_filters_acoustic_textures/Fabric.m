function [FabricCoeffa, FabricCoeffb] = Fabric(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.04 0.05 0.11 0.18 0.3 0.35 1];

% Syntax for function
[FabricCoeffa, FabricCoeffb] = yulewalk(fo,NormFreqVal,amps);
