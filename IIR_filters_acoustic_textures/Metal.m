function [MetalCoeffa, MetalCoeffb] = Metal(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.19 0.69 0.99 0.88 0.52 0.27 1];

% Syntax for function
[MetalCoeffa, MetalCoeffb] = yulewalk(fo,NormFreqVal,amps);
