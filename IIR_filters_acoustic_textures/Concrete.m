function [ConcreteCoeffa, ConcreteCoeffb] = Concrete(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.01 0.02 0.04 0.06 0.08 0.1 1];

% Syntax for function
[ConcreteCoeffa, ConcreteCoeffb] = yulewalk(fo,NormFreqVal,amps);