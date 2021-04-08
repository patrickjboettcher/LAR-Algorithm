function [PeopleCoeffa, PeopleCoeffb] = People(NormFreqVal,fo)

% Linear amplitudes for each frequency
amps = [0 0.25 0.35 0.42  0.46 0.5 0.5 1];

% Syntax for function
[PeopleCoeffa, PeopleCoeffb] = yulewalk(fo,NormFreqVal,amps);
