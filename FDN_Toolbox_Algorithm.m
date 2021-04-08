% Algorithm of a Feedback Delay Network for realistic artificial reverb
% Version 1: Using fdn_biquad.m by Piotr Majdak
% © Patrick J. Boettcher
%
% FDN Toolbox by Sebastian J. Schlecht
% https://github.com/SebastianJiroSchlecht/fdnToolbox
%
% SOFAspat by Piotr Majdak & Markus Noisternig as taken from SOFA API for
% MATLAB and OCTAVE
% https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)

% Reseting MATLAB before running this code
clear;clc

%% SET SAMPLE RATE AND NR OF DELAY LINES

% Set sample rate
fs = 48000;
Ts = 1/fs;

mix = 1;
N = 8;

numInput = 1;
numOutput = 1;
inputGain = ones(N,numInput);
outputGain = ones(numOutput,N);
direct = zeros(numOutput,numInput);
impulseResponseLength = fs*2;

%% SET TEMPERATURE AND PHYSICAL PROPERTIES

% set temperature filter order and filter gain
Ttemp = 20; % room temperature (set by user)
hAir = 30; % water vapor in percent
fo = 2;

% Initialize physical properties
Tspeed = 331.4; % speed of sound
Tcoeff = 0.607; % approximation coefficient

% calculate speed of sound in room
speed = Tspeed + (Tcoeff*Ttemp); % 337.47 at 20° C
% calculate speed of sound per sample (Equation 4.7)
speedperframe = (speed / 6000)*fs;

% initialize parameters for delay time calculation from meters
widthroom = 10;
lengthroom = 10;
heightroom = 3;

%% SET DELAY TIMES

% calculate distance from virtual loudspeaker
dl1256 = sqrt((lengthroom/2)^2 + (widthroom/2)^2); % Equation 4.8
dl3478 = sqrt((widthroom/2)^2 + (heightroom/2)^2); % Equation 4.9

dlfactor = [dl1256 dl1256*2 dl3478*3 dl3478*5 dl1256*7 dl1256*11 dl3478*13 dl3478*15];

dl = zeros(N,1);

for ii = 1:N
    dl(ii) = dlfactor(ii);
end

maxVal = max(dl);
dmax = maxVal * 1.46; % maximum delay time
dmaxMS = dmax / speedperframe; % maxmium delay time in milliseconds
maxDelay = ceil(dmaxMS*fs); % calculate max delay time in samples

dltimeMS = zeros(N,1);

for ii=1:N
    dltimeMS(ii) = dl(ii)/speedperframe;
end

for ii=1:N
    dl(ii) = fix(dltimeMS(ii)*fs);
end

diir = zeros(N, length(direct));

for ii=1:N
    diir(ii) = diir(ii) * (1/dlfactor(ii));
end

%% FEEDBACK MATRIX

feedbackMatrix = [0.707165350796018,-0.0816338697839416,0.193318349784906,0.153265374902657,0.138970421707016,-0.202460042687440,-0.585948215762172,0.169566092957095;-0.321315305732734,-0.302295414270738,0.417114728980825,0.469554834589229,0.224528913357556,0.546081140235545,-0.240240097384041,-0.0676340738705417;0.122769719226056,0.623134214681095,0.354365773346232,0.0876027438949645,0.101763605245239,0.232121846324512,0.305194047773877,0.553174669329243;0.429859900404027,-0.0709979469132607,0.446030688929616,0.116574804398857,-0.199131454942079,0.00373263379240895,0.496235104784466,-0.558327994057485;0.212162252817637,0.274340184332953,-0.584782057698353,0.598286581265253,-0.272926667087820,0.293686156317055,-0.0241588453145167,-0.135954329039721;0.384039579803479,-0.453053134686926,-0.296747341649147,-0.296832877975925,0.241629044174692,0.517959591793744,0.297741062954191,0.236159356753872;-0.0621615337546289,-0.476119844485083,0.0585116620778208,0.369304405925896,-0.448888215991392,-0.300705920227837,0.255673305489149,0.521864885308096;-0.0214108931669107,-0.0375053321658477,-0.175134769068075,0.390284154694317,0.738844117662729,-0.402433954450194,0.323966770705284,-0.0484027742024311];

%% FDN TOOLBOX ALGORYTHM

delays = dl';

% absorption filters
centerFrequencies = [ 63, 125, 250, 500, 1000, 2000, 4000, 8000]; % Hz
T60frequency = [1, centerFrequencies fs];
targetT60 = [2; 2; 2.2; 2.3; 2.1; 1.5; 1.1; 0.8; 0.7; 0.7];  % seconds

zAbsorption = zSOS(absorptionGEQ(targetT60, delays, fs),'isDiagonal',true);

% power correction filter
targetPower = [5; 5; 5; 3; 2; 1; -1; -3; -5; -5];  % dB
powerCorrectionSOS = designGEQ(targetPower);
outputFilters = zSOS(permute(powerCorrectionSOS,[3 4 1 2]) .* outputGain);

% compute impulse response
irTimeDomain = dss2impz(impulseResponseLength, delays, feedbackMatrix, inputGain, outputFilters, direct, 'absorptionFilters', zAbsorption);

% compute poles/zeros
[res, pol, directTerm, isConjugatePolePair,metaData] = dss2pr(delays, feedbackMatrix, inputGain, outputFilters, direct, 'absorptionFilters', zAbsorption);
irResPol = pr2impz(res, pol, directTerm, isConjugatePolePair, impulseResponseLength);

difference = irTimeDomain - irResPol;
fprintf('Maximum devation betwen time-domain and pole-residues is %f\n', permute(max(abs(difference),[],1),[2 3 1]));

% compute reverberation times 
[reverberationTimeEarly, reverberationTimeLate, F0, powerSpectrum, edr] = reverberationTime(irTimeDomain, fs);
targetPower = targetPower - mean(targetPower);
powerSpectrum = powerSpectrum - mean(powerSpectrum);

sound(irTimeDomain,fs);

%% plot
figure(1); hold on; grid on;
t = 1:size(irTimeDomain,1);
plot( t, difference(1:end) );
plot( t, irTimeDomain - 2 );
plot( t, irResPol - 4 );
legend('Difference', 'TimeDomain', 'Res Pol')

 
figure(2); hold on; grid on;
plot(T60frequency,targetT60);
plot(F0,reverberationTimeLate);
plot(F0,reverberationTimeEarly);
plot(rad2hertz(angle(pol),fs),slope2RT60(mag2db(abs(pol)), fs),'x');
set(gca,'XScale','log');
xlim([50 fs/2]);
xlabel('Frequency [Hz]')
ylabel('Reverberation Time [s]')
legend({'Target Curve','T60 Late','T60 Early','Poles'})

figure(3); hold on; grid on;
plot(T60frequency, targetPower);
plot(F0,powerSpectrum);
set(gca,'XScale','log');
xlim([50 fs/2]);
legend({'Target','Estimation'})
xlabel('Frequency [Hz]')
ylabel('Power Spectrum [dB]')


%% Test: Impulse Response Accuracy
assert(isAlmostZero(difference,'tol',10^0)) % TODO bad due to extra filters

%% Test: Reverberation Time Accuracy, 5% threshold
testF0 = linspace(65,8000,20);
testTarget = interp1(T60frequency,targetT60, testF0);
testEsti = interp1(F0, reverberationTimeLate, testF0);
assert(all(abs(testEsti ./ testTarget - 1) < 0.1)) % 10% error

%% Test: Power Spectrum Accuracy
testF0 = linspace(65,8000,20);
testTarget = interp1(T60frequency,targetPower, testF0);
testEsti = interp1(F0, powerSpectrum, testF0);
assert(all(abs(testEsti - testTarget) < 3)) % 3dB error

%% PLOT THE IMPULSE RESPONSE

% Plot the impulse response
faxsum=0:1/fs:length(irTimeDomain)*1/fs-1/fs;
  % plot the summed impulse response
figure(3);
plot(faxsum,irTimeDomain);
xlabel('Time (s)');
line ([0 4], [-60 -60],'Color','red','LineStyle','--');
title('The total (summed) impulse response of the FDN');
  % plot the summed impulse response (dB)
figure(5);
spect(fft(irTimeDomain),fs,'b');
title('amplitude spectrum of the total impulse response');

% Export the impulse response
fnOut = 'impulse_response_10_10_3_FDNToolbox.wav';
audiowrite(fnOut,irTimeDomain,fs);
