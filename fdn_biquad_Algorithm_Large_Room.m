% Algorithm of a Feedback Delay Network for realistic artificial reverb
% Version 1: Using fdn_biquad.m by Piotr Majdak
% © Patrick J. Boettcher
%
% SOFAspat by Piotr Majdak & Markus Noisternig as taken from SOFA API for
% MATLAB and OCTAVE
% https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)

% Reseting MATLAB before running this code
clear;clc
% -------------------------------------
% SET SAMPLE RATE AND NR OF DELAY LINES

% Set sample rate
Fs = 48000;
Ts = 1/Fs;

N = 8;

%% SET NORMALIZED FREQUENCIES

NrNormFreq = 6;
NormFreq = [125 250 500 1000 2000 4000];

% Calculate normalized frequencies for acoustic filters
fn = zeros(1,6);

for ii=1:NrNormFreq
    fn(ii) = 2 * (NormFreq(ii)/Fs);
end

% Normalized frequencies
NormFreqVal = [0 fn(1,1) fn(1,2) fn(1,3) fn(1,4) fn(1,5) fn(1,6) 1]; 

% Calculate normalized frequencies for air absorption
fn = zeros(1,3);
NrFreqAirAbs = 3;
fnCalc = [1000 2000 4000];

for ii =1:NrFreqAirAbs
    fn(ii) = 2 *(fnCalc(ii)/Fs);
end

% Normalized frequencies for Air Absorption
NormFreqAir = [0 fn(1,1) fn(1,2) fn(1,3) 1];

%% SET SOFA SETTINGS AND VIRTUAL LOUDSPEAKER POSITIONS

% SOFA settings
database='ari'; HRTFfilename='hrtf_nh4.sofa';

azi = [45 -45 45 -45 135 -135 135 -135];
ele = [-45 -45 45 45 -45 -45 45 45];

% Load the HRTF
fullfn=fullfile(SOFAdbPath, 'database', database, HRTFfilename);
disp(['Loading ' fullfn]); Obj=SOFAload(fullfn);

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
speedperframe = (speed / 6000)*Fs;

% initialize parameters for delay time calculation from meters
widthroom = 50;
lengthroom = 50;
heightroom = 20;

% set Variables for calculating Room Modes
RmAmp = 1;
RmOrder = 2;
RmBW = 50;

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
maxDelay = ceil(dmaxMS*Fs); % calculate max delay time in samples

dltimeMS = zeros(N,1);

for ii=1:N
    dltimeMS(ii) = dl(ii)/speedperframe;
end

for ii=1:N
    dl(ii) = fix(dltimeMS(ii)*Fs);
end

%% SET IMPULSE SIGNAL

% Delta Impulse
delta = repmat([1; zeros(2*Fs,1)],1,N);

% Sine signal
sine = audioread('440Hz_Sine.wav');
sinesig = repmat(sine,N);

in = delta;

%% SET AIR ABSORPTION FILTER

Aa = zeros(1,2+1);
Ab = Aa;

for ii=1:N
    [Aa(ii,:), Ab(ii,:)] = Atmosphere(NrFreqAirAbs,fnCalc,NormFreqAir,fo,hAir,Ttemp,dlfactor(ii));
end

% Array transpose
diir = zeros(length(in),N);

for ii=1:N
   diir(:,ii) = filter(Aa(ii,:), Ab(ii,:), in(:,ii));
end

%% SET FEEDBACK MATRIX

% Set feedback values for the feedback matrix
FeedbackMatrix = [0.707165350796018,-0.0816338697839416,0.193318349784906,0.153265374902657,0.138970421707016,-0.202460042687440,-0.585948215762172,0.169566092957095;-0.321315305732734,-0.302295414270738,0.417114728980825,0.469554834589229,0.224528913357556,0.546081140235545,-0.240240097384041,-0.0676340738705417;0.122769719226056,0.623134214681095,0.354365773346232,0.0876027438949645,0.101763605245239,0.232121846324512,0.305194047773877,0.553174669329243;0.429859900404027,-0.0709979469132607,0.446030688929616,0.116574804398857,-0.199131454942079,0.00373263379240895,0.496235104784466,-0.558327994057485;0.212162252817637,0.274340184332953,-0.584782057698353,0.598286581265253,-0.272926667087820,0.293686156317055,-0.0241588453145167,-0.135954329039721;0.384039579803479,-0.453053134686926,-0.296747341649147,-0.296832877975925,0.241629044174692,0.517959591793744,0.297741062954191,0.236159356753872;-0.0621615337546289,-0.476119844485083,0.0585116620778208,0.369304405925896,-0.448888215991392,-0.300705920227837,0.255673305489149,0.521864885308096;-0.0214108931669107,-0.0375053321658477,-0.175134769068075,0.390284154694317,0.738844117662729,-0.402433954450194,0.323966770705284,-0.0484027742024311];


%% CALCULATE ACOUSTIC TEXTURES

% Small test room
% fun = {'Concrete', 'Concrete', 'Concrete', 'Glass', 'Glass', 'Wood', 'Wood', 'Fabric'};
% Large test room
fun = {'Concrete', 'Concrete', 'Brick', 'Brick', 'Wood', 'Wood', 'Plaster', 'Glass'};

Ncoff = 2; % Number of filter coefficents

ha = zeros(N, Ncoff+1);
hb = ha;

for ii=1:N
    [ha(ii,:), hb(ii,:)] = feval(fun{ii}, NormFreqVal, fo);
end

fb = zeros(1,N);
inDL = zeros(1,N);

% Calculation inside the Feedback Matrix
outrefl = fdn_biquad(diir, FeedbackMatrix, dl, ha', hb');

aziout = zeros(1,N);
eleout = aziout;

out = zeros(1,2);

% Pan through SOFA
[out1,aziout1,eleout1,idx1]=SOFAspat(outrefl(:,1),Obj,azi(:,1),ele(:,1)); disp('Binaural signal 1 rendered');
[out2,aziout2,eleout2,idx2]=SOFAspat(outrefl(:,2),Obj,azi(:,2),ele(:,2)); disp('Binaural signal 2 rendered');
[out3,aziout3,eleout3,idx3]=SOFAspat(outrefl(:,3),Obj,azi(:,3),ele(:,3)); disp('Binaural signal 3 rendered');
[out4,aziout4,eleout4,idx4]=SOFAspat(outrefl(:,4),Obj,azi(:,4),ele(:,4)); disp('Binaural signal 4 rendered');
[out5,aziout5,eleout5,idx5]=SOFAspat(outrefl(:,5),Obj,azi(:,5),ele(:,5)); disp('Binaural signal 5 rendered');
[out6,aziout6,eleout6,idx6]=SOFAspat(outrefl(:,6),Obj,azi(:,6),ele(:,6)); disp('Binaural signal 6 rendered');
[out7,aziout7,eleout7,idx7]=SOFAspat(outrefl(:,7),Obj,azi(:,7),ele(:,7)); disp('Binaural signal 7 rendered');
[out8,aziout8,eleout8,idx8]=SOFAspat(outrefl(:,8),Obj,azi(:,8),ele(:,8)); disp('Binaural signal 8 rendered');

% Sum the output
outall = 0.25*(out1 + out2 + out3 + out4 + out5 + out6 + out7 + out8);

outreflAll = sum(outrefl,2);

% Apply room mode filters
RoomModesBiquad = RoomModes(lengthroom,widthroom,heightroom,speed,Fs,RmAmp,RmOrder,RmBW);
RoomModesSig = RoomModesBiquad(outall);

% Sum the Output for measurements
outSum = sum(RoomModesSig);

%%  CALCULATE CORRELATION COEFFICIENTS

CorrCoeff = zeros(N,N);

for ii=1:N
    Corrvalue1 = corrcoef(outrefl(:,1), outrefl(:,ii));
    Corrvalue2 = corrcoef(outrefl(:,2), outrefl(:,ii));
    Corrvalue3 = corrcoef(outrefl(:,3), outrefl(:,ii));
    Corrvalue4 = corrcoef(outrefl(:,4), outrefl(:,ii));
    Corrvalue5 = corrcoef(outrefl(:,5), outrefl(:,ii));
    Corrvalue6 = corrcoef(outrefl(:,6), outrefl(:,ii));
    Corrvalue7 = corrcoef(outrefl(:,7), outrefl(:,ii));
    Corrvalue8 = corrcoef(outrefl(:,8), outrefl(:,ii));
    CorrCoeff(ii,1) = Corrvalue1(2,1);
    CorrCoeff(ii,2) = Corrvalue2(2,1);
    CorrCoeff(ii,3) = Corrvalue3(2,1);
    CorrCoeff(ii,4) = Corrvalue4(2,1);
    CorrCoeff(ii,5) = Corrvalue5(2,1);
    CorrCoeff(ii,6) = Corrvalue6(2,1);
    CorrCoeff(ii,7) = Corrvalue7(2,1);
    CorrCoeff(ii,8) = Corrvalue8(2,1);
end

% Correlation Coefficient of the two channels
CorrCoeffSum = corrcoef(RoomModesSig(:,1), RoomModesSig(:,2));

% Plot the correlation coefficients of all delay lines
figure(1);
CorrPlot = pcolor(CorrCoeff);

% Plot the cross-correlation coefficient of the two channels
[XCor, lags] = xcorr(RoomModesSig(:,1), RoomModesSig(:,2),30,'normalized');
figure(2);
stem(lags,XCor);

%% CALCULATE MAXIMUM AMPLITUDE VALUES

ampvalue = zeros(N,1);

for ii=1:N
    ampvalue(ii) = max(outrefl(:,ii));
end

%% PLAY THE IMPULSE RESPONSE

% sound(RoomModesSigN,Fs);

if ~exist('dontplay','var'); 
  p=audioplayer(RoomModesSig, Obj.Data.SamplingRate);
  play(p); 
end

%% PLOT THE IMPULSE RESPONSE

% Plot the impulse response
fax=0:1/Fs:length(outrefl)*1/Fs-1/Fs;
faxsum=0:1/Fs:length(outreflAll)*1/Fs-1/Fs;

figure(3);
plot(fax,outrefl);
xlabel('Time (s)');
title('Impulse responses of each delay line'); 
  % plot the summed impulse response
figure(4);
plot(faxsum,etc(outreflAll));
xlabel('Time (s)');
line ([0 4], [-60 -60],'Color','red','LineStyle','--');
title('The total (summed) impulse response of the FDN (dB)');
  % plot the amplitude spectrum of the total impulse response
figure(5);
spect(fft(sum(outall,2)),Fs,'a');
title('amplitude spectrum of the total impulse response');
  % plot the amplitude spectra of each delay to show the filter effect
figure(6);
col='rgbkmc';
for ii=1:N
    spect(fft(outrefl(dl(ii):dl(ii)+200,ii)),Fs,'a',col(mod(ii,length(col))+1));
end
set(gca,'XScale','log');
axis([0 Fs/2 -40 10]);
title('Amp spectrum of the first reflection of each delay line');

%% PROOF THAT SOURCES DO NOT MOVE

% Plot only the first Loudspeaker
time = (1:length(aziout1))/Obj.Data.SamplingRate;

% figure
% subplot(2,1,1);
% plot(time,aziout1); % plot azimuthal trajectory
% ylabel('Azimuth (deg)');
% title('SOFAspat: Trajectory');

% subplot(2,1,2);
% plot(time,eleout1); % plot elevational trajectory
% ylabel('Elevation (deg)');
% xlabel('Time (s)');

% Calculate vector length
% Nplot = length(outSumN);

% Plot the reverberation time in MATLAB
% t = [0:Nplot-1]*Ts;
% plot(t,20*log10(abs(outSumN)));
% line ([0 4], [-60 -60],'Color','red','LineStyle','--');
% axis([0 4 -80 0]);
% title('RT60');

% Export the impulse response
fnOut = 'impulse_response_50_50_20.wav';
audiowrite(fnOut,outreflAll,Fs);