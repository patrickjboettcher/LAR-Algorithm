% Algorithm of a Feedback Delay Network for realistic artificial reverb
% Version 2: Using modDelay.m by Eric Tarr
% © Patrick J. Boettcher
%
% modDelay.m by Eric Tarr as published in:
% "Hack Audio: An Introduction to Computer Programming and Digital Signal
% Processing in MATLAB" © 2019 Taylor & Francis.

% SOFAspat by Piotr Majdak & Markus Noisternig as taken from SOFA API for
% MATLAB and OCTAVE
% https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)

% Reseting MATLAB before running this code
clear;clc

%% SETTINGS
% Set sample rate
Fs = 48000;
Ts = 1/Fs;

Nr = 8;

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

dl = zeros(Nr,1);

for ii = 1:Nr
    dl(ii) = dlfactor(ii);
end

maxVal = max(dl);
dmax = maxVal * 1.46; % maximum delay time
dmaxMS = dmax / speedperframe; % maxmium delay time in milliseconds
maxDelay = ceil(dmaxMS*Fs); % calculate max delay time in samples

dltimeMS = zeros(Nr,1);

for ii=1:Nr
    dltimeMS(ii) = dl(ii)/speedperframe;
end

for ii=1:Nr
    dl(ii) = fix(dltimeMS(ii)*Fs);
end

%% SET IMPULSE SIGNAL

% Delta Impulse
delta = repmat([1; zeros(2*Fs,1)],1,Nr);

% Sine signal
sine = audioread('440Hz_Sine.wav');
sinesig = repmat(sine,Nr);

in = delta;

%% SET AIR ABSORPTION FILTER

Aa = zeros(1,2+1);
Ab = Aa;

for ii=1:Nr
    [Aa(ii,:), Ab(ii,:)] = Atmosphere(NrFreqAirAbs,fnCalc,NormFreqAir,fo,hAir,Ttemp,dlfactor(ii));
end

% Array transpose
diir = zeros(length(in),Nr);

for ii=1:Nr
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

ha = zeros(Nr, Ncoff+1);
hb = ha;

for ii=1:Nr
    [ha(ii,:), hb(ii,:)] = feval(fun{ii}, NormFreqVal, fo);
end

%% FEEDBACK MATRIX

% Feed forward gains
b = [0 0 0 0 0 0 0 0];

% Intialize feedback holding variables
fb1 = 0; fb2 = 0; fb3 = 0; fb4 = 0; fb5 = 0; fb6 = 0; fb7 = 0; fb8 = 0;

% Gain to control reverb time
g = .67;

% feedbackmatrix = randomOrthogonal(Nr);

feedbackmatrix = [0.707165350796018,-0.0816338697839416,0.193318349784906,0.153265374902657,0.138970421707016,-0.202460042687440,-0.585948215762172,0.169566092957095;-0.321315305732734,-0.302295414270738,0.417114728980825,0.469554834589229,0.224528913357556,0.546081140235545,-0.240240097384041,-0.0676340738705417;0.122769719226056,0.623134214681095,0.354365773346232,0.0876027438949645,0.101763605245239,0.232121846324512,0.305194047773877,0.553174669329243;0.429859900404027,-0.0709979469132607,0.446030688929616,0.116574804398857,-0.199131454942079,0.00373263379240895,0.496235104784466,-0.558327994057485;0.212162252817637,0.274340184332953,-0.584782057698353,0.598286581265253,-0.272926667087820,0.293686156317055,-0.0241588453145167,-0.135954329039721;0.384039579803479,-0.453053134686926,-0.296747341649147,-0.296832877975925,0.241629044174692,0.517959591793744,0.297741062954191,0.236159356753872;-0.0621615337546289,-0.476119844485083,0.0585116620778208,0.369304405925896,-0.448888215991392,-0.300705920227837,0.255673305489149,0.521864885308096;-0.0214108931669107,-0.0375053321658477,-0.175134769068075,0.390284154694317,0.738844117662729,-0.402433954450194,0.323966770705284,-0.0484027742024311];

% Assign Feedback Matrix
g11 = feedbackmatrix(1,1); g12 = feedbackmatrix(1,2); g13 = feedbackmatrix(1,3); g14 = feedbackmatrix(1,4);
g15 = feedbackmatrix(1,5); g16 = feedbackmatrix(1,6); g17 = feedbackmatrix(1,7); g18 = feedbackmatrix(1,8);
g21 = feedbackmatrix(2,1); g22 = feedbackmatrix(2,2); g23 = feedbackmatrix(2,3); g24 = feedbackmatrix(2,4);
g25 = feedbackmatrix(2,5); g26 = feedbackmatrix(2,6); g27 = feedbackmatrix(2,7); g28 = feedbackmatrix(2,8);
g31 = feedbackmatrix(3,1); g32 = feedbackmatrix(3,2); g33 = feedbackmatrix(3,3); g34 = feedbackmatrix(3,4);
g35 = feedbackmatrix(3,5); g36 = feedbackmatrix(3,6); g37 = feedbackmatrix(3,7); g38 = feedbackmatrix(3,8);
g41 = feedbackmatrix(4,1); g42 = feedbackmatrix(4,2); g43 = feedbackmatrix(4,3); g44 = feedbackmatrix(4,4);
g45 = feedbackmatrix(4,5); g46 = feedbackmatrix(4,6); g47 = feedbackmatrix(4,7); g48 = feedbackmatrix(4,8);
g51 = feedbackmatrix(5,1); g52 = feedbackmatrix(5,2); g53 = feedbackmatrix(5,3); g54 = feedbackmatrix(5,4);
g55 = feedbackmatrix(5,5); g56 = feedbackmatrix(5,6); g57 = feedbackmatrix(5,7); g58 = feedbackmatrix(5,8);
g61 = feedbackmatrix(6,1); g62 = feedbackmatrix(6,2); g63 = feedbackmatrix(6,3); g64 = feedbackmatrix(6,4);
g65 = feedbackmatrix(6,5); g66 = feedbackmatrix(6,6); g67 = feedbackmatrix(6,7); g68 = feedbackmatrix(6,8);
g71 = feedbackmatrix(7,1); g72 = feedbackmatrix(7,2); g73 = feedbackmatrix(7,3); g74 = feedbackmatrix(7,4);
g75 = feedbackmatrix(7,5); g76 = feedbackmatrix(7,6); g77 = feedbackmatrix(7,7); g78 = feedbackmatrix(7,8);
g81 = feedbackmatrix(8,1); g82 = feedbackmatrix(8,2); g83 = feedbackmatrix(8,3); g84 = feedbackmatrix(8,4);
g85 = feedbackmatrix(8,5); g86 = feedbackmatrix(8,6); g87 = feedbackmatrix(8,7); g88 = feedbackmatrix(8,8);

%% BUFFERS & LFOS

buffer = zeros(maxDelay,Nr);

% for ii=1:Nr
%    buffer(ii)= zeros(maxDelay,ii);
% end

% Initialize buffers
buffer1 = zeros(maxDelay,1); buffer2 = zeros(maxDelay,1);
buffer3 = zeros(maxDelay,1); buffer4 = zeros(maxDelay,1);
buffer5 = zeros(maxDelay,1); buffer6 = zeros(maxDelay,1);
buffer7 = zeros(maxDelay,1); buffer8 = zeros(maxDelay,1);

rate = [0.25 0.38 0.4 0.51 0.63 0.75 0.87 0.95];
amp = [5 5 5 5 5 5 5 5];

% Initialize output signal
N = length(in);
out = zeros(N,1);

for n = 1:N
    
    inDL1 = filter(ha(1,:), hb(1,:), diir(n,1)) + fb1;
    inDL2 = filter(ha(2,:), hb(2,:), diir(n,1)) + fb2;
    inDL3 = filter(ha(3,:), hb(3,:), diir(n,1)) + fb3;
    inDL4 = filter(ha(4,:), hb(4,:), diir(n,1)) + fb4;
    inDL5 = filter(ha(5,:), hb(5,:), diir(n,1)) + fb5;
    inDL6 = filter(ha(6,:), hb(6,:), diir(n,1)) + fb6;
    inDL7 = filter(ha(7,:), hb(7,:), diir(n,1)) + fb7;
    inDL8 = filter(ha(8,:), hb(8,:), diir(n,1)) + fb8;

    % Eight Parallel Delay Lines
    [outDL1,buffer1] = modDelay(inDL1,buffer1,Fs,n,dl(1),amp(1),rate(1));
    [outDL2,buffer2] = modDelay(inDL2,buffer2,Fs,n,dl(2),amp(2),rate(2));
    [outDL3,buffer3] = modDelay(inDL3,buffer3,Fs,n,dl(3),amp(3),rate(3));
    [outDL4,buffer4] = modDelay(inDL4,buffer4,Fs,n,dl(4),amp(4),rate(4));
    [outDL5,buffer5] = modDelay(inDL5,buffer5,Fs,n,dl(5),amp(5),rate(5));
    [outDL6,buffer6] = modDelay(inDL6,buffer6,Fs,n,dl(6),amp(6),rate(6));
    [outDL7,buffer7] = modDelay(inDL7,buffer7,Fs,n,dl(7),amp(7),rate(7));
    [outDL8,buffer8] = modDelay(inDL8,buffer8,Fs,n,dl(8),amp(8),rate(8));
    
     % Set output paths
     outrefl1(n,1) = outDL1; outrefl2(n,1) = outDL2; outrefl3(n,1) = outDL3; outrefl4(n,1) = outDL4;
     outrefl5(n,1) = outDL5; outrefl6(n,1) = outDL6; outrefl7(n,1) = outDL7; outrefl8(n,1) = outDL8;
    
    % Calculate crossover feedback
    fb1 = g*(g11*outDL1 + g21*outDL2 + g31*outDL3 + g41*outDL4 + g51*outDL5 + g61*outDL6 + g71*outDL7 + g81*outDL8);
    fb2 = g*(g12*outDL1 + g22*outDL2 + g32*outDL3 + g42*outDL4 + g52*outDL5 + g62*outDL6 + g71*outDL7 + g81*outDL8);
    fb3 = g*(g13*outDL1 + g23*outDL2 + g33*outDL3 + g43*outDL4 + g53*outDL5 + g63*outDL6 + g71*outDL7 + g81*outDL8);
    fb4 = g*(g14*outDL1 + g24*outDL2 + g34*outDL3 + g44*outDL4 + g54*outDL5 + g64*outDL6 + g71*outDL7 + g81*outDL8);
    fb5 = g*(g15*outDL1 + g25*outDL2 + g35*outDL3 + g45*outDL4 + g55*outDL5 + g65*outDL6 + g71*outDL7 + g81*outDL8);
    fb6 = g*(g16*outDL1 + g26*outDL2 + g36*outDL3 + g46*outDL4 + g56*outDL5 + g66*outDL6 + g71*outDL7 + g81*outDL8);
    fb7 = g*(g16*outDL1 + g26*outDL2 + g36*outDL3 + g46*outDL4 + g56*outDL5 + g66*outDL6 + g71*outDL7 + g81*outDL8);
    fb8 = g*(g16*outDL1 + g26*outDL2 + g36*outDL3 + g46*outDL4 + g56*outDL5 + g66*outDL6 + g71*outDL7 + g81*outDL8);

end

%% SOFA

% Pan through SOFA
[out1,azi1,ele1,idx]=SOFAspat(outrefl1,Obj,azi(1),ele(1)); disp('Binaural signal 1 rendered');
[out2,azi2,ele2,idx]=SOFAspat(outrefl2,Obj,azi(2),ele(2)); disp('Binaural signal 2 rendered');
[out3,azi3,ele3,idx]=SOFAspat(outrefl3,Obj,azi(3),ele(3)); disp('Binaural signal 3 rendered');
[out4,azi4,ele4,idx]=SOFAspat(outrefl4,Obj,azi(4),ele(4)); disp('Binaural signal 4 rendered');
[out5,azi5,ele5,idx]=SOFAspat(outrefl5,Obj,azi(5),ele(5)); disp('Binaural signal 5 rendered');
[out6,azi6,ele6,idx]=SOFAspat(outrefl6,Obj,azi(6),ele(6)); disp('Binaural signal 6 rendered');
[out7,azi7,ele7,idx]=SOFAspat(outrefl7,Obj,azi(7),ele(7)); disp('Binaural signal 7 rendered');
[out8,azi8,ele8,idx]=SOFAspat(outrefl8,Obj,azi(8),ele(8)); disp('Binaural signal 8 rendered');

% Sum the output
outall = 0.25*(out1 + out2 + out3 + out4 + out5 + out6 + out7 + out8);

%% ROOMMODES FILTER

% Sum the output
% outall = 0.25*(outrefl(:,1) + outrefl(:,2) + outrefl(:,3) + outrefl(:,4) + outrefl(:,5) + outrefl(:,6) + outrefl(:,7) + outrefl(:,8));

% Apply room mode filters
RoomModesBiquad = RoomModes(lengthroom,widthroom,heightroom,speed,Fs,RmAmp,RmOrder,RmBW);
RoomModesSig = RoomModesBiquad(outall);

% Sum the Output for measurements
outSum = sum(RoomModesSig);

%%  CALCULATE CORRELATION COEFFICIENTS

outrefl = [outrefl1, outrefl2, outrefl3, outrefl4, outrefl5, outrefl6, outrefl7, outrefl8];

outreflAll = sum(outrefl,2);

CorrCoeff = zeros(Nr,Nr);

for ii=1:Nr
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
for ii=1:Nr
    spect(fft(outrefl(dl(ii):dl(ii)+200,ii)),Fs,'a',col(mod(ii,length(col))+1));
end
set(gca,'XScale','log');
axis([0 Fs/2 -40 10]);
title('Amp spectrum of the first reflection of each delay line');

%% SAVE FILE

% Export the impulse response
fnOut = 'impulse_response_50_50_20_outrefl8.wav';
audiowrite(fnOut,outrefl(:,8),Fs);