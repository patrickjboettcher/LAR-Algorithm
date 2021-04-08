function [RoomModesBiquad] = RoomModes(lengthroom,widthroom,heightroom,speed,Fs,RmAmp,RmOrder,RmBW)

NrRoomModes = 12;
rmfreq = zeros(NrRoomModes,1);

for ii=1:NrRoomModes
    rmfreq(ii) = (2 * ((speed/2*(sqrt(((ii)/lengthroom)^2+((ii)/widthroom)^2+((ii)/heightroom)^2)))/Fs))/(Fs/2);
end

% calculate axial, oblique and tangential room modes (Equation 2.7)
% rm1freq = (2 * ((speed/2*(sqrt((1/lengthroom)^2+(1/widthroom)^2+(1/heightroom)^2)))/Fs))/(Fs/2);
% rm2freq = (2 * ((speed/2*(sqrt((2/lengthroom)^2+(2/widthroom)^2+(2/heightroom)^2)))/Fs))/(Fs/2);
% rm3freq = (2 * ((speed/2*(sqrt((3/lengthroom)^2+(3/widthroom)^2+(3/heightroom)^2)))/Fs))/(Fs/2);
% rm4freq = (2 * ((speed/2*(sqrt((4/lengthroom)^2+(4/widthroom)^2+(4/heightroom)^2)))/Fs))/(Fs/2);
% rm5freq = (2 * ((speed/2*(sqrt((5/lengthroom)^2+(5/widthroom)^2+(5/heightroom)^2)))/Fs))/(Fs/2);
% rm6freq = (2 * ((speed/2*(sqrt((6/lengthroom)^2+(6/widthroom)^2+(6/heightroom)^2)))/Fs))/(Fs/2);
% rm7freq = (2 * ((speed/2*(sqrt((7/lengthroom)^2+(7/widthroom)^2+(7/heightroom)^2)))/Fs))/(Fs/2);
% rm8freq = (2 * ((speed/2*(sqrt((8/lengthroom)^2+(8/widthroom)^2+(8/heightroom)^2)))/Fs))/(Fs/2);
% rm9freq = (2 * ((speed/2*(sqrt((9/lengthroom)^2+(9/widthroom)^2+(9/heightroom)^2)))/Fs))/(Fs/2);
% rm10freq = (2 * ((speed/2*(sqrt((10/lengthroom)^2+(10/widthroom)^2+(10/heightroom)^2)))/Fs))/(Fs/2);
% rm11freq = (2 * ((speed/2*(sqrt((11/lengthroom)^2+(11/widthroom)^2+(11/heightroom)^2)))/Fs))/(Fs/2);
% rm12freq = (2 * ((speed/2*(sqrt((12/lengthroom)^2+(12/widthroom)^2+(12/heightroom)^2)))/Fs))/(Fs/2);

RmAmpVect = ((12:-1:RmAmp)/NrRoomModes)';

rma = zeros(3,NrRoomModes);
rmb = zeros(2,NrRoomModes);

for ii=1:NrRoomModes
   [rma(:,ii), rmb(:,ii)] = designParamEQ(RmOrder,RmAmpVect(ii,:),rmfreq(ii,:),RmBW);
end

% RoomModesBiquadRm = zeros(NrRoomModes);

% 'SOSMatrix' = [1 0.3 0.4 1 0.1 0.2];

RoomModesBiquadRm1 = dsp.BiquadFilter('SOSMatrix',[rma(:,1).',[1,rmb(:,1).']]);
RoomModesBiquadRm2 = dsp.BiquadFilter('SOSMatrix',[rma(:,2).',[1,rmb(:,2).']]);
RoomModesBiquadRm3 = dsp.BiquadFilter('SOSMatrix',[rma(:,3).',[1,rmb(:,3).']]);
RoomModesBiquadRm4 = dsp.BiquadFilter('SOSMatrix',[rma(:,4).',[1,rmb(:,4).']]);
RoomModesBiquadRm5 = dsp.BiquadFilter('SOSMatrix',[rma(:,5).',[1,rmb(:,5).']]);
RoomModesBiquadRm6 = dsp.BiquadFilter('SOSMatrix',[rma(:,6).',[1,rmb(:,6).']]);
RoomModesBiquadRm7 = dsp.BiquadFilter('SOSMatrix',[rma(:,7).',[1,rmb(:,7).']]);
RoomModesBiquadRm8 = dsp.BiquadFilter('SOSMatrix',[rma(:,8).',[1,rmb(:,8).']]);
RoomModesBiquadRm9 = dsp.BiquadFilter('SOSMatrix',[rma(:,9).',[1,rmb(:,9).']]);
RoomModesBiquadRm10 = dsp.BiquadFilter('SOSMatrix',[rma(:,10).',[1,rmb(:,10).']]);
RoomModesBiquadRm11 = dsp.BiquadFilter('SOSMatrix',[rma(:,11).',[1,rmb(:,11).']]);
RoomModesBiquadRm12 = dsp.BiquadFilter('SOSMatrix',[rma(:,12).',[1,rmb(:,12).']]);

% Cascade all Biquad filters
[RoomModesBiquad] = dsp.FilterCascade(RoomModesBiquadRm1, RoomModesBiquadRm2, RoomModesBiquadRm3, RoomModesBiquadRm4, RoomModesBiquadRm5, RoomModesBiquadRm6, RoomModesBiquadRm7, RoomModesBiquadRm8, RoomModesBiquadRm9, RoomModesBiquadRm10, RoomModesBiquadRm11, RoomModesBiquadRm12);