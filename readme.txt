LAR-Algorithm (Late Artificial Reverberation Algorithm)
© Patrick J. Boettcher

This folder contains three versions of the LAR Algorithm.

All versions require the MATLAB functions in the folder IIR_filters_acoustic_textures and the MATLAB function RoomModes.m.

All versions require the SOFA API for headphone binauralisation, which can be found here:
https://github.com/sofacoustics/API_MO

The versions are the following:

1. A version using modDelay.m by Eric Tarr as FDN core (modDelay_Algorithm.m)

This version requires modDelay.m, which modulates the reverberation time for preventing from building up at any frequency.
Just like fdn_biquad_Algorithm.m it makes use of a FDN core with a biquad filter for the acoustic textures in the feedback path.

modDelay.m is taken from "Hack Audio: An Introduction to Computer Programming and Digital Signal Processing in MATLAB" Copyright 2019 Taylor & Francis.

2. A version using fdn_biquad.m (fdn_biquad_Algorithm.m)

This version makes use of a FDN core with a biquad filter for the acoustic textures in the feedback path.
The delay time is not modulated through a LFO here.

3. A version where a part of the LAR algorithm has been included in the FDN Toolbox (FDN_Toolbox_Algorithm.m)
The FDN Toolbox by Sebastion J. Schlecht is required to run this version, which can be found here:
https://github.com/SebastianJiroSchlecht/fdnToolbox

Only the delay times will be calculated here, as the impulse response is calculated inside the FDN core.
Please note that the FDN Toolbax can't process delay times over 3000 samples.