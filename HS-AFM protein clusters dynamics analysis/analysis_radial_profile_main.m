%% radial profile analysis main: preparation
%%% This script take a mask file, data_mask, as the input for the HS-AFM
%%% protein cluster dynamics analysis. The data_mask file should have all
%%% pixels corresponding to the patch area equivalent to 255 and pixels
%%% characterizing the membrane area equivalent to 0.
%%% MIJI is recommanded for visualizing the analysis (MIJ) but not require.

%%%% input file
data_mask;    % The input mask of the protein cluster pixels

%%%% parameters
res = 0.5;   % resolution, unit: nm/pixel


close all
% MIJ.closeAllWindows();   % MIJI
%% process mask
%%%% create mask
BW = double(data_mask(:, :, :) == 255);   % mask file

%%%% denoise mask
% MIJ.createImage(BW)   % MIJI, display mask

[d1, d2, d3] = size(BW);    % find dimension of mask data
BW = double(imgaussfilt(BW) > 0);   % apply a gaussian filter to the mask to remove noise

%%%% apply walking average (3-frame) to the mask
BWa = zeros(d1, d2, d3-2);
for i = 1:d3-2
    BWa(:, :, i) = sum(BW(:, :, i:i+2), 3);
end
BWa = double(BWa > 2);
BW = BWa;

BW = double(bwmorph3(BW, "fill"));   % remove gap pixels within the mask
% MIJ.createImage(BW);   % MIJI, display mask


%% find patch of interest
%%%% find the patch of interest: should be the largest patch (most pixels)
CC = bwconncomp(BW, 26);    % find pixel connectivity groups
CC = CC.PixelIdxList;
Npixel = zeros(numel(CC), 1);    % find group sizes
for i = 1:numel(CC)
    CCi = CC{i};
    Npixel(i) = numel(CCi);
end
[~, idc] = max(Npixel);    % 'idc' should be the largest patch
CC = CC(idc);    % pixel ids in the largest patch
CC = cell2mat(CC);   % covert pixel ids into array format

%%% find the pixel coordinates in the largest patch
sz = size(BW);   
[XX, YY, ZZ] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));

BW2 = BW.*0;    % final mask file
XX = XX(CC);
YY = YY(CC);
ZZ = ZZ(CC);
for i = 1:numel(CC)
    BW2(XX(i), YY(i), ZZ(i)) = 1;
end

% MIJ.createImage(BW2);   % MIJI, display final mask file

%% run radial profile analysis core analysis
%%%% BW2 should be the mask file undergoing the core analysis
run("analysis_radial_profile_core.m");