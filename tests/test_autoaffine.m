%% test_autoaffine.m
% Test script for affine transformation estimation on a series of hologram images.
%
% This script loads a variable `referenceImages`, which is a complex double matrix.
% The third dimension of `referenceImages` corresponds to individual images in the sequence.
% These images are preprocessed, a region of interest (ROI) is selected, and
% an affine transformation is estimated between the first image and each subsequent image.
%
% Requirements:
%   - MATLAB with Image Processing Toolbox
%   - tests/holograms_for_testing.mat containing `referenceImages`
%
% Usage:
%   Simply run this script after placing the MAT-file in the specified folder.

%% Initialization
clc; clear; close all;

% Load test holograms (variable: referenceImages)
load("tests\holograms_for_testing.mat");

images = preprocessImages(referenceImages);

roiBox = getROIBox(images{1});

% Optionally: crop all images to the selected ROI
for ind = 1:length(images)
    images{ind} = imcrop(images{ind}, roiBox);
end

%% Estimation
% Initialize affine estimator with the first (cropped) image
estimator = getAffineEstimator(images{1}, "KAZE");
estimator.show_result = true;
estimator.params.Threshold = 0.001;
% Run estimation on subsequent images

tforms = {};
for ind = 2:length(images)
    tform = estimator.getAffineMatrix(images{ind});
    % (use tform as needed)
    tforms{ind} = tform;
end

%% Display results
ind = 2;

outputView = imref2d(size(images{1}));
warpedImage = imwarp(images{ind}, tforms{ind}, 'OutputView', outputView);

imshowpair(images{1}, warpedImage, 'falsecolor');
colormap('jet');
title('Result of affine transformation');

function roiBox = getROIBox(image)
    % Display the first image and let the user draw an ROI
    h = figure;
    imshow(image, []);
    title('Select ROI on the first image');
    % For newer MATLAB versions, use drawrectangle
    hRect = drawrectangle('StripeColor','r');
    % Wait until the user double-clicks or finishes drawing
    wait(hRect);
    % Get the ROI position [x, y, width, height] and round to integer
    roiBox = round(hRect.Position);
    close(h)
end

function images_out = preprocessImages(images)
    % Convert complex-valued hologram data into gradient magnitude images
    images_out = {};
    for ind = 1 : size(images, 3)
        images_out{ind} = (angle(images(:,:,ind)) + pi)/(2*pi);
        images_out{ind} = imgradient(images_out{ind}, 'Sobel');
    end
end
