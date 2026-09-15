% Improved full script for pO2 classification and visualization
clear; close all; clc;

% Folder containing images and .mat files
folderPath = 'D:\Partners HealthCare Dropbox\Nancy Ruiz Uribe\Harvard_Postdoc\2_TissuePO2\pO2_results\rTg4510_TissuePO2\selected for tangle analysis\z200';

% Get list of TIFF images
imageFiles = dir(fullfile(folderPath, '*.tif'));

% Initialize arrays to collect per-image mean pO2 values by classification
all_po2_positive = [];
all_po2_negative = [];

% Parameters for neighborhood mask
radius = 15;
[xx, yy] = meshgrid(-radius:radius, -radius:radius);
circularMask = (xx.^2 + yy.^2) <= radius^2;

for imgIdx = 1:length(imageFiles)
    % Load image
    imagePath = fullfile(folderPath, imageFiles(imgIdx).name);
    refImage = imread(imagePath);

    % Base name for associated .mat files
    [~, baseName, ~] = fileparts(imagePath);

    % Load associated .mat files (imported points & pO2 values)
    importedMat = load(fullfile(folderPath, ['imported_' baseName '.mat']));
    po2Mat = load(fullfile(folderPath, ['pO2_imported_' baseName '.mat']));

    % Convert to grayscale if RGB or RGBA
    if size(refImage,3) == 4
        refImage = refImage(:,:,1:3); % Drop alpha
    end
    if size(refImage,3) == 3
        grayImage = rgb2gray(refImage);
    else
        grayImage = refImage;
    end

    % Normalize image to [0,1]
    normImage = rescale(double(grayImage));

    % Adaptive threshold (Otsu)
    level = graythresh(normImage);

    % Binary mask
    binaryMask = imbinarize(normImage, level);

    % Extract point coordinates (scaled to image size and rounded)
    pointsX = round([importedMat.processed.xml.Sequence.PVFLIMPointScan.PVFLIMPointScanElement.ImagingPoints.PVGalvoPointElement.Point.X] * size(grayImage,2));
    pointsY = round([importedMat.processed.xml.Sequence.PVFLIMPointScan.PVFLIMPointScanElement.ImagingPoints.PVGalvoPointElement.Point.Y] * size(grayImage,1));
    
    % Clip points to image boundaries
    pointsX = max(min(pointsX, size(grayImage,2)), 1);
    pointsY = max(min(pointsY, size(grayImage,1)), 1);

    % Extract pO2 values (vectorized)
    pO2Values = reshape(po2Mat.pO2.pO2Value, [], 1);

    numPoints = length(pointsX);
    intensityValues = zeros(numPoints,1);

    % Compute local mean intensity around each point using circular mask
    for k = 1:numPoints
        x = pointsX(k);
        y = pointsY(k);

        % Bounding box indices clipped to image size
        x_min = max(1, x - radius);
        x_max = min(size(normImage, 2), x + radius);
        y_min = max(1, y - radius);
        y_max = min(size(normImage, 1), y + radius);

        % Extract image patch and corresponding mask patch
        patch = normImage(y_min:y_max, x_min:x_max);
        mask_patch = circularMask( (y_min - (y - radius) + 1):(y_max - (y - radius) + 1), ...
                                   (x_min - (x - radius) + 1):(x_max - (x - radius) + 1) );

        % Calculate mean intensity within circular mask
        intensityValues(k) = mean(patch(mask_patch));
    end

    % Classify points by intensity threshold
    isPositive = intensityValues >= level;

    po2_positive = pO2Values(isPositive);
    po2_negative = pO2Values(~isPositive);

    % Calculate means safely (handle empty sets)
    meanPositive = NaN; meanNegative = NaN;
    if ~isempty(po2_positive)
        meanPositive = mean(po2_positive);
    end
    if ~isempty(po2_negative)
        meanNegative = mean(po2_negative);
    end

    % Append per-image means
    all_po2_positive = [all_po2_positive; meanPositive];
    all_po2_negative = [all_po2_negative; meanNegative];

    %% Visualization of current image with classification overlay

    figure('Name', ['Image: ' baseName], 'NumberTitle', 'off');
    imshow(binaryMask);
    hold on;
    axis image;
    set(gca, 'YDir', 'reverse');

    % Scatter plot points colored by classification
    scatter(pointsX(isPositive), pointsY(isPositive), 30, 'r', 'filled');
    scatter(pointsX(~isPositive), pointsY(~isPositive), 30, 'b', 'filled');

    % Overlay white circles around all points
    theta = linspace(0, 2*pi, 50);
    for k = 1:numPoints
        x_circle = pointsX(k) + radius*cos(theta);
        y_circle = pointsY(k) + radius*sin(theta);
        plot(x_circle, y_circle, 'w-', 'LineWidth', 1.5);
    end
    title(sprintf('Classification Overlay: %s', baseName));
    legend({'Positive (Red)', 'Negative (Blue)'}, 'Location', 'bestoutside');
    hold off;
end

%% Aggregate plot: boxplot with scatter overlay for all images

% Remove NaN values (images without points classified as positive or negative)
validPos = ~isnan(all_po2_positive);
validNeg = ~isnan(all_po2_negative);

all_po2_positive = all_po2_positive(validPos);
all_po2_negative = all_po2_negative(validNeg);

all_po2 = [all_po2_positive; all_po2_negative];
groupLabels = [repmat({'Tangle-Positive'}, length(all_po2_positive), 1); 
               repmat({'Tangle-Negative'}, length(all_po2_negative), 1)];

figure('Name', 'pO2 Classification Summary', 'NumberTitle', 'off');
boxplot(all_po2, groupLabels, 'Colors', ['r' 'b']);
ylabel('Average pO2 Value');
title('Average pO2 by Tangle Classification Across Images');
hold on;

g = double(categorical(groupLabels));
scatter(g + 0.15*(rand(size(g))-0.5), all_po2, 40, 'filled', 'MarkerFaceAlpha', 0.6);

grid on;

%% Statistical Test: Wilcoxon Signed-Rank Test

% Ensure vectors are same length and paired correctly (remove unmatched NaNs)
minLen = min(length(all_po2_positive), length(all_po2_negative));
[~, sortPosIdx] = sort(all_po2_positive);
[~, sortNegIdx] = sort(all_po2_negative);

posData = all_po2_positive(sortPosIdx(1:minLen));
negData = all_po2_negative(sortNegIdx(1:minLen));

[p, h, stats] = signrank(posData, negData);

fprintf('Wilcoxon Signed-Rank Test Results:\n');
fprintf('----------------------------------\n');
fprintf('p-value       = %.4f\n', p);
fprintf('Test decision = %d (1 = reject null, 0 = fail to reject)\n', h);
fprintf('Signed rank statistic (W) = %d\n', stats.signedrank);

