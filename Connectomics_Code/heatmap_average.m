function [ap_mean, dv_density, centroid] = heatmap_average(thisMap, alpha, gwin)
%% HEATMAP_AVERAGE
% Blend transparency to white, trim upper region, keep ONLY largest
% connected region, then compute A-P mean and D-V density; smooth/normalize.
% Also overlays the largest region's boundary and centroid on the image.

% --- Blend alpha over white (handles transparency-as-black)
alpha   = double(alpha) / 255;
thisMap = im2double(thisMap);
thisMap = thisMap .* alpha + 1 .* (1 - alpha);

% --- Resize and take one channel
I = imresize(thisMap, 0.5, "box");
I = I(:,:,1);

% --- Mask background & remove upper region
blankThresh = 0.9;
I(I >= blankThresh) = NaN;
I(1:550, :) = NaN;

% --- Keep only the largest connected component
% --- Clean up sparse twigs before region selection
mask = ~isnan(I);

% 1) Remove small isolated bits
px_thresh = 600;
mask = bwareaopen(mask, px_thresh); % remove components <300 px (tune this)

% 2) Smooth jagged edges and fill small holes
mask = imclose(mask, strel('disk', 10));   % close gaps ≤5 px
mask = imfill(mask, 'holes');             % fill small holes inside regions

% 3) Optional: erode a bit to trim thin branches
mask = imerode(mask, strel('disk', 10));   % higher radius = more aggressive
centroid = [NaN NaN];  % default

if any(mask(:))
    CC = bwconncomp(mask, 8);
    areas = cellfun(@numel, CC.PixelIdxList);
    [~, idxMax] = max(areas);

    keep = false(size(mask));
    keep(CC.PixelIdxList{idxMax}) = true;

    % Zero out other components in the intensity map
    I(~keep) = NaN;

    % Centroid and boundary of the largest component
    stats = regionprops(keep, 'Centroid');
    centroid = stats.Centroid;              % [x y]
    B = bwboundaries(keep);                 % after masking, should be one region
    boundaryXY = fliplr(B{1});              % (x,y) from (col,row)
else
    keep = false(size(mask));
    boundaryXY = zeros(0,2);
end

% --- A-P mean (across rows) and D-V density (per row)
ap_mean = median(I, 1, 'omitnan');
if all(isnan(ap_mean))
    ap_mean(:) = 0;
else
    ap_mean(isnan(ap_mean)) = max(ap_mean, [], 'omitnan');
end
dv_density = sum(~isnan(I), 2);

% --- Smooth and normalize
ap_mean    = 1 - normalize(smoothdata(ap_mean,    "gaussian", gwin), "range");
dv_density =      normalize(smoothdata(dv_density, "gaussian", gwin), "range");

% --- Visualize with boundary + centroid overlay
% figure('Color','w'); 
% imshow(I, [], 'InitialMagnification', 'fit'); hold on
% if ~isempty(boundaryXY)
%     plot(boundaryXY(:,1), boundaryXY(:,2), '-', 'LineWidth', 1.5);          % border
% end
% if all(~isnan(centroid))
%     plot(centroid(1), centroid(2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5); % centroid
%     % optional: draw a small circle around centroid (visual cue)
%     viscircles(centroid, 8, 'LineWidth', 0.75);  % requires IPT; remove if undesired
% end
% title('Largest component (border) with centroid'); hold off
end
