% SAR_Pnt_Detector.m
% Author: John Summerfield
% Loads a SICD (complex) SAR image from an NTF using NGA's MATLAB SAR toolbox.
% Builds a memory-safe quicklook, checks corner edge lengths (no Mapping Toolbox),
% and reconciles Row/Col image-plane spans with ground-projected.
% ============================================================
%  NOTE ON AXIS ORIENTATION AND METADATA (SICD vs. MATLAB)
% ============================================================
%
%  According to the official SICD specification (NGA STDI-0002):
%     • Row direction  (first dimension)  →  Slant Range
%     • Column direction (second dimension) →  Azimuth (Cross-Range)
%
%  Example (SICD definition):
%     X(row, col) = complex pixel at (range_index, azimuth_index)
%
%  Python tools such as SARpy follow this convention literally:
%     Rows  = Range direction
%     Cols  = Azimuth direction
%     → data.shape = (N_range, N_azimuth)
%
%  MATLAB’s NGA SAR Toolbox (read_complex_data) transposes the array
%  intentionally when loading SICD/NTF data:
%     Rows  = Azimuth direction
%     Cols  = Range direction
%     → size(complex_data) = [N_azimuth, N_range]
%
%  This is done so that:
%     imagesc(x_range, y_az, abs(complex_data))
%     displays Range on the X-axis and Azimuth on the Y-axis
%     Note: Range on the X-axis and cross-range on the Y-axis is common 
%     for oldschool Engineers.  Logically things farther away are in
%     increased El angle, so Range should be in the Y-axis in my openion.   
%
%  Therefore:
%     meta_data.Grid.Col.* → corresponds to Range (horizontal axis)
%     meta_data.Grid.Row.* → corresponds to Azimuth (vertical axis)
%
%  Be careful when comparing to SICD files loaded in Python/SARpy!
%  In MATLAB, the data is effectively transposed relative to SICD spec.
%
% ============================================================


%% --- Setup ---
%Make sure NGA SAR Matlab toolbox is in Path
addpath(genpath('/home/jsummerfield/Matlab_files/Global_functions/MATLAB_SAR-master'));

close all;
clear all;

dynamic_range = 40;
dynamic_range2 = 60;

outputpath = './Images';
if exist(outputpath)==0,
    mkdir(outputpath);
end
%Alaska Airfield
NTF_filename = '/home/jsummerfield/Data/CAPELLA_C15_SP_SICD_HH_20250221063844_20250221063910.ntf';

% NTF_filename = '/home/jsummerfield/Data/CAPELLA_C15_SP_SICD_HH_20241105141453_20241105141500.ntf';

%% --- Read complex SICD ---
[complex_data, meta_data] = read_complex_data(NTF_filename);

% NGA SCID convention: rows = RANGE (slant), cols = AZIMUTH (cross-range)
% But when using NGA's Matlab reader, data is transposed.  The meta data 
% will still use Row to describe range data and Col for cross-range.  Here 
% in MATLAB, row is going to be Az/cross-range data and Col will be range.
Naz = size(complex_data, 1);   % azimuth samples 
Nr  = size(complex_data, 2);   % range samples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sanity check, spotlight (SL) mode SAR imagesare 5km x 5km on the ground.
%I verified by looking at the image corners and measuring the distance
%on the ground.
%Note this is not what I was expecting 
%  Nr*meta_data.Grid.Row.SS/cos(meta_data.SCPCOA.GrazeAng*pi/180) is 5km
%  Naz*meta_data.Grid.Col.SS is 5km
rngSS = meta_data.Grid.Row.SS/cos(meta_data.SCPCOA.GrazeAng*pi/180);
azSS  = meta_data.Grid.Col.SS;


% rng_res = meta_data.Grid.Row.ImpRespWid;
%.5*physconst('lightspeed')/meta_data.RadarCollection.Waveform.WFParameters.TxRFBandwidth
% is  0.2498 m.  This is Slant Range Res,GrazeAng was 69.4863 deg(At least
% for the NTF file I used to generate this code). Ground plane range res is 
% .7539. The Az (Cross-range) res that is in the meta data was 0.0451 (wow).
% rng_res/az_res is 16.7250.  I assume this was setup to do Multi-look Speckle 
% supression.  With this setup, I could average 16 single-look SAR images 
% resulting in a Multi-look SAR image with (aprrox) square pixels and supressed 
% speckle. FYI Speckle is noise-like but not truly noise (in theory if I 
% flew the exact same flight path and pulsed at the exact same positions 
% along the path, I should get the same speckle as long as the scene is 
% identical).

%From waveform bandwidth
%.5*physconst('lightspeed')/meta_data.RadarCollection.Waveform.WFParameters.TxRFBandwidth/cosd(meta_data.SCPCOA.GrazeAng)
%= 0.7129
%From meta data meta_data.Grid.Row.ImpRespWid/cosd(meta_data.SCPCOA.GrazeAng) = 0.7539
% Reasonable, you can lose res with envelope shape.  As long as the
% meta data is in the right ballpark as the waveform bandwidth res 

rng_res = meta_data.Grid.Row.ImpRespWid/cosd(meta_data.SCPCOA.GrazeAng);
az_res =  meta_data.Grid.Col.ImpRespWid;

%Note: rng_res/az_res is 16.7250

%% --- Memory-safe quicklook (single) ---
mag = abs(single(complex_data));
mx  = prctile(mag(:), 99.9);
mag = 20*log10(mag./max(mx, eps('single')) + eps('single'));

% Plot with range on X, azimuth on Y
% x_range = (0:Nr-1)  * rngSS;
% y_az    = (0:Naz-1) * azSS;

% I want the center of the image to be (0,0)
%I'm assumeing that Nr and Naz are even
x_range = (-Nr/2:Nr/2-1)  * rngSS;
y_az    = (-Naz/2:Naz/2-1) * azSS;


% figure; imagesc(x_grd_range, y_az, mag.'); axis xy equal tight;
%az fliped  - looks like Google Earth
fig_sar = figure; imagesc(x_range, y_az, mag); axis xy equal tight;
colormap gray; caxis([-dynamic_range 0]); colorbar;
ylabel('Ground range [m]'); xlabel('Azimuth (cross-range) [m]');
title('Capella SICD quicklook dB');
exportgraphics(fig_sar,[outputpath '/SAR_image.png'],'Resolution',300);

%


fil_sar = speckle_filter(single(complex_data), 'radius', 2, 'ENL', 1, 'filter', 'Lee');
fil_sar = fil_sar/max(abs(fil_sar(:)));

fig_Speckle = figure; imagesc(x_range, y_az, 20*log10(abs(fil_sar) + eps('single')),[-dynamic_range2 0]); axis xy equal tight;
colormap gray; colorbar;
ylabel('Ground range [m]'); xlabel('Azimuth (cross-range) [m]');
title('Speckle Filtered SAR image dB');
exportgraphics(fig_Speckle,[outputpath '/Speckle_Filtered_SAR_image.png'],'Resolution',300);

% --- SVD point-target detection (rows=azimuth, cols=range) ---

lin = abs(fil_sar);  % use your filtered image
lin_scale = prctile(lin(:), 99.9);
lin = lin ./ max(lin_scale, eps('single'));

patch   = 64;
stride  = 16;
minSep  = 32;
Kmax    = 300;
pctThr  = 99.7;

nr = Naz;          % rows = az
nc = Nr;           % cols = range
r_starts = 1:stride:(nr - patch + 1);
c_starts = 1:stride:(nc - patch + 1);
R = numel(r_starts);
C = numel(c_starts);

svd_metric  = zeros(R, C, 'single');

% NEW: arrays to store brightest-pixel location per window
bright_row  = zeros(R, C, 'uint32');   % absolute row in full image
bright_col  = zeros(R, C, 'uint32');   % absolute col in full image

for ii = 1:R
    r0 = r_starts(ii);
    slab = lin(r0:r0+patch-1, :);   % rows for this strip
    for jj = 1:C
        c0 = c_starts(jj);
        P = slab(:, c0:c0+patch-1); % patch (az x rng)

        % 1) SVD metric
        s = svd(P, 'econ');
        svd_metric(ii, jj) = s(1) / sum(s);

        % 2) Brightest point *inside this window*
        [~, idxMax] = max(P(:));
        [rLocal, cLocal] = ind2sub([patch, patch], idxMax);

        % Convert to ABSOLUTE image coords
        bright_row(ii, jj) = uint32(r0 + rLocal - 1);
        bright_col(ii, jj) = uint32(c0 + cLocal - 1);
    end
end

% Threshold by high percentile (CFAR-like)
tau = prctile(svd_metric(:), pctThr);
cand = svd_metric >= tau;

% Gather candidates, sort by score
[ii_win, jj_win] = find(cand);
scores = svd_metric(cand);
[score_sorted, order] = sort(scores, 'descend');

ii_win = ii_win(order);
jj_win = jj_win(order);

% Use brightest-pixel coords instead of window center
row_c = zeros(numel(ii_win),1);
col_c = zeros(numel(ii_win),1);
for k = 1:numel(ii_win)
    row_c(k) = double(bright_row(ii_win(k), jj_win(k)));
    col_c(k) = double(bright_col(ii_win(k), jj_win(k)));
end

keep = false(size(row_c));
for k = 1:numel(row_c)
    if ~keep(k)
        keep(k) = true;
        dr = row_c - row_c(k);
        dc = col_c - col_c(k);
        near = (dr.^2 + dc.^2) <= (minSep^2);
        near(k) = false;
        keep(near) = false;
    end
    if nnz(keep) >= Kmax, break; end
end

row_pk   = row_c(keep);
col_pk   = col_c(keep);
score_pk = score_sorted(keep);

%% --- Overlay detections on Speckle Filtered SAR image ---

% Convert pixel detections to meters (your centered axes)
x_m = ((col_pk - 1) - Nr/2)  .* rngSS;   % range (horizontal, x)
y_m = ((row_pk - 1) - Naz/2) .* azSS;    % azimuth (vertical, y)

% Reuse the speckle figure if it still exists
figure(fig_Speckle);
hold on;

% Plot detection markers
plot(x_m, y_m, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);

% Annotate with scores
for k = 1:numel(x_m)
    text(x_m(k), y_m(k), sprintf('  %.3f', score_pk(k)), ...
        'Color', 'y', 'FontSize', 8, 'VerticalAlignment', 'middle');
end

title(sprintf('Speckle Filtered SAR + SVD detections (patch=%d, stride=%d)', patch, stride));

% Save updated overlay
exportgraphics(gcf, [outputpath '/Speckle_Filtered_SAR_with_SVD_Detects.png'], 'Resolution', 300);


%% --- Chip the top 16 detections (length m × length m) and save as PNG ---

% How many chips to show
Nchips = min(16, numel(score_pk));
chip_length_m = 15;

% Window size in pixels that corresponds to chip_length_m in each dimension
win_rng_px = max(1, round(chip_length_m / rngSS));   % cols (range)
win_az_px  = max(1, round(chip_length_m / azSS));    % rows (az)

% Ensure odd sizes
if mod(win_rng_px,2) == 0, win_rng_px = win_rng_px + 1; end
if mod(win_az_px,2)  == 0, win_az_px  = win_az_px  + 1; end

half_rng = floor(win_rng_px/2);
half_az  = floor(win_az_px/2);

% Work from RAW magnitude to avoid speckle filter bias
lin_raw_mag = abs(single(complex_data));
[rows, cols] = size(lin_raw_mag);

% pad (replicate) so edge chips work
topPad    = repmat(lin_raw_mag(1, :),  half_az, 1);
bottomPad = repmat(lin_raw_mag(end,:), half_az, 1);
temp      = [topPad; lin_raw_mag; bottomPad];

leftPad  = repmat(temp(:,1), 1, half_rng);
rightPad = repmat(temp(:,end), 1, half_rng);
pad_mag  = [leftPad, temp, rightPad];

% figure + tiled layout FIRST
chipsFig = figure('Position',[100 100 900 900]);
t = tiledlayout(chipsFig, 4, 4, 'Padding','compact', 'TileSpacing','compact');

% preallocate metric arrays
IPR_rng_HPBW = zeros(Nchips,1);
IPR_rng_PSLR = zeros(Nchips,1);
IPR_rng_ISLR = zeros(Nchips,1);
IPR_az_HPBW  = zeros(Nchips,1);
IPR_az_PSLR  = zeros(Nchips,1);
IPR_az_ISLR  = zeros(Nchips,1);

for ii = 1:Nchips
    % center in original image
    r = row_pk(ii);
    c = col_pk(ii);

    % shift into padded coords
    r0 = r + half_az;
    c0 = c + half_rng;

    % extract chip
    r1 = r0 - half_az; r2 = r0 + half_az;
    c1 = c0 - half_rng; c2 = c0 + half_rng;
    chip_lin = pad_mag(r1:r2, c1:c2);

    % axes in meters, centered
    x_m_axis = linspace(-chip_length_m/2, chip_length_m/2, win_rng_px);  % range
    y_m_axis = linspace(-chip_length_m/2, chip_length_m/2, win_az_px);   % az

    % ---- IPR estimation ----
    [peak_lin, idx_peak] = max(chip_lin(:)); %#ok<ASGLU>
    [peak_row, peak_col] = ind2sub(size(chip_lin), idx_peak);

    rng_cut = chip_lin(peak_row, :);
    az_cut  = chip_lin(:, peak_col).';

    rng_metrics = ipr_metrics_1d(rng_cut, x_m_axis, 5);
    az_metrics  = ipr_metrics_1d(az_cut,  y_m_axis, 5);

    IPR_rng_HPBW(ii) = rng_metrics.HPBW_m;
    IPR_rng_PSLR(ii) = rng_metrics.PSLR_dB;
    IPR_rng_ISLR(ii) = rng_metrics.ISLR_dB;
    IPR_az_HPBW(ii)  = az_metrics.HPBW_m;
    IPR_az_PSLR(ii)  = az_metrics.PSLR_dB;
    IPR_az_ISLR(ii)  = az_metrics.ISLR_dB;

    % dB view
    chip_db = 20*log10(chip_lin ./ (max(chip_lin(:)) + eps('single')) + eps('single'));

    % === PLOT ===
    ax = nexttile(t); %#ok<LAXES>
    imagesc(ax, x_m_axis, y_m_axis, chip_db);
    axis(ax, 'xy', 'equal', 'tight');
    colormap(ax, gray);
    caxis(ax, [-dynamic_range2 0]);

    xlabel(ax, 'Range [m]',   'FontSize', 6);
    ylabel(ax, 'Azimuth [m]', 'FontSize', 6);

    title(ax, { ...
        sprintf('#%d  score=%.3f', ii, score_pk(ii)), ...
        sprintf('(x=%.1f m, y=%.1f m)', x_m(ii), y_m(ii)) ...
    }, 'FontSize', 6, 'FontWeight', 'normal');
end

sgtitle(t, sprintf('SVD Point Detections – %dm × %dm Chips', 10, 10), ...
    'FontSize', 12, 'FontWeight', 'bold');

exportgraphics(chipsFig, [outputpath '/svd_top_16_chips_10m.png'], 'Resolution', 200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Table with metrics

ChipIdx = (1:Nchips).';
IPR_Table = table(ChipIdx, ...
    score_pk(1:Nchips), ...
    IPR_rng_HPBW, IPR_rng_PSLR, IPR_rng_ISLR, ...
    IPR_az_HPBW,  IPR_az_PSLR,  IPR_az_ISLR, ...
    'VariableNames', { ...
        'Chip', 'SVDscore', ...
        'Rng_HPBW_m', 'Rng_PSLR_dB', 'Rng_ISLR_dB', ...
        'Az_HPBW_m',  'Az_PSLR_dB',  'Az_ISLR_dB' ...
    });

disp(IPR_Table);

