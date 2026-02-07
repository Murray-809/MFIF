function [R1_all, R2_all, wmap_all, alpha_all] = corrcoef_group_fusion_open()
close all;

% ---------- parameters ----------
input_path   = '';
output_path  = '';

GROUP_SIZE   = 10;
td           = 0.9;
r_corr_iter1 = 9;
r_corr_iter2 = 5;
gauss_sigma  = 2;
gauss_fsz    = 5;
gf_r         = 5;
gf_eps       = 0.3;
T_low        = 0.3;
T_high       = 1 - T_low;
area_ratio   = 0.1;
R_SUB        = 0.9;

if ~exist(output_path, 'dir'), mkdir(output_path); end

R1_all   = {};
R2_all   = {};
wmap_all = {};
alpha_all= {};

image_files = dir(fullfile(input_path, '*.bmp'));
if isempty(image_files), error('No .bmp files found in input_path'); end

file_numbers = zeros(length(image_files), 1);
for i = 1:length(image_files)
    [~, name, ~] = fileparts(image_files(i).name);
    nums = regexp(name, '\d+', 'match');
    if ~isempty(nums), file_numbers(i) = str2double(nums{1}); else, file_numbers(i) = i; end
end
[~, sort_idx] = sort(file_numbers);
image_files = image_files(sort_idx);

num_images = length(image_files);
all_images_orig = cell(num_images, 1);
all_images_gray = cell(num_images, 1);

for i = 1:num_images
    img_orig = imread(fullfile(input_path, image_files(i).name));
    all_images_orig{i} = img_orig;
    if size(img_orig,3) == 3
        all_images_gray{i} = single(rgb2gray(img_orig));
    else
        all_images_gray{i} = single(img_orig);
    end
end

num_groups = ceil(num_images / GROUP_SIZE);
group_fused_images = cell(num_groups, 1);

for g = 1:num_groups
    start_idx = (g-1) * GROUP_SIZE + 1;
    end_idx   = min(g * GROUP_SIZE, num_images);
    group_fused_images{g} = group_fusion(all_images_gray(start_idx:end_idx), ...
        all_images_orig(start_idx:end_idx));
end

current_result = group_fused_images{1};

for g = 2:num_groups
    [current_result, R1, R2, wmap, alpha] = pairwise_fusion_closed_form( ...
        current_result, group_fused_images{g}, ...
        td, r_corr_iter1, r_corr_iter2, gauss_sigma, gauss_fsz, ...
        gf_r, gf_eps, T_low, T_high, area_ratio, R_SUB);

    R1_all{end+1}   = R1;
    R2_all{end+1}   = R2;
    wmap_all{end+1} = wmap;
    alpha_all{end+1}= alpha;
end

imwrite(current_result, fullfile(output_path, 'final_fused_result.tif'));

end

function [fused_result, R1_out, R2_out, wmap_out, alpha_out] = pairwise_fusion_closed_form( ...
    I1_orig, I2_orig, td, r1, r2, gsigma, gfsz, gfr, gfeps, T_low, T_high, area_ratio, R_SUB)

if size(I1_orig,3) == 3
    I1 = single(rgb2gray(I1_orig));
else
    I1 = single(I1_orig);
end
if size(I2_orig,3) == 3
    I2 = single(rgb2gray(I2_orig));
else
    I2 = single(I2_orig);
end

D = (I1 + I2) * 0.5;

num_iterations = 2;
final_alpha = [];
gmap = [];

for iter = 1:num_iterations
    r_corr = (iter == 1) * r1 + (iter == 2) * r2;

    R1_raw = fast_correlation_integral(D, I1, r_corr, 1e-10);
    R2_raw = fast_correlation_integral(D, I2, r_corr, 1e-10);

    R1_s = imgaussfilt(R1_raw, gsigma, 'FilterSize', gfsz);
    R2_s = imgaussfilt(R2_raw, gsigma, 'FilterSize', gfsz);

    R1_t = R1_s .* (R1_s > td);
    R2_t = R2_s .* (R2_s > td);
    wmap = single(R1_t > R2_t);

    wmap(I2 < 5 & I1 >= 5) = 1;
    wmap(I1 < 5 & I2 >= 5) = 0;

    guide = (I1 / 255);
    gmap  = guidedfilter_LKH(guide, wmap, gfr, gfeps);
    gmap  = max(min(gmap,1),0);

    trimap = 0.5 * ones(size(gmap), 'single');
    trimap(gmap > T_high) = 1;
    trimap(gmap < T_low)  = 0;

    consts_map  = (trimap ~= 0.5);
    consts_vals = double(trimap == 1);

    I_matte_guide = I1_orig;
    if size(I_matte_guide,3) == 1
        I_matte_guide = repmat(I_matte_guide, [1 1 3]);
    end
    alpha = close_matte_w(I_matte_guide, consts_map, consts_vals);

    area = ceil(area_ratio * size(I1,1) * size(I1,2));
    final_alpha = double(bwareaopen(alpha >= 0.5, area));

    fused_gray = double(I1).*final_alpha + double(I2).*(1-final_alpha);
    final_alpha = guidedfilter_LKH(single(fused_gray/255), single(final_alpha), gfr, gfeps);
    final_alpha = double(max(min(final_alpha,1),0));

    D = single(double(I1).*final_alpha + double(I2).*(1-final_alpha));

    if iter == num_iterations
        R1_out = R1_raw - R_SUB;
        R2_out = R2_raw - R_SUB;
        wmap_out = wmap;
        alpha_out = final_alpha;
    end
end

alpha_up = imresize(single(alpha_out), [size(I1_orig,1), size(I1_orig,2)], 'nearest');
fused_result = uint8(double(I1_orig) .* repmat(double(alpha_up), [1,1,3]) + ...
                     double(I2_orig) .* repmat(double(1-alpha_up), [1,1,3]));
end

function fused_image = group_fusion(group_gray, group_orig)
num_imgs = numel(group_gray);
avg_img = zeros(size(group_gray{1}), 'single');
for i = 1:num_imgs
    avg_img = avg_img + group_gray{i};
end
avg_img = avg_img / num_imgs;

[h,w] = size(avg_img);
corr_maps = zeros(h,w,num_imgs,'single');
for i = 1:num_imgs
    corr_maps(:,:,i) = fast_correlation_integral(avg_img, group_gray{i}, 5, 1e-10);
end
[~, max_idx] = max(corr_maps, [], 3);

[oh, ow, ~] = size(group_orig{1});
max_idx_up = imresize(max_idx, [oh, ow], 'nearest');

fused_image = zeros(oh, ow, 3, 'uint8');
for i = 1:num_imgs
    mask = repmat(max_idx_up == i, [1 1 3]);
    fused_image(mask) = group_orig{i}(mask);
end
end

function corr_map = fast_correlation_integral(img1, img2, r, eps_ecc)
img1 = single(img1);
img2 = single(img2);

[h,w] = size(img1);
N  = boxfilter_fast(ones(h,w,'single'), r);

m1 = boxfilter_fast(img1, r) ./ N;
m2 = boxfilter_fast(img2, r) ./ N;

cov12 = boxfilter_fast(img1 .* img2, r) ./ N - m1 .* m2;

v1 = max(0, boxfilter_fast(img1 .* img1, r) ./ N - m1 .* m1);
v2 = max(0, boxfilter_fast(img2 .* img2, r) ./ N - m2 .* m2);

corr_map = cov12 ./ (sqrt(v1).*sqrt(v2) + eps_ecc);
corr_map = max(-1, min(1, corr_map));
end

function q = guidedfilter_LKH(I, p, r, eps)
I = single(I);
p = single(p);

[hei, wid] = size(I);
N = boxfilter_fast(ones(hei, wid,'single'), r);

mI  = boxfilter_fast(I, r) ./ N;
mp  = boxfilter_fast(p, r) ./ N;
mIp = boxfilter_fast(I .* p, r) ./ N;
covIp = mIp - mI .* mp;

mII = boxfilter_fast(I .* I, r) ./ N;
varI = mII - mI .* mI;

a = covIp ./ (varI + eps);
b = mp - a .* mI;

q = (boxfilter_fast(a, r) ./ N) .* I + (boxfilter_fast(b, r) ./ N);
end

function imDst = boxfilter_fast(imSrc, r)
if exist('imboxfilt', 'file') == 2
    imDst = imboxfilt(imSrc, 2*r+1, 'Padding', 'replicate');
else
    [hei, wid] = size(imSrc);
    imDst = zeros(size(imSrc), 'like', imSrc);
    imCum = cumsum(imSrc, 1);
    imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
    imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
    imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);
    imCum = cumsum(imDst, 2);
    imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
    imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
    imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
end
end

function alpha = close_matte_w(I, consts_map, consts_vals)
I = double(I) / 255;
alpha = solveAlpha_w(I, consts_map, consts_vals);
end

function alpha = solveAlpha_w(I, consts_map, consts_vals)
[h, w, ~] = size(I);
img_size = w * h;

A = getLaplacian1(I, consts_map);
D = spdiags(consts_map(:), 0, img_size, img_size);

lambda = 100;

b = lambda * consts_map(:) .* consts_vals(:);
x = (A + lambda * D) \ b;

alpha = max(min(reshape(x, h, w), 1), 0);
end

function A = getLaplacian1(I, consts)
epsilon = 1e-7;
win_size = 1;
neb_size = (win_size * 2 + 1)^2;

[h, w, c] = size(I);
img_size = w * h;

consts = imerode(consts, ones(win_size * 2 + 1));
indsM = reshape(1:img_size, h, w);

row_inds = [];
col_inds = [];
vals = [];

for j = 1 + win_size : w - win_size
    for i = 1 + win_size : h - win_size
        if consts(i, j), continue; end

        win_inds = indsM(i - win_size : i + win_size, j - win_size : j + win_size);
        win_inds = win_inds(:);

        winI = I(i - win_size : i + win_size, j - win_size : j + win_size, :);
        winI = reshape(winI, neb_size, c);

        win_mu = mean(winI, 1)';

        Sigma = (winI' * winI) / neb_size - win_mu * win_mu' + (epsilon / neb_size) * eye(c);
        win_var = Sigma \ eye(c);

        winI = winI - repmat(win_mu', neb_size, 1);
        tvals = (1 + winI * win_var * winI') / neb_size;

        row_inds = [row_inds; reshape(repmat(win_inds, 1, neb_size), neb_size^2, 1)];
        col_inds = [col_inds; reshape(repmat(win_inds', neb_size, 1), neb_size^2, 1)];
        vals = [vals; tvals(:)];
    end
end

A = sparse(row_inds, col_inds, vals, img_size, img_size);
A = spdiags(sum(A, 2), 0, img_size, img_size) - A;
end
