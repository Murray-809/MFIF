function results = alignment()

    clc; clear; close all;
    config.image_dir  = '';
    config.output_dir = '';

    config.levels    = 3;
    config.noi       = 50;
    config.transform = 'affine';

    if ~exist(config.output_dir, 'dir')
        mkdir(config.output_dir);
    end


    image_paths = get_image_files(config.image_dir);
    n_images = numel(image_paths);
    if n_images < 2
        error('Need at least 2 images in config.image_dir');
    end

    ref_idx = n_images;
    ref_img = imread(image_paths{ref_idx});
    ref_img = force_rgb_uint8(ref_img);

    % Save reference as-is (as 00X.tif)
    imwrite(ref_img, fullfile(config.output_dir, sprintf('%03d.tif', ref_idx)));

    % Only keep next_registered in memory
    next_registered = ref_img;

    % Store only metrics (not images)
    registration_results = cell(n_images, 1);
    registration_results{ref_idx} = struct( ...
        'success', true, 'correlation', 1.0, 'improvement', 0.0, 'time', 0.0, 'note', 'Reference');

    for i = (ref_idx - 1):-1:1
        moving = imread(image_paths{i});
        moving = force_rgb_uint8(moving);

        t0 = tic;
        [reg_img, success, correlation, improvement] = ecc_alignment_fullres( ...
            moving, next_registered, config);
        elapsed = toc(t0);

        registration_results{i} = struct( ...
            'success', success, ...
            'correlation', correlation, ...
            'improvement', improvement, ...
            'time', elapsed);

        imwrite(reg_img, fullfile(config.output_dir, sprintf('%03d.tif', i)));


        next_registered = reg_img;

        clear moving reg_img;
    end

    results = struct( ...
        'registration_results', {registration_results}, ...
        'reference_idx', ref_idx, ...
        'image_paths', {image_paths}, ...
        'config', config);
end

function [result_img, success, correlation, improvement] = ecc_alignment_fullres(img1_color, img2_color, config)
    
    img1_gray = rgb2gray(img1_color);
    img2_gray = rgb2gray(img2_color);

    
    img1_d = double(img1_color);
    img2_d = double(img2_color);
    mse_before = mean((img1_d(:) - img2_d(:)).^2);

    try
        
        [ecc_results, warp] = ecc( ...
            img1_gray, img2_gray, ...
            config.levels, config.noi, config.transform);

        correlation = ecc_results(1, end).rho;

        
        [h, w, c] = size(img2_color);
        result_img = zeros(h, w, c, 'uint8');

        for ch = 1:c
            warped_ch = spatial_interp( ...
                double(img1_color(:,:,ch)), ...
                warp, ...
                'linear', ...
                config.transform, ...
                1:w, ...
                1:h);
            result_img(:,:,ch) = uint8(warped_ch);
        end

        
        result_d  = double(result_img);
        mse_after = mean((result_d(:) - img2_d(:)).^2);

        if mse_before > 0
            improvement = (mse_before - mse_after) / mse_before * 100;
        else
            improvement = 0;
        end

        
        if strcmp(config.transform, 'affine') || strcmp(config.transform, 'euclidean')
            translation = hypot(warp(1,3), warp(2,3));
            det_val = det(warp(1:2, 1:2));
            is_reasonable = (translation < 300.0) && (det_val > 0.5) && (det_val < 2.0);
        else
            is_reasonable = true;
        end

        success = (correlation > 0.1) && is_reasonable;

    catch
        result_img  = img1_color;
        success     = false;
        correlation = 0;
        improvement = 0;
    end
end


function img = force_rgb_uint8(img)
    if ~isa(img, 'uint8')
        img = im2uint8(img);
    end
    if size(img,3) == 1
        img = repmat(img, [1 1 3]);
    end
end

function image_paths = get_image_files(image_dir)
    patterns = {'*.bmp','*.jpg','*.jpeg','*.png','*.tiff','*.tif'};
    all_files = {};
    for i = 1:numel(patterns)
        f1 = dir(fullfile(image_dir, patterns{i}));
        for j = 1:numel(f1)
            all_files{end+1} = fullfile(image_dir, f1(j).name); 
        end
        f2 = dir(fullfile(image_dir, upper(patterns{i})));
        for j = 1:numel(f2)
            all_files{end+1} = fullfile(image_dir, f2(j).name);
        end
    end
    all_files = unique(all_files);
    image_paths = natural_sort(all_files);
end

function sorted_list = natural_sort(file_list)
    [~, names, ~] = cellfun(@fileparts, file_list, 'UniformOutput', false);
    nums = zeros(numel(names), 1);
    for i = 1:numel(names)
        tokens = regexp(names{i}, '\d+', 'match');
        if ~isempty(tokens)
            nums(i) = str2double(tokens{1});
        else
            nums(i) = i;
        end
    end
    [~, idx] = sort(nums);
    sorted_list = file_list(idx);
end
