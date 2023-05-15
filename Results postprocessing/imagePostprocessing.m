% File, folder and variable definition

% Define the output folders
output_bmp_folder = fullfile('output_bmp');
output_mat_folder = fullfile('output_mat');
input_folder = 'input_img';

% Create output folders if they don't exist
if ~exist(output_bmp_folder, 'dir')
    mkdir(output_bmp_folder);
end

if ~exist(output_mat_folder, 'dir')
    mkdir(output_mat_folder);
end


% Initialize the file names
dmd_files = {'LG10_dmd.tif', 'LG11_dmd.tif', 'LG20_dmd.tif','LG21_dmd.tif', ...
             'HG10_dmd.tif', 'HG11_dmd.tif', 'HG20_dmd.tif', 'HG21_dmd.tif'};

film_files = {'LG10_film.tif', 'LG11_film.tif', 'LG20_film.tif', 'LG21_film.tif', ...
              'HG10_film.tif', 'HG11_film.tif', 'HG20_film.tif', 'HG21_film.tif'};

% Combine the two arrays into one
all_files = [dmd_files, film_files];

% Define the resize factor and Gaussian filter variables
resize_factor = 0.85; % Change this to modify the resize factor for _dmd images
gaussian_filter_sigma = 6; % Change this to modify the Gaussian filter effect

% Start of Image Post-processing

% Iterate over all_files and process each image
for i = 1:length(all_files)
    % Read the image
    img_path = fullfile(input_folder, all_files{i});
    img = imread(img_path);
    
    % Convert images to double
    if ~isa(img, 'double')
        img = im2double(img);
    end

    % Select the first plane 
    if ndims(img) > 2
        img = img(:, :, 1);
    end
    
    % Crop the image to 1:1 aspect ratio
    [height, width] = size(img);
    crop_size = min(height, width);
    x_start = floor((width - crop_size) / 2);
    y_start = floor((height - crop_size) / 2);
    cropped_img = imcrop(img, [x_start, y_start, crop_size-1, crop_size-1]);
   
    % Only for images from the DMD
    if ismember(all_files{i}, dmd_files)
        % Resize the image
        resized_img = imresize(cropped_img, resize_factor);
        
        % Get the dimensions of the original and resized images
        [original_height, original_width] = size(cropped_img);
        [resized_height, resized_width] = size(resized_img);
        
        % Calculate the padding dimensions
        pad_height = original_height - resized_height;
        pad_width = original_width - resized_width;
        
        % Calculate symmetric padding values for each side
        top_pad = floor(pad_height / 2);
        bottom_pad = pad_height - top_pad;
        left_pad = floor(pad_width / 2);
        right_pad = pad_width - left_pad;
        
        % Add zero-padding symmetrically to match the original scale
        padded_img = padarray(resized_img, [top_pad, left_pad], 'pre');
        padded_img = padarray(padded_img, [bottom_pad, right_pad], 'post');
        
        % Return in new variable
        scaled_img = padded_img;
    else
        % Return in new variable
        scaled_img = cropped_img;
    end

    % Continue Post-processing (Gaussian, save .bmp, save .mat, display)

    % Apply the Fourier Gaussian filter
    filtered_img = imgaussfilt(scaled_img, gaussian_filter_sigma);

    % Save the image as a .bmp file in the output_bmp_folder
    [filepath, name, ext] = fileparts(all_files{i});
    output_bmp_path = fullfile(output_bmp_folder, [name, '.bmp']);
    imwrite(filtered_img, output_bmp_path);
    
    % Save the image data as a .mat file in the output_mat_folder
    output_mat_path = fullfile(output_mat_folder, [name, '.mat']);
    save(output_mat_path, 'filtered_img');

    % Display the images before and after postprocessing
    subplot_index = mod(i-1, 4) + 1;
    
    % If we're starting a new batch of 4, create a new figure
    if subplot_index == 1
        figure;
    end
    
    % Display the original cropped image without scaling up
    subplot(2, 4, subplot_index);
    imshow(cropped_img);
    title(sprintf('Original Image %d', i));
    
    % Display the filtered image without scaling up
    subplot(2, 4, subplot_index + 4);
    imshow(filtered_img);
    title(sprintf('Filtered Image %d', i));
    
end

close all;
