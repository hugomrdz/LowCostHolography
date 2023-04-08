% Ensure the output_figures folder exists
output_fig_folder = 'output_figures';
if ~exist(output_fig_folder, 'dir')
    mkdir(output_fig_folder);
end

% File arrays
dmd_files = {'LG10_dmd.mat', 'LG11_dmd.mat', 'LG20_dmd.mat','LG21_dmd.mat', ...
             'HG10_dmd.mat', 'HG11_dmd.mat', 'HG20_dmd.mat', 'HG21_dmd.mat'};

film_files = {'LG10_film.mat', 'LG11_film.mat', 'LG20_film.mat', 'LG21_film.mat', ...
              'HG10_film.mat', 'HG11_film.mat', 'HG20_film.mat', 'HG21_film.mat'};

num_files = {'LG10_num.mat', 'LG11_num.mat', 'LG20_num.mat', 'LG21_num.mat', ...
              'HG10_num.mat', 'HG11_num.mat', 'HG20_num.mat', 'HG21_num.mat'};

n_files = length(dmd_files);

for i = 1:n_files

    % Load images from .mat files
    dmd_mat = load(fullfile('output_mat', dmd_files{i}));
    film_mat = load(fullfile('output_mat', film_files{i}));
    num_mat = load(fullfile('output_mat', num_files{i}));

    dmd_data = dmd_mat.filtered_img;  
    film_data = film_mat.filtered_img;
    num_data = num_mat.numericData;

    % Calculate the center row for each image
    dmd_center_row = floor(size(dmd_data, 1) / 2);
    film_center_row = floor(size(film_data, 1) / 2);
    num_center_row = floor(size(num_data, 1) / 2);
    
    % Extract the center row from each image
    dmd_center_line = dmd_data(dmd_center_row, :);
    film_center_line = film_data(film_center_row, :);
    num_center_line = num_data(num_center_row, :);

    % Normalize between 0 and 1
    dmd_center_line = (dmd_center_line - min(dmd_center_line)) / (max(dmd_center_line) - min(dmd_center_line));
    film_center_line = (film_center_line - min(film_center_line)) / (max(film_center_line) - min(film_center_line));
    num_center_line = (num_center_line - min(num_center_line)) / (max(num_center_line) - min(num_center_line));
    
    % Interpolate the shorter centerline to match the length of the longer one
    if length(dmd_center_line) > length(film_center_line)
        film_center_line = interp1(1:length(film_center_line), film_center_line, linspace(1, length(film_center_line), length(dmd_center_line)));
    elseif length(dmd_center_line) < length(film_center_line)
        dmd_center_line = interp1(1:length(dmd_center_line), dmd_center_line, linspace(1, length(dmd_center_line), length(film_center_line)));
    end
    
    % Align the DMD and Film centerlines
    [aligned_dmd, aligned_film] = alignsignals(dmd_center_line, film_center_line);
    
    % Align the DMD and Numerical centerlines
    [aligned_dmd, aligned_num] = alignsignals(aligned_dmd, num_center_line);
    
    % Pad the shorter signals with zeros, if necessary
    if length(aligned_film) < length(aligned_dmd)
        aligned_film = [aligned_film, zeros(1, length(aligned_dmd) - length(aligned_film))];
    elseif length(aligned_film) > length(aligned_dmd)
        aligned_dmd = [aligned_dmd, zeros(1, length(aligned_film) - length(aligned_dmd))];
    end
    
    if length(aligned_num) < length(aligned_dmd)
        aligned_num = [aligned_num, zeros(1, length(aligned_dmd) - length(aligned_num))];
    elseif length(aligned_num) > length(aligned_dmd)
        aligned_dmd = [aligned_dmd, zeros(1, length(aligned_num) - length(aligned_dmd))];
        aligned_film = [aligned_film, zeros(1, length(aligned_num) - length(aligned_film))];
    end
    
    % Calculate the mean squared error
    mse = mean((aligned_dmd - aligned_film).^2);
    
    % First figure with legends
    figure;
    plot(aligned_dmd, 'DisplayName', 'DMD');
    hold on;
    plot(aligned_film, 'DisplayName', 'Film');
    plot(aligned_num, 'k--', 'DisplayName', 'Numerical (reference)');
    title(sprintf('Center Line Plots for %s (MSE: %.4f)', dmd_files{i}(1:end-8), mse), 'Interpreter', 'latex');

    % Calculate x-axis positions in meters
    pixel_size = 4.65e-6; % Convert the pixel size to meters (4.65µm x 4.65µm)
    x_positions_meters = (0:length(aligned_dmd)-1) * pixel_size;
    
    % Update x-axis labels
    xlabel('Distance (m)', 'Interpreter', 'latex');
    xtick_positions = round(linspace(1, length(aligned_dmd), 5)); % Calculate and round tick positions
    xticks(xtick_positions); % Set tick positions
    
    xticklabels(arrayfun(@(x) create_tick_label(x), x_positions_meters(xtick_positions), 'UniformOutput', false));
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Update y-axis labels
    ylabel('Normalized Intensity', 'Interpreter', 'latex');
    ytick_positions = linspace(0, 1, 5); % Calculate tick positions
    yticks(ytick_positions); % Set tick positions
    yticklabels(arrayfun(@(y) sprintf('%.2f', y), ytick_positions, 'UniformOutput', false));
    set(gca, 'TickLabelInterpreter', 'latex');

    lgd = legend('DMD', 'Film', 'Numerical', 'Location', 'northeast');
    set(lgd, 'Interpreter', 'latex');

    xlim([1, length(aligned_dmd)]); % Adjust x-axis limits
    ylim([0, 1]); % Adjust y-axis limits
    box on;
    %pbaspect([1 1 1]); % Set aspect ratio manually
    hold off;

    filename_jpeg = sprintf('%s_centerlines.jpg', dmd_files{i}(1:end-8));
    filepath_jpeg = fullfile(output_fig_folder, filename_jpeg); 
    saveas(gcf, filepath_jpeg);

    % Second figure without legends
    figure;
    plot(aligned_dmd);
    hold on;
    plot(aligned_film);
    plot(aligned_num, 'k--');
    title('');
    x_ticks = linspace(1, max([length(aligned_dmd), length(aligned_film), length(aligned_num)]), 5); % Generate 5 equally spaced x-tick positions
    y_ticks = linspace(0, 1, 5); % Generate 5 equally spaced y-tick positions
    set(gca, 'XTick', x_ticks, 'YTick', y_ticks, 'XTickLabel', [], 'YTickLabel', []); % Set ticks without labels
    xlim([1, max([length(aligned_dmd), length(aligned_film), length(aligned_num)])]); % Adjust x-axis limits
    ylim([0, 1]); % Adjust y-axis limits
    pbaspect([1 1 1]); % Set aspect ratio manually
    box on;
    hold off;

    filename_svg = sprintf('%s_centerlines.svg', dmd_files{i}(1:end-8));
    filepath_svg = fullfile(output_fig_folder, filename_svg); 
    saveas(gcf, filepath_svg);
    
end

close all;

% Custom function to handle zero value
function label = create_tick_label(x)
    if x == 0
        label = '0';
    else
        label = sprintf('$%.2f\\times10^{%d}$', x / 10^floor(log10(x)), floor(log10(x)));
    end
end

