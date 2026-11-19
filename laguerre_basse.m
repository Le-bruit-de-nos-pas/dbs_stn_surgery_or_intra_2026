folder = 'C:\Users\paulo\Desktop\laguerre_basse';

amplitudePercentile = 99.99;

sides = {'Left', 'Right'};
prefixes = {'lt1d', 'rt1d'};

if ~isfolder(folder)
    error('Recording folder does not exist: %s', folder);
end

for s = 1:numel(sides)

    % Channel numbers and electrode positions differ between sides.
    if strcmp(sides{s}, 'Left')
        positions = {'Central', 'Medial', 'Anterior', 'Lateral'};
        channels = {'CSPK_01___Central', ...
                    'CSPK_02___Medial', ...
                    'CSPK_03___Anterior', ...
                    'CSPK_04___Lateral'};
    else
        positions = {'Central', 'Lateral', 'Posterior'};
        channels = {'CSPK_01___Central', ...
                    'CSPK_02___Lateral', ...
                    'CSPK_03___Posterior'};
    end

    nChannels = numel(channels);

    % Select f001 recordings and extract signed distance from target.
    files = dir(fullfile(folder, [prefixes{s} '*f001.mat']));
    depths = nan(numel(files), 1);

    for j = 1:numel(files)
        token = regexpi(files(j).name, ...
            '^[lr]t1d(-?\d+(?:\.\d+)?)f001\.mat$', 'tokens', 'once');

        if ~isempty(token)
            depths(j) = str2double(token{1});
        end
    end

    valid = isfinite(depths);
    files = files(valid);
    [depths, order] = sort(depths(valid), 'descend');
    files = files(order);

    if isempty(files)
        warning('No matching %s f001 recordings found in %s.', sides{s}, folder);
        continue
    end

    fprintf('\n%s: %d recordings, %d electrode positions\n', ...
        sides{s}, numel(files), nChannels);

    % Create one panel per electrode position.
    fig = figure('Color', 'w', ...
        'Name', [sides{s} ' SPK by depth'], ...
        'Position', [50 80 1500 800]);

    layout = tiledlayout(fig, nChannels, 1, ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    title(layout, [sides{s} ' SPK recordings by distance from target']);
    xlabel(layout, 'Cumulative recorded time (s)');
    ylabel(layout, 'Amplitude (exported units)');

    ax = gobjects(nChannels, 1);

    for k = 1:nChannels
        ax(k) = nexttile(layout);
        hold(ax(k), 'on');
        ylabel(ax(k), positions{k});
        box(ax(k), 'off');
    end

    edges = zeros(numel(files) + 1, 1);
    amplitudes = cell(numel(files), nChannels);

    % Plot each recording separately so lines do not bridge depth changes.
    for j = 1:numel(files)
        D = load(fullfile(files(j).folder, files(j).name), 'CSPK_*');

        available = find(isfield(D, channels));
        available = available(~cellfun(@(ch) isempty(D.(ch)), channels(available)));

        if isempty(available)
            error('No nonempty expected SPK channels in %s.', files(j).name);
        end

        starts = nan(1, nChannels);
        for k = available(:)'
            starts(k) = D.([channels{k} '_TimeBegin']);
        end

        firstStart = min(starts(available));
        duration = 0;

        for k = available(:)'
            ch = channels{k};
            y = double(D.(ch));
            fs = 1000 * D.([ch '_KHz']);

            offset = starts(k) - firstStart;
            t = edges(j) + offset + (0:numel(y)-1) / fs;

            plot(ax(k), t, y, ...
                'Color', [0.12 0.32 0.46], 'LineWidth', 0.4);

            amplitudes{j, k} = abs(y(:));
            duration = max(duration, offset + numel(y) / fs);
        end

        edges(j + 1) = edges(j) + duration;
    end

    % Set an independent y-limit for each electrode across all its depths.
    for k = 1:nChannels
        channelAmplitudes = vertcat(amplitudes{:, k});

        if isempty(channelAmplitudes)
            title(ax(k), 'No SPK samples');
            fprintf('%s %-9s | no samples\n', sides{s}, positions{k});
            continue
        end

        pLimit = prctile(channelAmplitudes, amplitudePercentile);
        peak = max(channelAmplitudes);
        pctZero = 100 * mean(channelAmplitudes == 0);

        if peak == 0
            title(ax(k), 'All recorded SPK samples are zero');
            yticks(ax(k), []);
            yLimit = 1;
        elseif pLimit <= 1 && peak > pLimit
            % Preserve sparse spikes when the percentile is near zero.
            yLimit = peak;
        else
            yLimit = pLimit;
        end

        ylim(ax(k), [-yLimit, yLimit]);

        fprintf(['%s %-9s | %.2f%% zeros | P%.1f = %.2f | ' ...
                 'peak = %.0f | displayed axis +/-%.2f\n'], ...
            sides{s}, positions{k}, pctZero, ...
            amplitudePercentile, pLimit, peak, yLimit);

        % Shade alternating depths. Peach marks recordings past target.
        for j = 1:numel(depths)
            if depths(j) < 0
                fillColor = [1.00 0.90 0.84];
            elseif mod(j, 2) == 1
                fillColor = [0.92 0.95 0.97];
            else
                continue
            end

            band = patch(ax(k), ...
                edges([j j+1 j+1 j]), ...
                [-yLimit -yLimit yLimit yLimit], ...
                fillColor, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');

            uistack(band, 'bottom');
        end
    end

    clear amplitudes channelAmplitudes

    % Synchronize time only. Each electrode retains its own y-axis.
    linkaxes(ax, 'x');
    xlim(ax(1), [0 edges(end)]);

    centers = (edges(1:end-1) + edges(2:end)) / 2;
    ax(1).XAxisLocation = 'top';
    xticks(ax(1), centers);
    xticklabels(ax(1), compose('%.3f', depths));
    xtickangle(ax(1), 90);
    xlabel(ax(1), 'Distance from target (mm)');

    % Keep time labels visible on the last panel on each side.
    for k = 2:nChannels-1
        ax(k).XTickLabel = [];
    end

    drawnow;

        % Export the complete figure for this side.
    exportFolder = fullfile(folder, 'SPK_figures');
    if ~isfolder(exportFolder)
        mkdir(exportFolder);
    end

    baseName = fullfile(exportFolder, [sides{s} '_SPK_by_depth']);

    % High-resolution image for slides or documents.
    exportgraphics(layout, [baseName '.png'], ...
        'Resolution', 300, 'BackgroundColor', 'white');

    % PDF with the full layout. Rasterize the dense waveform for file size.
    exportgraphics(layout, [baseName '.pdf'], ...
        'ContentType', 'image', 'Resolution', 300, ...
        'BackgroundColor', 'white');

    % Editable MATLAB figure: reopen to zoom or change axes and labels.
    savefig(fig, [baseName '.fig']);

    fprintf('Exported %s to %s\n', sides{s}, exportFolder);
end