folder = 'C:\Users\paulo\Desktop\veyries_vincent';
amplitudePercentile = 99.9;
exportFolder = fullfile(folder, 'SPK_figures');

if ~isfolder(folder)
    error('Recording folder does not exist: %s', folder);
end

if ~isfolder(exportFolder)
    mkdir(exportFolder);
end

sides = {'Left', 'Right'};
sideLetters = {'L', 'R'};

listing = dir(folder);
listing = listing(~[listing.isdir]);

for s = 1:numel(sides)

    % Match LT3D..., RT1D..., or another track number. Use F001 only.
    pattern = ['^' sideLetters{s} 'T(\d+)D(-?\d+(?:\.\d+)?)F001\.mat$'];
    tokens = regexpi({listing.name}, pattern, 'tokens', 'once');
    matched = ~cellfun(@isempty, tokens);

    files = listing(matched);
    tokens = tokens(matched);

    if isempty(files)
        warning('No %s F001 recordings found in %s.', sides{s}, folder);
        continue
    end

    tracks = cellfun(@(x) str2double(x{1}), tokens);
    depths = cellfun(@(x) str2double(x{2}), tokens);

    % Keep separate trajectories in separate figures.
    for track = unique(tracks(:))'

        onTrack = tracks == track;
        trackFiles = files(onTrack);
        trackDepths = depths(onTrack);

        [trackDepths, order] = sort(trackDepths(:), 'descend');
        trackFiles = trackFiles(order);
        nFiles = numel(trackFiles);

        % Discover continuous SPK waveforms across all depths.
        % Timing metadata and exact MATLAB variable names are kept per file.
        waveforms = cell(nFiles, 1);
        rateVars = cell(nFiles, 1);
        startVars = cell(nFiles, 1);

        for j = 1:nFiles
            filePath = fullfile(trackFiles(j).folder, trackFiles(j).name);
            vars = whos('-file', filePath);
            names = {vars.name};

            waveforms{j} = {};
            rateVars{j} = {};
            startVars{j} = {};

            for v = 1:numel(vars)
                name = vars(v).name;
                hasSPKName = ~isempty(regexpi(name, ...
                    '^(CSPK|SPK|CSpike|Spike)', 'once'));
                isVector = numel(vars(v).size) == 2 && any(vars(v).size == 1);
                hasSamples = prod(vars(v).size) > 1;

                if ~hasSPKName || ~isVector || ~hasSamples
                    continue
                end

                rateIdx = find(strcmpi(names, [name '_KHz']), 1);
                startIdx = find(strcmpi(names, [name '_TimeBegin']), 1);

                if isempty(rateIdx) || isempty(startIdx)
                    warning('Skipping %s in %s: missing KHz or TimeBegin.', ...
                        name, trackFiles(j).name);
                    continue
                end

                waveforms{j}{end+1} = name;
                rateVars{j}{end+1} = names{rateIdx};
                startVars{j}{end+1} = names{startIdx};
            end
        end

        channels = unique([waveforms{:}]);
        nChannels = numel(channels);

        if nChannels == 0
            warning('No continuous SPK waveforms identified for %s track %d.', ...
                sides{s}, track);
            fprintf('Variables in %s:\n', trackFiles(1).name);
            whos('-file', fullfile(trackFiles(1).folder, trackFiles(1).name))
            continue
        end

        % Use electrode positions as labels when present in the variable name.
        positions = regexprep(channels, ...
            '(?i)^(?:CSPK|SPK|CSpike|Spike)[_ ]*\d+[_ ]+(?=[A-Za-z])', '');
        positions = strrep(positions, '_', ' ');

        fprintf('\n%s | track %d | %d F001 recordings\n', sides{s}, track, nFiles);
        disp(table(channels(:), positions(:), ...
            'VariableNames', {'MATLAB_variable', 'Panel_label'}));

        fig = figure('Color', 'w', ...
            'Name', sprintf('%s track %d SPK by depth', sides{s}, track), ...
            'Position', [50 80 1500 800]);

        layout = tiledlayout(fig, nChannels, 1, ...
            'TileSpacing', 'compact', 'Padding', 'loose');

        title(layout, sprintf('%s track %d | SPK recordings by distance from target', ...
            sides{s}, track));
        xlabel(layout, 'Cumulative recorded time (s)');
        ylabel(layout, 'Amplitude (exported units)');

        ax = gobjects(nChannels, 1);
        for k = 1:nChannels
            ax(k) = nexttile(layout);
            hold(ax(k), 'on');
            ylabel(ax(k), positions{k}, 'Interpreter', 'none');
            box(ax(k), 'off');
        end

        edges = zeros(nFiles + 1, 1);
        amplitudes = cell(nFiles, nChannels);

        % Plot each depth at its full duration.
        for j = 1:nFiles
            edges(j+1) = edges(j);
            if isempty(waveforms{j})
                fprintf('Skipping %s: no SPK samples.\n', trackFiles(j).name);
                continue
            end

            filePath = fullfile(trackFiles(j).folder, trackFiles(j).name);
            needed = [waveforms{j}, rateVars{j}, startVars{j}];
            D = load(filePath, needed{:});

            starts = cellfun(@(name) double(D.(name)), startVars{j});
            ratesHz = cellfun(@(name) 1000 * double(D.(name)), rateVars{j});
            valid = isfinite(starts) & isfinite(ratesHz) & ratesHz > 0;

            if any(~valid)
                warning('Invalid timing metadata in %s for: %s', ...
                    trackFiles(j).name, strjoin(waveforms{j}(~valid), ', '));
            end
            if ~any(valid)
                continue
            end

            firstStart = min(starts(valid));
            duration = 0;

            for v = find(valid(:))'
                name = waveforms{j}{v};
                k = find(strcmp(channels, name), 1);

                y = double(D.(name));
                y = y(:);
                fs = ratesHz(v);
                offset = starts(v) - firstStart;
                t = edges(j) + offset + (0:numel(y)-1)' / fs;

                plot(ax(k), t, y, ...
                    'Color', [0.12 0.32 0.46], 'LineWidth', 0.4);

                amplitudes{j, k} = abs(y(isfinite(y)));
                duration = max(duration, offset + numel(y) / fs);
            end

            edges(j+1) = edges(j) + duration;
        end

        recorded = diff(edges) > 0;
        if ~any(recorded)
            warning('No plottable recordings for %s track %d.', sides{s}, track);
            close(fig);
            continue
        end

        % One y-limit per electrode, calculated across all recorded depths.
        for k = 1:nChannels
            values = vertcat(amplitudes{:, k});

            if isempty(values)
                title(ax(k), 'No finite SPK samples');
                continue
            end

            pLimit = prctile(values, amplitudePercentile);
            peak = max(values);

            if peak == 0
                yLimit = 1;
                title(ax(k), 'All recorded samples are zero');
                yticks(ax(k), []);
            elseif pLimit <= 1 && peak > pLimit
                yLimit = peak;
            else
                yLimit = pLimit;
            end

            ylim(ax(k), [-yLimit, yLimit]);
            fprintf('%s %-12s | P%.1f = %.2f | peak = %.0f | axis +/-%.2f\n', ...
                sides{s}, positions{k}, amplitudePercentile, pLimit, peak, yLimit);

            for j = find(recorded(:))'
                if trackDepths(j) < 0
                    fillColor = [1.00 0.90 0.84];
                elseif mod(j, 2) == 1
                    fillColor = [0.92 0.95 0.97];
                else
                    continue
                end

                band = patch(ax(k), edges([j j+1 j+1 j]), ...
                    [-yLimit -yLimit yLimit yLimit], fillColor, ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
                uistack(band, 'bottom');
            end
        end

        clear amplitudes values D

        linkaxes(ax, 'x');
        xlim(ax(1), [0 edges(end)]);
        centers = (edges(1:end-1) + edges(2:end)) / 2;

        ax(1).XAxisLocation = 'top';
        xticks(ax(1), centers(recorded));
        xticklabels(ax(1), compose('%.3f', trackDepths(recorded)));
        xtickangle(ax(1), 90);
        xlabel(ax(1), 'Distance from target (mm)');

        for k = 2:nChannels-1
            ax(k).XTickLabel = [];
        end

        % A single-channel figure needs a separate bottom axis for time.
        if nChannels == 1
            timeAx = axes(layout, 'Color', 'none', 'YAxisLocation', 'right', ...
                'YTick', [], 'Box', 'off');
            timeAx.Layout.Tile = 1;
            timeAx.YColor = 'none';
            xlim(timeAx, [0 edges(end)]);
            linkaxes([ax(1), timeAx], 'x');
        end

        drawnow;

        baseName = fullfile(exportFolder, ...
            sprintf('%s_Track%d_SPK_by_depth', sides{s}, track));

        exportgraphics(layout, [baseName '.png'], ...
            'Resolution', 300, 'BackgroundColor', 'white');
        exportgraphics(layout, [baseName '.pdf'], ...
            'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
        savefig(fig, [baseName '.fig']);

        fprintf('Exported %s track %d to %s\n', sides{s}, track, exportFolder);
    end
end