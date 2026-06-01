function topo_scatter(ch_names, values, clim, ax, ttl)
% TOPO_SCATTER  Simple scalp map at standard 10-20 positions.
%   Draws a head outline (circle + nose + ears), places colour-coded markers
%   at each found electrode, and fills the scalp with a scatteredInterpolant
%   background (if >= 4 channels are located in the lookup table).
%   Falls back to a bar chart when no channel names match.
%
%   ch_names   cell array of strings
%   values     numeric vector, same length as ch_names
%   clim       [cmin cmax] (pass [] for auto)
%   ax         axes handle (pass [] for current axes)
%   ttl        title string (optional)

    if nargin < 4 || isempty(ax),  ax  = gca; end
    if nargin < 5,                  ttl = '';  end

    % ── Standard 10-20 positions (x = right, y = front, unit circle) ─────────
    raw = {
        'Fp1',[-0.31, 0.87]; 'Fp2',[ 0.31, 0.87];
        'AF3',[-0.54, 0.67]; 'AF4',[ 0.54, 0.67];
        'F7', [-0.71, 0.49]; 'F3', [-0.39, 0.50];
        'Fz', [ 0.00, 0.50]; 'F4', [ 0.39, 0.50]; 'F8',[ 0.71, 0.49];
        'FC5',[-0.54, 0.30]; 'FC1',[-0.19, 0.30];
        'FCz',[ 0.00, 0.30]; 'FC2',[ 0.19, 0.30]; 'FC6',[ 0.54, 0.30];
        'T7', [-0.87, 0.00]; 'C3', [-0.44, 0.00];
        'Cz', [ 0.00, 0.00]; 'C4', [ 0.44, 0.00]; 'T8',[ 0.87, 0.00];
        'T3', [-0.87, 0.00]; 'T4', [ 0.87, 0.00];   % aliases
        'CP5',[-0.54,-0.30]; 'CP1',[-0.19,-0.30];
        'CPz',[ 0.00,-0.30]; 'CP2',[ 0.19,-0.30]; 'CP6',[ 0.54,-0.30];
        'TP7',[-0.54,-0.30]; 'TP8',[ 0.54,-0.30];
        'P7', [-0.71,-0.49]; 'P3', [-0.39,-0.50];
        'Pz', [ 0.00,-0.50]; 'P4', [ 0.39,-0.50]; 'P8',[ 0.71,-0.49];
        'T5', [-0.71,-0.49]; 'T6', [ 0.71,-0.49];   % aliases
        'PO3',[-0.35,-0.72]; 'PO4',[ 0.35,-0.72];
        'O1', [-0.31,-0.87]; 'Oz', [ 0.00,-0.87]; 'O2',[ 0.31,-0.87];
    };
    lookup = containers.Map(upper(raw(:,1)), raw(:,2));

    % Resolve positions (case-insensitive)
    n = numel(ch_names);
    xy = NaN(n, 2);
    for i = 1:n
        k = upper(strtrim(ch_names{i}));
        if lookup.isKey(k), xy(i,:) = lookup(k); end
    end
    found = ~any(isnan(xy), 2);

    if ~any(found)
        bar(ax, values);
        set(ax, 'XTickLabel', ch_names, 'XTickLabelRotation', 45);
        title(ax, ttl, 'Interpreter', 'none', 'FontSize', 8);
        return;
    end

    xf = xy(found,1);   yf = xy(found,2);   vf = values(found(:));
    cla(ax); hold(ax,'on'); axis(ax,'equal','off');

    % ── Background interpolation ─────────────────────────────────────────────
    [GX, GY] = meshgrid(linspace(-1.05,1.05,100));
    mask = GX.^2 + GY.^2 <= 1.0^2;
    if sum(found) >= 4
        try
            F  = scatteredInterpolant(xf, yf, vf(:), 'natural','none');
            GZ = F(GX, GY);
            GZ(~mask) = NaN;
            pcolor(ax, GX, GY, GZ);
            shading(ax,'interp');
            colormap(ax, 'jet');
        catch
        end
    end
    if ~isempty(clim) && numel(clim)==2 && diff(clim)>0
        caxis(ax, clim);
    end

    % ── Head outline ─────────────────────────────────────────────────────────
    th = linspace(0, 2*pi, 300);
    plot(ax, cos(th), sin(th), 'k-', 'LineWidth',1.5);
    % Nose
    plot(ax, [-0.08, 0, 0.08], [1.00, 1.13, 1.00], 'k-','LineWidth',1.5);
    % Ears
    ear_th = linspace(pi/2-0.4, pi/2+0.4, 25);
    ear_r  = 1.08;
    plot(ax, -ear_r*sin(ear_th), ear_r*cos(ear_th)-0.05, 'k-','LineWidth',1.5);
    plot(ax,  ear_r*sin(ear_th), ear_r*cos(ear_th)-0.05, 'k-','LineWidth',1.5);

    % ── Channel markers + labels ──────────────────────────────────────────────
    scatter(ax, xf, yf, 80, vf(:), 'filled', 'MarkerEdgeColor','k','LineWidth',0.5);
    found_idx = find(found);
    for i = 1:numel(found_idx)
        text(ax, xf(i), yf(i)+0.10, strtrim(ch_names{found_idx(i)}), ...
             'FontSize',5.5, 'HorizontalAlignment','center', 'Color','k');
    end

    set(ax,'XLim',[-1.3,1.3],'YLim',[-1.3,1.3]);
    colorbar(ax,'FontSize',7);
    if ~isempty(ttl)
        title(ax, ttl, 'FontSize',8, 'Interpreter','none');
    end
end
