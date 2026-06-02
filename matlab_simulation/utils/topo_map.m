function topo_map(ch_names, values, clim, ax, ttl, show_cbar, bg_zero)
% TOPO_MAP  Publication-quality scalp topography at standard 10-20 positions.
%   scatteredInterpolant on 200×200 grid, RdBu diverging colormap, contour at 0,
%   head outline with nose and ears. Each call uses its own clim (pass [] for auto).
%
%   ch_names   cell array of strings
%   values     numeric vector (same length as ch_names)
%   clim       [lo hi] — pass [] for auto (symmetric about 0 when data spans 0)
%   ax         axes handle — pass [] for gca
%   ttl        title string (optional)
%   show_cbar  logical, default true
%   bg_zero    logical (default false): if true, all known 10-20 channels not
%              in ch_names are added with value 0, giving smooth full-head maps
%              when only a subset of electrodes was selected.

    if nargin < 3, clim = []; end
    if nargin < 4 || isempty(ax), ax = gca; end
    if nargin < 5, ttl = ''; end
    if nargin < 6, show_cbar = true; end
    if nargin < 7, bg_zero = false; end

    pos = {
        'Fp1',[-0.31, 0.87]; 'Fp2',[ 0.31, 0.87];
        'AF3',[-0.54, 0.67]; 'AF4',[ 0.54, 0.67];
        'F7', [-0.71, 0.49]; 'F3', [-0.39, 0.50];
        'Fz', [ 0.00, 0.50]; 'F4', [ 0.39, 0.50]; 'F8',[ 0.71, 0.49];
        'FC5',[-0.54, 0.30]; 'FC1',[-0.19, 0.30];
        'FCz',[ 0.00, 0.30]; 'FC2',[ 0.19, 0.30]; 'FC6',[ 0.54, 0.30];
        'T7', [-0.87, 0.00]; 'C3', [-0.44, 0.00];
        'Cz', [ 0.00, 0.00]; 'C4', [ 0.44, 0.00]; 'T8',[ 0.87, 0.00];
        'T3', [-0.87, 0.00]; 'T4', [ 0.87, 0.00];
        'CP5',[-0.54,-0.30]; 'CP1',[-0.19,-0.30];
        'CPz',[ 0.00,-0.30]; 'CP2',[ 0.19,-0.30]; 'CP6',[ 0.54,-0.30];
        'TP7',[-0.54,-0.30]; 'TP8',[ 0.54,-0.30];
        'P7', [-0.71,-0.49]; 'P3', [-0.39,-0.50];
        'Pz', [ 0.00,-0.50]; 'P4', [ 0.39,-0.50]; 'P8',[ 0.71,-0.49];
        'T5', [-0.71,-0.49]; 'T6', [ 0.71,-0.49];
        'PO3',[-0.35,-0.72]; 'PO4',[ 0.35,-0.72];
        'O1', [-0.31,-0.87]; 'Oz', [ 0.00,-0.87]; 'O2',[ 0.31,-0.87];
    };
    lut = containers.Map(upper(pos(:,1)), pos(:,2));

    % Remember how many channels are "real" before zero-padding
    n_orig = numel(ch_names);

    % Fill all known electrodes with 0 so interpolation covers the full scalp
    if bg_zero
        known_keys  = lut.keys();
        input_upper = upper(cellfun(@strtrim, ch_names(:)', 'UniformOutput', false));
        for ki = 1:numel(known_keys)
            if ~any(strcmp(input_upper, known_keys{ki}))
                ch_names{end+1} = known_keys{ki}; %#ok<AGROW>
                values(end+1)   = 0;              %#ok<AGROW>
            end
        end
    end

    n  = numel(ch_names);
    xy = NaN(n,2);
    for i = 1:n
        k = upper(strtrim(ch_names{i}));
        if lut.isKey(k), xy(i,:) = lut(k); end
    end
    ok = ~any(isnan(xy),2);

    if ~any(ok)
        bar(ax, values);
        set(ax,'XTickLabel',ch_names,'XTickLabelRotation',45);
        if ~isempty(ttl), title(ax,ttl,'Interpreter','none','FontSize',8); end
        return
    end

    xf = xy(ok,1); yf = xy(ok,2); vf = double(values(ok));
    cla(ax); hold(ax,'on'); axis(ax,'equal','off');

    res = 200;
    [GX,GY] = meshgrid(linspace(-1.1,1.1,res));
    mask_head = GX.^2 + GY.^2 <= 1.0^2;
    GZ = NaN(res);
    if sum(ok) >= 3
        try
            F = scatteredInterpolant(xf(:), yf(:), vf(:), 'natural', 'linear');
            GZ = F(GX, GY);
            GZ(~mask_head) = NaN;
        catch; end
    end

    if isempty(clim)
        gf = GZ(isfinite(GZ));
        if isempty(gf), gf = vf(:); end
        vlo = min(gf); vhi = max(gf);
        if vlo < 0 && vhi > 0
            vm = max(abs([vlo vhi])); clim = [-vm vm];
        elseif vlo >= 0
            clim = [0 max(vhi, eps)];
        else
            clim = [min(vlo,-eps) 0];
        end
    end

    pcolor(ax, GX, GY, GZ); shading(ax,'interp');
    colormap(ax, topo_rdbu(256));
    caxis(ax, clim);
    if show_cbar
        colorbar(ax, 'FontSize',6, 'TickLabelInterpreter','none');
    end

    try; contour(ax, GX, GY, GZ, [0 0], 'k-', 'LineWidth', 0.5); catch; end

    th = linspace(0, 2*pi, 300);
    plot(ax, cos(th), sin(th), 'k-', 'LineWidth', 1.5);
    plot(ax, [-0.09,0,0.09], [0.99,1.13,0.99], 'k-', 'LineWidth', 1.5);
    ear_x = [0.98, 1.04, 1.06, 1.06, 1.04, 0.98];
    ear_y = [0.14, 0.12, 0.04, -0.04, -0.12, -0.14];
    plot(ax, -ear_x, ear_y, 'k-', 'LineWidth', 1.5);
    plot(ax,  ear_x, ear_y, 'k-', 'LineWidth', 1.5);

    scatter(ax, xf, yf, 18, vf, 'filled', 'MarkerEdgeColor','none');

    % Labels only for the original (non-zero-padded) channels
    orig_mask = find(ok) <= n_orig;
    fi_orig   = find(ok); fi_orig = fi_orig(orig_mask);
    for i = 1:numel(fi_orig)
        text(ax, xy(fi_orig(i),1), xy(fi_orig(i),2)+0.11, ...
             strtrim(ch_names{fi_orig(i)}), ...
             'FontSize',5.5,'HorizontalAlignment','center','Color','k', ...
             'FontWeight','bold');
    end
    set(ax, 'XLim',[-1.4,1.4], 'YLim',[-1.4,1.4]);
    if ~isempty(ttl)
        title(ax, ttl, 'FontSize', 7.5, 'Interpreter', 'none');
    end
end


function cm = topo_rdbu(n)
    if nargin < 1, n = 256; end
    h = floor(n/2);
    r = [linspace(0.17,1.00,h), linspace(1.00,0.70,h)];
    g = [linspace(0.51,1.00,h), linspace(1.00,0.09,h)];
    b = [linspace(0.73,1.00,h), linspace(1.00,0.07,h)];
    cm = [r(:), g(:), b(:)];
    if size(cm,1) > n, cm = cm(1:n,:); end
end
