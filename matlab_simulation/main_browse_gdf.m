%% MAIN_BROWSE_GDF  Interactive full-recording GDF browser.
%
%  Three synchronised panels (shared x-axis, time in seconds):
%    1. EEG montage — all channels with vertical offset (µV)
%    2. sLDA classifier outputs — MI (green), CVSA (orange), fused (purple)
%    3. Leaky-integrator buffer output — per class + threshold lines
%
%  Colour-coded vertical lines mark every GDF event.
%  Artifact frames are shaded in light red across all panels.
%
%  Controls
%    Slider         scroll time
%    ← / →          scroll by half-window
%    ↑ / ↓          EEG scale ÷1.4 / ×1.4
%    Win (s)        change window width
%    Scale (µV)     change per-channel spacing

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, ...
        fullfile(this_dir,'io'), fullfile(this_dir,'processing'), ...
        fullfile(this_dir,'artifacts'), fullfile(this_dir,'classifier'), ...
        fullfile(this_dir,'integrator'), fullfile(this_dir,'utils'));

HIT_CODE  = 897;  MISS_CODE = 898;
TO_CODE   = 899;  CF_CODE   = 781;

% ── File picker ───────────────────────────────────────────────────────────
def_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(def_dir), def_dir = fileparts(this_dir); end
[gdf_name, gdf_dir] = uigetfile({'*.gdf','GDF (*.gdf)'},'Select GDF', def_dir);
if isequal(gdf_name,0), return; end
gdf_path   = fullfile(gdf_dir, gdf_name);
[~,basename] = fileparts(gdf_path);
fprintf('[browser] %s\n', basename);

% ── Load data ─────────────────────────────────────────────────────────────
[signal, header, ~] = load_gdf(gdf_path);
[params,  ~]        = load_params_yaml(gdf_path);

paradigm   = params.integrator.paradigm;
fs         = double(params.acquisition.samplerate);
framerate  = double(params.acquisition.framerate);
chunk_size = round(fs / framerate);
if abs(fs - header.SampleRate) > 1e-3
    fs = header.SampleRate;  chunk_size = round(fs / framerate);
end
bufsize_proc = double(params.RingBufferCfg.params.size);
bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

do_car_mi = true;
if isfield(params,'processing_fbcsp_mi')
    do_car_mi = logical(params.processing_fbcsp_mi.do_car);
    nchannels = double(params.processing_fbcsp_mi.nchannels);
end
do_car_cvsa = true;
if isfield(params,'processing_fbcsp_cvsa')
    do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car);
    nchannels   = double(params.processing_fbcsp_cvsa.nchannels);
end
signal = signal(:, 1:nchannels);

use_mi   = ismember(paradigm,{'mi','hybrid'});
use_cvsa = ismember(paradigm,{'cvsa','hybrid'});
csp_mi=[]; slda_mi=[]; csp_cvsa=[]; slda_cvsa=[];
if use_mi,   csp_mi   = load_csp(params,'mi');   slda_mi   = load_slda(params,'mi');   end
if use_cvsa, csp_cvsa = load_csp(params,'cvsa'); slda_cvsa = load_slda(params,'cvsa'); end

% ── Pipeline ──────────────────────────────────────────────────────────────
fprintf('[browser] Running pipeline ...\n');
proc_base = struct('samplerate',fs,'chunk_size',chunk_size, ...
                   'bufsize',bufsize_proc,'filter_order',4);
p_mi_al=[]; header_mi=[];  feat_mi=[];
p_cv_al=[]; header_cv=[];  feat_cv=[];
if use_mi
    cfg = proc_base; cfg.do_car = do_car_mi; cfg.eog_names = eog_names;
    [feat_mi, header_mi] = apply_processing(signal, header, csp_mi, cfg);
    p_mi_al = apply_slda(feat_mi, slda_mi, csp_mi.bands);
end
if use_cvsa
    cfg = proc_base; cfg.do_car = do_car_cvsa; cfg.eog_names = eog_names;
    [feat_cv, header_cv] = apply_processing(signal, header, csp_cvsa, cfg);
    p_cv_al = apply_slda(feat_cv, slda_cvsa, csp_cvsa.bands);
end

art_cfg = params.ArtifactCfg.params;
art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
[art_flags,~] = detect_artifacts(signal, header, art_cfg, ...
    struct('samplerate',fs,'chunk_size',chunk_size,'bufsize_artifact',bufsize_art));

int_cfg = params.integrator;
if ~isfield(int_cfg,'increment'),            int_cfg.increment = 1; end
if ~isfield(int_cfg,'thresholds_rejection'), int_cfg.thresholds_rejection = []; end
if ~isfield(int_cfg,'cvsa_influence'),       int_cfg.cvsa_influence = 2.5; end
if ~isfield(int_cfg,'thresholds') || isempty(int_cfg.thresholds)
    int_cfg.thresholds = params.training_node.thresholds;
end
header_chunks = header_mi; if ~use_mi, header_chunks = header_cv; end
header_chunks.framerate = framerate;

trials    = integrate_signal(p_mi_al, p_cv_al, art_flags, header_chunks, int_cfg, paradigm);
n_trials  = numel(trials);
classes   = to_vec(int_cfg.classes);
n_cls     = numel(classes);
thresholds = to_vec(int_cfg.thresholds);
p_rest    = to_vec(int_cfg.init_val); p_rest = p_rest(1);

% ── Global time-series ────────────────────────────────────────────────────
n_samples = size(signal,1);
n_chunks  = max([size(p_mi_al,1), size(p_cv_al,1), numel(art_flags), 1]);
t_sig     = (0:n_samples-1) / fs;
t_chunk   = ((0:n_chunks-1) + 0.5) * chunk_size / fs;
t_total   = n_samples / fs;

p_mi_g    = NaN(n_chunks, n_cls);
p_cv_g    = NaN(n_chunks, n_cls);
p_fused_g = NaN(n_chunks, n_cls);
p_integ_g = NaN(n_chunks, n_cls);
art_g     = false(n_chunks,1);
naf = min(numel(art_flags), n_chunks);
art_g(1:naf) = art_flags(1:naf);

if ~isempty(p_mi_al)
    nr = min(size(p_mi_al,1),n_chunks);
    p_mi_g(1:nr,:) = p_mi_al(1:nr,:);
end
if ~isempty(p_cv_al)
    nr = min(size(p_cv_al,1),n_chunks);
    p_cv_g(1:nr,:) = p_cv_al(1:nr,:);
end
for ti = 1:n_trials
    np = trials(ti).n_pre; sc = trials(ti).start_chunk;
    % Reset frame: one chunk before CF starts (shows the 0.5 reset publish)
    if sc-1 >= 1 && sc-1 <= n_chunks
        p_integ_g(sc-1,:) = trials(ti).integrated(1,:);   % = init_val
    end
    % CF frames
    for j = 1:trials(ti).n_cf
        ci = sc+j-1; k = j+np;
        if ci>=1 && ci<=n_chunks && k<=size(trials(ti).raw,1)
            p_fused_g(ci,:) = trials(ti).raw(k,:);
            p_integ_g(ci,:) = trials(ti).integrated(k,:);
        end
    end
end

% ── Event colour map ──────────────────────────────────────────────────────
ev_code_list  = [CF_CODE, HIT_CODE, MISS_CODE, TO_CODE, classes(:)'];
ev_color_list = [[0.20 0.50 0.85]; [0.10 0.70 0.10]; [0.85 0.15 0.15]; ...
                 [0.92 0.82 0.00]; repmat([0.55 0.55 0.55],n_cls,1)];
ev_lbl_list   = [{'CF'},{'HIT'},{'MISS'},{'TIMEOUT'}, ...
                 arrayfun(@(ci)sprintf('Cue%d',ci), 1:n_cls, 'UniformOutput', false)];

fprintf('[browser] %.1f s  |  %d trials  |  %s\n', t_total, n_trials, upper(paradigm));

% ── Build GUI ─────────────────────────────────────────────────────────────
WIN_S = min(15.0, t_total);

fig = figure('Name', sprintf('GDF Browser  —  %s  [%s]', basename, upper(paradigm)), ...
             'Color','w','NumberTitle','off', ...
             'Position',[20 80 1560 780], ...
             'KeyPressFcn', @on_key);

ml=0.06; mr=0.01; mb=0.11; mt=0.04; gap=0.018;
h_int  = 0.38;
h_prob = 1 - mb - mt - h_int - gap;
ax_prob = axes('Parent',fig,'Position',[ml, mb+h_int+gap, 1-ml-mr, h_prob]);
ax_int  = axes('Parent',fig,'Position',[ml, mb,           1-ml-mr, h_int]);

slider = uicontrol('Style','slider','Units','normalized', ...
    'Position',[ml,0.005,1-ml-mr,0.040], ...
    'Min',0,'Max',max(1e-3,t_total-WIN_S),'Value',0, ...
    'SliderStep',[min(1,WIN_S/max(1,t_total-WIN_S))/20, ...
                  min(1,WIN_S/max(1,t_total-WIN_S))], ...
    'Callback',@on_slider);

uicontrol('Style','text','Units','normalized','BackgroundColor','w', ...
    'Position',[0.002,0.88,0.050,0.022],'String','Win (s):','FontSize',9, ...
    'HorizontalAlignment','right');
win_edit = uicontrol('Style','edit','Units','normalized', ...
    'Position',[0.002,0.855,0.050,0.027],'String',sprintf('%.0f',WIN_S), ...
    'FontSize',9,'Callback',@on_win_edit);

% ── State struct ───────────────────────────────────────────────────────────
S = struct();
S.t0 = 0;  S.win = WIN_S;
S.t_total = t_total;  S.fs = fs;  S.cs = chunk_size;
S.t_chunk = t_chunk;
S.p_mi_g  = p_mi_g;   S.p_cv_g  = p_cv_g;
S.p_fused_g = p_fused_g;  S.p_integ_g = p_integ_g;
S.art_g   = art_g;
S.ev_pos  = double(header.EVENT.POS) / fs;
S.ev_typ  = double(header.EVENT.TYP);
S.ev_code_list  = ev_code_list;
S.ev_color_list = ev_color_list;
S.thresholds = thresholds;  S.p_rest = p_rest;  S.n_cls = n_cls;
S.paradigm = paradigm;  S.use_mi = use_mi;  S.use_cvsa = use_cvsa;
S.ax_prob = ax_prob;  S.ax_int = ax_int;
S.slider = slider;  S.win_edit = win_edit;
S.basename = basename;

setappdata(fig, 'ev_lbl_list', ev_lbl_list);
guidata(fig, S);
draw_view(fig);

% ── Local functions (callable from callbacks via guidata/getappdata) ───────

function draw_view(fig_)
S_  = guidata(fig_);
ev_lbl_list_ = getappdata(fig_, 'ev_lbl_list');
t0_ = S_.t0;  t1_ = t0_ + S_.win;
fs_ = S_.fs;  cs_ = S_.cs;

c0 = max(1, floor(t0_*fs_/cs_)+1);
c1 = min(numel(S_.t_chunk), ceil(t1_*fs_/cs_));
if c0>c1, c1=c0; end

tc   = S_.t_chunk(c0:c1);
nf   = numel(tc);
dt_c = cs_ / fs_;

art_w = false(nf,1);
if nf>0 && c1<=numel(S_.art_g), art_w = S_.art_g(c0:c1); end
t_art = tc(art_w);

ev_in = S_.ev_pos >= t0_ & S_.ev_pos <= t1_;
ev_t  = S_.ev_pos(ev_in);
ev_ty = S_.ev_typ(ev_in);

% ── Probabilities  (P(class 1) view only) ─────────────────────────────────
ax = S_.ax_prob; cla(ax,'reset'); hold(ax,'on');
shade_ax(ax, t_art, dt_c, -0.04, 1.04);
c_mi=[0.15 0.68 0.30]; c_cv=[0.90 0.48 0.10]; c_fus=[0.55 0.25 0.72];
if S_.use_mi && c1>=c0
    yi = S_.p_mi_g(c0:c1,1);
    if ~all(isnan(yi))
        plot(ax,tc,yi,'Color',[c_mi,0.75],'LineWidth',1.8,'DisplayName','MI');
    end
end
if S_.use_cvsa && c1>=c0
    yi = S_.p_cv_g(c0:c1,1);
    if ~all(isnan(yi))
        plot(ax,tc,yi,'Color',[c_cv,0.75],'LineWidth',1.8,'DisplayName','CVSA');
    end
end
if strcmp(S_.paradigm,'hybrid') && c1>=c0
    yi = S_.p_fused_g(c0:c1,1);
    if ~all(isnan(yi))
        plot(ax,tc,yi,'Color',[c_fus,0.95],'LineWidth',2.4,'DisplayName','Fused');
    end
end
plot(ax,[t0_,t1_],[0.5,0.5],'--k','LineWidth',0.8,'HandleVisibility','off');
ev_lines(ax, ev_t, ev_ty, S_, ev_lbl_list_, -0.04, 1.04, true);
set(ax,'XLim',[t0_,t1_],'YLim',[-0.04,1.04], ...
    'Box','off','TickDir','out','FontSize',8,'XTickLabel',[]);
ylabel(ax,'P(cls 1)','FontSize',8);
title(ax, sprintf('[%s]  %s   t = %.1f – %.1f s  |  <- -> scroll', ...
      upper(S_.paradigm), S_.basename, t0_, t1_), ...
      'FontSize',8,'Interpreter','none');
legend(ax,'Location','northwest','FontSize',9);
grid(ax,'on');

% ── Integrator  (P(class 1) view, thresholds in P(cls1) space) ────────────
ax = S_.ax_int; cla(ax,'reset'); hold(ax,'on');
shade_ax(ax, t_art, dt_c, -0.04, 1.04);
yi = S_.p_integ_g(c0:c1,1);
if ~all(isnan(yi))
    plot(ax,tc,yi,'Color',[0.10 0.60 0.70],'LineWidth',2.2,'DisplayName','Integrator');
end
% Class-1 threshold (upper bound: buf[1] >= thr[1])
thr1 = S_.thresholds(1);
if isfinite(thr1)
    plot(ax,[t0_,t1_],[thr1,thr1],'--','Color',[0.15 0.50 0.80],'LineWidth',1.6, ...
         'DisplayName',sprintf('thr cls1 = %.2f',thr1));
end
% Class-2 threshold in P(cls1) space: buf[1] <= 1 - thr[2]
if S_.n_cls >= 2 && isfinite(S_.thresholds(2))
    thr2_view = 1 - S_.thresholds(2);
    plot(ax,[t0_,t1_],[thr2_view,thr2_view],'--','Color',[0.80 0.30 0.10],'LineWidth',1.6, ...
         'DisplayName',sprintf('1-thr cls2 = %.2f',thr2_view));
end
plot(ax,[t0_,t1_],[S_.p_rest,S_.p_rest],':k','LineWidth',1.0,'HandleVisibility','off');
ev_lines(ax, ev_t, ev_ty, S_, ev_lbl_list_, -0.04, 1.04, true);
set(ax,'XLim',[t0_,t1_],'YLim',[-0.04,1.04], ...
    'Box','off','TickDir','out','FontSize',8);
xlabel(ax,'time (s)','FontSize',9);
ylabel(ax,'buf[cls 1]','FontSize',8);
legend(ax,'Location','northwest','FontSize',9);
grid(ax,'on');

linkaxes([S_.ax_prob, S_.ax_int],'x');
drawnow limitrate;
end

function shade_ax(ax_, t_a, dt, ylo, yhi)
for ti = 1:numel(t_a)
    patch(ax_, t_a(ti)+[-dt/2,dt/2,dt/2,-dt/2], [ylo,ylo,yhi,yhi], ...
          [1 0.75 0.75],'EdgeColor','none','FaceAlpha',0.40,'HandleVisibility','off');
end
end

function ev_lines(ax_, et, ety, S__, ev_lbl_, ylo_, yhi_, show_lbl)
for ei = 1:numel(et)
    typ = ety(ei);
    idx = find(S__.ev_code_list == typ, 1);
    if ~isempty(idx)
        col = S__.ev_color_list(idx,:);
        lbl = ev_lbl_{idx};
    else
        col = [0.60 0.60 0.60]; lbl = num2str(typ);
    end
    plot(ax_,[et(ei),et(ei)],[ylo_,yhi_],'-','Color',[col,0.70],'LineWidth',1.3, ...
         'HandleVisibility','off');
    if show_lbl
        text(ax_, et(ei), yhi_*0.91+ylo_*0.09, lbl, ...
             'FontSize',6,'Color',col,'Rotation',90, ...
             'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
end
end

function on_slider(src,~)
fig_ = ancestor(src,'figure');
S_ = guidata(fig_);
S_.t0 = max(0, min(src.Value, S_.t_total - S_.win));
guidata(fig_, S_);
draw_view(fig_);
end

function on_win_edit(src,~)
fig_ = ancestor(src,'figure');
S_ = guidata(fig_);
v = str2double(src.String);
if isfinite(v) && v > 0.5
    S_.win = v;
    S_.slider.Max = max(S_.slider.Min+1e-3, S_.t_total - S_.win);
    S_.slider.Value = min(S_.slider.Value, S_.slider.Max);
    S_.t0 = S_.slider.Value;
    guidata(fig_, S_);
    draw_view(fig_);
end
end

function on_key(src,ev)
S_ = guidata(src);
step = S_.win * 0.45;
changed = true;
switch ev.Key
    case 'rightarrow', S_.t0 = min(S_.t0+step, S_.t_total-S_.win);
    case 'leftarrow',  S_.t0 = max(S_.t0-step, 0);
    otherwise,         changed = false;
end
if changed
    S_.slider.Value = max(S_.slider.Min, min(S_.t0, S_.slider.Max));
    guidata(src, S_);
    draw_view(src);
end
end
