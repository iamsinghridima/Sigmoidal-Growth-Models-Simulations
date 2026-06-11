% ===================================================================
%  SIGMOIDAL GROWTH MODELS
%  Logistic · Gompertz · Baranyi-Roberts · Richard's · Von Bertalanffy
%
%  Produces:
%    1. Console table  (t, N(t), N(t+1), dN/dt)  for t = 0..30
%    2. Figure matching the Python output:
%         - Top panel   : all 5 N(t) overlaid
%         - Middle rows : individual N(t) + dN/dt twin-axis panels
%         - Bottom-right: dN/dt comparison panel
%

clear; clc; close all;

%% ── PARAMETERS ─────────────────────────────────────────────────────
% Logistic
L_K  = 1000;  L_r  = 1.0;  L_N0  = 5;

% Gompertz
G_K  = 1000;  G_r  = 1.0;  G_N0  = 5;

% Baranyi-Roberts
BR_Nmax = 1000;  BR_mu = 1.0;  BR_N0 = 5;  BR_h0 = 5;

% Richard's (theta=1 reduces to Logistic)
RI_K = 1000;  RI_r = 1.0;  RI_N0 = 5;  RI_theta = 1.0;

% Von Bertalanffy
VB_Ninf = 1000;  VB_K = 1.0;  VB_N0 = 5;

% Time vectors
t_cont     = linspace(0, 30, 2000);   % continuous curve
t_discrete = 0:1:30;                  % integer points for scatter/table


%% ── COMPUTE CURVES ──────────────────────────────────────────────────
N_log  = logistic(t_cont, L_K, L_r, L_N0);
N_gomp = gompertz(t_cont, G_K, G_r, G_N0);
N_br   = baranyi(t_cont,  BR_Nmax, BR_mu, BR_N0, BR_h0);
N_ri   = richards(t_cont, RI_K, RI_r, RI_N0, RI_theta);
N_vb   = von_bertalanffy(t_cont, VB_Ninf, VB_K, VB_N0);

dN_log  = logistic_dNdt(N_log,  L_K, L_r);
dN_gomp = gompertz_dNdt(N_gomp, G_K, G_r);
dN_br   = baranyi_dNdt(t_cont,  N_br, BR_Nmax, BR_mu, BR_h0);
dN_ri   = richards_dNdt(N_ri,   RI_K, RI_r, RI_theta);
dN_vb   = von_bertalanffy_dNdt(N_vb, VB_Ninf, VB_K);

% Discrete points for scatter markers
Nd_log  = logistic(t_discrete, L_K, L_r, L_N0);
Nd_gomp = gompertz(t_discrete, G_K, G_r, G_N0);
Nd_br   = baranyi(t_discrete,  BR_Nmax, BR_mu, BR_N0, BR_h0);
Nd_ri   = richards(t_discrete, RI_K, RI_r, RI_N0, RI_theta);
Nd_vb   = von_bertalanffy(t_discrete, VB_Ninf, VB_K, VB_N0);


%% ── PRINT SIMULATION TABLES ─────────────────────────────────────────
models = {
    'LOGISTIC',        sprintf('K=%d  r=%.1f  N(0)=%.1f', L_K,  L_r,  L_N0),  ...
        @(t) logistic(t, L_K, L_r, L_N0),  ...
        @(t,N) logistic_dNdt(N, L_K, L_r);

    'GOMPERTZ',        sprintf('K=%d  r=%.1f  N(0)=%.1f', G_K,  G_r,  G_N0),  ...
        @(t) gompertz(t, G_K, G_r, G_N0),  ...
        @(t,N) gompertz_dNdt(N, G_K, G_r);

    'BARANYI-ROBERTS', sprintf('N_max=%d  mu=%.1f  N(0)=%.1f  h0=%.1f', BR_Nmax, BR_mu, BR_N0, BR_h0), ...
        @(t) baranyi(t, BR_Nmax, BR_mu, BR_N0, BR_h0), ...
        @(t,N) baranyi_dNdt(t, N, BR_Nmax, BR_mu, BR_h0);

    'RICHARD''S',      sprintf('K=%d  r=%.1f  N(0)=%.1f  theta=%.1f', RI_K, RI_r, RI_N0, RI_theta), ...
        @(t) richards(t, RI_K, RI_r, RI_N0, RI_theta), ...
        @(t,N) richards_dNdt(N, RI_K, RI_r, RI_theta);

    'VON BERTALANFFY', sprintf('N_inf=%d  K=%.1f  N(0)=%.1f', VB_Ninf, VB_K, VB_N0), ...
        @(t) von_bertalanffy(t, VB_Ninf, VB_K, VB_N0), ...
        @(t,N) von_bertalanffy_dNdt(N, VB_Ninf, VB_K);
};

for m = 1:size(models,1)
    name    = models{m,1};
    params  = models{m,2};
    Nt_fn   = models{m,3};
    dNdt_fn = models{m,4};

    fprintf('%s\n', repmat('=',1,72));
    fprintf('%s  |  %s\n', name, params);
    fprintf('%s\n', repmat('=',1,72));
    fprintf('%4s  %14s  %14s  %14s\n', 't', 'N(t)', 'N(t+1)', 'dN/dt');
    fprintf('%s\n', repmat('-',1,54));

    for ti = t_discrete(1:end-1)
        Ni   = Nt_fn(ti);
        Ni1  = Nt_fn(ti + 1);
        dNi  = dNdt_fn(ti, Ni);
        fprintf('%4.0f  %14.6f  %14.6f  %14.6f\n', ti, Ni, Ni1, dNi);
    end
    fprintf('\n');
end


%% ── COLOURS (matching Python) ───────────────────────────────────────
C_log  = [76  155 232]/255;
C_gomp = [78  203 113]/255;
C_br   = [245 166  35]/255;
C_ri   = [232  76 139]/255;
C_vb   = [181 123 238]/255;

BG    = [15  17  23]/255;
PANEL = [26  29  39]/255;
GRID  = [42  45  58]/255;
TXT   = [224 224 224]/255;


%% ── FIGURE ──────────────────────────────────────────────────────────
fig = figure('Color', BG, 'Position', [50 50 1200 1350]);

% ── Panel 1: all 5 N(t) overlaid ─────────────────────────────────
ax0 = axes('Parent', fig, 'Position', [0.06 0.78 0.90 0.17]);
hold(ax0, 'on');
plot(ax0, t_cont, N_log,  'Color', C_log,  'LineWidth', 2, 'DisplayName', 'Logistic');
plot(ax0, t_cont, N_gomp, 'Color', C_gomp, 'LineWidth', 2, 'DisplayName', 'Gompertz');
plot(ax0, t_cont, N_br,   'Color', C_br,   'LineWidth', 2, 'DisplayName', 'Baranyi-Roberts');
plot(ax0, t_cont, N_ri,   'Color', C_ri,   'LineWidth', 2, 'DisplayName', 'Richard''s');
plot(ax0, t_cont, N_vb,   'Color', C_vb,   'LineWidth', 2, 'DisplayName', 'Von Bertalanffy');
scatter(ax0, t_discrete, Nd_log,  18, C_log,  'filled', 'MarkerFaceAlpha', 0.8);
scatter(ax0, t_discrete, Nd_gomp, 18, C_gomp, 'filled', 'MarkerFaceAlpha', 0.8);
scatter(ax0, t_discrete, Nd_br,   18, C_br,   'filled', 'MarkerFaceAlpha', 0.8);
scatter(ax0, t_discrete, Nd_ri,   18, C_ri,   'filled', 'MarkerFaceAlpha', 0.8);
scatter(ax0, t_discrete, Nd_vb,   18, C_vb,   'filled', 'MarkerFaceAlpha', 0.8);
yline(ax0, 1000, 'w:', 'LineWidth', 0.7, 'Alpha', 0.35, 'DisplayName', 'Asymptote=1000');
style_ax(ax0, "All Five Sigmoidal Growth Models — N(t)  [N(0)=5 for all models]", 'Time (t)', 'N(t)', PANEL, GRID, TXT);
legend(ax0, {'Logistic','Gompertz','Baranyi-Roberts','Richard''s','Von Bertalanffy','Asymptote=1000'}, ...
    'TextColor', TXT, 'Color', PANEL, 'FontSize', 8.5, 'NumColumns', 3, 'Location', 'southeast');
hold(ax0,'off');

% ── Panels 2-6 + 7: layout positions ─────────────────────────────
positions = [
    0.06 0.535 0.41 0.20;
    0.55 0.535 0.41 0.20;
    0.06 0.295 0.41 0.20;
    0.55 0.295 0.41 0.20;
    0.06 0.055 0.41 0.20;
];

labels   = {'Logistic', 'Gompertz', 'Baranyi-Roberts', 'Richard''s', 'Von Bertalanffy'};
N_cells  = {N_log, N_gomp, N_br, N_ri, N_vb};
dN_cells = {dN_log, dN_gomp, dN_br, dN_ri, dN_vb};
Nd_cells = {Nd_log, Nd_gomp, Nd_br, Nd_ri, Nd_vb};
colors   = {C_log, C_gomp, C_br, C_ri, C_vb};

for p = 1:5
    col  = colors{p};
    Nc   = N_cells{p};
    dNc  = dN_cells{p};
    Ndc  = Nd_cells{p};
    lbl  = labels{p};

    ax = axes('Parent', fig, 'Position', positions(p,:)); %#ok<LAXES>
    hold(ax, 'on');
    fill(ax, [t_cont fliplr(t_cont)], [Nc zeros(size(Nc))], col, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none');
    plot(ax, t_cont, Nc, 'Color', col, 'LineWidth', 2, 'DisplayName', 'N(t)');
    scatter(ax, t_discrete, Ndc, 18, col, 'filled', 'MarkerFaceAlpha', 0.85);
    yline(ax, 1000, 'w:', 'LineWidth', 0.6, 'Alpha', 0.3);
    style_ax(ax, [lbl ' — N(t) & dN/dt'], 'Time (t)', 'N(t)', PANEL, GRID, TXT);

    ax2 = axes('Parent', fig, 'Position', positions(p,:), ...
                'Color','none', 'YAxisLocation','right', ...
                'XTick',[], 'YColor', col, 'XColor', TXT); %#ok<LAXES>
    hold(ax2, 'on');
    plot(ax2, t_cont, dNc, 'Color', col, 'LineWidth', 1.3, 'LineStyle', '--', 'DisplayName', 'dN/dt');
    ylabel(ax2, 'dN/dt', 'Color', col, 'FontSize', 7.5);
    ax2.GridColor = GRID;
    ax2.FontSize  = 7.5;

    % legend
    lh = legend([ax.Children(end-1); ax2.Children(end)], {'N(t)','dN/dt'}, ...
        'TextColor', TXT, 'Color', PANEL, 'FontSize', 7.5);
    try
        lh.BoxFace.ColorType = 'truecoloralpha';
        lh.BoxFace.ColorData = uint8([26 29 39 200]');
    catch; end

    hold(ax,'off'); hold(ax2,'off');
end

% ── Panel 7: dN/dt comparison ─────────────────────────────────────
ax_dN = axes('Parent', fig, 'Position', [0.55 0.055 0.41 0.20]);
hold(ax_dN, 'on');
plot(ax_dN, t_cont, dN_log,  'Color', C_log,  'LineWidth', 1.8, 'DisplayName', 'Logistic');
plot(ax_dN, t_cont, dN_gomp, 'Color', C_gomp, 'LineWidth', 1.8, 'DisplayName', 'Gompertz');
plot(ax_dN, t_cont, dN_br,   'Color', C_br,   'LineWidth', 1.8, 'DisplayName', 'Baranyi-Roberts');
plot(ax_dN, t_cont, dN_ri,   'Color', C_ri,   'LineWidth', 1.8, 'DisplayName', 'Richard''s');
plot(ax_dN, t_cont, dN_vb,   'Color', C_vb,   'LineWidth', 1.8, 'DisplayName', 'Von Bertalanffy');
style_ax(ax_dN, 'Growth Rate Comparison — dN/dt', 'Time (t)', 'dN/dt', PANEL, GRID, TXT);
legend(ax_dN, 'TextColor', TXT, 'Color', PANEL, 'FontSize', 7.5);
hold(ax_dN, 'off');

% Super-title
annotation(fig, 'textbox', [0 0.965 1 0.035], ...
    'String', {'Sigmoidal Growth Models  ·  Logistic · Gompertz · Baranyi-Roberts · Richard''s · Von Bertalanffy', ...
               'Parameters: K=1000, r=1.0, N(0)=5 for all models'}, ...
    'Color', 'white', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'BackgroundColor', 'none');

% Save
saveas(fig, 'sigmoidal_growth_models_plot.png');
fprintf('\nFigure saved: sigmoidal_growth_models_plot.png\n');


%% ════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  
% ════════════════════════════════════════════════════════════════════

function style_ax(ax, ttl, xlbl, ylbl, panel_c, grid_c, txt_c)
    ax.Color         = panel_c;
    ax.XColor        = txt_c;
    ax.YColor        = txt_c;
    ax.GridColor     = grid_c;
    ax.GridAlpha     = 0.7;
    ax.GridLineStyle = '--';
    ax.LineWidth     = 0.8;
    grid(ax, 'on');
    title(ax,  ttl,  'Color', txt_c, 'FontSize', 9,   'FontWeight', 'bold');
    xlabel(ax, xlbl, 'Color', txt_c, 'FontSize', 8.5);
    ylabel(ax, ylbl, 'Color', txt_c, 'FontSize', 8.5);
end

% ── Logistic ────────────────────────────────────────────────────────
function N = logistic(t, K, r, N0)
    N = K ./ (1 + ((K - N0) ./ N0) .* exp(-r .* t));
end
function dN = logistic_dNdt(N, K, r)
    dN = r .* N .* (1 - N ./ K);
end

% ── Gompertz ────────────────────────────────────────────────────────
function N = gompertz(t, K, r, N0)
    N = K .* exp(-log(K ./ N0) .* exp(-r .* t));
end
function dN = gompertz_dNdt(N, K, r)
    dN = r .* N .* log(K ./ N);
end

% ── Baranyi-Roberts ─────────────────────────────────────────────────
function A = baranyi_A(t, mu, h0)
    A = t + (1./mu) .* log(exp(-mu.*t) + exp(-h0) - exp(-mu.*t - h0));
end
function N = baranyi(t, N_max, mu, N0, h0)
    At  = baranyi_A(t, mu, h0);
    eAt = exp(mu .* At);
    N   = (N_max .* N0 .* eAt) ./ (N_max - N0 + N0 .* eAt);
end
function dN = baranyi_dNdt(t, N, N_max, mu, h0)
    num   = exp(-mu .* t);
    den   = exp(-mu .* t) + exp(-h0) - exp(-mu .* t - h0);
    alpha = num ./ den;
    dN    = alpha .* mu .* (1 - N ./ N_max) .* N;
end

% ── Richard's ───────────────────────────────────────────────────────
function N = richards(t, K, r, N0, theta)
    a = (K ./ N0).^theta - 1;
    N = K .* (1 + a .* exp(-r .* theta .* t)).^(-1./theta);
end
function dN = richards_dNdt(N, K, r, theta)
    dN = r .* N .* (1 - (N ./ K).^theta);
end

% ── Von Bertalanffy ─────────────────────────────────────────────────
function N = von_bertalanffy(t, N_inf, K, N0)
    N = N_inf - (N_inf - N0) .* exp(-K .* t);
end
function dN = von_bertalanffy_dNdt(N, N_inf, K)
    dN = K .* (N_inf - N);
end
