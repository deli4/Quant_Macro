function plot_VFI_results(gr, par, v0, v5, c0, c5, k0, k5, n0, n5, i0, i5, outdir)
% plot_VFI_results — tiled 2x2 policies + separate value figure, light theme, consistent styling.

    % Ensure output directory
    if ~exist(outdir,'dir'); mkdir(outdir); end

    % Indices for low/high productivity
    [~, j_low]  = min(gr.a);
    [~, j_high] = max(gr.a);

    % Visual settings
    lw       = 1.4;
    col_phi0 = [0 0 .8];    % blue for φ=0
    col_phi5 = [.85 0 0];   % red  for φ=5
    ls_lowA  = '-';         % solid for low A
    ls_highA = '--';        % dashed for high A
    lg = {'Low A, \phi=0','Low A, \phi=5','High A, \phi=0','High A, \phi=5'};

% -------------------- POLICIES: 2x2 TILE (LIGHT THEME) --------------------
fig1 = figure('Color','w');
set(fig1,'InvertHardcopy','off');

% Global graphic defaults for black text
set(fig1,'DefaultAxesXColor','k','DefaultAxesYColor','k', ...
         'DefaultAxesZColor','k','DefaultTextColor','k');

% Colors and styles
lw       = 1.4;
col_phi0 = [0 0 .8];    % blue
col_phi5 = [.85 0 0];   % red
ls_lowA  = '-';         % solid
ls_highA = '--';        % dashed
lg = {
    '\phi = 0, Low A', ...
    '\phi = 0, High A', ...
    '\phi = 5, Low A', ...
    '\phi = 5, High A'
};


t = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');
title(t,'Policy functions: Low A vs High A, \phi=0 vs \phi=5', ...
      'FontWeight','bold','Color','k');

% Consumption
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
plot(ax1, gr.k, c0(j_low,:),   'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_lowA);
plot(ax1, gr.k, c0(j_high,:),  'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_highA);
plot(ax1, gr.k, c5(j_low,:),   'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_lowA);
plot(ax1, gr.k, c5(j_high,:),  'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_highA);
title(ax1,'C(K,A)'); xlabel(ax1,'K'); ylabel(ax1,'Policy level');

% Capital
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
plot(ax2, gr.k, k0(j_low,:),   'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_lowA);
plot(ax2, gr.k, k0(j_high,:),  'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_highA);
plot(ax2, gr.k, k5(j_low,:),   'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_lowA);
plot(ax2, gr.k, k5(j_high,:),  'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_highA);
title(ax2,'K'''); xlabel(ax2,'K'); ylabel(ax2,'Policy level');

% Labor
ax3 = nexttile; hold(ax3,'on'); box(ax3,'on');
plot(ax3, gr.k, n0(j_low,:),   'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_lowA);
plot(ax3, gr.k, n0(j_high,:),  'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_highA);
plot(ax3, gr.k, n5(j_low,:),   'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_lowA);
plot(ax3, gr.k, n5(j_high,:),  'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_highA);
title(ax3,'N(K,A)'); xlabel(ax3,'K'); ylabel(ax3,'Policy level');

% Investment
ax4 = nexttile; hold(ax4,'on'); box(ax4,'on');
plot(ax4, gr.k, i0(j_low,:),   'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_lowA);
plot(ax4, gr.k, i0(j_high,:),  'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_highA);
plot(ax4, gr.k, i5(j_low,:),   'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_lowA);
plot(ax4, gr.k, i5(j_high,:),  'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_highA);
title(ax4,'I(K,A)'); xlabel(ax4,'K'); ylabel(ax4,'Policy level');

% Legend in white box, below grid
lgd = legend(ax1,lg,'Location','best');
try, lgd.Layout.Tile = 'south'; end
set(lgd,'Color','w','EdgeColor','k','TextColor','k');

exportgraphics(fig1, fullfile(outdir,'policy_functions.png'), 'Resolution',300);


    % -------------------- VALUE FUNCTION (LIGHT THEME) --------------------
    fig2 = figure('Color','w'); set(fig2,'InvertHardcopy','off');
    hold on; box on;
    plot(gr.k, v0(j_low,:),   'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_lowA);
    plot(gr.k, v0(j_high,:),  'LineWidth',lw,'Color',col_phi0,'LineStyle',ls_highA);
    plot(gr.k, v5(j_low,:),   'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_lowA);
    plot(gr.k, v5(j_high,:),  'LineWidth',lw,'Color',col_phi5,'LineStyle',ls_highA);
    xlabel('K'); ylabel('V(K,A)');
    title('Value function: low A vs high A, \phi=0 vs \phi=5');
% Legend: order by phi first, black text, white box
lg_v = {
    '\phi = 0, Low A', ...
    '\phi = 0, High A', ...
    '\phi = 5, Low A', ...
    '\phi = 5, High A'
};

lgd2 = legend(lg_v,'Location','best');
set(lgd2,'Color','w','EdgeColor','k','TextColor','k');



    exportgraphics(fig2, fullfile(outdir,'value_function.png'), 'Resolution',300);
end
