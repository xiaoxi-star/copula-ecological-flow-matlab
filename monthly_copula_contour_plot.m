Q = []; P = []

N_plot = 100;
nLevel = 20;

Qstar_all = zeros(1,12);
Pstar_all = zeros(1,12);
fmax_all  = zeros(1,12);

figure('Position', [50, 100, 1000, 1600]);
set(gcf,'Color','w');
t = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:12

    useNew = (m == 7) || (m == 8);

    q = Q(:, m);
    p = P(:, m);

    [pdQ, uQ] = fitQ(q, Dist(1, m));
    [pdP, uP, pdfP] = fitP(p, Dist(2, m));

    uv = [uQ, uP];
    theta = fitCopula(uv, Dist(3, m));
    printCopulaInfo(theta, Dist(3, m), m);

    pdfQ = @(x) arrayfun(@(xx) pdf(pdQ,xx), x);

    q_grid = linspace(min(q), max(q), N_plot);
    p_grid = linspace(min(p), max(p), N_plot);
    [Qmat, Pmat] = meshgrid(q_grid, p_grid);

    Fmat = zeros(size(Qmat));
    for i = 1:N_plot
        for j = 1:N_plot
            Fmat(i,j) = jointDensity( ...
                Qmat(i,j), Pmat(i,j), ...
                pdfQ, pdfP, pdQ, pdP, ...
                Dist(1,m), Dist(2,m), theta, Dist(3,m));
        end
    end

    [fmax, ~] = max(Fmat(:));

    rel_threshold = 0.9;
    c = rel_threshold * fmax;

    if useNew
        P_threshold = 0.90;

        mask_region = (Fmat >= c);
        mask_valid  = mask_region & (Pmat >= P_threshold);

        if ~any(mask_valid(:))
            error('Month %d: 高概率区域内没有满足 P≥%.2f 的点', m, P_threshold);
        end

        Qstar = max(Qmat(mask_valid));

        tol = 1e-12 * max(1, abs(Qstar));
        idx_star = mask_valid & (abs(Qmat - Qstar) <= tol);
        Pstar = max(Pmat(idx_star));

    else
        C = contourc(q_grid, p_grid, Fmat, [c c]);

        idx = 1;
        Q_candidates = [];
        P_candidates = [];

        while idx < size(C, 2)
            num_points = C(2, idx);
            if num_points == 0
                break;
            end
            q_line = C(1, idx+1 : idx+num_points);
            p_line = C(2, idx+1 : idx+num_points);
            
            Q_candidates = [Q_candidates, q_line];
            P_candidates = [P_candidates, p_line];
            
            idx = idx + num_points + 1;
        end

        if isempty(Q_candidates)
            error('未找到等值线，请降低 rel_threshold（如改为 0.95 或 0.90）');
        end

        [Qstar, idx_max] = max(Q_candidates);
        Pstar = P_candidates(idx_max);
    end

    Qstar_all(m) = Qstar;
    Pstar_all(m) = Pstar;
    fmax_all(m)  = c;

    nexttile
    hold on

    if useNew
        v = linspace(min(Fmat(:)), c, nLevel);
        contour(Qmat, Pmat, Fmat, v, 'LineWidth', 1.0);
    else
        v = linspace(min(Fmat(:)), fmax, nLevel);
        if ~any(abs(v - c) < 1e-10)
            v = sort([v, c]);
        end
        contour(Qmat, Pmat, Fmat, v, 'LineWidth', 1.2);
    end

    scatter(Qstar, Pstar, 100, 'r', 'filled');

    plot([Qstar Qstar], [min(p_grid) Pstar], '--r', 'LineWidth', 1.5);
    plot([min(q_grid) Qstar], [Pstar Pstar], '--r', 'LineWidth', 1.5);

    xline(min(q_grid), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xline(max(q_grid), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(min(p_grid), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(max(p_grid), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

    text(Qstar, Pstar, sprintf('  (%.2f, %.0f%%)', Qstar, 100*Pstar), ...
        'FontSize', 9, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    xlabel('$x^{EQ}\,(m^3/s)$', 'Interpreter','latex', 'FontSize', 11, 'FontName', 'Times New Roman');
    ylabel('$x^{GR}$', 'Interpreter','latex', 'FontSize', 11, 'FontName', 'Times New Roman');
    title(sprintf('Month %d', m), 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    box on
    grid on
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'FontName', 'Times New Roman');
    yticklabels(compose('%g%%', 100*yticks));
end
