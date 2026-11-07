function [pd, u] = fitQ(x, id)
% ·µ»Ø pd (distribution object) Óë u = CDF(x)
switch id
    case 1 % P3 -> Gamma
        ph = gamfit(x);
        pd = makedist('Gamma','a',ph(1),'b',ph(2));
    case 2 % GEV
        ph = gevfit(x);
        pd = makedist('GeneralizedExtremeValue','k',ph(1),'sigma',ph(2),'mu',ph(3));
    case 3 % Gumbel (ExtremeValue)
        ph = evfit(x);
        pd = makedist('ExtremeValue','mu',ph(1),'sigma',ph(2));
    case 4 % Lognormal
        ph = lognfit(x);
        pd = makedist('Lognormal','mu',ph(1),'sigma',ph(2));
    otherwise
        error('Unknown q distribution id %d', id);
end
u = cdf(pd, x);
% clip to (eps,1-eps)
u = min(max(u,1e-10),1-1e-10);
end