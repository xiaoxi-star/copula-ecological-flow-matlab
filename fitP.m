function [pd, u, pdfP] = fitP(x, id)
% 对保证率 P 做拟合
eps1 = 1e-6;
x = min(max(x, eps1), 1-eps1);

switch id
    case 1 % Logit-Normal
        z = log(x./(1-x));
        [mu,sig] = normfit(z);
        pd = makedist('Normal','mu',mu,'sigma',sig);
        u = cdf(pd, log(x./(1-x)));
        u = min(max(u,1e-10),1-1e-10);
        pdfP = @(y) arrayfun(@(yy) pdf(pd, log(yy/(1-yy))) ./ (yy.*(1-yy)), y);
        
case 2 % Kumaraswamy
    xs = sort(x);
    n = numel(xs);
    empirical = (1:n)'/(n+1);
    
    % 对数变换：par(1)=log(a), par(2)=log(b)，保证 a,b > 0
    obj = @(par) sum(( (1 - (1 - xs.^exp(par(1))).^exp(par(2))) - empirical ).^2);
    
    par0 = [0 0];          % 对应 a=1, b=1
    options = optimset('Display','off');
    par = fminsearch(obj, par0, options);
    
    a = exp(par(1));       % 还原参数
    b = exp(par(2));
    
    pd = struct('a', a, 'b', b);
    u = 1 - (1 - x.^a).^b;
    u = min(max(u, 1e-10), 1-1e-10);
    pdfP = @(y) arrayfun(@(yy) a*b .* yy.^(a-1) .* (1 - yy.^a).^(b-1), y);
        
    case 3 % Beta 分布
        [phat, ~] = betafit(x);
        alpha = phat(1); beta = phat(2);
        pd = makedist('Beta','a',alpha,'b',beta);
        u = cdf(pd, x);
        u = min(max(u,1e-10),1-1e-10);
        pdfP = @(y) pdf(pd, y);
        
    case 4 % Unit-Logistic
        z = log(x./(1-x));
        mu = mean(z);
        sig = fminsearch(@(s) -sum(log(pdf(makedist('Logistic','mu',mu,'sigma',s), z))), std(z)*0.5);
        pd = struct('mu', mu, 'sigma', sig);
        u = arrayfun(@(zz) cdf(makedist('Logistic','mu',mu,'sigma',sig), zz), z);
        u = min(max(u,1e-10),1-1e-10);
        pdfP = @(y) arrayfun(@(yy) pdf(makedist('Logistic','mu',mu,'sigma',sig), ...
            log(yy/(1-yy))) ./ (yy.*(1-yy)), y);
        
    case 5 % Unit-Weibull
        z = -log(x);
        pd_wbl = fitdist(z, 'Weibull');
        alpha = pd_wbl.ParameterValues(1);
        beta = pd_wbl.ParameterValues(2);
        pd = struct('alpha', alpha, 'beta', beta);
        u = exp(-(z/beta).^alpha);
        u = min(max(u,1e-10),1-1e-10);
        pdfP = @(y) arrayfun(@(yy) (alpha/beta) * ((-log(yy))/beta)^(alpha-1) ...
            .* exp(-((-log(yy))/beta)^alpha) ./ yy, y);
        
    otherwise
        error('Unknown p distribution id %d', id);
end
end