function val = jointDensity(qv, pv, pdfQ, pdfP, pdQ, pdP, qId, pId, theta, cId)
    uq = uCDF_Q(qv, pdQ, qId);
    up = uCDF_P(pv, pdP, pId);
    cden = copDen(uq, up, theta, cId);
    fq = pdfQ(qv);
    fp = pdfP(pv);
    val = fq .* fp .* cden;
    val = max(val, 0);
    val(~isfinite(val)) = 0;
end