function [Em, Ep, SINm, SINp, COSm, COSp] = EXP(xmz, xpz, DT)
argm = xmz.^2/(4*DT);
argp = xpz.^2/(4*DT);

SINm = sin(argm);
SINp = sin(argp);
COSm = cos(argm);
COSp = cos(argp);
Em = exp(-1i*argm);
Ep = exp(-1i*argp);
end

