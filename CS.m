function [Sm, Sp, Cm, Cp] = CS(xmz, xpz, DT)
frargm =  xmz/sqrt(2*pi*DT);
frargp =  xpz/sqrt(2*pi*DT);

Sm = fresnels(frargm);
Sp = fresnels(frargp);
Cm = fresnelc(frargm);
Cp = fresnelc(frargp);
end

