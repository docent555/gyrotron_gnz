clear;

Nz = 1000;
DT = 0.025;
Zout = 47;
DZ = Zout/Nz;
ZAxis = 0:Zout/Nz:Zout;
f = @(z) cos(25*2*pi*z/Zout) + 1i*sin(25*2*pi*z/Zout);
fg = @(z) f(z).*G(z,Zout/2,DT);

fs = f(ZAxis);
fre = real(fs);
fim = imag(fs);
[reb,rec,red] = spline(Nz+1,ZAxis,fre);
[imb,imc,imd] = spline(Nz+1,ZAxis,fim);
B = reb + 1i*imb;
C = rec + 1i*imc;
D = red + 1i*imd;
S1 = @(z) seval_cmplx(z, Nz+1, ZAxis, fre, fim, reb, rec, red, imb, imc, imd);


for i=1:Nz+1
    y1(i) = fg(ZAxis(i));
    y2(i) = S1(ZAxis(i))*G(ZAxis(i),Zout/2,DT);
end

timerVal = tic;
[xmz, xpz] = coord(Nz,DZ);
[Sm, Sp, Cm, Cp] = CS(xmz, xpz, DT);
[Em, Ep] = EXP(xmz, xpz, DT);

j = 1:Nz;
[I1m, I1p, I2m, I2p, I3m, I3p, I4m, I4p] = ...
    I_for_lineq(Sm, Sp, Cm, Cp, Em, Ep, DZ, DT, Zout, Nz, j);

Ii = integral(fg,0,Zout)

% for i=1:Nz
%    S1(i) = 
%       + fs(i)*(I1p(i) + I1m(i)) ...
%       + B(i)*(I2p(i) + I2m(i)) ...
%       + C(i)*(I3p(i) + I3m(i)) ...
%       + D(i)*(I4p(i) + I4m(i));
% end
% Is = Is*1/2*sqrt(1i/(pi*DT))

plot(ZAxis,imag(y1),ZAxis,imag(y2))

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

