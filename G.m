function val = G(x, z, t)
argp = -1i*(x + z).^2/(4*t);
argm = -1i*(x - z).^2/(4*t);

val = 1/2*sqrt(1i/(pi*t))*(exp(argm) - exp(argp));
end

