function [xmz, xpz] = coord(Nz,DZ)
j=0:Nz;
% i=(0:Nz)';
% xmz = DZ*(j - i);
% xpz = DZ*(j + i);
xmz = DZ*(j - Nz);
xpz = DZ*(j + Nz);
end

