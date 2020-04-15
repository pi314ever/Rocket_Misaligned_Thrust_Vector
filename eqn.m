function [dx] = eqn(t,x,al)
%eqn Equation of motions for rocket with misaligned engine
%   Rocket trajectory on flat earth, constant gravity, and  misaligned
%   engine with respect to rocket vertical axis.
%   INDEX   x    y    z    psi  theta phi
%   f(t)    1    3    5    7    9     11
%   f'(t)   2    4    6    8    10    12
global G m l Ix Iz T me  
dx = ones(12,1); % Initialize output
% Trivial setup
for ii = 1:6
    dx(2*ii-1) = x(2*ii);
end
w1 = x(10);
w2 = x(8).*sin(x(9));
w3 = x(8).*cos(x(9))+x(12);
wz = x(8).*cos(x(9));
R = sqrt(x(1).^2+x(3).^2+x(5).^2);

dx(2) = T/m*sin(x(9)).*sin(x(7))*cos(al)-G*me*x(1)./R.^3-T/m*sin(al)*(cos(x(11)).*cos(x(7))-sin(x(11)).*cos(x(9)).*sin(x(7)));
dx(4) = -T/m*sin(x(9)).*cos(x(7))*cos(al)-G*me*x(3)./R.^3-T/m*sin(al)*(cos(x(11)).*sin(x(7))+sin(x(11)).*cos(x(9)).*cos(x(7)));
dx(6) = T/m*cos(x(9))*cos(al)-G*me*x(5)./R.^3-T/m*sin(al)*sin(x(11)).*sin(x(9));
dx(8) = (Iz/Ix*w3.*w1-w1.*wz+l/2/Ix*T*sin(al)*cos(x(11))-wz.*x(10))./sin(x(9));
dx(10) = w2.*wz-Iz/Ix*w2.*w3-l/2/Ix*T*sin(al)*sin(x(11));
dx(12) = w2.*x(10)-dx(8).*cos(x(9));

end

