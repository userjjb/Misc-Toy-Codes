function Jackle29=Jackle29(Gi,Gii, g, rho, gamma, k, h, df, f)
omega = 2*pi*(f+(1i*df)); %complex freuqnecy

Nu = ((Gii/omega)-1i*(Gi/omega))/rho;   %kinematic viscocity

zeta = (1-(1-1i*omega)./(Nu * k^2)).^0.5;

z122 = (1+zeta.^2).^2;

coshkh = cosh(k*h);
coshzkh = cosh(zeta*k*h);
sinhkh = sinh(k*h);
sinhzkh = sinh(zeta*k*h);

N = -4*zeta.*(1+zeta.^2)+zeta.*(4+z122).*coshkh.*coshzkh-(4 * zeta.^2+z122).*sinhkh.*sinhzkh;
D = zeta.*sinhkh.*coshzkh-sinhzkh.*coshkh;

Jackle29 = 1./(gamma*k^2 + g*rho +rho.* Nu.^2 .*k^3.*(N./D));
end