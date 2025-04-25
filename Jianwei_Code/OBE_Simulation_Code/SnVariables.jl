const lam=303.4e-9;
const gam = 2 * pi * 31.8e6;
const normalizedBohrMag = 0.044; #\mu_{B}/\hbar\Gamma in units 1/Gauss
const mass = (118) * 1.67e-27;#mass of Sn

const kA = 2 * pi / lam; #wavenumber
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2