function [molecule] = def_molecule(moleculeName)

if moleculeName == "Sn"
    molecule.kA = 2*pi/303e-9;
    molecule.gam = 2*pi*31.8e6;
    molecule.mass = 118*1.67e-27;
end

if moleculeName == "Sr"
    molecule.kA = 2*pi/461e-9;
    molecule.gam = 2*pi*31.99e6;
    molecule.mass = 88*1.67e-27;
end

if moleculeName == "Ag"
    molecule.kA = 2*pi/328e-9;
    molecule.gam = 2*pi*23.6e6;
    molecule.mass = 109*1.67e-27;
end

end