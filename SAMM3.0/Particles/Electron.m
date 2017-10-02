classdef Electron
    % Electron
    
    properties (Constant)
        
        charge = -PhysicalConstants.PositronCharge;
        mass   =  PhysicalConstants.ElectronMass;
        g      =  PhysicalConstants.PositrongFactor;
        
        Cq     = 55*PhysicalConstants.PlanckConstant/32/sqrt(3)/2/pi/PhysicalConstants.SpeedOfLight/PhysicalConstants.ElectronMass;
        Cgamma = PhysicalConstants.PositronCharge^2/3/PhysicalConstants.VacuumPermittivity/(PhysicalConstants.ElectronMass*PhysicalConstants.SpeedOfLight^2)^4;

    end % properties
    
end % classdef