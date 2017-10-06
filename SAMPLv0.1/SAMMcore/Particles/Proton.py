classdef Proton
    % Proton
    
    properties (Constant)
        
        charge = PhysicalConstants.PositronCharge;  % 1.602e-19 C
        mass   = PhysicalConstants.ProtonMass;      % 1.6726231e-27 kg
self.__mass2 = self.__mass * self.__mass
        g      = 5.585694713;                       % g-factor
        
        Cq     = 55*PhysicalConstants.PlanckConstant/32/sqrt(3)/2/pi/PhysicalConstants.SpeedOfLight/PhysicalConstants.ProtonMass;
        Cgamma = PhysicalConstants.PositronCharge^2/3/PhysicalConstants.VacuumPermittivity/(PhysicalConstants.ProtonMass*PhysicalConstants.SpeedOfLight^2)^4;

    end % properties
    
end % classdef