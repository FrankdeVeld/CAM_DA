function [primary,secondary] = generateInitShort(ind)
mu     = 398600.4418;    % [m^3/s^2]
load("conjunctions_leo.mat");
% load("train.mat");
B = table2array(data); clear data;

%% Primary
x0 = toColumn(B(ind,3:8));
primary = cartesian2kepler(x0,mu);
primary.x0 = x0;
primary.C0 = [[B(ind,9) B(ind,12)  B(ind,13);
               B(ind,12)  B(ind, 10) B(ind,14);
               B(ind,13)  B(ind, 14)  B(ind,11)] zeros(3,3); 
                                            zeros(3,6)];
primary.T      = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = B(ind,2)/2;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.ctrlMax = 1e-7;

%% Secondary
x0s                = toColumn(B(ind,15:20)); % [km] [km/s] Secondary initial state in ECI
secondary.x0       = x0s;         
secondary.C0       = [[B(ind,21) B(ind,24) B(ind,25);
                    B(ind,24)  B(ind,22) B(ind,26);
                    B(ind,25)  B(ind,26) B(ind,23)] zeros(3,3); 
                                            zeros(3,6)];
secondary.HBR      = B(ind,2)/2 + primary.HBR;         % [km]
secondary.mass     = 100;          % [kg] mass
secondary.A_drag   = 1;            % [m^2] drag surface area
secondary.Cd       = 2.2;          % [-] shape coefficient for drag
secondary.A_srp    = 1;            % [m^2] SRP surface area
secondary.Cr       = 1.31;         % [-] shape coefficient for SRP


end