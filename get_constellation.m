function X = get_constellation(Nbps, type_mod, pwr)
%   X = get_constellation(Nbps, type_mod)
%   Generate all the constellation points using Matlab comm toolbox, using
%   the default Gray mapping.
% _____________________________________________________________________________
%	Inputs:
% 		Nbps:       scalar, number of bits per symbol
%       type_mod:   string, 'PSK' for phase shift keying and QAM for quadrature 
%                   amplitude modulation
%       pwr:        scalar, average power of the modulated symbols
%	Outputs:
%		X:			2 ^ Nbps-by-1 vector, the modulated constellations
% _____________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 12/18/2014
% Codename: Dunkirk
% _____________________________________________________________________________

Q = 2 ^ Nbps; % Number of constellation points

if strcmpi(type_mod, 'PSK')
    hmod = comm.PSKModulator(Q);
elseif strcmpi(type_mod, 'QAM')
    hmod = comm.RectangularQAMModulator(Q);
else
    error('type_mode must either be PSK or QAM')
end

X = step(hmod, (0 : Q - 1)');
pwr_tmp = X' * X / Q; % Average power of the constellation points
X = X * sqrt(pwr / pwr_tmp); % Normalization

end

