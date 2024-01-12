%% lueck_spectrum
% Generate the non-dimensional shear spectrum proposed by Lueck (2020).
%%
% <latex>\index{Functions!lueck_spectrum}</latex>
%
%%% Syntax
%   [phi, varargout] = lueck_spectrum( varargin )
%
% * [e] Rate of dissipation in units of watts per kg.
% * [nu] Kinematic viscosity in units of metre-squared per second. Default
%       value is 1e-6. 
% * [N] Number of spectral points. Default value is 1000.
% * [k] Wavenumber in cpm.
% * []
% * [phi] Lueck spectrum in units of per seconds-squared per cpm. The
%       length of phi is N, or the length of k.
% * [k] Wavenumber in cpm.
%
%%% Description
% This function generates a Lueck universal shear spectrum.  There are
% four basic forms of this function:
%
%    1) [phi,k] = lueck_spectrum( e, nu, N);
%    2)  phi    = lueck_spectrum( e, nu, k);
%    3) [phi,k] = lueck_spectrum( 0, N );
%    4)  phi    = lueck_spectrum( 0, k );
%
% Form 1: Return the Lueck spectrum for the dissipation rate $\texttt{e}$
% and viscosityn $\texttt{nu}$. The length of the returned spectrum
% $\texttt{phi}$ is $\texttt{N}$ and $\texttt{k}$ is the wavenumber in cpm.
% Default values are $\texttt{nu} = \num{1e-6}$ and $\texttt{N}$ = 1000. 
%
% Form 2: Same as form 1) except that the Lueck spectrum is evaluated at
% the wavenumbers given in the input vector $\texttt{k}$ [in cpm].
% $\texttt{k}$ must be a vector.
%
% Form 3: Return the non-dimensional Lueck spectrum (the G2 function in
% Oakey, 1982) of length $\texttt{N}$ points. The wavenumber is $k = k'/ks$
% [where $k'$ is in cpm, $k_s = (\epsilon/\nu^3)^{1/4}$ (see Oakey 1982)]
% and runs from $\num{1e-4}$ to 1.
%
% Form 4: Same as form 3) except that the non-dimensional spectrum is
% evaluated at the wavenumbers given in the input vector $\texttt{k}$ (in
% $k'/k_s$).
%
%%% Note
% For forms 1) and 2), the dissipation rate can be a vector, e.g.
% $\texttt{e}$ = [1e-7 1e-6 1e-5], in which case $\texttt{phi}$ is a matrix
% whose columns contain the dimensional Lueck spectra for the dissipation
% rates specified in $\texttt{e}$. 
%
% The form of the spectrum is computed from a fit by Lueck that is 
% documented in Internal Technical Note XXX, Spectral models of isotropic
% turbulence. 
%
%%% Examples
% Form 1:
%
%    >> [phi,k] = lueck_spectrum( 1e-7, 1.2e-6, 512 )
%    >> [phi,k] = lueck_spectrum( 1e-7, 1.2e-6 )
%    >> [phi,k] = lueck_spectrum( 1e-7 )
%
% Form 2:
%
%    >> phi = lueck_spectrum( 1e-7, 1.2e-6, logspace(-1,3,512) )
%
% Form 4:
%
%    >> phi = lueck_spectrum( 0, logspace(-3,0,512) )
%
%%% References:
%
% # Oakey, N. S., 1982: J. Phys. Ocean., 12, 256-271.
% # McMillan, J.M., A.E. Hay, R.G. Lueck and F. Wolk, 2016: Rates of
% dissipation of turbulent kinetic energy in a high Reynolds Number tidal
% channel. J. Atmos. and Oceanic. Techno., 33, 817-837.

% *Version History:*
%
% * 2020-03-04, RGL original version. 

function [phi,varargout] = lueck_spectrum(varargin)

% argument checking
narginchk(1,3);
[scaled,e,nu,N,k] = checkArgs(varargin,nargin);

if scaled  % forms 1) and 2)
   e = e(:)';
   Ne = length(e);
   ks = (e./nu.^3).^(1/4); % Kolmogorov wavenumber(s)
   ks = ks(ones(N,1),:);
   if isempty(k) % form 1)
      x = logspace(-4,0,N)';
      x = x(:, ones(1,Ne));
   else          % form 2)
      k = k(:,ones(1,Ne));
      x = k./ks;
   end
   
   y = (x/0.0152).^2;
   G2 = 8.048 * x.^(1/3) ./ (1 + (21.7*x).^3); % Lueck's improved fit slightly modified
   G2 = G2 .* (1 ./ (1 + (6.6*x).^(5/2))); % More modifications
   G2 = G2 .* (1 + 0.36 * y ./ ((y-1).^2 + 2*y)); % This catches the peaking of the spectrum
   k = x.*ks;
   e = e(ones(N,1),:);
   phi = e.^(3/4) * nu^(-1/4) .* G2;
   varargout{1} = k;
else    % forms 3) and 4)
   if isempty(k) % form 3)
      k = logspace(-4,0,N)';  % k = k_hat/k_s, as in Oakey 1982.
      varargout{1} = k;
   end
   x = k;
   y = (x/0.0152).^2;
   G2 = 8.048 * x.^(1/3) ./ (1 + (21.7*x).^3); % Lueck's improved fit slightly modified
   G2 = G2 .* (1 ./ (1 + (6.6*x).^(5/2))); % More modifications
   G2 = G2 .* (1 + 0.36 * y ./ ((y-1).^2 + 2*y)); % This catches the peaking of the spectrum
   phi = G2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [scaled,e,nu,N,k] = checkArgs(arg,narg);
% Helper function for Lueck.m

% the default values:
scaled = 0;
e  = 1e-6;
nu = 1e-6;
N  = 1000;
k  = [];


if any(arg{1} ~= 0) % it's form 1) or 2)
   scaled = 1;
   if narg == 3
      e = arg{1};
      nu = arg{2};
      N = arg{3};
   elseif narg == 2
      e = arg{1};
      nu = arg{2};
   elseif narg == 1
      e = arg{1};
   end
else                % it is form 3) or 4)
   scaled = 0;
   N = arg{2};
end

if length(N) > 1 % last argument is vector means it's a wavenumber vector
   if all(size(N)>1)
      error('Sorry, can''t have matrix wavenumber input.');
   else
      k = N(:);
      N = length(k);
   end
end
