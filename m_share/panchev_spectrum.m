%% panchev_spectrum
% Generate the non-dimensional shear spectrum proposed by Kanchev-Kesich (1969).
%%
% <latex>\index{Functions!panchev_spectrum}</latex>
%
%%% Syntax
%   [phi, varargout] = panchev_spectrum( varargin )
%
% * [e] Rate of dissipation in units of watts per kg.
% * [nu] Kinematic viscosity in units of metre-squared per second. Default
%       value is 1e-6. 
% * [N] Number of spectral points. Default value is 1000.
% * [k] Wavenumber in cpm.
% * []
% * [phi] Panchev-Kesich spectrum in units of per seconds-squared per cpm. The
%       length of phi is N, or the length of k.
% * [k] Wavenumber in cpm.
%
%%% Description
% This function generates a Panchev-Kesich universal shear spectrum.  There are
% four basic forms of this function:
%
%    1) [phi,k] = panchev_spectrum( e, nu, N);
%    2)  phi    = panchev_spectrum( e, nu, k);
%    3) [phi,k] = panchev_spectrum( 0, N );
%    4)  phi    = panchev_spectrum( 0, k );
%
% Form 1: Return the Panchev-Kesich spectrum for the dissipation rate $\texttt{e}$
% and viscosityn $\texttt{nu}$. The length of the returned spectrum
% $\texttt{phi}$ is $\texttt{N}$ and $\texttt{k}$ is the wavenumber in cpm.
% Default values are $\texttt{nu} = \num{1e-6}$ and $\texttt{N}$ = 1000. 
%
% Form 2: Same as form 1) except that the Panchev-Kesich spectrum is evaluated at
% the wavenumbers given in the input vector $\texttt{k}$ [in cpm].
% $\texttt{k}$ must be a vector.
%
% Form 3: Return the non-dimensional Panchev-Kesich spectrum (the G2 function in
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
% whose columns contain the dimensional Panchev-Kesich spectra for the dissipation
% rates specified in $\texttt{e}$. 
%
% The form of the spectrum is computed from a fit by Panchev-Kesich that is 
% documented in Internal Technical Note XXX, Spectral models of isotropic
% turbulence. 
%
%%% Examples
% Form 1:
%
%    >> [phi,k] = panchev_spectrum( 1e-7, 1.2e-6, 512 )
%    >> [phi,k] = panchev_spectrum( 1e-7, 1.2e-6 )
%    >> [phi,k] = panchev_spectrum( 1e-7 )
%
% Form 2:
%
%    >> phi = panchev_spectrum( 1e-7, 1.2e-6, logspace(-1,3,512) )
%
% Form 4:
%
%    >> phi = panchev_spectrum( 0, logspace(-3,0,512) )
%
%%% References:
%
% # Panchev, S., Kesich, D., 1969. Energy spepctrum of isotropic turbulence 
% at large wavenumbers. Comptes rendus de l’acade ́mie Bulgare des
% sciences 22, 627–630.
% # Roget, E., I. Lozovatsky, X. Sanchez, and M. Figueroa, 2006:
% Microstructure measurements in natural waters: Methodology and
% applications. Progress in Oceanography, 70, 126-148.

% *Version History:*
%
% * 2021-10-20, RGL original version. 

function [phi,varargout] = panchev_spectrum(varargin)

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
   
   G2 = 0.954 * (x*2*pi).^(0.372) .* exp(-5.824 .* (x*2*pi).^(1.495) ); % The Panchev-Kesich spectrum approximation by Roget et al, 2006.
   G2= 2*pi*G2;

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
   
   G2 = 0.954 * (x*2*pi).^(0.372) .* exp(-5.824 .* (x*2*pi).^(1.495) ); % The Panchev-Kesich spectrum approximation by Roget et al, 2006.
   G2 = 2*pi*G2;

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
