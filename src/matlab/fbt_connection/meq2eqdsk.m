function eqdsk=meq2eqdsk(L,LY,varargin)
% MEQ2EQDSK - Convert meq output structures into G-EQDSK format
%
%   eqdsk=meq2eqdsk(L,LY [,ParameterName,ParameterValue])
%
% INPUT
%   L   STRUCT  MEQ parameter structure
%   LY  STRUCT  MEQ result structure; must contain surface integrals
%               ('iterq'>1) and only a single time.
% 
% PARAMETERS
%   'fname'   CHAR    Name of the file to be written, if empty no file will be written
%   'facb0'   DOUBLE  Vacuum B-field multiplier (default: 1 except if
%                     tok=TCV then 0.996, see the following link for details:
%                     https://gitlab.epfl.ch/spc/tcv/tbx/meq/-/issues/41)
%   'cocos'   NUMERIC COCOS convention for output structure. By default COCOS 17
%                     as in the rest of the MEQ suite will be used. Only 11 and 17
%                     have been implemented so far.
%                     See https://spc.epfl.ch/cocos for more details.
%
% OUTPUT
%   eqdsk  STRUCT  Equilibrium structure containing all G-EQDSK data.
%
% EXAMPLE
%   [L,LX,LY]=liuqe(65402,1,'iterq',50);
%   eqdsk=meq2eqdsk(L,LY,'facB0',0.996);
%
% Based on 'meq2eqdskval' from the CHEASEgui MATLAB toolbox by Holger Reimerdes (SPC-EPFL).
% The main differences include:
%   - the plasma boundary is the one computed by RTCI and has L.noq points
%   - the q-profile is computed via cubic spline interpolation of LY.iqQ
%     up to rho=1
%
% Details on the G-EQDSK file format can be found at https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf
%
% [+MEQ MatlabEQuilibrium Toolbox+] Swiss Plasma Center EPFL Lausanne 2022. All rights reserved.

%% Parameters
facb0=1;
if strcmpi(L.P.tokamak,'tcv')
  facb0=0.996;
end
in = inputParser;
in.addParameter('fname','',@ischar);
in.addParameter('facb0',facb0,@isdouble);
in.addParameter('cocos',17,@isnumeric);
in.parse(varargin{:});
p = in.Results;

%% Checks
assert(LY.nA==1,'Plasma must be single-axis.');
assert(LY.PQ(end)==0,'Edge pressure must be zero.')
assert(isfield(LY,'iqQ'),'Surface integral must have been calculated (''iterq''>1).')
assert(numel(LY.t)==1,'LIUQE structure must only contain a single time.')

%% FBT legacy parameters
if strcmp(func2str(L.bfct),'bffbt')
  it = iround(L.P.t,LY.t);
  L.bfp = L.bfp(:,it);
end

%% Constants
mu0=4e-7*pi;

%% Conventions
psiNeq=linspace(0,1,L.nrx)'; % All profiles on equidistant-psiN grid

eqdsk.cocos = p.cocos;
switch p.cocos
  case 11,   facpsi = -1;
  case 17,   facpsi = +1; % MEQ internal convention (same as TCV)
  otherwise, error('meq2eqdsk:cocos','meq2eqdsk can only produce structure with COCOS 11 or 17');
end

%% Fill structure going through fields as they appear in G-EQDSK file format
% 1st line: 48 character, dummy integer, nr, nz
ss=sprintf('%s (MEQ) #%0d t=%0.4f %s',L.P.tokamak,LY.shot,LY.t,date);        
eqdsk.stitle=sprintf('%-48s',ss);
eqdsk.ind1=3;
eqdsk.nr=L.nrx;
eqdsk.nz=L.nzx;
        
% 2nd line: rboxlen, zboxlen, r0, rboxlft, zboxmid
eqdsk.rboxlen=L.rx(end)-L.rx(1);
eqdsk.zboxlen=L.zx(end)-L.zx(1); 
eqdsk.r0=L.P.r0;
eqdsk.rboxleft=L.rx(1);
eqdsk.zmid= 0.5*(L.zx(1)+L.zx(end));

% 3rd line: rmag, zmag, psimag, psiedge, B0
eqdsk.raxis=LY.rA;
eqdsk.zaxis=LY.zA;
eqdsk.psiaxis=(LY.FA-LY.FB)*facpsi;
eqdsk.psiedge=0.0*facpsi;
eqdsk.b0=LY.rBt*p.facb0/L.P.r0;

% 4th line: Ip, psiax1 (duplicated), psiax2 (duplicated), raxis1 (duplicated), raxis2 (duplicated)
eqdsk.ip=LY.Ip;

% 5th line: zaxis1, zaxis2, psi_sep, R_xpoint, Z_xpoint
if LY.lX
  eqdsk.rxpt=LY.rB;
  eqdsk.zxpt=LY.zB;
else
  eqdsk.rxpt=0.;
  eqdsk.zxpt=0.;  
end

% Calculate profiles on the EQDSK-'Psi'mesh; include B-field correction.
[PpPsi,TTpPsi,PPsi,TPsi] = meqprof(L.fPg,L.fTg,LY.ag,psiNeq,LY.FA,LY.FB,LY.rBt*p.facb0,L.bfct,L.bfp,L.idsx,L.smalldia);
% 6th entry: T(psiN) on nr equidistant psiN mesh
eqdsk.F=TPsi;
% 7th entry: p(psiN) on nr equidistant psiN mesh
eqdsk.p=PPsi;
% 8th entry: TT'(psiN) on nr equidistant psiN mesh (in MKSA)
eqdsk.FFprime=TTpPsi*facpsi;
% 9th entry: p'(psiN) on nr equidistant psiN mesh (in MKSA)
eqdsk.pprime=PpPsi*facpsi;

% 10th entry: psi(i,j)
eqdsk.psirz=reshape(LY.Fx'-LY.FB,L.nrx*L.nzx,1)*facpsi;

% 11th entry: q profile on nr equidistant psiN mesh
% use cubic spline interpolation with 0 derivative in the center vs 0 and not-a-knot bc at 1
% (Multiplication/Division by T is used to account for facB0).
[M,tau] = csdec(L.pQ,L.pQ,'w');
eqdsk.q = TPsi./bspsum(tau,M*(LY.iqQ./LY.iTQ),sqrt(psiNeq),0)*facpsi;
if LY.lX, eqdsk.q(end) = eqdsk.q(end)*Inf; end

eqdsk.psimesh=psiNeq;
eqdsk.rhopsi=sqrt(psiNeq);

% 12th entry: (R,Z) plasma boundary and wall position
% eqdsk.nbbound=L.noq;
eqdsk.nblim=L.G.nl;
% Plasma boundary
eqdsk.rplas=LY.rq(:,end);
eqdsk.zplas=LY.zq(:,end);
eqdsk.nbbound = length(eqdsk.rplas);
% Wall position
eqdsk.rlim=L.G.rl(:);
eqdsk.zlim=L.G.zl(:);

%% End of official part of eqdsk. Add same last lines as expeq file for info
% eqdsk.extralines{1}=sprintf('%18.8e   psiaxis, psi0/2pi = %18.8e', ...
%     (LY.FA-LY.FB)*facpsi, (LY.FA-LY.FB)*facpsi/(2*pi));
% eqdsk.extralines{end+1}=sprintf('%18.8e   r-magaxe', LY.rA);
% eqdsk.extralines{end+1}=sprintf('%18.8e   z-magaxe', LY.zA);
% 
% zrmax = max(eqdsk.rplas);
% zrmin = min(eqdsk.rplas);
% zzmax = max(eqdsk.zplas);
% zzmin = min(eqdsk.zplas);
% za      = 0.5*(zrmax-zrmin);
% zkappa  = 0.5*(zzmax-zzmin)./za;
% zaver   = 0.5*(zzmax+zzmin);
% eqdsk.extralines{end+1}=sprintf('%18.8e   Z0 (zaver)', zaver);
% eqdsk.extralines{end+1}=sprintf('%18.8e   R0 [m]', L.P.r0);
% eqdsk.extralines{end+1}=sprintf('%18.8e   B0 [T]', LY.rBt*p.facb0/L.P.r0);
% eqdsk.extralines{end+1}=sprintf('%18.8e   SIGN OF B0 IN EXPERIMENT', sign(LY.rBt));
% eqdsk.extralines{end+1}=sprintf('%18.8e   TOTAL CURRENT -> I-p [A]:  %18.8e', ...
%     mu0*LY.Ip./(LY.rBt*p.facb0), LY.Ip);
% eqdsk.extralines{end+1}=sprintf('%18.8e   SIGN OF IP IN EXPERIMENT', sign(LY.Ip));
% eqdsk.extralines{end+1}=sprintf('%18.8e   kappa', zkappa);
%  
% eqdsk.extralines{end+1}=sprintf('%18.8e   q_0        zrhopol(1)  %18.8e', ...
%     eqdsk.q(1), 0);
% eqdsk.extralines{end+1}=sprintf('%18.8e   q_edge', eqdsk.q(end));
% eqdsk.extralines{end+1}=sprintf('%18.8e   beta_pol', LY.bp);
% eqdsk.extralines{end+1}=sprintf('%18.8e   beta_tor', LY.bt);
% eqdsk.extralines{end+1}=sprintf('%18.8e   li', LY.li);      % Is li(3)?
% 
% eqdsk.extralines{end+1}=sprintf('%18.8e   number of X points', LY.nX);
% for ix=1:LY.nX
%   eqdsk.extralines{end+1}=sprintf('%18.8e  R of Xpt nb: %d  -> in [m]: %18.8e', ...
%     LY.rX(ix)/L.P.r0, ix, LY.rX(ix));
%   eqdsk.extralines{end+1}=sprintf('%18.8e  Z of Xpt nb: %d  -> in [m]: %18.8e', ...
%     LY.zX(ix)/L.P.r0, ix, LY.zX(ix));
% end
% 
% eqdsk.extralines{end+1}=sprintf('%18.8e   time', LY.t);
% eqdsk.extralines{end+1}=sprintf('  %5i   shot number', LY.shot);
% eqdsk.extralines{end+1}=sprintf('  %2d   COCOS',eqdsk.cocos);

% % Duplicate information
% eqdsk.psi=(LY.Fx'-LY.FB)*facpsi;
% eqdsk.rmesh=L.rx;
% eqdsk.zmesh=L.zx;

%%
% sign and scale corrections

eqdsk.ip = -eqdsk.ip;
eqdsk.b0 = -eqdsk.b0;
eqdsk.FFprime = -eqdsk.FFprime;
eqdsk.F = -eqdsk.F;
eqdsk.psirz = -eqdsk.psirz / (2*pi);
eqdsk.psiaxis = -eqdsk.psiaxis / (2*pi);
eqdsk.psiedge = -eqdsk.psiedge / (2*pi);

if isinf(eqdsk.q(end))
  eqdsk.q(end) = eqdsk.q(end-1);
end

%% Write file
if ~isempty(p.fname)
  meqweqdsk(eqdsk,p.fname);
end