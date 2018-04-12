% HW9.m
%
% Homework #9
% ECE 53922 -- 21ST CENTURY ELECTROMAGNETICS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;
clear all;

% UNITS
degrees = pi/180;

% OPEN FIGURE WINDOW WITH WHITE BACKGROUND
fig = figure('Color','w','Units','normalized','Position',[0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL PARAMETERS
a = 1;
Nxu = 512;
Nyu = Nxu;

% GRID PARAMETERS
Sx = 10;
Sy = Sx;
NRESLO = 10;
NRESHI = 10;

% SVL PARAMETERS
NP = 21;
NQ = NP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #1: BUILD A GRAYSCALE TRIANGLE UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL GRID
dxh = a/Nxu;
dyh = a/Nyu;

% CREATE AXES
xah = [0:Nxu-1]*dxh; xah = xah - mean(xah);
yah = [0:Nyu-1]*dyh; yah = yah - mean(yah);

% UNIT CELL MESHGRID
[YH,XH] = meshgrid(yah,xah);

% DEFINE TRIANGLE VERTICES
w = 0.9*a;
h = w*sqrt(3)/2;
v1 = [ 0 ; h/2 ];
v2 = [ -w/2 ; -h/2 ];
v3 = [ +w/2 ; -h/2 ];

% SECTION 1
p1 = v1;
p2 = v2;
D1 = (p2(2) - p1(2))*XH - (p2(1) - p1(1))*YH + p2(1)*p1(2) - p2(2)*p1(1);
D1 = -D1./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% SECTION 2
p1 = v1;
p2 = v3;
D2 = (p2(2) - p1(2))*XH - (p2(1) - p1(1))*YH + p2(1)*p1(2) - p2(2)*p1(1);
D2 = D2./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% SECTION 3
p1 = v2;
p2 = v3;
D3 = (p2(2) - p1(2))*XH - (p2(1) - p1(1))*YH + p2(1)*p1(2) - p2(2)*p1(1);
D3 = -D3./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% BUILD UNIT CELL
UC = min(D1,D2);
UC = min(UC,D3);
UC = UC - min(UC(:));
UC = UC / max(UC(:));

% DISPLAY UNIT CELL
ha = imagesc(xah,yah,UC');
ha = get(ha,'Parent');
set(ha,'YDir','normal');
axis equal tight;
colorbar;
colormap('Gray');
title('UNIT CELL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #2: CALCULATE GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOW RESOLUTION GRID
dx = a/NRESLO;
dy = a/NRESLO;

% CALCULATE GRID
Nx = ceil(Sx/dx);
dx = Sx/Nx;
Ny = ceil(Sy/dy);
dy = Sy/Ny;

% CALCULATE AXES
xa = [0 : Nx-1]*dx; xa = xa - mean(xa);
ya = [0 : Ny-1]*dy; ya = ya - mean(ya);

% CREATE A LOW-RES MESGRID
[Y,X] = meshgrid(ya,xa);

% HIGH RESOLUTION GRID
Kmax = (2*pi/a) * [floor(NP/2) ; floor(NQ/2)];
amin = 2*pi/norm(Kmax);
dx2 = amin/NRESHI;
dy2 = amin/NRESHI;

% CALCULATE GRID
Nx2 = ceil(Sx/dx2);
Ny2 = ceil(Sy/dy2);

% CALCULATE AXES
xa2 = linspace(xa(1),xa(Nx),Nx2);
ya2 = linspace(ya(1),ya(Nx),Ny2);
dx2 = xa2(2) - xa2(1);
dy2 = ya2(2) - ya2(1);

% CREATE A HI-RES MESHGRID
[Y2,X2] = meshgrid(ya2,xa2);

% REPORT VALUES
display(['Nx  = ' num2str(Nx)]);
display(['Ny  = ' num2str(Ny)]);
display(['dx  = ' num2str(dx)]);
display(['dy  = ' num2str(dy)]);
display(['Nx2 = ' num2str(Nx2)]);
display(['Ny2 = ' num2str(Ny2)]);
display(['dx2 = ' num2str(dx2)]);
display(['dy2 = ' num2str(dy2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #3: GENERATE INPUT MAPS FOR SPATIAL VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERIOD
RSQ = X.^2 + Y.^2;
r   = 3.33*a;
PER = RSQ <= r^2;
PER = a*(1-PER) + 0.5*a*PER;

% ANGLE
x = [0 : Nx-1]*dx;
y = [0 : Ny-1]*dy;
[Yy,Xx] = meshgrid(y,x);
THETA = atan2(Xx,Yy);
THETA = flipud(fliplr(flipud(THETA)));
THETA = THETA/degrees;

% THRESHOLD
thresh = linspace(0.15,0.35,Nx2*Ny2)';
THRESH = ones(Nx2,Ny2);
THRESH = thresh.*THRESH(:);
THRESH = reshape(THRESH,Nx2,Ny2);
THRESH = THRESH';

% VISUALIZE
figure('Color','w','Units','normalized','Position',[0 0 1 1]);
subplot(131);
imagesc(xa,ya,PER');
colorbar;
axis equal tight off
title('PER');

subplot(132);
imagesc(xa,ya,THETA');
colorbar;
axis equal tight off
title('THETA');

subplot(133);
imagesc(xa2,ya2,THRESH');
colorbar;
axis equal tight off
title('THRESH');

gth_dat = linspace(0,1,100);
ff_dat  = 0*gth_dat;
for n = 1 : length(gth_dat)
   % GENERATE BINARY UNIT CELL
   UCB = UC > gth_dat(n);
   
   % FILL FRACTION
   ff_dat(n) = sum(UCB(:))/(Nxu*Nyu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #4: GENERATE LIST OF PLANAR GRATINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DECOMPOSE UNIT CELL
ER = fftshift(fft2(UC))/(Nxu*Nyu);

% DEFINE np and nq
np = NP;
nq = np;

% TRUNCATE EXPANSION
p0 = 1 + floor(Nxu/2);
q0 = 1 + floor(Nyu/2);
p1 = p0 - floor(np/2);
p2 = p0 + floor(np/2);
q1 = q0 - floor(nq/2);
q2 = q0 + floor(nq/2);
ERF = ER(p1:p2,q1:q2);

% CALCULATE FOURIER AXES
pa = [-floor(np/2) : +floor(np/2)];
qa = [-floor(nq/2) : +floor(nq/2)];

% CALCULATE GRATING VECTOR EXPANSION
KX = 2*pi*pa/a;
KY = 2*pi*qa/a;
[KY,KX] = meshgrid(KY,KX);

% DISPLAY TRUNCATED UNIT CELL
figure('Color','w','Units','normalized','Position',[0 0 1 1]);
imagesc(pa,qa,abs(ERF)');
colorbar;
hold on;
quiver(pa,qa,KX',KY','Color','w');
hold off;
axis equal tight off
title('PLANAR GRATING EXPANSION');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #5: GENERATE SPATIALLY VARIANT LATTICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE SWITCH VARIABLES
LATTICE = 1;

% DEFINE INPUTS
switch LATTICE
    % Lattice 1: Uniform Lattice
    case 1
        PER     = a*ones(Nx2,Ny2);
        THETA   = 0*degrees*ones(Nx2,Ny2);
        THRESH  = 0.3*ones(Nx2,Ny2);     
end