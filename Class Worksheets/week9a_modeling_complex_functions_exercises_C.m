%% OPT 212: Modeling Complex Functions
% Category summary: Complex numbers, amplitude-phase form, phase visualization, 
% branch cuts, complex optical fields, plane waves, and superposition.
% 
% 

clear
close all
clc
format short
set(groot,'DefaultAxesFontSize',14,'DefaultTextFontSize',14, ...
    'DefaultLegendFontSize',14,'DefaultColorbarFontSize',14);
%% Quick Review:  Elementary Complex Functions
% A complex variable may be written in Cartesian form as
% 
% $$z=x+iy.$$
% 
% A complex function is equivalently a pair of real functions of two real variables:
% 
% $$f(z)=u(x,y)+iv(x,y).$$
% 
% Polynomial and rational functions of complex variables behave much as their 
% real counterparts do. The familiar Taylor series also extend naturally:
% 
% $$e^z=1+z+\frac{z^2}{2!}+\frac{z^3}{3!}+\cdots,$$
% 
% $$e^{iz}=1+iz-\frac{z^2}{2!}-\frac{iz^3}{3!}+\cdots,$$
% 
% $$e^{-iz}=1-iz-\frac{z^2}{2!}+\frac{iz^3}{3!}+\cdots.$$
% 
% From these series,
% 
% $$\sin z=\frac{e^{iz}-e^{-iz}}{2i},\qquad \cos z=\frac{e^{iz}+e^{-iz}}{2}.$$
% 
% The hyperbolic functions are similarly tied to exponentials:
% 
% $$\sinh z=\frac{e^z-e^{-z}}{2}=i\sin(-iz),$$
% 
% $$\cosh z=\frac{e^z+e^{-z}}{2}=\cos(-iz).$$
%% Multiplying and dividing complex numbers
% Let's consider two complex numbers:  $z_1=x_1+iy_1$  and $z_2=x_2+iy_2$.   
% 
% The addition or subtraction operators can just be written as the sum (or difference) 
% of the real and imaginary parts $z_1 + z_2=(x_1+ x_2) +i(y_1+y_2)$.  
% 
% Multiplication is like the product of two binomials:   $z_1 z_2=(x_1+ iy_1) 
% (x_2+ iy_2)=(x_1x_2-y_1y_2)+i(x_1y_2 + x_2y_1)$.   
% 
% MatLab has built-in options that automatically carry out this functions:

z1 = 1-2i ;
z2 = 2+1i ; 
zzplus = z1 + z2;
zzminus = z1 - z2;
disp(['zzplus = ' num2str(zzplus)  '    zzminus = ' num2str(zzminus)])
zzprod = z1*z2   ;
zzquot = z1/z2   ;
disp(['zzprod = ' num2str(zzprod)  '    zzquot = ' num2str(zzquot)])
%% 
% The sum, difference and products can be verified by inspection.  But what 
% about the quotient?   The ordinary way to work out the quotion is by turning 
% the denominator real: $\frac{z_1}{ z_2}=\frac{z_1z_2^*}{ z_2z_2^*}=\frac{(x_1+ 
% iy_1)(x_2- iy_2)}{ (x_2.^2 + y_2^2)}$.
% 
% But there is an even easier way to visualize products and quotients.    
%% Amplitude-Phase Form in the complex plane
% The complex plane uses the real part as the horizontal coordinate and the 
% imaginary part as the vertical coordinate. The same number may be written in 
% amplitude-phase form:
% 
% $$z=re^{i\phi}=x+iy.$$
% 
% The coordinate transformations are
% 
% $$x=r\cos\phi,\qquad y=r\sin\phi,$$
% 
% $$r=|z|=\sqrt{x^2+y^2},\qquad \phi=\tan^{-1}(y/x).$$
% 
% The complex conjugate is
% 
% $$z^*=x-iy=re^{-i\phi}.$$
% 
% The amplitude-phase form greatly simplifies multiplication and division. If 
% $z_1=r_1e^{i\phi_1}$ and $z_2=r_2e^{i\phi_2}$, then
% 
% $$z_1z_2=r_1r_2e^{i(\phi_1+\phi_2)},$$
% 
% and
% 
% $$\frac{z_1}{z_2}=\frac{r_1}{r_2}e^{i(\phi_1-\phi_2)}.$$
% 
% Thus multiplication adds phases, while division subtracts them.
% 
% There are built in functions for the magnitude and phase in MatLab:  

z1_abs = abs(z1)
z2_abs = abs(z2)
z1_phase = angle(z1)
z2_phase = angle(z2)

product_check = z1_abs*z2_abs*exp(1i*(z1_phase+z2_phase))
quotient_check = z1_abs/z2_abs*exp(1i*(z1_phase-z2_phase))
%% Visualizing Complex Functions
% 
% 
% Complex functions require more than one surface plot because the output has 
% both real and imaginary parts. For example, if $f(z)=z^2$ and $z=x+iy$,
% 
% $$f(z)=(x+iy)^2=x^2-y^2+i\,2xy.$$
% 
% Thus
% 
% $$\it{Re}\{f(z)\}=x^2-y^2,\qquad \it{Im}\{f(z)\}=2xy.$$
% 
% In amplitude-phase form, with $z=re^{i\phi}$,
% 
% $$f(z)=r^2e^{i2\phi}.$$
% 
% The magnitude and phase are therefore
% 
% $$|f(z)|=r^2=x^2+y^2,\qquad \arg f(z)=2\phi \pmod{2\pi}.$$
% 
% 
% 
% The results below should look something like this: 
% 
% 

% Set up a meshgrid of real values of x & y
range = linspace(-3,3,201);
[x,y]=meshgrid(range);
% Define an array of complex values z:
zz = x + 1i.*y; 

% Define functions


figure
subplot(1,2,1)
surf(x,y,fmag); shading interp
title('Amplitude')
subplot(1,2,2)
surf(x,y,fphase); shading interp
title('Phase')
%% 
% Note that MATLAB's |angle| function reports phase wrapped into the interval 
% $[-\pi,\pi]$. This wrapping creates jumps in the plotted phase even when the 
% underlying algebraic phase is smooth.  The |unwrap| function removes jumps of 
% size $2\pi$ along a chosen direction when a continuous phase surface is desired. 
% However, sometimes there are singularities that don't allow full unwrapping.   


%% 
% 
% 
% Zeros of a complex function appear as points where the magnitude vanishes 
% and the phase becomes undefined. For example,
% 
% $$f(z)=(z-z_1)(z-z_0)^2$$
% 
% has a simple zero at $z_1$ and a double zero at $z_0$. Around a zero, the 
% phase winds by a multiple of $2\pi$, with the multiplicity of the zero controlling 
% the number of phase cycles.

% Define new functions
z1= 1;
z0=2i;
fz = @(z) (z-z1).*(z-z0).^3;
fmag = abs(fz(zz));
fphase = angle(fz(zz));

% 
figure
surf(x,y,fmag); shading interp
colorbar
figure
imagesc(range,range,fphase); axis image
colorbar
%% 
% Sometimes it is useful to visualize both on the same plot.  

figure
surf(x,y,fmag,fphase); shading interp
colorbar('vert')
%% 
% Here we are using the phase value to determine the color.  
%% 
%% 
%% 
%% Multivalued Functions and Branch Cuts
% The amplitude-phase representation makes logarithms and roots possible for 
% complex numbers, but it also reveals why they are multivalued. Since
% 
% $$z=re^{i\phi}=re^{i(\phi+2m\pi)},\qquad m\in\mathbb{Z},$$
% 
% where m is any integer, the square root becomes
% 
% $$z^{1/2}=r^{1/2}e^{i(\phi+2m\pi)/2} =r^{1/2}e^{i\phi/2}e^{im\pi}.$$
% 
% If m=0,2,4 etc.  the square root will have a positive  value and if m=1,3,5 
% etc, it will have a negative value.   

z1 = 2 - 1i;
sqrt(z1)


%% The complex logarithm
% Since adding $2m\pi$ to a phase leaves $z=re^{i\phi}$ unchanged, logarithms  
% are multivalued:
% 
% $$\log z=\log r+i(\phi+2m\pi).$$
% 
% 

nsheets = 5;
% initialize the independent variable for multiple sheets
zabs = [];
zphase = [];
xx = [];
yy = []; 
for ij = 1:nsheets
    zabs = [zabs abs(zz) abs(zz).*NaN];
    zphase = [zphase angle(zz)+(ij-1)*2*pi abs(zz).*NaN];
    xx = [xx x abs(zz).*NaN];
    yy = [yy y abs(zz).*NaN];
end

logz = log(zabs)+1i*(zphase);
figure; hold on
surf(xx,yy,abs(logz)); shading interp
colormap(jet)
colorbar
hold off

view([-30.29 35.38])
%% 


%% Plane Waves in Two and Three Dimensions
% Adapted from the original PDF discussion, Part 3, pages 18-26. Figures are 
% omitted here; the runnable MATLAB cells below regenerate the computational figures 
% used in the notes.
% 
% The multidimensional wave equation uses the Laplacian:
% 
% $$\nabla^2 A(\mathbf{r},t) =\frac{n^2}{c^2}\frac{\partial^2 A(\mathbf{r},t)}{\partial 
% t^2},$$
% 
% where
% 
% $$\nabla^2= \frac{\partial^2}{\partial x^2} +\frac{\partial^2}{\partial y^2} 
% +\frac{\partial^2}{\partial z^2}.$$
% 
% For a two-dimensional plane wave,
% 
% $$A(y,z,t)=A_0\exp[i(k_zz+k_yy-\omega t)].$$
% 
% The wave equation is satisfied when
% 
% $$k_z^2+k_y^2=n^2k_0^2,\qquad k_0=\frac{2\pi}{\lambda_0}.$$
% 
% In three dimensions the wave vector is
% 
% $$\mathbf{k}=k_x\hat{x}+k_y\hat{y}+k_z\hat{z},$$
% 
% with magnitude
% 
% $$|\mathbf{k}|=\frac{2\pi n}{\lambda_0}.$$
% 
% Using a unit propagation direction $\hat{\alpha}$,
% 
% $$\mathbf{k}=\frac{2\pi n}{\lambda_0}\hat{\alpha}.$$
% 
% Surfaces of constant phase satisfy
% 
% $$\mathbf{k}\cdot\mathbf{r}-\omega t=2m\pi.$$
% 
% At fixed time these are planes perpendicular to $\mathbf{k}$ and separated 
% by one wavelength in the direction of propagation.
%% 


%% 


%% Superposition of Plane Waves
% Adapted from the original PDF discussion, Part 3, pages 27-37. Figures are 
% omitted here; the runnable MATLAB cells below regenerate the computational figures 
% used in the notes.
% 
% The wave equation is linear. If $f$ and $g$ are solutions, then
% 
% $$\nabla^2 f=\frac{n^2}{c^2}\frac{\partial^2 f}{\partial t^2}, \qquad \nabla^2 
% g=\frac{n^2}{c^2}\frac{\partial^2 g}{\partial t^2}.$$
% 
% For constants $\alpha$ and $\beta$,
% 
% $$\nabla^2(\alpha f+\beta g) =\alpha\nabla^2 f+\beta\nabla^2 g,$$
% 
% and
% 
% $$\frac{\partial^2}{\partial t^2}(\alpha f+\beta g) =\alpha\frac{\partial^2 
% f}{\partial t^2} +\beta\frac{\partial^2 g}{\partial t^2}.$$
% 
% Therefore $\alpha f+\beta g$ is also a solution.
% 
% A general superposition of plane waves is
% 
% $$A(\mathbf{r},t)= \sum_{m=1}^{N}a_m \exp[i(\mathbf{k}_m\cdot\mathbf{r}-\omega_m 
% t)].$$
% 
% The coefficients $a_m$ are complex numbers that carry relative amplitude and 
% phase. For a monochromatic superposition, every wave has the same $\omega$:
% 
% $$A(\mathbf{r},t)=e^{-i\omega t} \sum_{m=1}^{N}a_m\exp(i\mathbf{k}_m\cdot\mathbf{r}).$$
% 
% The common time factor disappears from the intensity,
% 
% $$|A(\mathbf{r},t)|^2= \left|\sum_{m=1}^{N}a_m\exp(i\mathbf{k}_m\cdot\mathbf{r})\right|^2,$$
% 
% so the irradiance pattern is time independent even though the field oscillates. 
% If the frequencies differ, no common factor can be removed and the field becomes 
% a temporal pulse or beat pattern.
%% Superposition of Waves
% The wave equation is linear, so a linear combination of solutions is also 
% a solution:
% 
% $$\alpha f(\mathbf{r},t)+\beta g(\mathbf{r},t).$$
% 
% For a collection of plane waves,
% 
% $$A(\mathbf{r},t)= \sum_{m=1}^N a_m\exp[i(\mathbf{k}_m\cdot\mathbf{r}-\omega_m 
% t)].$$

wave_int_snapshot(0.2,2,0,'surface')
wave_int_snapshot(0.7,64,0,'image')
%% Intensity of a Monochromatic Superposition
% If all plane waves have the same frequency, the common factor $e^{-i\omega 
% t}$ disappears when computing $|A|^2$. The intensity is therefore time-independent 
% even while the field oscillates.

[field0,zplot,yplot] = wave_int_field(0.7,64,0);
[field1,~,~] = wave_int_field(0.7,64,0.25);
intensity_difference = max(abs(abs(field0(:)).^2 - abs(field1(:)).^2))

figure('Name','Monochromatic wave-superposition intensity')
imagesc(zplot(1,:),yplot(:,1),abs(field0).^2)
axis image
colormap(parula)
colorbar
xlabel('z (microns)')
ylabel('y (microns)')

%% Temporal Superposition
% A pulse can be built from plane waves with the same direction but different 
% frequencies:
% 
% $$A(z,t)=\sum_m a_m\exp[i(k_mz-\omega_mt)].$$
% 
% Different frequencies imply different wave-vector magnitudes.

freq_super_snapshot(1e14,64,0)
freq_super_snapshot(1e14,64,0.4)

subplot(2,1,1)
view([353.07 38.92])
%% Local Functions


function wave_int_snapshot(half_angle,points,cycle_fraction,displayMode)
%WAVE_INT_SNAPSHOT Display a snapshot of equal-frequency plane-wave interference.

[ff,zplot,yplot] = wave_int_field(half_angle,points,cycle_fraction);

figure('Name',sprintf('Plane-wave superposition: %.2f rad, %d waves',half_angle,points))
if strcmp(displayMode,'image')
    imagesc(zplot(1,:),yplot(:,1),real(ff))
    axis image
    colormap(parula)
    colorbar
    xlabel('z (microns)')
    ylabel('y (microns)')
    title('Real part of wave superposition');
else
    surf(zplot,yplot,real(ff),'EdgeColor','none')
    shading interp
    xlabel('z (microns)')
    ylabel('y (microns)')
    zlabel('Displacement')
    title('Real part of wave superposition');
end
end

function [ff,zplot,yplot] = wave_int_field(half_angle,points,cycle_fraction)
%WAVE_INT_FIELD Compute equal-frequency plane-wave superposition.

n = 1.5;
thetaz = linspace(-half_angle,half_angle,points);
wavl0 = 0.5e-6;
k = 2*pi*n/wavl0;
kz = k.*cos(thetaz);
ky = k.*sin(thetaz);
c = 3e8;
freq = c/wavl0;
omega = 2*pi*freq;
period = 1/freq;
amplitude = 2;
zrange = 10*wavl0;
yrange = 10*wavl0;
zpoints = 128;
ypoints = 128;
t = cycle_fraction*period;

zz = linspace(-zrange/2,zrange/2,zpoints);
yy = linspace(-yrange/2,yrange/2,ypoints);
[z,y] = meshgrid(zz,yy);
ff = zeros(size(z));

for ij = 1:length(thetaz)
    ff = ff + amplitude*exp(1i.*(kz(ij).*z + ky(ij).*y - omega.*t));
end

ff = ff./max(abs(ff),[],'all');
zplot = z*1e6;
yplot = y*1e6;
end

function freq_super_snapshot(spectral_width,points,cycle_fraction)
%FREQ_SUPER_SNAPSHOT Display one frame of temporal frequency superposition.

[ff,zplot,yplot] = freq_super_field(spectral_width,points,cycle_fraction);

figure('Name',sprintf('Temporal superposition: %.2g Hz, %d frequencies',spectral_width,points))
subplot(2,1,1)
surf(zplot,yplot,real(ff),'EdgeColor','none')
shading interp
zlim([-1 1])
xlabel('z (microns)')
ylabel('y (microns)')
zlabel('Displacement')
title('Real field');

subplot(2,1,2)
imagesc(zplot(1,:),yplot(:,1),abs(ff).^2)
axis image
colormap(parula)
colorbar
xlabel('z (microns)')
ylabel('y (microns)')
title('Intensity');
end

function [ff,zplot,yplot] = freq_super_field(spectral_width,points,cycle_fraction)
%FREQ_SUPER_FIELD Compute superposition of waves with different frequencies.

n = 1.5;
wavl0 = 0.5e-6;
c = 3e8;
thetaz = 0;
center_freq = c/wavl0;
freq = linspace(center_freq-spectral_width/2,center_freq+spectral_width/2,points);
omega = 2*pi*freq;
k = omega*n/c;
kz = k.*cos(thetaz);
ky = k.*sin(thetaz);
period = 1/center_freq;
amplitude = 1;
zrange = 20*wavl0;
yrange = 5*wavl0;
zpoints = 256;
ypoints = 32;
t = (cycle_fraction - 0.5)*5*period;

zz = linspace(-zrange/2,zrange/2,zpoints);
yy = linspace(-yrange/2,yrange/2,ypoints);
[z,y] = meshgrid(zz,yy);
ff = zeros(size(z));

for ij = 1:length(omega)
    ff = ff + amplitude*exp(1i.*(kz(ij).*z + ky(ij).*y - omega(ij).*t));
end

ff = ff./length(omega);
zplot = z*1e6;
yplot = y*1e6;
end