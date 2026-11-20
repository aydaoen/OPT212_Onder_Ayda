%% OPT 212: Complex Functions in wave optics
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
%% Big picture of Fourier Optics
% 
% 
% A fourier transform connects the pupil fields to the fields in the focal plane.  
% Numerically, this becomes a fast Fourier transform (FFT).   
%% Fourier Decomposition and Numerical Fourier Transforms
% Functions of finite value can be represented, over a chosen interval, as superpositions 
% of plane-wave-like functions. A periodic function may be represented by a discrete 
% plane-wave expansion:
% 
% $$f(\mathbf{r})=\sum_{m=1}^{N}a_m \exp(i\mathbf{k}_m\cdot\mathbf{r}).$$
% 
% For a one-dimensional periodic function of period $\Lambda$, the fundamental 
% spatial frequency is
% 
% $$f_0=\frac{1}{\Lambda}.$$
% 
% The allowed wave vectors satisfy
% 
% $$k_m\Lambda=2m\pi,\qquad k_m=2\pi m f_0.$$
% 
% The Fourier series is therefore
% 
% $$g(x)=\sum_{m=-N}^{N}a_m\exp[i2\pi(mf_0)x].$$
% 
% Orthogonality over one period gives
% 
% $$\int_{-\Lambda/2}^{\Lambda/2} e^{i2\pi(mf_0)x}e^{-i2\pi(nf_0)x}\,dx =\Lambda\delta_{mn}.$$
% 
% The coefficients are
% 
% $$a_n=\frac{1}{\Lambda}\int_{-\Lambda/2}^{\Lambda/2} g(x)e^{-i2\pi(nf_0)x}\,dx.$$
% 
% As $\Lambda\to\infty$, the spacing between allowed frequencies approaches 
% zero and the Fourier series becomes the continuous transform pair
% 
% $$g(x)=\int_{-\infty}^{\infty}G(f_x)e^{i2\pi f_xx}\,df_x,$$
% 
% $$G(f_x)=\int_{-\infty}^{\infty}g(x)e^{-i2\pi f_xx}\,dx.$$
% 
% The rectangle function is
% 
% $$\it{rect}(x)= \begin{cases} 1, & |x|\leq 1/2,\\ 0, & |x|>1/2. \end{cases}$$
% 
% Its transform is
% 
% $$\mathcal{F}\{\it{rect}(x)\} =\frac{\sin(\pi f_x)}{\pi f_x} =\it{sinc}(f_x).$$
% 
% For sampled data, a direct numerical Fourier transform has the form
% 
% $$G(f_j)\approx \Delta x\sum_k g(x_k)e^{-i2\pi f_jx_k}.$$
% 
% The FFT evaluates the same kind of discrete transform much more efficiently 
% at a regular grid of frequencies.
%% Fourier Decomposition
% A periodic function can be represented as
% 
% $$g(x)=\sum_{m=-N}^{N}a_m\exp[i2\pi(mf_0)x],$$
% 
% where $f_0=1/\Lambda$ is the fundamental spatial frequency. The orthogonality 
% relationship is
% 
% $$\int_{-\Lambda/2}^{\Lambda/2} e^{i2\pi(mf_0)x}e^{-i2\pi(nf_0)x}\,dx =\Lambda\delta_{mn}.$$
% 
% Therefore,
% 
% $$a_n=\frac{1}{\Lambda}\int_{-\Lambda/2}^{\Lambda/2} g(x)e^{-i2\pi(nf_0)x}\,dx.$$
%% The Continuous Fourier Transform
% In the limit $\Lambda\to\infty$, the frequency spacing approaches zero and 
% the Fourier series becomes the transform pair
% 
% $$g(x)=\int_{-\infty}^{\infty} G(f_x)e^{i2\pi f_xx}\,df_x,$$
% 
% $$G(f_x)=\int_{-\infty}^{\infty} g(x)e^{-i2\pi f_xx}\,dx.$$
%% *The Continuous Fourier Transform and Fourier Optics*
% The same Fourier formalism extends to two dimensions by allowing both the 
% position and frequency variables to be vectors. The transform pair is
% 
% $$g(\mathbf{r})= \int\!\!\int G(\mathbf{f}) e^{i2\pi\mathbf{f}\cdot\mathbf{r}}\,d^2f,$$
% 
% $$G(\mathbf{f})= \int\!\!\int g(\mathbf{r}) e^{-i2\pi\mathbf{f}\cdot\mathbf{r}}\,d^2r.$$
% 
% In coordinates, $\mathbf{r}=x\hat{x}+y\hat{y}$ and $\mathbf{f}=f_x\hat{x}+f_y\hat{y}$, 
% so
% 
% $$\mathbf{f}\cdot\mathbf{r}=f_xx+f_yy.$$
% 
% For free-space propagation to the far field, the field is proportional to 
% the Fourier transform of the aperture field:
% 
% $$E(x,y,z)= \int\!\!\int E_\mathrm{aperture}(u,v) \exp\left[-i2\pi\left(xu 
% +yv\right)\right]\,du\,dv.$$
% 
% where $x$ and  $y$ are in units of $\frac{\lambda}{NA}$, $NA$ is the numerical 
% aperture of the system, and $u$ & $v$ are in pupil coordinates.  
% 
% This is what we will simulate.  The last bit of information is how the discrete 
% FFT relates to the continuous Fourier transform in units.  The basic relation 
% is $\Delta x \Delta u = \frac{1}{N}$ ,  where $N$ is the number of points in 
% the window.   Alternately,  $\Delta x  = \frac{1}{N\Delta u} = \frac{1}{W}$ 
% , where $W$ is the window size.   
%% Step one: Simulating the point spread function of a perfect optical system
%% 
% # Define the number of points.
% # Establish a range at least 10 times the diameter of the pupil (the pupil 
% is a unit circle).
% # Create a 2d meshgrid over that range.
% # Use a logical variable to identify points inside the pupil (with radius 
% < 1).
% # Set points outside the pupil to 0, inside to 1.
% # Calculate the Fourier transform of the pupil.  (Amplitude squared = Irradiance)
% # Define the units for the axes.
% # Create an image to visualize the focus spot of the perfect system.  \


Npts = 1024 ;  % points in the sampling of the pupil
W = 20  ;  % window size
range = linspace(-W/2,W/2,Npts) ; % The pupil diameter is 1
% Create a meshgrid
[u,v] = meshgrid(range);
%% Set the pupil field Puv to zero
Puv = zeros(size(u))  ; 
%% Define a radius and set points inside the pupil  1 
r2 = sqrt(u.^2+v.^2)   ;  % square of the radius
lvp = r2<=1   ;  % points inside the pupil 
% 
Puv(lvp) = 1; 

figure
imagesc(range,range,Puv);
axis image
colormap("pink")
xlim([-W/5 W/5])
ylim([-W/5 W/5])
title('Pupil Image')
%% 
% Your pupil image should look like this:  
% 
% 

% Compute the Fourier transform of the pupil  
Eimage = fftshift(fft2(fftshift(Puv)))     ;
% figure
% Compute the irradiance
II = abs(Eimage).^2         ;
% 
% % Calculate the coordinates for the focal spot
dx = 1/W;
xrange = dx.*(-(Npts-1)/2:Npts/2);
% 
% % display an image of the focal spot (point spread function)
imagesc(xrange,xrange,sqrt(II)); axis image
colormap("pink")
% % zoom in to see the point spread function
xlim([-3 3])
ylim([-3 3])
title('Point Spread Function')
%% 
% Your PSF should look like this: 
% 
% 

%%
% Display two interfering beams
figure
Puv(lvp)=0;
u1 = 20;
v1 = 20;
Puv(Npts/2 + v1,Npts/2 + u1)=1;
Puv(Npts/2 - v1,Npts/2 - u1)=1;
% 
% Now calculate the diffraction pattern (the Fourier Transform)
Eimage = fftshift(fft2(fftshift(Puv)));
II = abs(Eimage).^2;

imagesc(xrange,xrange,II); shading interp
axis image
colormap("pink")
colorbar
% % % zoom in to see the fringes
xlim([-5 5])
ylim([-5 5])
%% 
% Your interference pattern should look something like this:
% 
% 

% now create a random number of interfering beams
NN = 10;
up = 2*rand(NN)-1;
vp = 2*rand(NN)-1;
ruv = sqrt(up.^2+vp.^2);
lp = ruv<=1;
up = up(lp);
vp = vp(lp);
ra = 0.1;
for ij=1:NN
    lap = sqrt((u-up(ij)).^2+(v-vp(ij)).^2) <= ra;
    Puv(lap) = exp(1i*2*pi*rand(1));
   
end

% Visualize the phase in the pupil 
figure
imagesc(range,range,angle(Puv))
shading interp
axis image
% % colormap(jet)
colorbar
xlim([-W/5 W/5])
ylim([-W/5 W/5])


% Find the Fourier transform of the pupil and display the pattern 
Eimage = fftshift(fft2(fftshift(Puv)));
II = abs(Eimage).^2;

figure
imagesc(xrange,xrange,sqrt(II)); shading interp
axis image
colormap("pink")
% 
xlim([-10 10])
ylim([-10 10])
% Now display the phase

%% 
%