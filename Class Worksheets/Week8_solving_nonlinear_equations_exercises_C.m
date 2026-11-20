%% OPT 212: Solving Nonlinear Equations: Newton's method
% Solutions
% 
% These exercises accompany the Solving Nonlinear Equations active script. They 
% focus on bracketing, bisection, Newton iteration, and ray intersection with 
% a curved surface.
% 
% 

clear
close all
clc
format compact
set(groot,'DefaultAxesFontSize',14,'DefaultTextFontSize',14, ...
    'DefaultLegendFontSize',14,'DefaultColorbarFontSize',14);
%% Example -- simple brute force search for a polynomial
% For the brute force search, we simply sample the search region and look for 
% zero crossings

format long e
x = linspace(0,2,10003);
dx = x(2) - x(1); %distance between values of x
y = @(x) 2.*x.^3 + x.^2 - pi.*x - 2; 
yy = y(x);
signtest = yy(1:end-1).*yy(2:end);
lz = signtest<=0; 
xz = x(lz) + dx/2;
y(xz);
figure
plot(x,y(x),'-r','LineWidth',1.2); hold on
plot(xz,zeros(size(xz)),'x','MarkerSize',12)
plot([min(x) max(x)],[0 0],'--b')
hold off
ylabel('function value')
xlabel('x')
legend('function','roots')

%% 
%% 
%% Exercise 1: Brute force search for a Fringe Level
% A normalized two-beam interference signal can be written
% 
% $$I(z)={1\over 2}\left[1+\cos\left({2\pi z\over p}+\phi\right)\right].$$
% 
% Tasks:
% 
% 1. Use a fringe period |p = 3.0 mm|, phase |0.40 rad|, and threshold irradiance 
% |0.80| in normalized arbitrary units. 
% 
% 2. Plot the function.
% 
% 3. Locate the points at which the Irradiance crosses the threshold.
% 
% Your result should look something like this:
% 
% 

% Insert code here
% Define parameters
period_mm = 3;
phase = 0.4;
threshold = 0.8;
fringe = @(z) (1/2)*(1+cos((2*pi*z./period_mm)+phase));
f = @(z) fringe(z) - threshold;
% Choose sampling points
z = 1;
% Sample the function
ff = f(z) 
% Perform logical test for threshold
signtest = ff(1:end-1).*ff(2:end);
lz = signtest<=0;
z0 = z(lz) + dx/2;  % identify the zeros

figure
plot(z,fringe(z),'b-','LineWidth',1.4)
hold on
yline(threshold,'r--','LineWidth',1.2)
plot(z0,fringe(z0),'ko','MarkerFaceColor','k')
hold off
grid on
xlabel('z (mm)')
ylabel('normalized irradiance')
title('Brute force search for fringe thresholds');
%% 
% 
% 
% *Newton's method* uses a local Taylor expansion. If $x_n$ is the current estimate, 
% then
% 
% $$g(x_n+h)\approx g(x_n)+h\,g'(x_n).$$
% 
% Setting this linear approximation equal to zero gives
% 
% $$h=-\frac{g(x_n)}{g'(x_n)},$$
% 
% and therefore the iteration
% 
% $$x_{n+1}=x_n-\frac{g(x_n)}{g'(x_n)}.$$
% 
% Newton's method can converge very rapidly near a well-behaved root, but the 
% notes emphasize that it is not automatically robust. A poor starting point or 
% a misleading local slope can send the iteration far from the intended solution. 
% A practical strategy is to use a brute-force search to obtain starting estimates, 
% then refine them with Newton's method.
%% Exercise 2: Newton's Method for the Grating Equation
% A transmission grating obeys
% 
% $$\Lambda(\sin\theta_m-\sin\theta_i)=m\lambda,$$
% 
% where $\Lambda$represents the grating period, $\theta_i$ represents the angle 
% of incidence, $\theta_m$ the angle of diffraction, $\lambda$ the wavelength 
% and $m$ the 'order' of the diffraction.
% 
% There is, of course, an exact solution.  You can get there by solving the 
% linear equation for $\sin\theta_m$.
% 
% $$\sin\theta_m=\frac{m\lambda}{\Lambda} + \sin\theta_i$$
% 
% and then solving for the arcsin (inverse sin) using the asin function.   
% 
% $$\theta_m=\sin^{-1}\left(\frac{m\lambda}{\Lambda} + \sin\theta_i\right)$$
% 
% This has a built-in function in MatLab but we will use Newton's method and 
% check our accuracy.  To do so, we first define it so we are finding the roots: 
% 
% $$\frac{m\lambda}{\Lambda} + \sin\theta_i - \sin\theta_m= 0$$
% 
% Use wavelength |633 nm|, grating density |600 lines/mm|, incident angle |10 
% deg|, and order |m = 1|. Convert the wavelength from nm to mm before combining 
% it with the groove spacing.  Pay attention to the units of angle:  If you use 
% an angle of degrees you will need to use the sind function in MatLab.  If you 
% use radians, you can use the sin function.   
% 
% Tasks:
% 
% 1. Convert grooves per mm to a period $\Lambda$.  
% 
% 2. Choose a starting value based on the small angle approximation.    
% 
% 3. Use Newton's method to solve for the angle of the first diffraction order. 
% 
% 4. Plot the residual magnitude versus iteration number.
% 
% The result should look like this:
% 
% 
% 
% 

% Put code in here
wl_nm = 633;
wl_mm = wl_nm*1e-6;
grooves_mm = 600;
Lambda = 1/grooves_mm;
theta_i = deg2rad(10);
%order
m = 1;

%parameter for dx calc
h = 1e-6;
g = @(theta) m*wl_mm/Lambda + sin(theta_i) - sin(theta);
dg = @(theta) (g(theta+h/2) - g(theta-h/2));

%starting value for theta
theta(1) = m*wl_mm/Lambda + sin(theta_i);
%exact value for theta
theta0 = asin(m*wl_mm/Lambda + sin(theta_i));

abs_error(1) = abs(theta(1)-theta0);
for ij=1:5
    theta(ij+1)=theta(ij) - g(theta(ij))/dg(theta(ij));
    change(ij+1) = theta(ij+1) - theta(ij);
    abs_error(ij+1) = abs(theta(ij+1)-theta0);
end

disp([theta(:) change(:) abs_error(:)])

%%