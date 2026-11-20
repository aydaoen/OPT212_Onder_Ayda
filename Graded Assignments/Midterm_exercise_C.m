%% OPT 212 midterm exercise:  Solutions of the slab waveguide equation
% A symmetric slab waveguide has the following geometry:  
% 
% 
% 
% For a transverse electric (TE) mode, the field points along y, out of the 
% page and parallel to the interfaces. The core occupies -d/2 < x < d/2 and the 
% propagation is along the z direction.  The effective index of a mode is a solution 
% to the following equation that has no closed form and must be solved numerically:  
% 
% $$g(n_{\rm eff})=\frac{(\beta^2-\alpha^2)\sin(\beta d)-2\beta\alpha\cos(\beta 
% d)}{\beta^2+\alpha^2}=0.$$
% 
% where  $\beta=\sqrt{n_{\rm core}^2-n_{\rm eff}^2},\qquad \alpha=\sqrt{n_{\rm 
% eff}^2-n_{\rm clad}^2}$ .
% 
% In this equation,  the thickness d is in units of wavelengths.  
% 
% Guided modes are those with an effective index between the cladding index 
% and the core index.   Do the following:
% 
% 1. Create a plot of the equation $g(n_{\rm eff})$ with $n_{\rm eff}$ ranging 
% from 1.45 to 1.5.  
% 
% 2. Draw a horizontal dashed line to mark where the function crosses zero.  
% 
% 3. Use the brute force method similar to what was done in the class example 
% to identify and mark where the zeros of the function are.  
% 
% 4. Write code to use those zeros as a starting value for Newton's method to 
% find the exact effective index of the mode.   
% 
% 

% constants and variable definitions
n_clad = 1.45;
n_core = 1.5;
d = 5; 

% beta and alpha written as a function of n_eff
Beta = @(n_eff) sqrt(n_core.^2 - n_eff.^2);
Alpha = @(n_eff) sqrt(n_eff.^2 - n_clad.^2);

% g writen as a function of n_eff
g = @(n_eff) ((Beta(n_eff).^2 - Alpha(n_eff).^2).*sin(Beta(n_eff).*d) - 2.*Beta(n_eff).*Alpha(n_eff).*cos(Beta(n_eff).*d))./(Beta(n_eff).^2+Alpha(n_eff).^2);

% sampling points 
x = linspace(1.45,1.5,10003);
dx = x(2) - x(1); % distance between two adjacent values of x
% sampling the function
gg = g(x);
% logical test to see where the function crossed zero
signtest = gg(1:end-1).*gg(2:end);
lz = signtest<=0; % if two signs are opposite it measn it crossed zero + the product will be negative
xz = x(lz) + dx/2;

% plotting
g(xz);
figure
plot(x,g(x),'-r','LineWidth',1.2); hold on
plot(xz,g(xz),'x','MarkerSize',12);
plot([min(x) max(x)],[0 0],'--b');
hold off
ylabel('g(n_(eff)) value')
xlabel('n_(eff)')
legend('function','roots')

% newton's method
% parameter for dx calc
h = 1e-6;
dg = @(n_eff) (g(n_eff+h/2) - g(n_eff-h/2))/h;

% starting value for n_eff
n_eff = zeros(6, length(xz)); % 6 rows, 2 columns
n_eff(1,:) = xz; % first row is equal to the two xz values found from brute force

% for the loop, i wanted it so that R1 would be equal to the values of xz i got from brute force, and the
% remaining 5 rows would be populated by the loop for newton

% loop
for ij=1:5
    n_eff(ij+1,:) = n_eff(ij,:) - g(n_eff(ij,:))./dg(n_eff(ij,:)); %newtons eqn
    change(ij+1,:) = n_eff(ij+1,:) - n_eff(ij,:); %filling the matrix
end

format long % new code learned!
disp(n_eff)