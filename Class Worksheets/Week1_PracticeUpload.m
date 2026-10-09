%% 
% Ayda Önder
% 
% This is an example of text, including equations. 
% 
% $$B=x+\phi_o$$
% 
% $$Y=sin(B)$$
% 
% $$I(\lambda)=I_0exp((\lambda-\lambda_o)/\Delta\lambda)^2)$$
% 
% $$I(\lambda)=I_0exp\left(\frac{(\lambda-\lambda_0)^2}{\Delta\lambda^2}\right)$$

x = linspace(0,2*pi,400)';
phase = linspace(0,pi,5);
Y = sin(x + phase);
plot(x,Y,'LineWidth',1.3)
%% 
% 
% 
% 
% 
%