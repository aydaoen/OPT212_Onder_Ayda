% Solving nonlinear eqns, brute force method
x = linspace(0,2,101);
y = @(x) 2.*x.^2 -pi.*x + 1;

figure
plot(x,y(x), '-r', 'LineWidth',2);

hold on
plot([-10 10], [0 0], '-b');

ff = abs(y(x));
plot(x, ff, '-m', 'LineWidth',2)
xmin = x(ff == min(ff))

y(xmin)
plot (xmin, zeros(size(xmin)),'x');

%% another approach

yy = y(x);
signtest = yy(1:end-1).*yy(2:end);
lz = signtest <=0;
xz = x(lz);

y(xz)

figure
plot (x, y(x), '-r'); hold on
plot (xz, zeros(size(xz)), 'x');