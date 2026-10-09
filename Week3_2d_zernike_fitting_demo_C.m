%[text] # OPT 212: 2d polynomial fitting class activity
%[text] This exercise is aimed at practicing 2d fitting using Zernike wavefront fitting as the example. 
clear
close all
clc
format compact
set(groot,'DefaultAxesFontSize',14,'DefaultTextFontSize',14, ...
    'DefaultLegendFontSize',11,'DefaultColorbarFontSize',14);
%[text] ## Exercise : Zernike Fitting on a Circular Pupil
%[text] A sampled pupil phase map can be fitted with a Zernike expansion. Use a square grid with `129` points in each direction, coordinates from `-1` to `+1`, and the circular pupil condition
%[text]{"align":"center"} $\\rho=\\sqrt{\\rho\_x^2+\\rho\_y^2}\\le 1.$
%[text] Use the following u Wyant-ordered basis:
%[text]{"align":"center"} $Z\_0=1,$ $Z\_1=\\rho\\cos\\theta=\\rho\_x,$ $Z\_2=\\rho\\sin\\theta=\\rho\_y,$ $Z\_3=2\\rho^2-1,$ $Z\_4=\\rho^2\\cos(2\\theta)=\\rho\_x^2-\\rho\_y^2,$ $Z\_5=\\rho^2\\sin(2\\theta)=2\\rho\_x\\rho\_y,$ $Z\_6=(3\\rho^3-2\\rho)\\cos\\theta=(3\\rho^2-2)\\rho\_x,$ $Z\_7=(3\\rho^3-2\\rho)\\sin\\theta=(3\\rho^2-2)\\rho\_y,$ $Z\_8=6\\rho^4-6\\rho^2+1,$ $Z\_9=\\rho^3\\cos(3\\theta)=\\rho\_x^3-3\\rho\_x\\rho\_y^2.$
%[text] 1\. Build the pupil grid and mask. 
%[text] 2\. Assemble a Zernike matrix with one row for each valid pupil sample and one column for each Zernike mode. 
%[text] 3\. Choose an example wavefront and find the Zernike coefficients that match the wavefront
%[text] $W\_{xy} = 4x + 5x^2 - y^2$
%[text] 4\. Invert the Zernike matrix to find the coefficients that match the wavefront.  
% Define a 2d meshgrid for x and y.  
[x,y] = meshgrid(linspace(-1,1,201));

% Turn x & y into column vectors
x = x(:);y = y(:);


% choose only values of x&y within the unit circle 
ro2 = x.^2 + y.^2;
lvr = ro2<1;
x = x(lvr);
y = y(lvr);
figure %[output:690eb36a]
plot(x,y,'or') %[output:690eb36a]
axis equal %[output:690eb36a]

% Define an aberration function with coefficients of waves
wxy = 4.*x + 5.*x.^2 -y.*2;

% Define the Zernike matrix
ZZ = Zernike(x,y);
size(ZZ) %[output:0bb56482]

% Invert to find the coefficients.  
coeffs = ZZ\wxy %[output:6debc09b]



%[text] 
%[text] 
%[text] 
%[text] 
%[text] 
%[text] 
%[text] 
% Defines the Zernike Polynomials for our calculations 
function Zout=Zernike(x,y)
ro2=x.^2+y.^2; 
ro2(ro2>1.00001)=nan;
ro=sqrt(ro2);
ro3=ro.^3;
ro4=ro.^4;
ro5=ro.^5;
ro6=ro.^6;
ro7=ro.^7;
ro8=ro.^8;
ro9=ro.^9; 
ro10 = ro.^10; 

theta=atan2(y,x);

theta=atan2(y,x);
Zer(:,1)=ro.*cos(theta);%x tilt
Zer(:,2)=ro.*sin(theta);
Zer(:,3)=-1+2.*ro2;
Zer(:,4)=ro2.*cos(2.*theta);
Zer(:,5)=ro2.*sin(2.*theta);
Zer(:,6)=ro.*(-2+3.*ro2).*cos(theta);
Zer(:,7)=ro.*(-2+3.*ro2).*sin(theta);
Zer(:,8)=1-6.*ro2+6.*ro4;
Zer(:,9)=ro3.*cos(3.*theta);
Zer(:,10)=ro3.*sin(3.*theta);
Zer(:,11)=ro2.*(-3+4.*ro2).*cos(2.*theta);
Zer(:,12)=ro2.*(-3+4.*ro2).*sin(2.*theta);
Zer(:,13)=ro.*(3-12.*ro2+10.*ro4).*cos(theta);
Zer(:,14)=ro.*(3-12.*ro2+10.*ro4).*sin(theta);
Zer(:,15)=-1+12.*ro2-30.*ro4+20.*ro6;
Zer(:,16)=ro4.*cos(4.*theta);
Zer(:,17)=ro4.*sin(4.*theta);
Zer(:,18)=ro3.*(-4+5.*ro2).*cos(3.*theta);

Zer(:,19) = (5 * ro5 - 4 * ro3) .* sin(3*theta);
Zer(:,20) = (15 * ro6 - 20 * ro4 + 6 * ro2) .* cos(2*theta);
Zer(:,21) = (15 * ro6 - 20 * ro4 + 6 * ro2) .* sin(2*theta);
% 
Zer(:,22) = (35 * ro7 - 60 * ro5 + 30 * ro3 - 4 * ro) .* cos(theta);
Zer(:,23) = (35 * ro7 - 60 * ro5 + 30 * ro3 - 4 * ro) .* sin(theta);
Zer(:,24) = (70 * ro.^8 - 140 * ro.^6 + 90 * ro.^4 - 20 * ro.^2 + 1);
Zer(:,25) = ro5 .* cos(5*theta);
Zer(:,26) = ro5 .* sin(5*theta);
Zer(:,27) = (6 * ro6 - 5 * ro4) .* cos(4*theta);
Zer(:,28) = (6 * ro6 - 5 * ro4) .* sin(4*theta);
Zer(:,29) = (21 * ro7 - 30 * ro5 + 10 * ro3) .* cos(3*theta);
Zer(:,30) = (21 * ro7 - 30 * ro5 + 10 * ro3) .* sin(3*theta);

Zer(:,31) = (56 * ro8 - 105 * ro6 + 60 * ro4 - 10 * ro2) .* cos(2*theta);
Zer(:,32) = (56 * ro8 - 105 * ro6 + 60 * ro4 - 10 * ro2) .* sin(2*theta);
Zer(:,33) = ro.*(5-60*ro2 + 210*ro4 - 280*ro6 + 126*ro8).* cos(theta);
Zer(:,34) = ro.*(5-60*ro2 + 210*ro4 - 280*ro6 + 126*ro8).* sin(theta);
Zer(:,35) = -1 + 30*ro2 - 210*ro4 + 560*ro6 -630*ro8 + 252*ro10;  
Zout = [ones(size(Zer(:,1))) Zer]; 
end
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:690eb36a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAgsAAAE7CAYAAABTzcL+AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnXuMVdXZ\/x8JAwFUKEQjclGMgmgj0qRgQghGbBMlobxvE\/NKDSEBo4KQcmlFQU0ICjVcRKpgBUobxKiJL9BW+KvBTIiCBjTEGybQ4iCN9qcDLzrCDMMva9ejw8w+5+zL2nuv9ezPSQgo6\/r5rj3ry7Oevc5F\/fv3Py98IAABCEAAAhCAQBUCF2EWWBsQgAAEIAABCNQigFlgfUAAAhCAAAQgUJMAZoEFAgEIQAACEIAAZoE1AAEIQAACEIBAcgJEFpKzoyYEIAABCECgFAQwC6WQmUlCAAIQgAAEkhPALCRnR00IQAACEIBAKQhgFkohM5OEAAQgAAEIJCeAWUjOjpoQgAAEIACBUhDALJRCZiYJAQhAAAIQSE4As5CcHTUhAAEIQAACpSCAWSiFzEwSAhCAAAQgkJwAZiE5O2pCAAIQgAAESkEAs1AKmZkkBCAAAQhAIDkBzEJydtSEAAQgAAEIlIIAZqEUMjNJCEAAAhCAQHICmIXk7KgJAQhAAAIQKAUBzEIpZGaSZSYwtL1dftvSIuNaW8X8uePnWLduF\/y\/zv9tykb5f2nKvNSzp+zt3l32NjSUWSbmDgGnCWAWnJaHwUEgGoEwQxC2gUdrrbhSHcds\/oyRKE4LeoZARwKYBdYDBDwiUCtK4NE0Ug0VE5EKH5UhkIgAZiERNipBIB8CmIP6nDEP9RlRAgJpCWAW0hKkPgQyIGBMwv+cOSMPtbRk0LreJs2xxYN9+uidIDODQEEEMAsh4AcOHCjmFx8I5Emg+\/HjMu3IkdBExDzH4WNfYYmaxjgcGTJEjgwd6uOUGLPHBE6cOCHml6YPZqGTmsYkLFmyRH7yk59o0pm5OEig4fhxGbBunfTav1\/Mn\/lkR6B10CA59V\/\/Jd+MHSstY8Zk1xEtQ0BEDhw4IMuWLVNlGDALnZa2MQnPPvusOqFtPcGjR4+WmTNnwich0CHnzgVHC9d8+ikGISHDtNWMcdhz9dW8rpkQJD8DaoOr8Jk9e3ZgGrR8MAtVzII2oW0t2IqZgk90oiQpRmeVd0mSI+MT52dAbWZa+WAWMAuxflqYY5o777xTXn\/9dVUhtlgQIhYmSTEiKEeKkRwZTQh+BmAWoq0U5aW0ukLlsjk3vUo04e4zZ5wbGwOqTsDcJPm7Xr24TZJFkpiA1j2EyAKRhcQPBRUvJMBxg54VwfGEHi3znglmIW\/iBfWnVeiCcJaiW6IIumU2xxMv9ehBtEG3zNZmp3UPIbJAZMHaQ1K2hshJKJfi5DSUS++ks8UsJCXnWT2tQnsmg\/PDJZrgvESZDJCchkywqmpU6x5CZIHIgqoHNevJYBKyJuxH+xxN+KFTEaPELBRBvYA+tQpdAEo1XZK4qEZK6xMhEdI6Uu8b1LqHEFkgsuD9w5nVBMhJyIqsznbJadCpa9xZYRbiEvO0vFahPZWjsGFjFApD73XHJqdh8qWXej0HBp+OgNY9hMgCkYV0T4bC2uQlKBQ1xymRBJkjbAe7wiw4KEoWQ9IqdBastLVJNEGbosXOh2OJYvkX1bvWPYTIApGFop4pp\/olmuCUHGoGQ5RBjZSRJ4JZiIzK74JahfZblexGTzQhO7a0\/AMBogzlWQ1a9xAiC0QWyvMUh8z0ty0t8lBLS6kZMPl8CGAY8uFcdC+YhaIVyKl\/rULnhM+bbogoeCOVqoFiGFTJGToZrXsIkQUiC\/qf3g4z5IKlUsnt5GS5yMlJWawNCrNgDaXbDWkV2m3q2Y+OSEL2jOkhPgEiDfGZuV5D6x5CZIHIguvPXurxYRRSI6SBDAlwkVOGcAtoGrNQAPQiutQqdBEsXemTJEZXlGAc1QgQYdCzNrTuIUQWiCzoeUpDZmKiCgebm1XPkcnpIDD5kktkb0ODjsmUeBaYhZKIr1Xoksj3\/TS5ZKlsivs\/XxIf\/dfQzEDrHkJkgciCjif0u1mQn6BKztJOhmMJf6XHLPirXayRaxU6FgSPC5Of4LF4DP0CAhgGPxeE1j2EyAKRBT+fSPIT1OjGRKoTII\/Bv9WBWfBPs0Qj1ip0IhieVOLowROhGGYiAkQYEmErrJLWPYTIApGFwh4qGx1jFGxQpA3XCXAXg+sK\/TA+zII\/WqUaqVahU0FxtDJGwVFhGFYmBDAMmWC13qjWPYTIApEF6w9LXg2SzJgXafpxhQBHEq4oUX0cmAX3NbIyQq1CW4HjUCNctuSQGAwlVwIkPeaKO3ZnWvcQIgtEFmI\/DEVW4OihSPr07QoBIgyuKNF1HJgFd7WxOjKtQluFVFBjGIWCwNOtkwTIYXBSFm5wdFMW+6PCLNhnaqtFchRskaQdLQSIMLinpNY9hGMIjiHce9pCRkSOghcyMcgCCJDDUAD0Gl1iFtzSI7PRaBU6M2AZN8zRQ8aAaV4FASIM7siodQ8hskBkwZ2nLGQkHD04LQ+Dc4gAhsENMTALbuiQ+Si0Cp05uAw64OghA6g0qZoARxLFy6t1DyGyQGSh+Kerygh2njol49ranB0fA4OAawSILhSvCGaheA1yGYFWoXOBZ7ETogoWYdJUqQgQXShWbq17CJEFIgvFPlkhvRujsOPUKTG\/84EABOIRMPcv\/K5XL9nb0BCvIqWtEMAsWMHofiNahXafvATm4H\/OnJGHWlp8GC5jhIDTBDiSKEYerXsIkQUiC8U8USG9mjcf7j5zhoiCM4owEN8JYBjyVxCzkD\/zQnrUKnQhMGN0WslRMCFUkhpjgKMoBKoQMM\/SkPZ2ebBPH44kclwlWvcQIgtEFnJ8jKp3Zd58MD\/YyFNwQg4GoYgA0YV8xcQs5Mu7sN60Cl0Y0Agd8+ZDBEgUgUAKAgP6909Rm6pxCGjdQ4gsEFmI8xxYLUtCo1WcNAaBqgT4hsr8FgdmIT\/WhfakVehCoVbpnIRGF1VhTFoJcByRj7Ja9xAiC0QW8nmCOvVCQmMh2Om05AS4sCn7BYBZyJ6xEz1oFdoJuB0GQUKja4ownjIQ4Dgie5W17iFEFogsZP\/0hEQVuKExd+x0CIGAAIYh24WAWbDId+DAgTJr1iwZO3asXHzxxUHLzc3N0tjYKFu2bJETJ07E7s20d88998hFF11Ute7p06dlyZIlsm\/fvqpltAodG2iGFfja6Qzh0jQEIhAgfyECpIRFtO4huUcWhg8fLk888YQMHjw4VIqmpiZZvHixHD58OJZUK1askAkTJtSsg1mIhTSTwrwmmQlWGoVAbALkL8RGFqkCZiESptqFunfvLuvWrZObb75Zzp49K7t27ZINGzYElebOnSsTJ06UHj16yP79+2XBggXSFvHriU27mzdvluuuu07efffdIDoR9mltbZWPPvpIvvnmGyILFvRM0gRfO52EGnUgYJ8A0QX7TE2LmAULXCdPnizz588Xs7lv375dVq5ceUGrCxculClTpsi3334b\/N3u3bsj9XrNNdcE5c3xxssvvyxPP\/10pHphhbQKnRiIpYrcqWAJJM1AwCIBDINFmN81pXUPyfUYonJU8MUXX4gxBp2PGswRhdn0L7vssiDqsHTp0khK3nbbbfLwww8HUYlnn31WXnnllUj1MAuJMcWuSJ5CbGRUgEAuBDAMdjFjFlLyvOKKK2Tt2rUydOhQee+99+T+++8PbdEcS4waNUqOHj0qDzzwgJw8ebJuz9OmTZMZM2ZIS0tLkMD4zjvv1K1TrYBWoRMDsVCRPAULEGkCAhkSIH\/BHlyte0hukYWRI0eKiSxcfvnlNaMGy5cvl1tvvTV4I8JEH44cOVJXxccee0zuuOMOOXbsWNC2+bM5kmhoaJCvv\/5a3nzzTXnuuecivWWhVei6EDMsQJ5ChnBpGgIWCBBdsACRYwg7EM1rksuWLQteldy2bVuQ6Bj2mTNnjkydOlWivLlQqV+JRpw7dy54dbJbt25dmv7qq6\/kqaeekj179tScEGbBjt4dWznY3My3SdrHSosQsEbgWLduMrpfP2vtlbkhrXtIbpGFrMzCkCFDZM2aNTJo0CBpb2+XQ4cOycaNG4OjCJMDYY4nbrnlliCf4fjx4\/LII4\/UfC2zIrRp4+DBg8GaN1GOJHc\/lPmBqcydIwhWAQT8IMBRRHKdTCTb\/DIf87s5Dp89e7YcOHAgeaOO1czNLEyaNCl4HbJXr15WIwsmv8EcQwwYMCC4bMnc0dD5lUvzBoZ5y8K8hfHaa691eQujoyYVs9Dx\/xnjsGnTJsek82M4HEH4oROjhABHEcnXgPlH6cyZMy9oALOQkGdWkYUowzFO75lnngkugqqXOFkxC+bIpBJNILIQhXLXMkQVknGjFgSKIkB0IRn5jpGF0aNHB8YBs5CMZXC1c1Y5C1GGZF7JHDduXN3ESa3nTVEY2S5DVME2UdqDQLYEiC6k56t1D8ntGKLj2xA7duwI3owI+yR5GyKKvJU7Hj7\/\/HNZtGiRfPjhh6HVtAodhZHtMiQ22iZKexDIlgCJjun5at1DcjMLffv2lfXr18uwYcPk7bffDq53DvskuWehX79+ctVVV8knn3xS9SrnSrsmyXHevHny6aefYhbSPxehLXBbY0ZgaRYCORAgupAOMmYhHb+gduVf95999llwh4LJH+j4iXoXQ8c6lVctzXdNVLu9sePNkLUuhDLtahXagnyRm+C2xsioKAgBJwlgGJLLonUPyS2yYNBH\/W4I8zbD6tWrZefOnXUVM1c9m9ch+\/TpE3yJlDEPHd+GMG9AmG+5HD9+fPDlVS+88IK8+OKLVdvVKnRdkJYKVJIaTTjT\/JkPBCDgHwHz\/D7Yp4\/sbWjwb\/AFj1jrHpKrWTBHEatWrZIbb7yx5rdOhm36lQiCWQcdL3UyZsC0OWbMGDl\/\/ry8\/\/778vzzzwf3LNx0001y7733BtECc1FTWLud15VWofN6fkxS45D2doxCXsDpBwIZESC6kAys1j0kV7Ng0JsjgSeffDK4RCns09TUFNyV0PlLpqqZBdPG1VdfLY8\/\/riMGDEiuMGx88eYiA8++EAeffTRupcraRU62bKPX4ukxvjMqAEBFwmQ7JhMFa17SO5mweA376TOmjUreJ3SXP9sPs3NzdLY2ChbtmwJ3dBrmQVTv3fv3jJ9+nS5\/fbbgwuazI2Nra2t8uWXXwZfh71169YulzWFLQWtQidb9vFqca9CPF6UhoDrBLh3Ib5CWveQQsxCfPz51dAqdB4EuVchD8r0AYH8CHAUEZ+11j0Es9BpLWgVOv6Sj1eDqEI8XpSGgC8EiC7EU0rrHoJZwCzEexKqlCaqYAUjjUDAOQJEF+JJglmIx8vb0lqFzlIQE1XYceoUb0BkCZm2IVAggQH9+xfYu19da91DiCwQWUj9JHIJU2qENAABpwkQXYguD2YhOiuvS2oVOitRuIQpK7K0CwF3CHBJU3QttO4hRBaILER\/CkJKcglTKnxUhoA3BIguRJMKsxCNk\/eltAqdlTBcwpQVWdqFgFsEuKQpmh5a9xAiC0QWoj0BIaV4XTIxOipCwEsCvEZZXzbMQn1GKkpoFToLcXhdMguqtAkBdwlwFFFfG617CJEFIgv1V3+VEhxBJEZHRQh4SYCjiPqyYRbqM1JRQqvQWYjz\/778MotmaRMCEHCYAHcu1BZH6x5CZIHIQqIfS+QrJMJGJQh4T4C8BcyC94vYxgS0ukIbbDq2Qb6CbaK0BwE\/CJC3gFnwY6VmPErMQjTA5CtE40QpCGgjQN4CZkHbmk40H8xCNGzkK0TjRCkIaCRA3kJ1VbXuIeQskLOQ6GcZZiERNipBQAUBzAJmQcVCTjMJra4wDZPOdUlutEmTtiDgHwGSHDEL\/q1ayyPGLNQHSnJjfUaUgIBmAiQ5YhY0r+9Ic8Ms1MZkogo7Tp0S8zsfCECgvAQ4igjXXuseQs5CJ721Cm3rR9pvW1rkoZYWW83RDgQg4CkBoguYBU+Xrp1hYxaqcyRXwc4aoxUIaCFA7kJXJbXuIUQWiCxE\/rlFrkJkVBSEQCkIEF3ALJRioYdNUqsrtCFo5SImcykLOQs2iNIGBPwlYH4OmM\/ofv38nUQGI9e6hxBZILIQ+XExdytgFCLjoiAESkGARMcLZcYslGLZi2gV2oZ8XMRkgyJtQEAXAcwCZkHXio44G8xCdVCYhYiLiGIQKBEBzAJmoUTL\/YepYhbCZedNiFI+DkwaAnUJ8EYEZqHuItFYALMQripvQmhc7cwJAukJ8EYEZiH9KvKwBcxCuGh8JbWHi5khQyAHAnxlNWYhh2XmXheYhXBNyFdwb60yIgi4QoC8Bf1H2bw62elpwyxgFlz5Acw4IOALAcwCZsGXtWptnJgFjiGsLSYagkAJCHAMwTFECZZ51yliFogslHLhM2kIpCBAZIHIQorl42dVzEJX3Xht0s+1zKghkBcBXp\/ELOS11pzpB7PQVQpem3RmeTIQCDhJgNcnMQtOLswsB4VZ6EqX1yazXHG0DQH\/CZC3gFnwfxXHnAFmoSswXpuMuYgoDoESEiBv4T+ia91DeHWy00OtVeg0P7swC2noURcC5SCAWcAslGOlfzdLzAKRhVIteCYLAUsEMAuYBUtLyY9mMAuYBT9WKqOEgFsEMAuYBbdWZMajwSx0BUyCY8aLjuYh4DkBEhx\/EFDrHkLOAjkLdX9MkbNQFxEFIFB6AkQWiCyU6iHQ6grTiIhZSEOPuhAoBwHMAmahHCv9u1liFshZKNWCZ7IQsEQAs4BZsLSUsm9m\/PjxMn36dLn22mulR48e0traKidOnJDt27fLq6++Km1tbXUHgVkgZ6HuIqEABCBwAQFyFshZ8OaRuOuuu+S+++6T3r17dxlze3u7NDY2ypIlS+oaBswCkQVvFj0DhYBDBIgsEFlwaDmGD8Vs8MuWLZMf\/ehH8u9\/\/1s2b94sf\/nLX+SGG26Q+fPny\/Dhw+Xs2bOyZcuW4FetD2YBs+D8gmeAEHCQAGYBs+DgsrxwSEuXLpWf\/exncvLkSVmxYoXs2bPn+wL9+vWTNWvWyPXXXy\/\/+Mc\/ZN68efKvf\/2r6pwwC5gF5xc8A4SAgwQwC5gFB5flD0MaNmyYrFy5Uq688kp58803g0hC58+vfvUruffee8UcR6xatUr+9re\/YRZiqMrbEDFgURQCJSWAWcAsOL30J02aJAsWLAgSGv\/0pz\/JCy+80GW8o0aNkuXLlwfHFDt27AiiD9U+RBa6kuFSJqcfAQYHgcIJkOD4gwRa9xDvL2WaNm2azJgxQ86dO1c1ajBkyJDgKGLQoEGyd+9eWbhwIWYhxo8XIgsxYFEUAiUlQGSByILTS3\/OnDkydepUOX36dPC2w759+0LHaxIbR4wYIR9\/\/HHweiWRheiyYhais6IkBMpKALOAWXB67WMWspcHs5A9Y3qAgO8EMAuYBafX8GOPPSZ33HGH9cjCxo0b5eDBg8HczcVO5ldZP+QslFV55g2BaATKnrMwcOBAMb\/Mx\/xuotyzZ8+WAwcORAPoQSnvcxayiix01M4Yh02bNnkgZzZDJLKQDVdahYAmAmWOLJi8uZkzZ14gJ2bBsdWdlVkwlzxVoglljyxgFhxb9AwHAg4SKLNZ6BhZGD16dGAcMAuOLdLK2xDmhkbzeuTf\/\/73LiPkbYh0omEW0vGjNgTKQKDMZqGjvrw66ehqv+222+Thhx8O7lkwRwV\/\/vOfu4yUexbSiYdZSMeP2hAoAwHMwn9Uxiw4uto73uD4xhtvyKJFi7qMNMpdDJVKWoVOIx8JjmnoURcC+gmUPcGRyIInazzqd0M0NTXJ3Llza77ZgFnoKjqRBU8eBIYJgQIJEFkgslDg8ovW9bhx42Tx4sU1v3Wyra1Ntm7dKn\/4wx9qNopZwCxEW3WUggAEOhLALGAWvHgi7r77brnvvvukZ8+eXcZrvkCqsbExePfVmIZaH8wCZsGLBc8gIeAYAcwCZsGxJVl9OOPHjw+ucr722muDhMfW1tbgyGH79u3y6quv1jUKpmXMAmbBmwXPQCHgEAHMAmbBoeWY\/VAwC10Zk+CY\/bqjBwj4TIAExx\/U07qHeH+Do+0HTKvQaTiR4JiGHnUhUA4CRBaILJRjpX83S8wCxxClWvBMFgKWCGAWMAuWlpIfzWAWMAt+rFRGCQG3CGAWMAturciMR4NZIGch4yVG8xBQR4CcBXIW1C3qehPCLBBZqLdG+HsIQKArASILRBZK9VxgFjALpVrwTBYClghgFjALlpaSH81gFjALfqxURgkBtwhgFjALbq3IjEeDWcAsZLzEaB4CKglgFjALKhd2tUlhFrqS4VKmUj0CTBYCsQmQ4PgDMq17CJcydXostAod++nvUOH3X38td585k6YJ6kIAAooJ7O3eXSZfeqniGUafmtY9BLOAWaj7FIxrbZWd\/\/d\/dctRAAIQKCeByZdcInsbGso5+ZLsIZiFkgid9inmyue0BKkPAb0EyFfgGELv6q4yM60hpLRCkreQliD1IaCTAPkKF+qqdQ8hskBkIdJPMCILkTBRCAKlJEBkgchC6Ra+VleYVkjMQlqC1IeAXgKYBcyC3tXNMUQsbTmGiIWLwhAoDQGOITiGKM1i7zhRIgvhsvP6ZCkfByYNgboEeG0Ss1B3kWgsgFkIV5XXJzWuduYEgfQEeG0Ss5B+FXnYAmahumjkLXi4oBkyBDImQL4CZiHjJeZm85gFzIKbK5NRQcBNApgFzIKbKzPjUWEWapsFk8w0tL09YxVoHgIQ8IUAZgGz4MtatTpOzEJ1nOaNiMoHw2B12dEYBLwjYP7hYD6j+\/XzbuxZDljrHsKlTJ1WjVahbTwcvBFhgyJtQEAPAd6E6Kql1j0Es4BZiPyTizciIqOiIARKQYA3ITALpVjoYZPU6gptCUp0wRZJ2oGA3wRe6tlTHuzTx+9JZDB6rXsIkQUiC7EfF25zjI2MChBQRYBbG6vLiVlQtdTLJ7RN+Ygu2KRJWxDwjwC5CuXbQ4gsEFmI\/ZOK3IXYyKgAAVUEyFXALKha0EkmozWElIRFrTrc5mibKO1BwB8C3K2AWfBntWY0UsxCNLCYhWicKAUBjQQwC5gFjes61pwwC9FwkeQYjROlIKCNAMmNtRXVuoeQs0DOQqKfZSQ5JsJGJQh4T4DkRsyC94vYxgS0ukIbbDq2QZKjbaK0BwE\/CJDciFnwY6VmPErMQnTA5C1EZ0VJCGghQL4CZkHLWk41D8xCdHzkLURnRUkIaCBAvkJ9FbXuIeQskLNQf\/VXKUHeQmJ0VISAlwTIV6gvG2ahPiMVJbQKnYU45C1kQZU2IeAuAfIV6mujdQ8hskBkof7qr1GCo4hU+KgMAW8IcAQRTSrMQjRO3pfSKnRWwpijCBNhGNrenlUXtAsBCDhAgCOIaCJo3UOILBBZiPYEVClljIIxDOaDYUiFksoQcJaAiSqYr6Pe29Dg7BhdGRhmwRUlMh6HVqGzxEaiY5Z0aRsCxRN4qWfPwCzwqU9A6x5CZIHIQv3VH6EEuQsRIFEEAh4SIFchnmiYhXi8vC2tVeisBSG6kDVh2odAMQTIVYjHXeseQmSByEK8J6FKaV6jtIKRRiDgHAFel4wnCWYhHq+apQcOHCizZs2SsWPHysUXXxyUbW5ulsbGRtmyZYucOHEidm+mvXvuuUcuuuiiqnVPnz4tS5YskX379lUto1Xo2EATVNh56pSMa2tLUJMqEICAiwSIKsRXResekntkYfjw4fLEE0\/I4MGDQ1VoamqSxYsXy+HDh2OptGLFCpkwYULNOpiFWEhjFya6EBsZFSDgNAGiCvHlwSzEZ9alRvfu3WXdunVy8803y9mzZ2XXrl2yYcOGoNzcuXNl4sSJ0qNHD9m\/f78sWLBA2iL+K9W0u3nzZrnuuuvk3XffDaITYZ\/W1lb56KOP5JtvviGyYEHPsCZIdMwILM1CIGcCJDYmA45ZSMbtglqTJ0+W+fPni9nct2\/fLitXrrzg7xcuXChTpkyRb7\/9Nvi73bt3R+r1mmuuCcqb442XX35Znn766Uj1wgppFToxkJgVuaQpJjCKQ8BRAhxBJBNG6x6S6zFE5ajgiy++EGMMOh81mCMKs+lfdtllQdRh6dKlkdS67bbb5OGHHw6iEs8++6y88sorkephFhJjqlqRS5rsM6VFCORNgEuYkhPHLCRnF9S84oorZO3atTJ06FB577335P777w9t0RxLjBo1So4ePSoPPPCAnDx5sm7P06ZNkxkzZkhLS0uQwPjOO+\/UrVOtgFahEwNJUJHXKBNAowoEHCLAJUzJxdC6h+QWWRg5cqSYyMLll19eM2qwfPlyufXWW4M3Ikz04ciRI3VVe+yxx+SOO+6QY8eOBW2bP5sjiYaGBvn666\/lzTfflOeeey7SWxZaha4L0XKB37a0yEMtLZZbpTkIQCBrAr\/r1Uue6tUr627Utq91D8nNLJjXJJctWxa8Krlt27Yg0THsM2fOHJk6dapEeXOhUr8SjTh37lzw6mS3bt26NP3VV1\/JU089JXv27Km5SLUKXcSTSbJjEdTpEwLJCZDUmJxdpabWPcR7szBkyBBZs2aNDBo0SNrb2+XQoUOycePG4CjC5ECY44lbbrklyGc4fvy4PPLIIzVfy9QqdPpHIH4LHEfEZ0YNCBRJgKTG9PS17iG5mYVJkyYFr0P26tXLamTB5DeYY4gBAwYEly2ZOxo6v3Jp3sAwb1mYtzBee+21Lm9hdFweFaGN4Th48GDwV+ZIJMlFUemXnd8tcO+C3\/ox+vIR4F6FZJqbY2\/zy3zM7yZ3bvbs2XLgwIFkDTpYKzezkOUxRD2uRrxnnnkmuAiqXuJkxSx0bNMYh02bNtXrhr8PIcCtjiwLCPhBgKhCcp1MBHvmzJkXNIBZCOE5bNiw4F\/rV155ZZe\/NZcvmdcZ\/\/nPf2aWsxBFYjO+cePG1U2crJgFk19RiSYQWYhCOLwM0YXk7KgJgTwJEFVITrtjZGH06NHl9DCzAAAOl0lEQVSBccAsJDQLJpeg8jbEjh07gj+HfZK8DRFF4sodD59\/\/rksWrRIPvzww9BqWs+bojDKqgyJjlmRpV0I2CFAYqMdjqYVrXtIbscQffv2lfXr14uJQrz99tvB9c5hnyT3LPTr10+uuuoq+eSTT6pe5Vxp1yQ5zps3Tz799FPMgr3no2ZLJDrmBJpuIJCQAEcQCcGFVMMsWGBZ+df9Z599FtyhYPIHOn6i3sXQsU7lVcvKcUfY7Y0db4asdSGUZldoQb7ETXAUkRgdFSGQCwGOIOxhxixYYBn1uyHM2wyrV6+WnTt31u3VXPVsXofs06dP8CVSxjx0fBvCvAFhvuVy\/PjxwZdXvfDCC\/Liiy9WbVer0HVBZlyA6ELGgGkeAgkJcFtjQnBVqmndQ3I7hjBczVHEqlWr5MYbb6z5rZNhm34lgmDa6XipkzEDps0xY8bI+fPn5f3335fnn38+uGfhpptuknvvvTc4QzIXNYW121lvrULbfRyStcatjsm4UQsCWRHgtkb7ZLXuIbmaBSOLORJ48skng0uUwj5NTU3BXQmdv2SqmlkwbVx99dXy+OOPy4gRI4IbHDt\/jIn44IMP5NFHH617X4JWoe0\/EslaJMKQjBu1IGCbAHkKton+pz2te0juZsHANK+ZzJo1S8zdC+b6Z\/Npbm6WxsZG2bJlS+iGXsssmPq9e\/eW6dOny+233x5c0GRubGxtbZUvv\/wy+DrsrVu3drmsKWypaBU6m8cifqvkL8RnRg0IZEGAPIUsqGIWsqHqYKuYhexFIbqQPWN6gEAtAuQpZLc+tO4hhUQWspMpfctahU5Pxm4L3L1glyetQSAqAYxCVFLJymndQzALndaDVqGTLfvsapnogjmSGNrenl0ntAwBCFxAgDyF7BeE1j0Es4BZyP7pCenBGAVjGD7t1k3GtbUVMgY6hUDZCJCnkL3imIXsGTvRg1ahnYDbaRDkLrioCmPSSoDjh3yU1bqHEFkgspDPE1SlF+5eKBQ\/nZeEAEYhP6ExC\/mxLrQnrUIXCrVO5yZvYcepU+QvuCwSY\/OWAF8Sla90WvcQIgtEFvJ9kqr0RsKjEzIwCIUESGrMV1TMQr68C+tNq9CFAY3YMQmPEUFRDAIRCZiIgvk82KeP7G1oiFiLYmkJaN1DiCwQWUj7bFirT8KjNZQ0BIGAALkK+S8EzEL+zAvpUavQhcBM0CkJjwmgUQUCIQT4kqhiloXWPYTIApGFYp6oGr2ahMffnz7N\/QvOKcOAfCBgjh84eihOKcxCcexz7Vmr0LlCtNAZXzhlASJNlJIAFy8VK7vWPYTIApGFYp+sGr3vPHWK6IKz6jAwFwnw5kPxqmAWitcglxFoFToXeJY7IbpgGSjNqSdAVKF4ibXuIUQWiCwU\/3TVGAFvSDgtD4NziABvPrghBmbBDR0yH4VWoTMHl2EHvCGRIVyaVkGANx\/ckVHrHkJkgciCO09ZjZFwJOGFTAyyAALmzQcTVeDjBgHMghs6ZD4KrUJnDi6HDjiSyAEyXXhFgKMH9+TSuocQWSCy4N7TRg6DV5ow2GIIYBSK4V6vV8xCPUJK\/l6r0ErkCaZBDoMmNZlLEgLkKCShlk8drXsIkQUiC\/k8QZZ7IYfBMlCa84YAOQpuS4VZcFsfa6PTKrQ1QA41RA6DQ2IwlFwIcPSQC+ZUnWjdQ4gsEFlI9WAUXRnDULQC9J8XAYxCXqTT9YNZSMfPm9pahfZGgAQDxTAkgEYVrwhgFPyRS+seQmSByII\/T2GNkZL0qEJGJhFCgGRGv5YFZsEvvRKPVqvQiYF4VJGkR4\/EYqiRCJDMGAmTU4W07iFEFogsOPWgpR0MRxJpCVLfFQIcPbiiRLxxYBbi8fK2tFahvRUkwcA5kkgAjSpOEeDowSk5Yg1G6x5CZIHIQqwHwZfCQ9vbg8ub7j5zxpchM86SEzjWrZvsbWiQl3r0CH7n4ycBzIKfusUetVahY4NQUoE8BiVClmAaky+5BJOgQGetewiRBSILCh7P2lMgj0G9xN5PkPwE7yX8fgKYBT1a1pyJVqFLIl\/VaWIYyr4C3J0\/RsFdbZKMTOseQmSByEKS58HLOiQ+eimb6kGTyKhPXsyCPk1DZ6RV6JLIV3eaJD7WRUSBjAmQyJgx4IKb17qHEFkgslDwo1Vc90QaimNf1p6JJOhXHrOgX+NghlqFLol8sadJLkNsZFRISIDchITgPKumdQ8hskBkwbNH0f5wiTDYZ0qLFxIgolCeFYFZKInWWoUuiXyJp2lyGX5\/+rSMa2tL3AYVIdCZwN7u3cUYBS5ZKs\/a0LqHEFkgslCepzjCTIkyRIBEkUgEiCZEwqSuEGZBnaThE9IqdEnkszJNogxWMJa2EaIJpZU+mLjWPYTIApGFcj\/ZNWZP8iNLIy4BkhjjEtNXHrOgT9PQGWkVuiTyWZ8mxxLWkaptkGMHtdLGmpjWPYTIApGFWA9CGQtzkVMZVY82Zy5YisapTKUwCyVRW6vQJZEv82liHDJH7EUH5riBr5L2QqrcB6l1DyGyQGQh94dJQ4ckQWpQMf4cSF6Mz6xsNTALJVFcq9AlkS\/3aZLTkDvywjokJ6Ew9F51rHUPIbJAZMGrB9HFwXI04aIq9sbEkYM9lmVoCbNQBpUVvyNbEvkKnSamoVD8VjsncdEqzlI1hlkoidxahS6JfM5Mk5wGZ6SINRByEmLhonAIAa17CMcQHEPEeuAHDhwod955p7z++uty4sSJWHXLWJicBn9UJychmlb8DKjNCbMQbR15X0qr0LaEgU98khxPxGeWVw2OG+KT5mcAZiH+qlFYo\/IgbNy4UQ4ePKhwhummZP5VsWTJEoFPMo5Dzp2Th1paZMRbbyVrgFqpCbQOGiRHhgwRk7jIt0HGx8nPgNrMKnxmz54tBw4ciA\/Y0RocQ3QSpiK0MQ18IJAlgYbjx2XAunVy6f\/+b5bd0LaIGIPQMmaMnPzv\/w5+5wOBLAkYk7Bs2TJVR7WYhZAVYwyD+cUHAnkS6H78uEw7ckTuPnMmz25V9tXxeKFl7FiVc2RS7hIw+VzacrowC+6uN0ZWcgIkR0ZbAMYYmI\/JDSFJMRozSkEgLgHMQlxilIdAjgRIjqwPmyTF+owoAYG0BDALaQlSHwI5EsA8iGAOclxwdAWB7whgFlgKEFBAQKOJwBQoWJhMQQ0BzIIaKZkIBMIJuGwkMASsWgj4QQCz4IdOjBICEIAABCBQGAHMQmHo6RgCEIAABCDgBwHMgh86MUoIQAACEIBAYQQwC4Whp2MIQAACEICAHwQwC37oxCghAAEIQAAChRHALBSGno4hAAEIQAACfhDALPihE6OEAAQgAAEIFEYAs1AYel0d\/\/znP5ff\/OY3cvz4cZk+fbquyTGb4IvVZs2aJWPHjpWLL744INLc3CyNjY2yZcsWdV+ag+S1CZhn3Px66623ZNGiReAqAQHMQglEznqK\/fr1kzVr1sj1118vH3\/8MWYha+A5tz98+HB54oknZPDgwaE9NzU1yeLFi+Xw4cM5j4zuiiBg1sOKFSsCA\/nGG29gFooQoYA+MQsFQNfUpTEKq1atkhtuuCGYFmZBk7oi3bt3l3Xr1snNN98sZ8+elV27dsmGDRuCSc6dO1cmTpwoPXr0kP3798uCBQukra1NFwBmcwEBYxSefPJJGTRoUPD\/MQvlWSCYhfJobX2mP\/7xj+WRRx6RYcOGfd82ZsE65kIbnDx5ssyfPz8wDdu3b5eVK1deMJ6FCxfKlClT5Ntvvw3+bvfu3YWOl86zI3DrrbfKvHnz5PLLL\/++E8xCdrxdaxmz4JoiHoynd+\/eMnv2bJk0aZL07NlTvvrqKzl\/\/rz079+fyIIH+sUZogk3T5gwQb744gsxxqDzUYP5l6YxCZdddlkQdVi6dGmc5inrAQFz3PDrX\/9abrnlliCK9Pnnnwe\/m6giZsEDAS0NEbNgCWSZmrnrrrsCs2D+tXno0KFgs1iyZImMGDECs6BoIVxxxRWydu1aGTp0qLz33nty\/\/33h87OHEuMGjVKjh49Kg888ICcPHlSEQWmMmfOHJk6dWpwDGUSGl9++eUgR+XKK6\/ELJRoeWAWSiS2ran+8pe\/lF\/84hfyyiuvyF\/\/+tegWZMRj1mwRdiNdkaOHBkkspmwc62owfLly8WEqE+cOBFEH44cOeLGBBiFFQLGAP70pz+VP\/7xj8HbL+bY0fwDAbNgBa83jWAWvJHK7YFiFtzWJ8nozGuSy5YtC16V3LZtW5DoGPap\/Mvz9OnTQYRp3759SbqjjicEMAueCGV5mJgFy0DL2hxmQZ\/ymAV9mtqYEWbBBkX\/2sAs+KeZkyPGLDgpS6pBmQRW8zpkr169iCykIqmrMmZBl55RZ4NZiEqKcjUJYBb0LRAiC\/o0tTEjzIINiv61gVnwT7PMRlzZ8MM6qHVmbcpjFjKTpbCGMQuFoXe6Y8yC0\/JkNjjMQmZo\/WsYs+CfZlmOuOPbEDt27AjejAj78DZEliq41zZmwT1N8hgRZiEPyiXog8iCPpH79u0r69evD16Ve\/vtt4PrncM+3LOgT\/taM8IslEvvymwxC+XU3fqsMQvWkTrRYOUGx88++yy4Q8FcvNTxE\/UuBicmwyCsEMAsWMHoXSOYBe8kc3PAmAU3dUk7qqjfDWG+QGr16tWyc+fOtF1S33ECmAXHBcpoeJiFjMCWrVnMgk7FzVGE+VbRG2+8sea3Tr777rtiLmfiWyd1roOOs8Is6Nc4bIaYhXLqbn3WmAXrSJ1psPPXEnceWFNTU\/BdAZ2\/ZMqZCTAQqwQwC1ZxetMYZsEbqdweKGbBbX3Sjs588+CsWbPEvE5prn82n+bm5uC7Aoz25nsh+JSDAGahHDp3niVmoZy6M2sIQAACEIBAZAKYhcioKAgBCEAAAhAoJwHMQjl1Z9YQgAAEIACByAQwC5FRURACEIAABCBQTgKYhXLqzqwhAAEIQAACkQn8f0YuMxfK1J2SAAAAAElFTkSuQmCC","height":210,"width":349}}
%---
%[output:0bb56482]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["31397","36"]]}}
%---
%[output:6debc09b]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"coeffs","rows":36,"type":"double","value":[["1.2500"],["4.0000"],["-2.0000"],["1.2500"],["2.5000"],["-0.0000"],["0.0000"],["0.0000"],["-0.0000"],["0.0000"],["0.0000"],["-0.0000"],["-0.0000"],["0.0000"],["0.0000"]]}}
%---
