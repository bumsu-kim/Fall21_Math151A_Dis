%% Prob 1
disp(' ====== Problem 1 ====== ');
% (a)
a = 1; b = 2;
f = @(x) x.*log(x);
n = 4;

% analytic form of the indefinite integral:
% .25*x^2*(2*log(x)-1)
I = 2*log(2) - 1 + .25;

disp_res('part (a)', a,b,f,n,I);

% (b)
a = 1; b = 3;
f = @(x) x./(x.^2+4);
n = 8;
% analytic form of the indefinite integral:
% .5*log(x^2+4)
I = .5*(log(13)-log(5));
disp_res('part (b)', a,b,f,n,I);

% (c)
a = 0; b = 3*pi/8; 
f = @(x) tan(x);
n = 8;
% analytic form of the indefinite integral:
% -log(cos(x))
I = log(cos(0)) - log(cos(3*pi/8));
disp_res('part (c)', a,b,f,n,I);

% (d)
a = -1; b = 1;
f = @(x) 1/sqrt(2*pi) * exp(-x.^2/2);
n = 8;
% no analytic solutions in terms of simple functions
% but see "error function (erf)" for more
I = 0.5*(erf(1/sqrt(2)) - erf(-1/sqrt(2)));
disp_res('part (d)', a,b,f,n,I);


%% Prob 3
fprintf('\n');
disp(' ====== Problem 3 ====== ');
x = [-.3, -.2, -.1, 0];
h = .1;
fx = [-.27652, -.25074, -.16134, 0];

% part a
dfx = zeros(size(x));
dfx(1) = 1/(2*h)*(-3*fx(1) + 4*fx(2) - fx(3));
for i=2:3
    dfx(i) = 1/(2*h)*(fx(i+1)-fx(i-1));
end
dfx(end) = 1/(2*h)*(3*fx(end) - 4*fx(end-1) + fx(end-2));

% part b
f = @(x) exp(2*x) - cos(2*x);
df =@(x) 2*exp(2*x) + 2*sin(2*x);

fprintf('actual    f''(x) = '); disp(df(x));
fprintf('estimated f''(x) = '); disp(dfx);
fprintf('    errors      = '); disp(df(x)-dfx);
xgrid = linspace(-.3, 0);
figure(); hold on;
plot(xgrid, f(xgrid), 'r');
window = linspace(-0.02, 0.02);
for i=1:4
    plot(x(i)+window, f(x(i)) + df(x(i))*window, 'k', 'LineWidth',2);
    plot(x(i)+window, f(x(i)) + dfx(i)*window, 'b', 'LineWidth',2);
end
plot(x, fx, 'r*');

%% Prob 4
fprintf('\n');
disp(' ====== Problem 4 ====== ');
f = @(x) 3*x.*exp(x) - cos(x);
d2f = @(x) 6*exp(x)+3*x*exp(x)+cos(x);
x0 = 1.3;
for k=1:11
    h = 10^(-k);
    x = x0 + [-1, 0, 1]*h;
    fx = f(x);
    d2f_num(k) = (fx(1)-2*fx(2)+fx(3))/h/h;
end

fprintf('actual f''''(1.3) = %e\n', d2f(1.3));
fprintf('approximations:\n\th\t\t\tf''''\t\t\t abs error\n');
for k=1:11
    fprintf('\t%.e\t%e\t%e\n',10^(-k), d2f_num(k), d2f_num(k)-d2f(1.3));
end

%% Prob 5
fprintf('\n');
disp(' ====== Problem 5 ====== ');
q_pts = [-sqrt(3/5), 0, sqrt(3/5)]; % for [-1,1]
weights = [5/9, 8/9, 5/9];

% part a
disp('(a)');
a = 0; b = 1;
f = @(x) x.^2.*exp(-x);
I = 2-5/exp(1); % ~ 0.1606
I_Gauss = (b-a)/2* sum(weights.*f( (b+a)/2 + (b-a)/2*q_pts ));
fprintf('\tI = %f\n', I);
fprintf('\tI_Gauss = %f\n', I_Gauss);
fprintf('\tError = %e\n', I-I_Gauss);

% part b
disp('(b)');
a = 0; b = pi/4;
f = @(x) x.^2.*sin(x);
I = -2+sqrt(2)- (pi-8)*pi/16/sqrt(2); % ~ 0.08876
I_Gauss = (b-a)/2* sum(weights.*f( (b+a)/2 + (b-a)/2*q_pts ));
fprintf('\tI       = %f\n', I);
fprintf('\tI_Gauss = %f\n', I_Gauss);
fprintf('\tError   = %e\n', I-I_Gauss);

%% Prob 6
a = 1; b = 2;
f = @(x) x.*log(x);
I = 2*log(2) - 1 + .25;
fprintf('\nNumerical Integration\n');
fprintf('\tI             = %f\n', (I));
fprintf('Errors\n');
fprintf('\tn\t\tTrapezoidal\t\tSimpson\t\t\tMidpoint\n');
for k=1:14
n = 2^(k+2);
h = (b-a)/n;
% Simpson
idx = 0:n;
x = a + idx*h;
evenidx = (2:2:(n-2)) + 1; % matlab idx starts from 0
oddidx  = (1:2:(n-1)) + 1; % matlab idx starts from 0
ssum = 2*sum(f(x(evenidx))) + 4*sum(f(x(oddidx)));
I_Simpson = h/3*(f(a) + ssum + f(b));

% Trapezoidal
tsum = 2*sum(f(x(2:end-1)));
I_Trap = h/2*(f(a) + tsum + f(b));

% midpoint
h_mid = (b-a)/(n+2);
x_mid = a + (1:2:(n+2)) * h_mid;
I_midpoint = 2*h_mid*sum(f(x_mid));
% fprintf(' n = %d\n', n);
fprintf('%8d\t%e\t%e\t%e\n', n, I-I_Trap, I-I_Simpson, I-I_midpoint);
end

%% Helper functions
function disp_res(str, a, b, f, n, I)
h = (b-a)/n;
% Simpson
idx = 0:n;
x = a + idx*h;
evenidx = (2:2:(n-2)) + 1; % matlab idx starts from 0
oddidx  = (1:2:(n-1)) + 1; % matlab idx starts from 0
ssum = 2*sum(f(x(evenidx))) + 4*sum(f(x(oddidx)));
I_Simpson = h/3*(f(a) + ssum + f(b));

% Trapezoidal
tsum = 2*sum(f(x(2:end-1)));
I_Trap = h/2*(f(a) + tsum + f(b));
fprintf('\nNumerical Integration (%s)\n', str);
fprintf('\tI             = %f\n', (I));
fprintf('\tI_Simpson     = %f\n', (I_Simpson));
fprintf('\tI_Trapezoidal = %f\n', (I_Trap));
fprintf('\n\terr_Simpson     = %f\n', (I_Simpson-I));
fprintf('\terr_Trapezoidal = %f\n', (I_Trap-I));
end