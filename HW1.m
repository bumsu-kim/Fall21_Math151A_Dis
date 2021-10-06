%% Problem 1
disp(' ');
disp('==== Problem 1 ====');
% f(x) = 1.01*exp(4*x) - 4.62*exp(3*x) - 3.11*exp(2*x) + 12.2*exp(x) - 1.99

d = 3; % 3-digit
r = true; % rounding

fval = f(1.53);
faval = fa(1.53, d, r);
fbval = fb(1.53, d, r);
fprintf('f val    = %4f\n\nf val(a) = %4f\nf val(b) = %4f\n', fval, faval, fbval);
fprintf('\nerror(a) = %4f\nerror(b) = %4f\n', fval-faval, fval-fbval);

% if you're more interested in the errors:
% x = linspace(1.50,1.56,1001)';
% y = f(x);
% ya = fa(x, d, r); erra = y-ya;
% yb = fb(x, d, r); errb = y-yb;
% 
% figure(); hold on;
% plot(x,y);
% plot(x,ya);
% plot(x,yb);
% legend('f', 'a', 'b', 'location','best');
% 
% figure(); hold on;
% plot(x, abs(erra));
% plot(x, abs(errb));
% legend('a err','b err', 'location', 'best');

%% Problem 2
disp(' ');
disp('==== Problem 2 ====');
fl4 = @(x) fl(x, 4, false); % 4-digit chopping

a = fl4(1/3);
b = fl4(123/4);
c = fl4(1/6);

% part a
sqrtpart = fl4(sqrt( fl4(b*b) - fl4(fl4(4*a)*c)));
numerator = fl4(b - sqrtpart);
x1a = fl4(numerator/ fl4(2*a));
fprintf('\n(a) x1      = %5f\n', x1a);

% part b
x1_actual = fl4(3/2*(123/4 - sqrt( 123*123/4/4 - 4/3/6)));
% or x1_actual = 0.005420;
fprintf('\n(b) abs err = %5f\n(b) rel err = %5f\n',...
    abs(x1_actual-x1a), abs( (x1_actual-x1a)/x1_actual ) );

% part c
x1c = fl4( fl4(2*c) / fl4( b + sqrtpart));
fprintf('\n(c) x1      = %5f\n', x1c);

% part d
fprintf('\n(d) abs err = %5f\n(d) rel err = %5f\n',...
    abs(x1_actual-x1c), abs( (x1_actual-x1c)/x1_actual ) );

% part e
fprintf('\n(e) Note that \nb = %5f, fl_sqrt(..) = %5f,\n\t and fl(b - fl_sqrt(..)) = %5f,\n',...
    b, sqrtpart, b-sqrtpart);
fprintf('while the actual sqrt(...) = %5f,\n\t and b-sqrt(...) = %5f\n',...
    sqrt(b*b-4*a*c), b-sqrt(b*b-4*a*c));

%% Problem 3
disp(' ');
disp('==== Problem 3 ====');

% part a
a1 = 4/5 + 1/3;
a2 = 1/3 + 3/11 - 3/20;

fprintf('\n (a)\n\t(i)  = %4f\n\t(ii) = %4f\n',a1,a2);

% part b
fl3 = @(x) fl(x, 3, false); % 3-digit chopping
b1 = fl3( fl3(4/5) + fl3(1/3));
b2 = fl3( fl3( fl3(1/3) + fl3(3/11)) - fl3(3/20) );
fprintf('\n (b)\n\t(i)  = %4f\n\t(ii) = %4f\n',b1,b2);

% part c
fl3 = @(x) fl(x, 3, true); % now rounding
c1 = fl3( fl3(4/5) + fl3(1/3));
c2 = fl3( fl3( fl3(1/3) + fl3(3/11)) - fl3(3/20) );
fprintf('\n (c)\n\t(i)  = %4f\n\t(ii) = %4f\n',c1,c2);

% part d
err_b1 = (a1-b1)/a1;
err_c1 = (a1-c1)/a1;
err_b2 = (a2-b2)/a2;
err_c2 = (a2-c2)/a2;
fprintf('\n (d) rel err\n\t(i)\n\t\tchopping = %4f\n\t\trounding = %4f\n',err_b1,err_c1);
fprintf('\n\t(ii)\n\t\tchopping = %4f\n\t\trounding = %4f\n',err_b2,err_c2);


%% Coding Problem
disp(' ');
disp('==== Coding Problem ====');

fcoding = @(x) x*x-3;
sol = sqrt(3); % desired sol

a = 0; f_a = fcoding(a);
b = 4; f_b = fcoding(b);
if (f_a*f_b > 0)
    disp('Wrong initial interval. Try with different a and b');
else
    c = (a+b)/2;
    itcnt = 0;
    rel_err = @(c) abs(sol-c)/sol;
    
    while( rel_err(c) > 5*10^-5 )
        itcnt = itcnt + 1;
        c = (a+b)/2;
        f_c = fcoding(c);
        if (f_a*f_c < 0)
            b = c;
            f_b = f_c;
        else
            a = c;
            f_a = f_c;
        end
    end
    fprintf( '\nsqrt(3) = %6f\napprox = %6f\nrel err = %.5g\n', ...
        sol, c, rel_err(c));
end


%% Helper functions

function y = f(x)
 y = 1.01*exp(4*x) - 4.62*exp(3*x) - 3.11*exp(2*x) + 12.2*exp(x) - 1.99;
%  y = fl3(y);
end

function y = fa(x, d, r)
t1 = mult( 1.01, pow(fl(exp(x), d, r), 4, d, r), d, r);
t2 = mult( -4.62, pow(fl(exp(x), d, r), 3, d, r), d, r);
t3 = mult( -3.11, pow(fl(exp(x), d, r), 2, d, r), d, r);
t4 = mult( 12.2, fl(exp(x), d, r), d, r);
t5 = fl(-1.99, d, r);

y = add_fl(add_fl(add_fl(add_fl(t1, t2, d, r), t3, d, r), t4, d, r), t5, d, r);
end

function y = fb(x, d, r)
y =  mult(fl(exp(x), d, r), 1.01, d, r);
y = add_fl(-4.62, y, d, r);
y = mult(fl(exp(x), d, r), y, d, r);
y = add_fl(-3.11, y, d, r);
y = mult(fl(exp(x), d, r), y, d, r);
y = add_fl(12.2, y, d, r);
y = mult(fl(exp(x), d, r), y, d, r);
y = add_fl(-1.99, y, d, r);
end

% Define d-digit rounding/chopping mantissa operations

% floating point approx
function y = fl(x, d, rounding)
% get sign and make x nonnegative
s = sign(x);
x = abs(x);
ex = zeros(size(x));
nonzeros = x>0;
N = 10^(d-1);
%get exp and mantissa
ex(nonzeros) = floor(log10(x(nonzeros)));
if rounding
    m = floor(x.*10.^(-ex)*N + 0.5)/N;
else
    m = floor(x.*10.^(-ex)*N)/N;
end
y = s.*m.*10.^ex;
end

% multiplication
function z = mult(x,y,d,r)
% truncation
x = fl(x, d, r);
y = fl(y, d, r);
z = fl(x.*y, d, r);
end


% integer power
function y = pow(x,n, d, r)
x = fl(x, d, r); y = ones(size(x));
if n>0
    for i=1:n
        y = mult(y,x, d, r);
    end
elseif n<0
    % not implemented because we don't need it
end
y = fl(y, d, r);
end

% addition
function z = add_fl(x,y, d, r)
% truncation
x = fl(x, d, r);
y = fl(y, d, r);
z = fl(x+y, d, r);
end