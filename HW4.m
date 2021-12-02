% For HW4

%% Problem 3
fprintf('\n\n======Problem 3======');
xnodes = [0, 10, 20, 30, 40, 50];
fnodes = [179323, 203302, 226542, 249633, 281422, 308746];

syms x % x as a symbol

L = Lag_int(x, xnodes, fnodes);
fprintf('P_5''s coefficients: ');
disp(sym2poly(L));
% Plotting
N = 61;
% % naive idea
% tic;
% xplot = linspace(0,60,N);
% yplot = zeros(N,1);
% for i=1:N
%     x = xplot(i);
%     yplot(i) = subs(L);
% end
% figure(); plot(xplot, yplot);
% t_naive = toc;
% fplot
tic;
figure(); fplot(L, [-10,60]);
title('US Population History, x=Year-1960')
t_fplot = toc;

fprintf('\nTIMING:\n\tnaive = %f\n\tfplot = %f\n', t_naive, t_fplot);
x = 2020-1960;
fprintf('In 2020, prediction = %d\n', subs(L));

%% Problem 4
fprintf('\n\n======Problem 4======\n');
xnodes = [0, .25, .5, .75];
fnodes = [1, 1.64872, 2.71828, 4.48169];
n = length(xnodes);

dd_1 = zeros(n-1, 1);
for i=1:(n-1)
    dd_1(i) = (fnodes(i)-fnodes(i+1))/(xnodes(i)-xnodes(i+1));
end

dd_2 = zeros(n-2, 1);
for i=1:(n-2)
    dd_2(i) = (dd_1(i)-dd_1(i+1))/(xnodes(i)-xnodes(i+2));
end


dd_3 = zeros(n-3, 1);
for i=1:(n-3)
    dd_3(i) = (dd_2(i)-dd_2(i+1))/(xnodes(i)-xnodes(i+3));
end

P1 = @(x) fnodes(1) + dd_1(1)*(x-xnodes(1));
P2 = @(x) P1(x) + dd_2(1)*(x-xnodes(1))*(x-xnodes(2));
P3 = @(x) P2(x) + dd_3(1)*(x-xnodes(1))*(x-xnodes(2))*(x-xnodes(3));
syms x;
disp(sym2poly(P1(x)));
disp(sym2poly(P2(x)));
disp(sym2poly(P3(x)));

fprintf('P1(0.43) = %f\nP2(0.43) = %f\nP3(0.43) = %f\n',P1(0.43), P2(0.43), P3(0.43));

%% Problem 5
fprintf('\n\n======Problem 5======');
xnodes = [-.1, 0, .2, .3];
fnodes = [5.3, 2, 3.19, 1];
n = length(xnodes);

dd_1 = zeros(n-1, 1);
for i=1:(n-1)
    dd_1(i) = (fnodes(i)-fnodes(i+1))/(xnodes(i)-xnodes(i+1));
end

dd_2 = zeros(n-2, 1);
for i=1:(n-2)
    dd_2(i) = (dd_1(i)-dd_1(i+1))/(xnodes(i)-xnodes(i+2));
end


dd_3 = zeros(n-3, 1);
for i=1:(n-3)
    dd_3(i) = (dd_2(i)-dd_2(i+1))/(xnodes(i)-xnodes(i+3));
end

P1 = @(x) fnodes(1) + dd_1(1).*(x-xnodes(1));
P2 = @(x) P1(x) + dd_2(1).*(x-xnodes(1)).*(x-xnodes(2));
P3 = @(x) P2(x) + dd_3(1).*(x-xnodes(1)).*(x-xnodes(2)).*(x-xnodes(3));

syms x;
P3x = expand(P3(x));
disp(P3x); % - (1670*x^3)/3 + (371*x^2)/2 - (533*x)/60 + 2
% now you can do the nested evaluation and appy 3-digit-chopping arithmetic

% adding a new point
xnodes(end+1) = .35; fnodes(end+1) = .973;
% need to do only one computation for each dd_i
i = n;   dd_1(i) = (fnodes(i)-fnodes(i+1))/(xnodes(i)-xnodes(i+1));
i = n-1; dd_2(i) = (dd_1(i)-dd_1(i+1))/(xnodes(i)-xnodes(i+2));
i = n-2; dd_3(i) = (dd_2(i)-dd_2(i+1))/(xnodes(i)-xnodes(i+3));
i = 1; dd_4(i) = (dd_3(i)-dd_3(i+1))/(xnodes(i)-xnodes(i+4));
P4 = @(x) P3(x) + dd_4(1).*(x-xnodes(1)).*(x-xnodes(2)).*(x-xnodes(3)).*(x-xnodes(4));

figure();
plot(xnodes, fnodes, '*');
hold on;
xplot = linspace(-0.2, 0.4);
plot(xplot, P1(xplot));
plot(xplot, P2(xplot));
plot(xplot, P3(xplot));
plot(xplot, P4(xplot));
legend('nodes', ...
    'P1', 'P2', 'P3', 'P4', ...
    'Location', 'best');
title('Polynomical Interpolation using Newton''s Divided Differences');

%% Helper functions
% Lagrange interpolation as symbolic calculations
function L = Lag_int(x, xnodes, fnodes)
n = length(xnodes);
L = sym(0);
for i=1:n
    L = L + Lagr(x, i, xnodes)*fnodes(i);
end
end

function P = Lagr(x, k, xnodes)
n = length(xnodes);
P = sym(1);
for i=1:n
    if i~=k
        P = P*(x-xnodes(i))/(xnodes(k)-xnodes(i));
    end
end
end