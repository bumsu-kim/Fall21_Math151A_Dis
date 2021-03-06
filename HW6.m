%% Part 1
fprintf('\nQuestion 1\n');
plot_on = true;

% Finite difference
f = @(x) exp(2*x);
df = @(x) 2*exp(2*x);

x = [1.1, 1.2, 1.3, 1.4]; % x(i+1)-x(i) is const
fx = [9.025013, 11.02318, 13.46374, 16.44465];

dfx = (fx(2:end) - fx(1:end-1))./ (x(2:end)-x(1:end-1)); % fwd diff
dfx(end+1) = (fx(end)-fx(end-1)) / (x(end)-x(end-1)); % bwd diff

fprintf('fw/bw FD results:\n');
fprintf('       x        f''\n');
disp([x',dfx']);

xgrid = linspace(1.1, 1.4);
fvals = f(xgrid);
dfvals = 2*fvals; % f' = 2*f

% plot
if plot_on
    figure(); hold on;
    plot(xgrid, fvals, 'r');
    plot(x,fx, 'ro');
    plot(xgrid, dfvals, 'b');
    plot(x,dfx, 'bo');
end
% Remark: central difference for the interior pts
h = x(2)-x(1);
dfx2 = zeros(size(fx));
dfx2(1) = [-3/2, 2, -1/2]*fx(1:3)'/h; % 2nd order accuracy
dfx2(2:end-1) = (fx(3:end) - fx(1:end-2))/2/h; % central diff
dfx2(end) = [1/2, -2, 3/2]*fx(end-2:end)'/h; % 2nd order accuracy

% plot
if plot_on
    plot(x, dfx2, 'mo');
    title('Estimated derivatives');
    legend('f','sampled pts','f''', 'fw/bw diff', '2nd order accurate method', 'Location', 'best' );
end

% Errors
fprintf('\nQuestion 2\n');
err = abs(df(x)-dfx);
err2 = abs(df(x)-dfx2);
errbd = h/2*(4*f(x(2:end)));
errbd(end+1) = errbd(end); % last elt(bw) = same error bd with the prev one
fprintf('\tFD errors (absolute):\n');
fprintf('       x     f''(est)   err(fw/bw) err bd    err(2nd order method) \n');
disp([x',dfx', err', errbd', err2']);
% plot
if plot_on
    figure(); hold on;
    plot(x, err);
    plot(x, err2);
    title('Absolute errors of finite difference');
    legend('fw/bw diff', '2nd order accurate method', 'Location', 'best' );
    % set(gca, 'YScale', 'log')
end

%% Part 2
% skipped
