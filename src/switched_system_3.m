clear;

% Part 1, 2

A1 = [-3 1; 1 2];
A2 = [1 -0.5; -3 -5];

eig_A1 = eig(A1);
eig_A2 = eig(A2);

disp('Eigenvalues of A1:');
disp(eig_A1);
disp('Eigenvalues of A2:');
disp(eig_A2);

% Part 3

w1_range = linspace(0, 1, 100); % Considering w1 values between 0 and 1

for w1 = w1_range
    w2 = 1 - w1;
    if w1 >= 0 && w2 >= 0 % Ensure that both w1 and w2 are non-negative
        A_ave = w1 * A1 + w2 * A2;

        eig_A_ave = eig(A_ave);

        % Checking if all eigenvalues have negative real part and none have zero real part
        if all(real(eig_A_ave) < 0) && all(real(eig_A_ave) ~= 0)
            disp('A stable combination of w1 and w2:');
            disp([w1, w2]);
            break;
        end
    end
end

% Part 4

T_range = linspace(0.01, 10, 100000); % Sample range of T values

for T = T_range
    A_d = expm(w2 * A2 * T) * expm(w1 * A1 * T);
    eig_A_d = eig(A_d);
    
    % Check if all eigenvalues lie strictly inside the unit circle
    if all(abs(eig_A_d) < 1)
        disp(['Stable T value: ', num2str(T)]);
        break;
    end
end

% Part 5

x0 = [4, 3]; % Sample initial condition
t_values = linspace(0, 150, 150000);
[x_values, sigma_values] = periodic_switching_response(A1, A2, w1, T, x0, t_values);

figure;
subplot(2, 1, 1);
plot(x_values(1, :), x_values(2, :));
xlabel('x1');
ylabel('x2');
title('State trajectory');
grid on;

subplot(2, 1, 2);
stairs(t_values, sigma_values, 'r');
xlabel('time');
ylabel('\sigma(t)');
title('Switching signal');
grid on;

% Function for the periodic switching response
function [x_values, sigma_values] = periodic_switching_response(A1, A2, w1, T, x0, t_values)
    x_values = zeros(2, length(t_values));
    sigma_values = zeros(1, length(t_values));
    x_values(:, 1) = x0;

    for i = 1:length(t_values)
        if mod(t_values(i), T) < w1 * T
            A_current = A1;
            sigma_values(i) = 1;
        else
            A_current = A2;
            sigma_values(i) = 2;
        end
        
        if i < length(t_values)
            dt = t_values(i+1) - t_values(i);
            % Using approximate discrete update for small dt
            x_values(:, i+1) = (eye(2) + A_current * dt) * x_values(:, i);
        end
    end
end
