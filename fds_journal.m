
% Finite Difference Solution of the Reynolds Equation for journal bearings
image_output_format = '-dpng';
image_suffix = '.png';
    %% 1) Domain Discretisation
    % Setup variables
        % These are always the same
         dtheta = 5;
         M = 360/dtheta + 1;     % Divisions X
         N = M;                  % Divisions Y
         D = 0.04;      % Bearing Diameter
         L = 0.02;      % Bearing Width Julias = 0.01
         c = 0.00005;   % Bearing Clearance
         epsilon = 0.7; % Eccentricity ratio
         RPM = 2000;    % Rotational speed
         eta = 0.1;     % Kinematic viscosity
         converge = 1e-3;

    %% 2) Bearing Geometry
    % Setup bearing information
    dx = D / (2*M); %(D/2) * 1/M
    dy = L*2 / (N-1);
    omega = 209.44;
    U = omega*D/2;
    
    fprintf('D        = %f\n', D)
    fprintf('L        = %f\n', L)
    fprintf('c        = %f\n', c)
    fprintf('epsilon  = %f\n', epsilon)
    fprintf('RPM      = %f\n', RPM)
    fprintf('omega    = %f\n', omega)
    fprintf('U        = %f\n', U)
    fprintf('eta      = %f\n', eta)
    fprintf('converge = %f\n', converge)
    
    fprintf('dx  = %f\n', dx)
    fprintf('dy  = %f\n', dy)

    %% 3) Setup Arrays
    % Setup the arrays nessessary for the problem

    y_values = zeros(1, N);
    for i = 1:N
        y_values(i) = L - dy * (i-1);
    end

    
    % On x axis
    theta_values = zeros(1, M);
    dtheta = 360 / (M-1);
    for i = 1:M
        theta_values(i) = (i-1) * dtheta;
    end
 
%% Theoretical Pressure Distribution Sommerfeld

P_sommerfeld = zeros(1, M);
r = (D/2);

for i = 2:M-1
    sin_value = sin(deg2rad(theta_values(i)));
    cos_value = cos(deg2rad(theta_values(i)));
    

    t_aa_ = (eta * U * r)/(c^2);
    t_ba_ = 6 * epsilon * sin_value * (2 + epsilon * cos_value);
    t_bb_ = (2 * epsilon^2) * (1 + epsilon * cos_value)^2;

    P_sommerfeld(i) = t_aa_ * (t_ba_ / t_bb_);
    
    if P_sommerfeld(i) < 0
        x_right_sommerfeld = i;
        break
    end
end

h = figure();
plot(theta_values(1:i-1), P_sommerfeld(1:i-1));
title('Theoretical Pressure Distribution (Sommerfeld)');
xlabel('Theta [deg]');
ylabel('Pressure [pa]');
print(h, image_output_format, ['pressure_dist_sommerfeld' image_suffix]);
close

%% Attitude angle

beta = atan((pi * sqrt(1 - epsilon^2))/(4 * epsilon));

%% Load Capacity

W_th_a = (eta * U * L^3)/(c^2) * ...
    ((epsilon * (pi^2 * (1 - epsilon^2) + 16 * epsilon^2)^.5)/ ...
    (4* (1-epsilon^2)^2));


%% Theoretical Pressure Distribution Ocvirk

P_ocvirk = zeros(N, M);
r = (D/2);

for i = 2:M-1
    for j = 1:N
        sin_value = sin(deg2rad(theta_values(i)));
        cos_value = cos(deg2rad(theta_values(i)));
        
        t_aa_ = (eta * U) / (r * c^2);
        t_bb_ = (L^2)/4 - y_values(j)^2;
        t_cc_ = (3 * epsilon * sin_value) / (1 + epsilon * cos_value)^3;
        P_ocvirk(j, i) = t_aa_ * t_bb_ * t_cc_;
    end
end

% Find the actual area

x_right = 0;
y_top = 0;
y_bottom = 0;

for j = 1:N
    if round(P_ocvirk(j, 2)) == 0
        if x_right == 0
            for i = 2:M
                if round(P_ocvirk(j+1, i)) == 0
                    x_right = i;
                    break
                end
            end
        end
        
        if y_top == 0
            y_top = j;
            continue
        end
        
        if y_bottom == 0
            y_bottom = j;
            break
        end
            
    end
end

width = x_right;
height = y_bottom - y_top;

P_ocvirk_reduced = zeros(height, width);

for j = y_top:y_bottom
    for i = 1:x_right
        P_ocvirk_reduced(i, j-y_top+1) = P_ocvirk(j, i);
    end
end

%% Make plots

h = figure('Name','Theoretical (Ocvirk) Pressure Distribution, Front View','NumberTitle','off');
surf(y_values(y_top:y_bottom), theta_values(1:x_right), P_ocvirk_reduced, 'FaceLighting','phong');
camorbit(37.5,-31);
title('Theoretical (Ocvirk) Pressure Distribution, Front View');
ylabel('Theta [deg]');
xlabel('Bearing Width [m]');
zlabel('Pressure [Pa]');
print(h, image_output_format, ['pressure_dist_front' image_suffix]);
close

h = figure('Name','Theoretical (Ocvirk) Pressure Distribution, Side View','NumberTitle','off');
surf(y_values(y_top:y_bottom), theta_values(1:x_right), P_ocvirk_reduced, 'FaceLighting','phong');
camorbit(127,-31);
title('Theoretical (Ocvirk) Pressure Distribution, Side View');
ylabel('Theta [deg]');
xlabel('Bearing Width [m]');
zlabel('Pressure [Pa]');
print(h, image_output_format, ['pressure_dist_right' image_suffix]);
close

h = figure('Name','Theoretical (Ocvirk) Pressure Distribution, Isometric View','NumberTitle','off');
surf(y_values(y_top:y_bottom), theta_values(1:x_right), P_ocvirk_reduced, 'FaceLighting','phong');
camorbit(77,15);
title('Theoretical (Ocvirk) Pressure Distribution, Isometric View');
ylabel('Theta [deg]');
xlabel('Bearing Width [m]');
zlabel('Pressure [Pa]');
print(h, image_output_format, ['pressure_dist_isometric' image_suffix]);
close

P_ocvirk_mid = zeros(1, x_right);

[p_max idx] = max(P_ocvirk_reduced);
[p_max idx] = max(p_max);

for i = 1:x_right
    P_ocvirk_mid(i) = P_ocvirk_reduced(i, idx);
end

h = figure('Name','Theoretical Pressure Distribution, Comparison','NumberTitle','off');
plot(theta_values(1:x_right_sommerfeld-1), P_sommerfeld(1:x_right_sommerfeld-1), ...
     theta_values(1:x_right), P_ocvirk_mid, '-.');
title('Theoretical Pressure Distribution, Comparison');
legend('Sommerfeld', 'Ocvirk', 'location', 'northwest')
xlabel('Theta [deg]');
ylabel('Pressure [Pa]');
print(h, image_output_format, ['pressure_dist_theoretical' image_suffix]);
close

% Theoretical Work
fprintf('\nTheoretical Load Capacity = %f\n', W_th_a)
fprintf('\nAttitude angle beta       = %f\n rad', beta)
