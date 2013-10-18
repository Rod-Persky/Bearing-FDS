function max_p = fds_plain(N, M, L, B, h0, dhdx, U, eta, converge)
% Finite Difference Solution of the Reynolds Equation

    %% 1) Domain Discretisation
    % Setup variables
        % These are always the same
%         clear
%         clc
%         M = 50;      % Divisions X
%         N = 50;      % Divisions Y
%         B = 0.10;    % Bearing length (Direction of sliding) (x) [m]
%         h0 = 0.0001; % Bearing minimum height
%         U = -20;     % Runner speed [m/s], see page 6
%         eta = 0.2;   % Dynamic viscosity [Pa.s]
%         converge = 1e-10;
%         
%         % These are different between tasks
%         L = 0.50;      % Bearing width (y) [m]
%         dhdx = 0.0012; % Bearing slope





    %% 2) Bearing Geometry
    % Setup bearing information
    K = (L * dhdx) / h0;	%Convergance ratio... %((L*k+h0)-h0)/h0

    % Modified from original because code is linear
    dx = B/M;       %B/M
    dy = L/N;       %L/N
    dx2 = dx^2;     %dx^2
    dy2 = dy^2;     %dy^2
    
    fprintf('N        = %f\n', N)
    fprintf('M        = %f\n', M)
    fprintf('B        = %f\n', B)
    fprintf('L        = %f\n', L)
    fprintf('h0       = %f\n', h0)
    fprintf('m = dhdx = %f\n', dhdx)
    fprintf('U        = %f\n', U)
    fprintf('eta      = %f\n', eta)
    fprintf('converge = %f\n', converge)
    fprintf('kr  = %f\n', K)
    fprintf('dx  = %f\n', dx)
    fprintf('dy  = %f\n', dy)
    fprintf('dx2 = %f\n', dx2)
    fprintf('dy2 = %f\n', dy2)
    
    %% 3) Setup Arrays
    % Setup the arrays nessessary for the problem

    P = zeros((N+1), (M+1));
    H = zeros(1, M+1);
    CONV = zeros(1, 3000);

    %% 4) Setup Boundary Conditions
    % Set boundary contitions to zero

    % Already done as the arrays were initialised to 0

    %% 5) Setup H
    % Setup the film thickness array

    for i = 1:M+1
        H(i) = dhdx * dx * (i-1) + h0;
    end


    %% 6) Calculate
    % * Sweep each node (ex, boundary nodes) using eq. 5
    % * Repeat and calculate the change in the max
    % * Repeat until the change in the max is under 1% (or less, whatever)


    %The i-1, i+1 point is the boundary conditions
    last_max = 0;
    change = 0;
    iterations = 0;

    while(converge < (1-change))
        this_max = 0;
        for i = 2:M  
            for j = 2:N
                a = (P(j, i+1) + P(j, i-1))/dx2;
                b = (P(j+1, i) + P(j-1, i))/dy2;
                c = ((6*U*eta)/(H(i)^3)) * dhdx;
                da = (3/H(i)) * dhdx;
                db =  (P(j,i+1) - P(j,i-1))/(2*dx);
                d = da * db;
                e = 2/dx2 + 2/dy2;

                P(j, i) = (a + b - c + d)/e;

                if this_max < P(j, i)
                    this_max = P(j, i);
                end
            end
        end

        change = last_max / this_max;
        iterations = iterations + 1;

        if mod(iterations, 10) == 0
            CONV(iterations/10) = (1-change)*100;
        end
        
        if mod(iterations, 200) == 0
            fprintf('At iteration %i with convergence = %f\n', iterations, (1-change)*100);
        end


        last_max = this_max;
    end

    %% 7) Validate
    % Note that the convergence ratio is K is not k bearing slope
    % instead K is (h1 - h0)/h0

    P_t = zeros(1, M+1);
    
    for i = 2:M
        kx_b = 1 + (K * (i-1) * dx)/B;
        
        a = (6 * U * eta * B)/(K * h0^2);
        b = 1/kx_b;
        c = (1+K)/(2+K);
        d = 1/(kx_b^2);
        e = 1/(2+K);

        P_t(i) = a * (-b + c * d + e);
    end
    
    P_th = -((6 * U * eta * B) / (h0^2)) * (K)/(4 * (K + 2)*(K + 1));

    fprintf('\nPressure Theoretical Max      = %f \n', P_th);
    fprintf('Pressure Theoretical from P_t = %f \n', max(P_t))
    fprintf('Pressure Numerical Max        = %f \n\n', max(max(P)));
    
    %% 8) Task 3
    disp('TASK 3 - Calculate load capacity ---------------------------------------')
    %Calculate load capacity W of the bearing for each operating condition.
    % Also calcuate the theoretical load capacity Wth and compare SLFs


    % Calculate numerical pressure
    W = sum(sum(P)) * (dx * dy);
    
    a = 6 * U * eta * B^2;
    b = K^2 * h0^2;
    c = -log(K + 1)
    d = (2 * K) / (K + 2)
    
    W_th = (a/b) * (c + d) * L;
    SLF = W/W_th;

    fprintf('Work Theoretical    = %f \n', W_th)
    fprintf('Work Numerical      = %f \n', W)
    fprintf('Side Leakage Factor = %f percent \n\n', SLF*100)
    
    %% Task 4
    disp('TASK 4 - Calculate lubricant flow rate ---------------------------------')
    % Lubricant flow rate,
    % Calculae lubricant flow rates from isdes, entry, and exit
    % and compare to see if sum of flow out equals flow in. Calculate
    % leakage ratio (out sides / in) for each configuration


    % Remember, column/row 1 and [N,M]+1 are boundary conditions

    Qtxl = 0;
    Qtxr = 0;
    Qtyr = 0;
    Qtyf = 0;

    for i = 1:M+1
            QLeft =  (P(i, 2) - P(i, 1))   / (dy);                       %Left,  dp/dy
            QRight = (P(i, N+1) - P(i, N)) / (dy);                       %Right, dp/dy

            QLeft  = -((H(i)^3)/(12 * eta)) * QLeft;
            QRight = -((H(i)^3)/(12 * eta)) * QRight;

            Qtxl = Qtxl + QLeft  * dy;
            Qtxr = Qtxr + QRight * dy;
    end

    for j = 1:N+1
            QRear  = (-3*P(2, j) + 4*P(3, j) - P(4, j)) / (2*dx);        %Rear dp/dx
            QFront = (3*P(M, j) - 4*P(M-1, j) + P(M-2, j)) / (2*dx);     %Front dp/dy

            QRear  = -((H(2)^3)/(12 * eta)) * QRear  + U * (H(2) / 2);
            QFront = -((H(M)^3)/(12 * eta)) * QFront + U * (H(M) / 2);

            Qtyr = Qtyr + QRear  * dy;
            Qtyf = Qtyf + QFront * dy;

    %Excel formula
    %print *, '=(-3*',P(2, j),'+4*',P(3, j),'-',P(4, j),')/(2*',dx,')*-((',h(1),'^3)/(12*',eta,'))+',U,'*(',h(1),' / 2)'
    end

    q_in = -Qtyf;
    q_out = Qtyr + Qtxl - Qtxr;
    Q_th = U * h0 * (K+1)/(K+2) * L;

    fprintf('Flow Numerical Front = %f\n', Qtyf);
    fprintf('Flow Numerical Rear  = %f\n', Qtyr);
    fprintf('Flow Numerical Left  = %f\n', Qtxl);
    sprintf('Flow Numerical Right = %f\n', Qtxr);

    fprintf('Flow In  = %f\n', q_in);
    fprintf('Flow Out = %f\n', q_out);
    fprintf('Flow Theoretical  = %f\n', Q_th);
    fprintf('Leakage Ratio     = %f percent \n\n', (Qtxl - Qtxr)/Qtyf*100);

    %% Task 5
    disp('TASK 5 - Calculate Friction Force --------------------------------------')
    % Calculate friction force acting on the pad and the friction coefficient
    % of the bearing for each operating condition. Compare this with theoretical
    % friction force

    % Calculate Force theoretical from given equation
    F_th = ((U * eta * B)/h0) * (6/(K+2) - (4 * log(K+1))/K) * L;
    %Print *, '=((',U,' * ',eta,' * ',B,')/',h0,') * (6/(',Kr,'+2) - (4 * ln(',Kr,'+1))/',Kr,') * ',L

    DU_DZ = zeros(M+1, N+1);
    % Calculate internal points using centeral difference method (3a)
    for i = 2:M
        for j = 1:N+1
            DU_DZ(i, j) = (2*H(i) - H(i))/(2 * eta) * (P(i+1, j)-P(i-1, j))/2*dx - U/H(i);
        end
    end

    % We need to use a different formula for i = 1 and i = N+1
    for j = 1:N+1
        DU_DZ(1, j)   = (2*H(i) - H(1))/(2 * eta)   * (P(2,   j)-P(1, j))/dx - U/H(1);
        DU_DZ(M+1, j) = (2*H(i) - H(M+1))/(2 * eta) * (P(M+1, j)-P(M, j))/dx - U/H(M+1);
    end

    F_numeric = sum(sum(DU_DZ)) * eta * dx * dy;

    fprintf('Force Theoretical = %f\n', F_th)
    fprintf('Force numeric     = %f\n', F_numeric)



%     if (0 < tfr)
%         if (0.0001 < abs(tfr - F_numeric) / F_numeric)
%             disp('FTH FAIL')
%             sprinf('Check f_numeric = %f\n', tfr)
%         end
%     end

    disp('Done...')
    
%% Setup axis values

    %used for plotting
    x_values = zeros(1, M+1);
    y_values = zeros(1, N+1);
    
    for i = 1:M+1
        x_values(i) = (i-1) * dx;
    end
    
    for i = 1:N+1
        y_values(i) = (i-1) * dy;
    end
    
    max_axis = max(L, B);
    [max_Pth idx] = max(P);
    max_Pth = max(max_Pth);
    max_P_Pth = max(max(P_t), max_Pth);
    
    midline_P = zeros(1, M+1);
    idx = max(idx);
    fprintf('max index is %d\n', idx)
    
    for i = 1:M+1
        midline_P(i) = P(idx, i);
    end
    

    %% Results
    % Calculate and Show Results
    
    image_output_format = '-dpng';
    image_suffix = '.png';

    h = figure('Name','Convergence','NumberTitle','off');
    semilogy(CONV);
    title('Convergence');
    xlabel('Iteration x10-1');
    ylabel('% Change');
    print(h, image_output_format, ['convergence' image_suffix]);
    close
    

    h = figure('Name','Pressure Distribution, Right Side View','NumberTitle','off');
    surf(x_values, y_values, P, 'FaceLighting','phong');
    axis([0 max_axis 0 max_axis 0 max_Pth])
    camorbit(37,-30);
    title('Pressure Distribution, Right Side View');
    xlabel('Bearing Length [m]');
    ylabel('Bearing Width [m]');
    zlabel('Pressure [Pa]');
    print(h, image_output_format, ['pressure_dist_right' image_suffix]);
    close
    
    h = figure('Name','Pressure Distribution, Front View','NumberTitle','off');
    surf(x_values, y_values, P, 'FaceLighting','phong');
    axis([0 max_axis 0 max_axis 0 max_Pth])
    camorbit(127,-30);
    title('Pressure Distribution, Front View');
    xlabel('Bearing Length [m]');
    ylabel('Bearing Width [m]');
    zlabel('Pressure [Pa]');
    print(h, image_output_format, ['pressure_dist_front' image_suffix]);
    close
    
    h = figure('Name','Pressure Distribution, Isometric View','NumberTitle','off');
    surf(x_values, y_values, P, 'FaceLighting','phong');
    axis([0 max_axis 0 max_axis 0 max_Pth])
    camorbit(77,15);
    title('Pressure Distribution, Isometric View');
    xlabel('Bearing Length [m]');
    ylabel('Bearing Width [m]');
    zlabel('Pressure [Pa]');
    print(h, image_output_format, ['pressure_dist_isometric' image_suffix]);
    close
    
    h = figure('Name','Pressure Distribution, Isometric View. Skewed aspect ratio','NumberTitle','off');
    surf(x_values, y_values, P, 'FaceLighting','phong');
    axis tight
    camorbit(77,15);
    title('Pressure Distribution, Isometric View. Skewed aspect ratio');
    xlabel('Bearing Length [m]');
    ylabel('Bearing Width [m]');
    zlabel('Pressure [Pa]');
    print(h, image_output_format, ['pressure_dist_isometric_skewed' image_suffix]);
    close
    
    h = figure('Name','Pressure Distribution, Theoretical','NumberTitle','off');
    plot(x_values, P_t, x_values, midline_P, '-.');
    ylim([0 max_P_Pth])
    title('Pressure Distribution, Comparison to Theoretical (Numerical is dashed)');
    print(h, image_output_format, ['pressure_dist_theoretical' image_suffix]);
    close

    save('FDS_Plain_Workspace.mat')


    fprintf('done\n')
end


