clear all; close all; clc

delta_t = 0.001;
t_start = 0;
t_end = 0.001;
t = t_start: delta_t: t_end;

V_dot_target_initial = -10;

% Add disturbance?
%d = 2*(0.5-rand());% Random # between -1 and 1
%d = 1;
d = 0;

% Alpha's are one type of 'gain' for the SMC
alpha0 = 1;
alpha1 = 1;
% Eta is the other 'gain' for the SMC
eta = 1.1;

% For crude FBL
Kp = -1;

% Saturation on the Lyap controller
u_max = 2;
u_min = -2;
apply_saturation = true;

count_V_incr = 0;

% Loop through a 3D grid of x1, x2, x3. Go at least through 3^0.5~1.74
delta = 0.1;
for (x1_IC = -1.74: delta: 1.74)
    for (x2_IC = -1.74: delta: 1.74)
        for (x3_IC = -1.74: delta: 1.74)
            
            % if on the sphere
            distance = x1_IC^2+x2_IC^2+x3_IC^2;
            if ( (2.9 <= distance) && (3.1 >= distance))
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % IC's and simulation parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                x_IC = [x1_IC x2_IC x3_IC];
                
                % Pre-allocate
                x_OL = zeros(length(t), 3);
                x_OL(1,:) = x_IC;
                y_OL = zeros(length(t), 1);
                y_OL(1) = x_IC(2);      % y = x2
                
                x_CL = zeros(length(t), 3);
                x_CL(1,:) = x_IC;
                y_CL = zeros(length(t),1 );
                y_CL(1) = x_IC(2);      % y = x2
                
                V = zeros(length(t), 1);
                                
                u_lyap = zeros(length(t),1 );
                u_robust = zeros(length(t),1 );
                
                
                for epoch = 2: length(t)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Put the system in normal form, i.e. calculate xi
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    r = 2; % Relative order of the system
                    
                    dh_dx = [0 1 0]; % 1x3
                    
                    f = [ -x_CL(epoch-1,1);
                        x_CL(epoch-1,3);
                        x_CL(epoch-1,1)*x_CL(epoch-1,3) ];
                    % 3x1
                    
                    Lf_h = dh_dx * f; % scalar
                    
                    dLf_h_dx = [0 0 1]; % 1x3
                    
                    Lf_2_h = dLf_h_dx*f; % scalar
                    
                    g = [(2+x_CL(epoch-1,3)^2)/(1+x_CL(epoch-1,3)^2); 0; 1]; % 3x1
                    
                    Lg_Lf_h = dLf_h_dx * g; % scalar, dLf_h/dx*g
                    
                    xi(1) = x_CL(epoch-1,2);  % xi(1) = h(x) = x2
                    xi(2) = Lf_h;
                    
                    V(1) = 0.5*(xi*xi');
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Calculate u_lyap with the switched Lyapunov algorithm
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Calculate dV1_dot_du
                    dV1_dot_du = xi(r)*Lg_Lf_h;
                    
                    % Calculate dV2_dot_du
                    dV2_dot_du = (xi(r)*(0.9+0.1*abs(xi(r)-1)) + 0.1*V(epoch-1)*sign( xi(r)-1 ) )*Lg_Lf_h;
                    
                    % Calculate V_dot_target
                    V_dot_target = (V(epoch-1)/V(1))^2*V_dot_target_initial;
                    
                    % Compare dV1_dot_du and dV2_dot_du to choose the CLF
                    dV_dot_du(epoch,:) = [dV1_dot_du dV2_dot_du];
                    
                    [M,I] = max(abs(dV_dot_du));
                    
                    % Calculate u_lyap with the CLF of choice
                    if ( I==1 ) % use V1
                        using_V1(epoch) = y_CL(epoch);
                        u_lyap(epoch) = (V_dot_target - xi(1)*Lf_h - xi(r)*Lf_2_h) /...
                            dV1_dot_du;
                    else %use V2
                        using_V2(epoch) = y_CL(epoch);
                        u_lyap(epoch) = (V_dot_target -...
                            xi(1)*(0.9+0.1*abs(xi(r)-1))*Lf_h-... % for xi(1)
                            (xi(r)*(0.9+0.1*abs(xi(r)-1) )+0.1*V(epoch-1)*sign(xi(r)-1))*Lf_2_h )/... % for xi(r)
                            dV2_dot_du;
                    end
                    
                    % Saturation
                    if (apply_saturation)
                        if (u_lyap(epoch) > u_max)
                            u_lyap(epoch) = u_max;
                        end
                        if (u_lyap(epoch) < u_min)
                            u_lyap(epoch) = u_min;
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Calculate u_robust with crude FBL, then proportional
                    % control
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    u_robust(epoch) = -x_CL(epoch-1,2)*x_CL(epoch-1,3) +Kp * xi(1);
                    
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     % Calculate u_robust with the SMC algorithm
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                     % Calculate omega (the surface)
%                     
%                     dxi1_dt = xi(2);    % xi(1)_dot = xi(2)
%                     %dxi2_dt = 0; % TO DO - filter to calculate this
%                     
%                     omega = +alpha0*xi(1)+alpha1*dxi1_dt;
%                     %omega(2) = -alpha0*xi(2)-alpha1*dxi2_dt;
%                     
%                     % Calculate u_robust
%                     u_robust(epoch) = -eta*sign(omega);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Apply (u = u_lyap + u_robust) to the system and simulate for one time step
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    xy = simulate_sys( x_CL(epoch-1,:), y_CL(epoch-1), u_robust(epoch), delta_t); %u_lyap(epoch)+d
                    x_CL(epoch,:) = xy(1:3);    % First 3 elements are x
                    y_CL(epoch) = xy(end);         % Final element is y
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Simulate the open-loop system
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xy = simulate_sys( x_OL(epoch-1,:), y_OL(epoch-1), 0, delta_t);
                    x_OL(epoch,:) = xy(1:3);    % First 3 elements are x
                    y_OL(epoch) = xy(end);         % Final element is y
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Update V
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    V(epoch) = 0.5*(xi*xi');
                    
                    % Check if V increased
                    if ( V(epoch)>V(epoch-1) )
                        count_V_incr = count_V_incr+1;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hold on
                axis equal
                set(gcf,'color','w');
                plot3(x_CL(:,1), x_CL(:,2), x_CL(:,3), '-o')
                xlabel('x_1')
                ylabel('x_2')
                zlabel('x_3')
                title('x: Closed Loop')
            end
        end
    end
end

disp('Count of increased Vs:')
disp(count_V_incr)
