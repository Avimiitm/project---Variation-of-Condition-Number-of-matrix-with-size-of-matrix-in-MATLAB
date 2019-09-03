clear all;
clc;
N=[10, 20, 50, 100, 200, 500,1000,2000];
for lp= 1:length(N)
    dt = (2.*pi)./N(lp); % discretize the boundary using N nodes
    
    % defining the function  f(t)=1/(2*pi*(|y(t)-z(t)|))
    % where, y(t) = integration(sqrt(9-5cos^2t)) and z(t) = 6/sqrt(4*cos^2t+ 9*sin^2t)
    % and dy = dt*sqrt(9-5cos^2t) and z1 = 3*cost and z2 = 2*sint
    f =@(t)((sqrt(9-5.*(cos(t.*(180./pi))).^2))./((2*pi).*sqrt((36./(4.*(cos(t.*(180./pi))).^2 + 9.*(sin(t.*(180./pi))).^2)))+ (2*ellipticE(t.*(180./pi), -5/4)).^2));
    
    % value of z^2 at x=(2,1) and at x(2,1) the value of t is = 0.6435
    z1 = 36 ./(4.*(cos(0.6435)).^2 + 9.*(sin(0.6435)).^2);
    % function defination for calculating the value of u(2,1)
    f1 =@(t)((sqrt(9-5.*(cos(t.*(180./pi))).^2))./((2*pi).* sqrt(z1+ (2*ellipticE(t.*(180./pi), -5/4)).^2)));
    
    % integration by Trapezoidal rule
    for k = 1:N(lp)
        if(k == 1)
            f_n(k) = 0.5.*dt.*f(k.*dt);
            f_n1(k) = 0.5.*dt.*f1(k.*dt);
        elseif(k == N(lp))
            f_n(k) = 0.5.*dt.*f(k.*dt);
            f_n1(k) = 0.5.*dt.*f1(k.*dt);
        else f_n(k) = dt.*f(k.*dt);
            f_n1(k) = dt.*f1(k.*dt);
        end
        g(k) = ((3.*cos(k.*dt.*(180./pi)))+(2.*sin(k.*dt.*180./pi))).^2;
    end
    f_nv{lp} =f_n ;
    g_v{lp} = g;
    B= repmat(f_nv{lp},N(lp),1);
    B_v{lp} = B;
    A = (0.5.* eye(N(lp)))+ B_v{lp};
    
    A_v{lp} = A; %the reqired matrix
    
    % LU decomposition (direct method)
    [l u] = lu(A_v{lp});
    l_v{lp} = l;
    u_v{lp}  =u;
    Y = l\(g_v{lp})';
    Y_v{lp} = Y;
    h = u_v{lp}\Y_v{lp};
    h_v{lp} = h;
    f_nv1{lp} = f_n1;
    u = f_nv1{lp}*h_v{lp};
    u_v{lp}= u;
    
    % iteration method 
    % Jacobi method
    n_j = length(g);
    h_j= zeros(n_j,1);
    h_jnew = zeros(n_j,1);
    h_j(:) = 0;
    it = 100;
    tol = 0.0000001;
    for iteration = 1 : it
        convergence = true;
        for i = 1 : n_j 
            Sum = 0;
            for j = 1 : n_j 
                if (j ~= i)
                    Sum = Sum + A(i,j)*h_j(j);
                end
            end
            h_jnew(i) = -1/A(i,i) * (Sum - g(i));
            if abs(h_jnew(i) - h_j(i)) > tol
                convergence = false;
            end
        end
        if convergence
            break
        end
        h_j = h_jnew;
    end
    
    h_jv{lp} = h_jnew;
    u_j = f_nv1{lp}*h_jv{lp};
    u_jv{lp} = u_j;
    
    
    c_n(lp) = cond(A_v{lp});
end

fprintf( 'The value of u at x=(2,1) for direct(LU decomposition) method are\n');
u_v
fprintf('\n');
fprintf('The value of u at x =(2,1) for iterative(Jacobi) method are');
u_jv

plot(N,c_n,'-o')
grid on
xlabel('Size(N)of matrix ');
ylabel('Condition number');
title('Variation of the condition number with size of matrix')
    