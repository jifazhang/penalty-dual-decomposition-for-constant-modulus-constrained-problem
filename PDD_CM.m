clc,clear;
close all;

N = 10;

A = randn(N,N) + 1i*rand(N,N);
A = (A*A');
A = (A + A')/2;
c = 0.9;
gamma=1e-3;
iter_max = 1e4;
x = zeros(N,iter_max+1);
x(:,1) = exp(1i*2*pi*rand(N,1));
y = zeros(N,iter_max+1);
y(:,1) = exp(1i*2*pi*rand(N,1));
u = zeros(N,iter_max+1);
u(:,1) = rand(N,1) + 1i * rand(N,1);

varrho = 0.1;
dis = zeros(iter_max+1,1);
for iter = 1:iter_max

    x(:,iter+1) =1/varrho *  ((2*A + 1/varrho * eye(N))\(y(:,iter) - varrho * u(:,iter)));
    % x(:,iter+1) =1/varrho * (2*A + 1/varrho * eye(N)) \  (y(:,iter) - u(:,iter));
    y(:,iter+1) = exp(1i * angle(x(:,iter+1) + varrho * u(:,iter)));
    % y(:,iter+1) = exp(1i * angle(x(:,iter+1) + u(:,iter)));
    eta = 0.99 * norm(x(:,iter) - y(:,iter),'inf');

    if (norm(x(:,iter + 1) - y(:,iter+1),'inf') <= eta)
        u(:,iter+1) = u(:,iter) + 1/varrho * (x(:,iter + 1) - y(:,iter+1));
        % u(:,iter+1) = u(:,iter) +   (x(:,iter + 1) - y(:,iter+1));
        
    else
        u(:,iter+1) = u(:,iter);
        varrho = varrho * c;


    end
    dis(iter + 1) = norm(x(:,iter + 1) - y(:,iter+1),'inf') ;
    if (dis(iter + 1) <= gamma)
        break;

    end
    disp(iter)
end
 
figure
plot(dis(2:iter + 1),'b-','linewidth',1);
box on
grid on
xlabel('Iterations')
ylabel('Dual gap')
abs(x(:,iter + 1))