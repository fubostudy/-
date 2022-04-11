clear size

    gradient = [];
    threshold = 1e-18;
    error = [];
    cx = [];
    [m, n] = size(A);
    x = zeros(n, 1);
    cnt = 1;
    x1_list=[];
    while true
        r = abs((A*x)'*(A*x - b)) / norm(A*x);
        error = [error, r];      
        cx = [cx, norm(A*x-b)^2];
        p = A'*(A*x - b);
        gradient = [gradient, p];
        alpha = (norm(p)^2) / norm(A*p)^2;
        x_iteration = x - alpha*p;
        x1_list=[x1_list;x];
        delta = abs((norm(A*x_iteration-b)^2 - norm(A*x-b)^2) / norm(A*x-b)^2);
        if delta <= threshold
            x = x_iteration;
            break
        else
            x = x_iteration;
        end
        cnt = cnt + 1;
    end
    
    output = reshape(x1_list,n,cnt);
    output = output';
    plot(output, '-b')
 
    i = 1:1:cnt;
    figure(1);
    plot(i, cx, '-b');
    disp(error(cnt));

    disp('cnt=');
    disp(cnt);
    disp('cx = ')
    disp(cx(1,cnt));
    disp('cx_GD - cx_LS = ')
    disp(cx(1,cnt) - cx1);
    disp('误差Proj(Ax,Ax-b) = ')
    disp(error(1,cnt));
    


