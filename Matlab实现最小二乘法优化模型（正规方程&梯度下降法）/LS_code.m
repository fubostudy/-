clear size

    [m,n] = size(A);
    if n==rank(A)
        disp(['rank=',num2str(rank(A))])
        [Q, R] = householder(A);
        d = Q'*b;
        x = backsub(R, d);
        cx1 = norm(A*x-b)^2;
        disp(cx1)
        error = abs((A*x)'*(A*x - b)) / norm(A*x);
        disp(error);
    end

  function [Q, R] = householder(A)
    [m,n]=size(A);
    Q=zeros(m,n);
    R=zeros(n,n);

    R=A;
    H=eye(m);

    for k=1:n
        y=R(k:m,k);
        e=[1;zeros(m-k,1)];
        w=y+y(1,1)/abs(y(1,1))*norm(y)*e;
        v=w/norm(w);
        I=eye(m);
        I(k:m,k:m)=I(k:m,k:m)-2*v*v';
        h=I;
        R=h*R;
        H=h*H;
%       R 主元为0，退出
        if abs(R(k,k) - 0) < 1e-18
            disp("error！it is linear dependent！");
            break;
        end
    end
    R = R(1:n, :);
    Q = H';
    Q = Q(:, 1:n);
  end
  
function x = backsub(R, d)
    [~, n] = size(R);
    x = d;
    for i = n:-1:1
        x(i) = x(i)/R(i,i);
        x(1:i-1) = x(1:i-1)-x(i)*R(1:i-1,i);
    end
end