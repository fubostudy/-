function x = backsub(R, b)
[rm, rn] = size(R);
[bm, ~] = size(b);
if rm ~= bm
    fprintf("dimension is not the same!\n");
end
x(rn) = b(rn)/R(rn,rn); 
for j = rn-1:-1:1
x(j) = b(j) - sum((R(j, j+1:end).* x(j+1:end)));
x(j) = x(j)/R(j, j);
end
x = x(:);
end