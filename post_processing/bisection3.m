function p = bisection3(a,b,tol,N,func)
fa = func(a);
i = 0;
while i <= N;
    p = a+(b-a)/2;
    fp = func(p);
    if abs(fp) < tol
        break;
    end
    i = i+1;
    if fa * fp > 0
        a = p;
        fa = fp;
    else
        b = p;
    end
end
if i == N;
    disp('Method failed to find a solution after max number of iterations')
end