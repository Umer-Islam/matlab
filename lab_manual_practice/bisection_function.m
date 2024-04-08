function [my_root, iterations]  = bisection_function(f,a,b,min_step)
f_a = evalin ('base','fa');
iterations = 0;
step_size = abs(a-b)/2;
x = (a+b)/2;

while step_size > min_step
    f_x = f(x);
    if f_x ==0
    iterations = iterations +1;
    break
    end

    if f_a*f_b <0
        b =x;
    else 
        a = x;
    end
    x = (a+b)/2;
    step_size = abs(a-b)/2;
    iterations = iterations +1;
end
my_root = x;