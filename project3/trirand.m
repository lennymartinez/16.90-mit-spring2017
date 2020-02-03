
function x = trirand(a, b, c)
% draw a sample from a triangular distribution with min=a, mode=b, max=c
a1 = a;
b1 = c;
c1 = b;

z = rand;

area = (c1-a1)/(b1-a1);

if z < area
    x = a1 + sqrt(z*(b1-a1)*(c1-a1));
elseif z > area
    x = b1 - sqrt((1-z)*(b1-a1)*(b1-c1));
else
    x = c1;
end

end