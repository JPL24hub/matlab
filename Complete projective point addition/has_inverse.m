%return 1 if has inverse, p must be a prime number
function[yes] = has_inverse(a, p)
yes = 0;
gcd_a_p = gcd(vpi(a),vpi(p));

if gcd_a_p==1
    if mod(vpi(a),vpi(2))~=0
        yes=1;
    end
end
end