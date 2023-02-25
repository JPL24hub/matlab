% x has inverse if
% x is odd, and
% x is relatively prime to prime_number_p, i.e., gcd(x, prime_number_p) = 1.
clc
clear
%find prime
%get prime https://www.numberempire.com/primenumbers.php
%check prime https://bigprimes.org/primality-test
p = 205; %find immediate prime number greater than this number

a=vpi(p)-vpi(10);

while has_inverse(a, p)~=1
    a = vpi(a) + vpi(1);
end
p
a
a_inverse = binary_inversion_algorithm(a, p)
