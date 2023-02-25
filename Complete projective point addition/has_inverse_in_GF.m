function has_inverse = has_inverse_in_GF(x, P)
    % Check if x is odd
    if mod(x, 2) == 0
        has_inverse = false;
        return
    end
    
    % Check if x is relatively prime to 2^n
    if gcd(x, P) == 1
        has_inverse = true;
    else
        has_inverse = false;
    end
end
