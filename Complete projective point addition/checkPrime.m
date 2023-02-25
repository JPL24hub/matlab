function isPrime = checkPrime(number)
% This function checks if a number is prime or not
% Input: number - the number to check for primality
% Output: isPrime - a boolean indicating if the number is prime or not
    
    % Check if the number is divisible by any number from 2 to sqrt(number)
 
    i=vpi(vpi(number)-10000);
    while i<sqrt(vpi(number))
        if mod(vpi(number), vpi(i)) == 0
            isPrime = false;
            return;
        end
        i=vpi(i) + 1;
    end
    
    % If the loop completes without finding a factor, the number is prime
    isPrime = true;
end
