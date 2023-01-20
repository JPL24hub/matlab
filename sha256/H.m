function[c] = H(i)

    k={'6a09e667','0xbb67ae85','0x3c6ef372','0xa54ff53a','0x510e527f','0x9b05688c','0x1f83d9ab','0x5be0cd19'};
    c = hexToBinaryVector(k(i));
    if (length(c) < 32)
        z = zeros(1,32);
        z(1,(33-length(c)):32) = c;
        c=z;
    end
end