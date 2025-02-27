function decimalVal = bin2dec_custom(binaryStr)
    [numRows, numBits] = size(binaryStr);
    powers = pow2(numBits-1:-1:0);
    decimalVal = binaryStr * powers.';
    decimalVal = reshape(decimalVal, numRows, []);
 

end
