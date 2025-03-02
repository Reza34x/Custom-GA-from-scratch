function binaryMat = dec2bin_custom(decimalVals, numBits)
    if nargin < 2
        numBits = floor(log2(max(decimalVals))) + 1;
    end

    binaryMat = rem(floor(decimalVals(:) * pow2(1 - numBits:0)), 2);
end
