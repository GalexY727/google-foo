/**
     * Given an array, xs, of integers,
     * return the largest possible number that can be made
     * by multiplying a subset of ints in xs
     * @param xs the array of integers
     * @return   the largest multiple of xs' subset
     */
    public static String solution(int[] xs) {
        // If there is only one solar panel in the array
        // then we need not multiply any values
        if (xs.length == 1) return Integer.toString(xs[0]);

        /* Explanation
         * The largest multiple is the product of all positive numbers
         * and the two largest negative numbers.
         * 
         * We should only do one loop to preserve memory, 
         * so we need to find the two largest negative numbers
         * if there are two
         */

        // The product of our array multiplications
        // It should be a BigInteger since we could be multiplying
        // 1000^n
        BigInteger total = new BigInteger("1");
        
        // Make a variable to track the largest negative number
        // This will be used to negate the sign change if there are an odd number of negatives
        int largestNegative = Integer.MIN_VALUE;
        
        // Make a variable to track the amount of zeros in the array
        // If there are only zeros in the array, we must return "0"
        int zeroCount = 0;
        
        // Use an enhanced for loop to capture each value
        for (int value : xs) 
        {
            
            // Track the amount of zeros in the array
            // using zeroCount, but continue;
            // to reduce time complexity
            if (value == 0) {
                zeroCount++;
                continue;
            }
            // If we are working with a negative value
            if (value < 0)
            {
                
                // Check if it is the largest encountered
                // to use in division if the end product is negative
                if (value > largestNegative) 
                {
                    largestNegative = value;
                }
            }
            // Finally, multiply the total by the current array value
            total = total.multiply(new BigInteger(value + ""));
        }   
        
        // If there are only zeros in the array, 
        // Then we must return "0"
        if (zeroCount == xs.length) 
        {
            return "0";
        }

        // If our final result is currently negative,
        // then divide it by the largest negative value
        // to make it positive, and thus as large as it can be.
        if (total.compareTo(new BigInteger("0")) < 0 && xs.length - zeroCount > 1) 
        {
            total = total.divide(new BigInteger(largestNegative + ""));
        }
        // If it is STILL less than zero
        // then the array must look something along the lines of
        // {0, -5, 0}
        if (total.compareTo(new BigInteger("0")) < 0) {
            return "0";
        }   
        // Return the total in string format
        // since it can be very large.
        return total.toString();
    }