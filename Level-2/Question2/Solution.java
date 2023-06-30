public class Solution {
    public static String[] solution(String[] l) {
        /*
         * Explanation:
         * We need to sort the array, but not normally
         * We need to sort it by major, minor, and revision
         * Given that no value is less than a value of 0,
         * we can treat it as if it were -1
         */
        // Use the quick sort algorithm to sort the array
        // Do it first by major, then minor, then reivison
        // This will first sort the majors
        quickSort(l, 0, l.length - 1, 0);

        // Then sort the minors and revisions at the same time
        subVersionSort(l, 0,  l.length - 1, 1);

        return l;
    }

    private static void subVersionSort(String[] l, int startIndex, int finalIndex, int digit) {
        // Look at the element and see when it increments
        int endIndex = startIndex;
        for (int currentIndex = startIndex; currentIndex < finalIndex; currentIndex++) {
            // Split the string by periods, to get our major, minor, or revision -> "currentDigit"

            String[] currentDigitArray = l[currentIndex].split("\\.");
            String[] nextDigitArray = l[currentIndex + 1].split("\\.");
            // If the current digit does not exist, then we can assume it to be -1
            int currentIndexDigit = -1;
            int nextIndexDigit = -1;

            // If digits exist, get them out of the string array
            // Recall that digit is on the same level as the length of the array
            // Because we are looking at one level above the current digit
            if (currentDigitArray.length >= digit) 
                currentIndexDigit = Integer.parseInt(currentDigitArray[digit-1]);
            if (nextDigitArray.length >= digit)
                nextIndexDigit = Integer.parseInt(nextDigitArray[digit-1]);
            
            // If the current digit is less than the next digit, then we have found the end of the scope
            if ((currentIndexDigit < nextIndexDigit)) {
                endIndex = currentIndex;
                // If there is only one value in the scope, 
                // then we can skip it
                quickSort(l, startIndex, endIndex, digit);
                if (endIndex - startIndex > 0 && digit < 2) {
                    subVersionSort(l, startIndex, endIndex, digit + 1);
                }
                startIndex = currentIndex + 1;
            }
        }
        // We need to run one final sort to make sure
        // that there are no edge cases on the end of the array
        quickSort(l, startIndex, finalIndex, digit);
        if (l.length - 1 - startIndex > 0 && digit < 2) {
            subVersionSort(l, startIndex, finalIndex, digit + 1);
        }

    }

    /**
     * Swaps the elements in the array, first element as pivot
     * @param array The array to utilize in swapping
     * @param first the first element to swap
     * @param second the second element to swap
     */
    private static void swapElements(String[] array, int first, int second) 
    {
        String temp = array[first];
        array[first] = array[second];
        array[second] = temp;
    }

	/**
     * Find the partition for the quick sort algorithm
     * @param arr The array to partition
     * @param startIndex The starting index
     * @param endIndex The ending index
     * @return The partition index
     */
	private static int partition(String[] arr, int startIndex, int endIndex, int currentDigit)
	{
        // just use the last element as the pivot
        // if the digit does not exist, then we can treat it as if it were -1
        int pivotValue = -1;
        String[] pivotDigits = arr[endIndex].split("\\.");
        // If the digit exists in the pivot value array, assign the special pivot value
        if (!(pivotDigits.length - 1 < currentDigit)) {
            pivotValue = Integer.parseInt(pivotDigits[currentDigit]);
        }
        // loop through the array and swap elements
        // if the current element is less than the pivot
        int partitionIndex = startIndex;
        
        for (int currentIndex = startIndex; currentIndex < endIndex; currentIndex++) 
        {
            // Get the amount of digits in the current element
            String[] digitArray = arr[currentIndex].split("\\.");
            // If the digit does not exist, then we need can treat it as if it were -1
            if (digitArray.length - 1 < currentDigit) 
            {
                swapElements(arr, partitionIndex++, currentIndex);
                continue;
            }
            // Split the string by periods, to get our major, minor, or revision -> "currentDigit"
            int currentIndexDigit = Integer.parseInt(digitArray[currentDigit]);

            if (currentIndexDigit < pivotValue) 
            {
                swapElements(arr, partitionIndex++, currentIndex);
            }
        }
        swapElements(arr, partitionIndex, endIndex);
        return partitionIndex;
	}

	/**
     * Sorts the array using the quick sort algorithm
     * @param l The list to sort
     * @param low The starting index
     * @param high The ending index
     */
	private static void quickSort(String[] l, int low, int high, int digit)
	{
		//check if values are given correctly and then move ahead
		if (low < high) {

			int partitionIndex = partition(l, low, high, digit);

			// Separately sort elements before
			// partition and after partition
			quickSort(l, low, partitionIndex - 1, digit);
			quickSort(l, partitionIndex + 1, high, digit);
		}
	}
    public static void main(String[] args) {
        
        for (int i = 0; i < 10; i++) {
            String[] test = {"1.11", "2.0.0", "1.2", "2", "0.1", "1.2.1", "1.1.1", "2.0"};
            String[] result = solution(test);
            
            for (String s : result) {
                System.out.print(s + ", ");
            }
            System.out.println();
        }
    }
}