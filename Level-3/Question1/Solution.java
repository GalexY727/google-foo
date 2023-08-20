import java.math.BigInteger;
import java.util.ArrayList;

public class Solution {
    private static int firstIndex = -1;
    public static int[] solution(int[][] m) {

        // Using the Absorbing Markov Chain algorithm, we first must
        // re-order the matrix into standard form, where the absorbing
        // states are at the bottom right of the matrix.
        // In this case the absorbing states are the terminal states

        // Find the terminal states in the matrix, and re-order columns to reflect the
        // change
        ArrayList<Integer> terminalStates = findTerminalStates(m);
        
        // If the first state is terminal, we're already done.
        if (terminalStates.get(0) == 0) {
            int[] finalArray = new int[terminalStates.size() + 1];
            finalArray[0] = 1;
            finalArray[terminalStates.size()] = 1;
            return finalArray;
        }

        // Re-order the matrix in standard form
        fraction[][] standardizedMatrix = toStandardForm(m, terminalStates);

        // Which requires the standardized matrix and the identity matrix
        fraction[][] F = getFundamentalMatrix(standardizedMatrix, terminalStates.size());

        // Compute the limiting matrix
        // i.e. F * R
        return finalArray(multiplyMatrices(F, getR(standardizedMatrix, terminalStates.size())));
    }

    /**
     * Converts the final matrix into the final array
     * @param FR the final matrix
     * @return an array of integers representing fractions 
     *   with the last index being the common denominator
     */
    private static int[] finalArray(fraction[][] FR) {
        // The final array is the row of the FR matrix that corresponds to the s0
        int[] finalArray = new int[FR[firstIndex].length+1];

        finalArray[FR[firstIndex].length] = 1;

        BigInteger lcm = FR[firstIndex][0].getDenominator();
        for (int i = 0; i < FR[firstIndex].length; i++) {
            lcm = lcm(lcm, FR[firstIndex][i].getDenominator());
        }

        finalArray[FR[firstIndex].length] = lcm.intValue();

        for (int i = 0; i < FR[firstIndex].length; i++) {
            finalArray[i] = FR[firstIndex][i].getNumerator().multiply(lcm.divide(FR[firstIndex][i].getDenominator())).intValue();
        }

        return finalArray;
    }

    private static BigInteger lcm(BigInteger a, BigInteger b) {
        return a.multiply(b.divide(gcd(a, b)));
    }

    private static BigInteger gcd(BigInteger a, BigInteger b) {
        if (b.equals(BigInteger.valueOf(0))) {
            return a;
        }
        return gcd(b, a.mod(b));
    }

    /**
     * Finds all termainal states in a matrix
     * A terminal state is one that has no exits
     * 
     * @param matrix the matrix to search
     */
    private static ArrayList<Integer> findTerminalStates(int[][] matrix) {

        ArrayList<Integer> terminalStates = new ArrayList<Integer>();

        // terminal state = state that has no exits
        for (int i = 0; i < matrix.length; i++) {
            boolean isTerminal = true;
            for (double num : matrix[i]) {
                if (num != 0) {
                    isTerminal = false;
                    if (firstIndex == -1) {
                        firstIndex = i;
                    }
                    break;
                }
            }
            if (isTerminal) {
                terminalStates.add(i);
            }
        }

        return terminalStates;

    }

    /**
     * Re-orders the matrix into standard form
     * Where terminal states (identity matricies) are top left and
     * zeros are top right
     * with the bottom left and bottom right being
     * the transition matricies
     * Since we are standardizing this matrix,
     * we will also be converting the matrix to fractions
     * for easier computation later on.
     * 
     * @param matrix the matrix to standardize
     * @param terminalStates the terminal states in the matrix
     * @return the computed matrix
     */
    private static fraction[][] toStandardForm(int[][] matrix, ArrayList<Integer> terminalStates) {

        // First we should re-order the rows
        // This is done by swapping the rows of the terminal states
        // with the highest non-terminal row
        for (int i = 0; i < terminalStates.size(); i++) {
            if (terminalStates.get(i) == i) {
                continue;
            }
            if (i == firstIndex) {
                firstIndex = terminalStates.get(i);
            }
            int[] temp = matrix[i];
            matrix[i] = matrix[terminalStates.get(i)];
            matrix[terminalStates.get(i)] = temp;
        }
        firstIndex -= terminalStates.size();

        // Now, the columns
        // This is done by swapping the columns of the terminal states
        // with the highest non-terminal column
        // this algorithm must be the same as the row swaps
        for (int i = 0; i < terminalStates.size(); i++) {
            if (terminalStates.get(i) == i) {
                continue;
            }

            for (int j = 0; j < matrix.length; j++) {
                int temp = matrix[j][i];
                matrix[j][i] = matrix[j][terminalStates.get(i)];
                matrix[j][terminalStates.get(i)] = temp;
            }
        }

        fraction[][] fractionMatrix = new fraction[matrix.length][matrix.length];
        // Now, we need to convert the matrix to fractions
        // We need to add all of the values in a row up to get our denominator
        for (int i = 0; i < matrix.length; i++) {
            BigInteger denominator = BigInteger.valueOf(0);
            for (int j = 0; j < matrix.length; j++) {
                denominator = denominator.add(BigInteger.valueOf(matrix[i][j]));
            }
            if (denominator.equals(BigInteger.valueOf(0))) {
                denominator = BigInteger.valueOf(1);
            }
            for (int j = 0; j < matrix.length; j++) {
                fractionMatrix[i][j] = new fraction(BigInteger.valueOf(matrix[i][j]), denominator);
            }
        }

        return fractionMatrix;
    }

    /**
     * Gets the identity matrix of a given size
     * 1 0 0
     * 0 1 0 <-- example of a 3x3 identity matrix
     * 0 0 1
     * @param size the size of the matrix
     * @return the identity matrix
     */
    private static fraction[][] getIdentityMatrix(int size) {
        fraction[][] identityMatrix = new fraction[size][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    identityMatrix[i][j] = new fraction(BigInteger.valueOf(1), BigInteger.valueOf(1));
                } else {
                    identityMatrix[i][j] = new fraction(BigInteger.valueOf(0), BigInteger.valueOf(1));
                }
            }
        }

        return identityMatrix;
    }

    /**
     * Gets the fundamental matrix of the full matrix
     * (I - Q)^-1
     * @param fullMatrix the full matrix
     * @param terminalStateCount the number of terminal states
     * @return
     */
    private static fraction[][] getFundamentalMatrix(fraction[][] fullMatrix, int terminalStateCount) {

        // Matrix F = (I - Q)^-1
        // Where I is the identity matrix
        // Q is the transition matrix
        // Q starts at fullMatrix(I.lenth, I.length),
        // and ends at fullMatrix(fullMatrix.length, fullMatrix.length)

        // F will temporarily be the Q, until F is fully computed.
        // Since they are the same size
        fraction[][] F = getQ(fullMatrix, terminalStateCount);

        // I - Q
        // Recall that Q if currently F.
        F = subtractMatrices(getIdentityMatrix(fullMatrix.length - terminalStateCount), F);

        return invert(F);
    }

    /**
     * Gets the Q submatrix of the full matrix
     * @param fullMatrix the full matrix
     * @param terminalStateCount the number of terminal states
     * @return
     */
    private static fraction[][] getQ(fraction[][] fullMatrix, int terminalStateCount) {

        fraction[][] Q = new fraction[fullMatrix.length - terminalStateCount][fullMatrix.length - terminalStateCount];

        for (int i = terminalStateCount; i < fullMatrix.length; i++) {
            if (fullMatrix.length - terminalStateCount >= 0)
                System.arraycopy(fullMatrix[i], terminalStateCount, Q[i - terminalStateCount], 0, fullMatrix.length - terminalStateCount);
        }

        return Q;
    }

    /**
     * Gets the R submatrix of the full matrix
     * @param fullMatrix the full matrix
     * @param terminalStateCount the number of terminal states
     * @return
     */
    private static fraction[][] getR(fraction[][] fullMatrix, int terminalStateCount) {

        fraction[][] R = new fraction[fullMatrix.length - terminalStateCount][terminalStateCount];

        for (int i = terminalStateCount; i < fullMatrix.length; i++) {
            System.arraycopy(fullMatrix[i], 0, R[i - terminalStateCount], 0, terminalStateCount);
        }
        return R;
    }

    private static fraction[][] subtractMatrices(fraction[][] matrix1, fraction[][] matrix2) {
        fraction[][] subtractedMatrix = new fraction[matrix2.length][matrix2[0].length];

        for (int i = 0; i < matrix2.length; i++) {
            for (int j = 0; j < matrix2[0].length; j++) {
                subtractedMatrix[i][j] = matrix1[i][j].subtract(matrix2[i][j]);
            }
        }

        return subtractedMatrix;
    }

    private static fraction[][] multiplyMatrices(fraction[][] matrix1, fraction[][] matrix2) {
        fraction[][] multipliedMatrix = new fraction[matrix1.length][matrix2[0].length];

        for (int i = 0; i < multipliedMatrix.length; i++) {
            for (int j = 0; j < multipliedMatrix[0].length; j++) {
                fraction sum = new fraction(BigInteger.valueOf(0), BigInteger.valueOf(1));
                for (int k = 0; k < matrix1[0].length; k++) {
                    sum = sum.add(matrix1[i][k].multiply(matrix2[k][j]));
                }
                multipliedMatrix[i][j] = sum;
            }
        }

        return multipliedMatrix;
    }

    public static class fraction {
        private BigInteger numerator = BigInteger.valueOf(0);
        private BigInteger denominator = BigInteger.valueOf(1);

        public fraction(BigInteger numerator, BigInteger denominator) {
            this.numerator = numerator;
            this.denominator = denominator;
            simplify();
        }

        public BigInteger getNumerator() {
            return this.numerator;
        }

        public BigInteger getDenominator() {
            return this.denominator;
        }

        private BigInteger gcd(BigInteger a, BigInteger b) {
            if (b.equals(BigInteger.valueOf(0))) {
                return a;
            }
            return gcd(b, a.mod(b));
        }

        public void simplify() {
            BigInteger gcd = gcd(numerator, denominator).abs();
            this.numerator = numerator.divide(gcd);
            this.denominator = denominator.divide(gcd);
        }

        public fraction divide(fraction b) {
            BigInteger numerator = this.getNumerator().multiply(b.getDenominator());
            BigInteger denominator = this.getDenominator().multiply(b.getNumerator());
            return new fraction(numerator, denominator);
        }

        public fraction multiply(fraction b) {
            BigInteger numerator = this.getNumerator().multiply(b.getNumerator());
            BigInteger denominator = this.getDenominator().multiply(b.getDenominator());
            return new fraction(numerator, denominator);
        }

        public fraction add(fraction b) {
            BigInteger numerator = this.getNumerator().multiply(b.getDenominator()).add(b.getNumerator().multiply(this.getDenominator()));
            BigInteger denominator = this.getDenominator().multiply(b.getDenominator());
            return new fraction(numerator, denominator);
        }

        public fraction subtract(fraction b) {
            BigInteger numerator = this.getNumerator().multiply(b.getDenominator()).subtract(b.getNumerator().multiply(this.getDenominator()));
            BigInteger denominator = this.getDenominator().multiply(b.getDenominator());
            return new fraction(numerator, denominator);
        }

        public fraction abs() {
            return new fraction(numerator.abs(), denominator);
        }

        public double getDecimal() {
            return numerator.doubleValue() / denominator.doubleValue();
        }
    }

    /**
     * Method to carry out the partial-pivoting Gaussian
     * @param a input matrix
     * @return the inverted matrix
     */
    public static fraction[][] invert(fraction[][] a)
    {
        int n = a.length;
        fraction[][] x = new fraction[n][n];
        fraction[][] b = new fraction[n][n];
        int[] index = new int[n];
        b = getIdentityMatrix(n);
    
        // Transform the matrix into an upper triangle
        gaussian(a, index);
    
        // Update the matrix b[i][j] with the ratios stored
        for (int i = 0; i < n-1; i++)
            for (int j = i+1; j < n; j++)
                for (int k = 0; k < n; k++) 
                    b[index[j]][k] = b[index[j]][k].subtract(a[index[j]][i].multiply(b[index[i]][k]));
    
        // Perform backward substitutions
        for (int i = 0; i < n; i++) 
        {
            x[n-1][i] = b[index[n-1]][i].divide(a[index[n-1]][n-1]);
            for (int j = n-2; j >= 0; j--) 
            {
                x[j][i] = b[index[j]][i];
                for (int k = j+1; k < n; k++) 
                {
                    x[j][i] = x[j][i].subtract(a[index[j]][k].multiply(x[k][i]));
                }
                x[j][i] = x[j][i].divide(a[index[j]][j]);
            }
        }
        return x;
    }
	    
    /**
     * Method to carry out the partial-pivoting Gaussian
     * @param a input matrix
     * @param index index
     */
    public static void gaussian(fraction[][] a, int[] index) 
    {
        int n = index.length;
        fraction[] c = new fraction[n];
    
        // Initialize the index
        for (int i = 0; i < n; i++) 
            index[i] = i;
    
        // Find the rescaling factors, one from each row
        for (int i = 0; i < n; i++) 
        {
            fraction c1 = new fraction(BigInteger.valueOf(0), BigInteger.valueOf(1));
            for (int j = 0; j < n; j++) 
            {
                fraction c0 = a[i][j].abs();
                if (c0.getDecimal() > c1.getDecimal()) c1 = c0;
            }
            c[i] = c1;
        }
        int k = 0;
        for (int j = 0; j < n-1; j++) 
        {
            double pi1 = 0;
            for (int i = j; i < n; i++) 
            {
                double pi0 = Math.abs(a[index[i]][j].getDecimal());
                pi0 /= c[index[i]].getDecimal();
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i = j+1; i < n; i++) 	
            {
                fraction pj = a[index[i]][j].divide(a[index[j]][j]);
                a[index[i]][j] = pj;
                for (int l = j+1; l < n; l++)
                    a[index[i]][l] = a[index[i]][l].subtract(pj.multiply(a[index[j]][l]));
            }
        }
    }
    public static void main(String[] args) {
        int[][] example10 = new int[][]{
                {5, 13, 0, 0, 7, 21, 1024},
                {1, 2, 3, 0, 8, 6, 0},
                {7, 5, 2, 1, 2048, 0, 9},
                {0, 0, 0, 0, 0, 0 , 0},
                {10, 9, 8, 228, 729, 630, 1},
                {15, 10000, 3, 0, 0, 0, 8},
                {0, 0, 0, 0, 0, 0 , 0}
        };
        int[] sol = solution(example10);
        for (int i : sol) {
            System.out.print(i + " ");
        }
    }
}