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

        // Re-order the matrix in standard form
        fraction[][] standardizedMatrix = toStandardForm(m, terminalStates);

        // Which requires the standardized matrix and the identity matrix
        fraction[][] F = getFundamentalMatrix(standardizedMatrix, getIdentityMatrix(terminalStates.size()));

        // Compute the limiting matrix
        // i.e. F * R

        return finalArray(multiplyMatrices(F, getR(standardizedMatrix, terminalStates.size())));
    }

    private static int[] finalArray(fraction[][] FR) {
        // The final array is the first row of the FR matrix
        // We need to find the lowest common multiple of the denominators

        int[] finalArray = new int[FR[firstIndex].length+1];

        finalArray[FR[firstIndex].length] = 1;
        for (int i = 0; i < FR[firstIndex].length; i++) {
            finalArray[FR[firstIndex].length] = lcm(finalArray[FR[firstIndex].length], FR[firstIndex][i].getDenominator());
        }

        for (int i = 0; i < FR[firstIndex].length; i++) {
            finalArray[i] = FR[firstIndex][i].getNumerator() * (finalArray[FR[firstIndex].length] / FR[firstIndex][i].getDenominator());
        }

        return finalArray;
    }

    private static int lcm(int a, int b) {
        return a * (b / gcd(a, b));
    }

    private static int gcd(int a, int b) {
        if (b == 0) {
            return a;
        }
        return gcd(b, a % b);
    }

    /**
     * Finds all termainal states in a matrix
     * A terminal state is one that has no exits
     * 
     * @param matrix
     */
    private static ArrayList<Integer> findTerminalStates(int[][] matrix) {

        ArrayList<Integer> terminalStates = new ArrayList<Integer>();

        // Iterate backwards so that the top of the stack is the
        // smallest terminal state
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
     * 
     * Since we are standardizing this matrix,
     * we will also be converting the matrix to fractions
     * for easier computation later on.
     * 
     * @param matrix
     * @param terminalStates
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
            matrix[terminalStates.get(i)][0] = terminalStates.get(i);
            matrix[i] = matrix[terminalStates.get(i)];
            matrix[terminalStates.get(i)] = temp;

        }
        firstIndex -= terminalStates.size();
        // Now, the columns
        // This is done by swapping the columns of the terminal states
        // with the highest non-terminal column
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
            int denominator = 0;
            for (int j = 0; j < matrix.length; j++) {
                denominator += matrix[i][j];
            }
            if (denominator == 0) {
                denominator = 1;
            }
            for (int j = 0; j < matrix.length; j++) {
                fractionMatrix[i][j] = new fraction(matrix[i][j], denominator);
            }
        }

        return fractionMatrix;
    }

    private static fraction[][] getIdentityMatrix(int size) {
        fraction[][] identityMatrix = new fraction[size][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    identityMatrix[i][j] = new fraction(1, 1);
                } else {
                    identityMatrix[i][j] = new fraction(0, 1);
                }
            }
        }

        return identityMatrix;
    }

    private static fraction[][] getFundamentalMatrix(fraction[][] fullMatrix, fraction[][] I) {

        // Matrix F = (I - Q)^-1
        // Where I is the identity matrix
        // Q is the transition matrix
        // Q starts at fullMatrix(I.lenth, I.length),
        // and ends at fullMatrix(fullMatrix.length, fullMatrix.length)

        // F will temporarily be the Q, until F is fully computed.
        // Since they are the same size
        fraction[][] F = getQ(fullMatrix, I.length);

        // I - Q
        // Recall that Q if currently F.
        F = subtractMatrices(I, F);

        return invert(F);
    }

    private static fraction[][] getQ(fraction[][] fullMatrix, int terminalStateCount) {

        fraction[][] Q = new fraction[fullMatrix.length - terminalStateCount][fullMatrix.length - terminalStateCount];

        for (int i = terminalStateCount; i < fullMatrix.length; i++) {
            for (int j = terminalStateCount; j < fullMatrix.length; j++) {
                Q[i - terminalStateCount][j - terminalStateCount] = fullMatrix[i][j];
            }
        }

        return Q;
    }

    private static fraction[][] getR(fraction[][] fullMatrix, int terminalStateCount) {

        fraction[][] R = new fraction[fullMatrix.length - terminalStateCount][terminalStateCount];

        for (int i = terminalStateCount; i < fullMatrix.length; i++) {
            for (int j = 0; j < terminalStateCount; j++) {
                R[i - terminalStateCount][j] = fullMatrix[i][j];
            }
        }
        return R;
    }

    private static fraction[][] subtractMatrices(fraction[][] matrix1, fraction[][] matrix2) {
        fraction[][] subtractedMatrix = new fraction[matrix2.length][matrix2.length];

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
                fraction sum = new fraction(0, 1);
                for (int k = 0; k < matrix1[0].length; k++) {
                    sum = sum.add(matrix1[i][k].multiply(matrix2[k][j]));
                }
                multipliedMatrix[i][j] = sum;
            }
        }

        return multipliedMatrix;
    }

    public static class fraction {
        int numerator = 0;
        int denominator = 1;

        public fraction(int numerator, int denominator) {
            this.numerator = numerator;
            this.denominator = denominator;
            simplify();
        }

        public int getNumerator() {
            return this.numerator;
        }

        public int getDenominator() {
            return this.denominator;
        }

        private int gcd(int a, int b) {
            if (b == 0) {
                return a;
            }
            return gcd(b, a % b);
        }

        public void simplify() {
            int gcd = Math.abs(gcd(numerator, denominator));
            this.numerator /= gcd;
            this.denominator /= gcd;
        }

        public void setNumerator(int numerator) {
            this.numerator = numerator;
        }

        public void setDenominator(int denominator) {
            this.denominator = denominator;
        }

        public void convertDenominator(int newDenominator) {
            int factor = newDenominator / denominator;
            numerator = numerator * factor;
            denominator = denominator * factor;
        }

        public fraction divide(fraction b) {
            int numerator = this.getNumerator() * b.getDenominator();
            int denominator = this.getDenominator() * b.getNumerator();
            return new fraction(numerator, denominator);
        }

        public fraction multiply(fraction b) {
            int numerator = this.getNumerator() * b.getNumerator();
            int denominator = this.getDenominator() * b.getDenominator();
            return new fraction(numerator, denominator);
        }

        public fraction add(fraction b) {
            int numerator = this.getNumerator() * b.getDenominator() + b.getNumerator() * this.getDenominator();
            int denominator = this.getDenominator() * b.getDenominator();
            return new fraction(numerator, denominator);
        }

        public fraction subtract(fraction b) {
            int numerator = this.getNumerator() * b.getDenominator() - b.getNumerator() * this.getDenominator();
            int denominator = this.getDenominator() * b.getDenominator();
            return new fraction(numerator, denominator);
        }

        public fraction inverse() {
            return new fraction(numerator, denominator);
        }

        public fraction abs() {
            return new fraction(Math.abs(numerator), denominator);
        }

        public double getDecimal() {
            return (double) numerator / (double) denominator;
        }
    }



    public static fraction[][] invert(fraction[][] a)
    {
        int n = a.length;
        fraction[][] x = new fraction[n][n];
        fraction[][] b = new fraction[n][n];
        int[] index = new int[n];
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < n; j++) 
                if (i == j)
                    b[i][i] = new fraction(1, 1);
                else
                    b[i][j] = new fraction(0, 1);
    
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
            fraction c1 = new fraction(0, 1);
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
        int[][] matrix1 = { { 0, 1, 0, 0, 0, 1 }, { 4, 0, 0, 3, 2, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
        int[][] matrix2 = {{0, 2, 1, 0, 0}, {0, 0, 0, 3, 4}, {0, 0, 0, 0, 0}, {0, 0, 0, 0,0}, {0, 0, 0, 0, 0}};
        solution(matrix2);
    }
}

