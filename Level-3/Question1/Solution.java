import java.util.ArrayList;

public class Solution {
    public static int[] solution(int[][] m) {
        // Using the Absorbing Markov Chain algorithm, we first must
        // re-order the matrix into standard form, where the absorbing
        // states are at the bottom right of the matrix. 
        // In this case the absorbing states are the terminal states
  
        // Find the terminal states in the matrix, and re-order columns to reflect the change
        ArrayList<Integer> terminalStates = findTerminalStates(m);

        // Re-order the matrix in standard form
        fraction[][] standardizedMatrix = toStandardForm(m, terminalStates);

        // Which requires the standardized matrix and the identity matrix
        fraction[][] F = getFundamentalMatrix(standardizedMatrix, getIdentityMatrix(terminalStates.size()));
  
        // Compute the limiting matrix 
        // i.e. F * R
        F = multiplyMatrices(F, getR(standardizedMatrix, terminalStates.size()));

        return null;
    }

    /**
     * Finds all termainal states in a matrix
     * A terminal state is one that has no exits
     * @param matrix
     */
    private static ArrayList<Integer> findTerminalStates(int[][] matrix) {

        ArrayList<Integer> terminalStates = new ArrayList<Integer>();

        // Iterate backwards so that the top of the stack is the 
        // smallest terminal state
        for (int i = 0; i < matrix.length; i++) {
            boolean isTerminal = true;
            for (int num : matrix[i]) {
                if (num != 0) {
                    isTerminal = false;
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

            int[] temp = matrix[i];
            matrix[i] = matrix[terminalStates.get(i)];
            matrix[terminalStates.get(i)] = temp;

        }

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

        // (I - Q)^-1
        return inverseMatrix(F);
    }

    private static fraction[][] getQ(fraction[][] fullMatrix, int terminalStateCount) {

        fraction[][] Q = new fraction[fullMatrix.length-terminalStateCount][fullMatrix.length-terminalStateCount];

        for (int i = terminalStateCount; i < fullMatrix.length; i++) {
            for (int j = terminalStateCount; j < fullMatrix.length; j++) {
                Q[i-terminalStateCount][j-terminalStateCount] = fullMatrix[i][j];
            }
        }

        return Q;
    }

    private static fraction[][] getR(fraction[][] fullMatrix, int terminalStateCount) {
  
        fraction[][] R = new fraction[fullMatrix.length-terminalStateCount][terminalStateCount];

        for (int i = terminalStateCount; i < fullMatrix.length; i++) {
            for (int j = 0; j < terminalStateCount; j++) {
                R[i-terminalStateCount][j] = fullMatrix[i][j];
            }
        }

        return R;
    }

    private static fraction[][] subtractMatrices(fraction[][] matrix1, fraction[][] matrix2) {
        fraction[][] subtractedMatrix = new fraction[matrix1.length][matrix1.length];

        for (int i = 0; i < matrix2.length; i++) {
            for (int j = 0; j < matrix2.length; j++) {
                subtractedMatrix[i][j] = new fraction(matrix1[i][j].getNumerator() - matrix2[i][j].getNumerator(), 1);
            }
        }

        return subtractedMatrix;
    }

    private static fraction[][] inverseMatrix(fraction[][] matrix) {
        // Using row operations, we can find the inverse of a matrix
        GaussElimination gaussElimination = new GaussElimination(matrix);
        return gaussElimination.inverse();
    }

    private static fraction[][] multiplyMatrices(fraction[][] matrix1, fraction[][] matrix2) {
        fraction[][] multipliedMatrix = new fraction[matrix1.length][matrix1.length];

        for (int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix1.length; j++) {

                int numerator = 0;
                int denominator = 0;

                for (int k = 0; k < matrix1.length; k++) {
                    numerator += matrix1[i][k].getNumerator() * matrix2[k][j].getNumerator();
                    denominator += matrix1[i][k].getDenominator() * matrix2[k][j].getDenominator();
                }

                multipliedMatrix[i][j] = new fraction(numerator, denominator);
            }
        }

        return multipliedMatrix;
    }


    public static class fraction {
        int numerator;
        int denominator;
        public fraction (int numerator, int denominator) {
            this.numerator = numerator;
            this.denominator = denominator;
        }

        public int getNumerator() {
            return this.numerator;
        }

        public int getDenominator() {
            return this.denominator;
        }

        private int gcd(int a, int b) {
            if(b == 0){
                return a;
            }
            return gcd(b, a%b);
        }

        public void simplify() {
            int gcd = gcd(numerator, denominator);
            numerator = numerator/gcd;
            denominator = denominator/gcd;
        }

        public void setNumerator(int numerator) {
            this.numerator = numerator;
        }

        public void setDenominator(int denominator) {
            this.denominator = denominator;
        }

        public void convertDenominator (int newDenominator) {
            int factor = newDenominator/denominator;
            numerator = numerator*factor;
            denominator = denominator*factor;
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

        public fraction inverse() {
            return new fraction(denominator, numerator);
        }
    }

    public static void main(String[] args) {
        int[][] matrix = {{0, 1, 0, 0, 0, 1}, 
                          {4, 0, 0, 3, 2, 0}, 
                          {0, 0, 0, 0, 0, 0}, 
                          {0, 0, 0, 0, 0, 0}, 
                          {0, 0, 0, 0, 0, 0}, 
                          {0, 0, 0, 0, 0, 0}};
        solution(matrix);
    }

    public static class GaussElimination {                                                            
        private fraction[][] m;
        int row;
        int column;

        public GaussElimination(fraction[][] A){
  

            row = A.length;
            column = A[0].length;

            m = new fraction[row][2 * column];
            for(int i = 0; i < row; i++){

                for(int j = 0; j < column; j++){
 
                    m[i][j] = A[i][j];
                }
            }
            for(int i = 0; i < row; i++){ 

                m[i][column + i] = new fraction(1, 1);

            }
 
        }


 
        /** ifZero checks if the element in position (c,c) is zero, if so swap rows until the element is not zero
             * @param c - the element to be checked.*/
        private void ifZero(int c){
  
            int i;
            boolean processing = m[c][c].getNumerator() == 0;
            while(processing){
                for(i = c+1; i < m.length; i++){
 
                    swapRows(m, c, i);
                    if (!(m[c][c].getNumerator() == 0)){
  
                        i = m.length;
                        processing = false;
 
                    }
                }
  
                if(i == m.length && (m[c][c].getNumerator() == 0)){
 
                    processing = false;

                }

            }
  
        }

        /** swapRows swaps the rows of the matrix m.
             * @param m - the matrix.
             * @param r1 - the first row.
             * @param r2 - the second row.*/
        private void swapRows(fraction[][] m, int r1, int r2){

                fraction[] temp = m[r1];
                m[r1] = m[r2];
                m[r2] = temp;

        }
  
        /** diagonalize manipulates with rows of the matrix m so that 
             all the elements outside the diagonal become zero.*/
        private void diagonalize(){
  
            for(int i = 0; i < m.length; i++){ 

                ifZero(i);  
                for(int j = 0; j < m.length; j++){ 
                    if (i != j && m[i][i].getNumerator() != 0){ 
                    fraction d = (m[j][i].divide(m[i][i]));

                    for(int a = 0; a < m[0].length; a++) {
                        fraction value = m[i][a].multiply(d);
                        m[j][a] = m[j][a].add(value);
                    }  
                    }

                }

            }
  
        }
  
        /** divideByDiagonal turns the elements of the diagonal into 1's.*/
        private void divideByDiagonal(){
  
            for(int i = 0; i < m.length; i++) {
                if(m[i][i].getNumerator() != 0) {
                    for (int a = 0; a < m.length; a++) {
                        m[i][a] = m[i][i].inverse();
                    }
                }
            }  
        }
  
        /** inverse calculates and returns the inverse of the matrix m using the upper methods.*/
        public fraction[][] inverse(){

            diagonalize();
            divideByDiagonal();

            //the inverse matrix
            fraction[][] inverse_m = new fraction[m.length][m[0].length]; 
            for(int i = 0; i < inverse_m .length; i++){
 
                for(int j = 0; j <inverse_m .length; j++){
 
                    inverse_m [i][j] = m[i][j + inverse_m .length];
  
                }
 
            }

            return inverse_m;
        }
    }
}                                                                                                                                                                                                                                                                       