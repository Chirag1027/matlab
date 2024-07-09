#include <stdio.h>
#include <stdlib.h>

// Function to add two matrices
void add(int **A, int **B, int **C, int size) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            C[i][j] = A[i][j] + B[i][j];
}

// Function to subtract two matrices
void subtract(int **A, int **B, int **C, int size) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            C[i][j] = A[i][j] - B[i][j];
}

// Function to implement Strassen's Algorithm
void strassen(int **A, int **B, int **C, int size) {
    if (size == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }

    int newSize = size / 2;
    int **A11, **A12, **A21, **A22;
    int **B11, **B12, **B21, **B22;
    int **C11, **C12, **C21, **C22;
    int **M1, **M2, **M3, **M4, **M5, **M6, **M7;
    int **tempA, **tempB;

    // Allocate memory for submatrices
    A11 = (int **)malloc(newSize * sizeof(int *));
    A12 = (int **)malloc(newSize * sizeof(int *));
    A21 = (int **)malloc(newSize * sizeof(int *));
    A22 = (int **)malloc(newSize * sizeof(int *));
    B11 = (int **)malloc(newSize * sizeof(int *));
    B12 = (int **)malloc(newSize * sizeof(int *));
    B21 = (int **)malloc(newSize * sizeof(int *));
    B22 = (int **)malloc(newSize * sizeof(int *));
    C11 = (int **)malloc(newSize * sizeof(int *));
    C12 = (int **)malloc(newSize * sizeof(int *));
    C21 = (int **)malloc(newSize * sizeof(int *));
    C22 = (int **)malloc(newSize * sizeof(int *));
    M1 = (int **)malloc(newSize * sizeof(int *));
    M2 = (int **)malloc(newSize * sizeof(int *));
    M3 = (int **)malloc(newSize * sizeof(int *));
    M4 = (int **)malloc(newSize * sizeof(int *));
    M5 = (int **)malloc(newSize * sizeof(int *));
    M6 = (int **)malloc(newSize * sizeof(int *));
    M7 = (int **)malloc(newSize * sizeof(int *));
    tempA = (int **)malloc(newSize * sizeof(int *));
    tempB = (int **)malloc(newSize * sizeof(int *));
    for (int i = 0; i < newSize; i++) {
        A11[i] = (int *)malloc(newSize * sizeof(int));
        A12[i] = (int *)malloc(newSize * sizeof(int));
        A21[i] = (int *)malloc(newSize * sizeof(int));
        A22[i] = (int *)malloc(newSize * sizeof(int));
        B11[i] = (int *)malloc(newSize * sizeof(int));
        B12[i] = (int *)malloc(newSize * sizeof(int));
        B21[i] = (int *)malloc(newSize * sizeof(int));
        B22[i] = (int *)malloc(newSize * sizeof(int));
        C11[i] = (int *)malloc(newSize * sizeof(int));
        C12[i] = (int *)malloc(newSize * sizeof(int));
        C21[i] = (int *)malloc(newSize * sizeof(int));
        C22[i] = (int *)malloc(newSize * sizeof(int));
        M1[i] = (int *)malloc(newSize * sizeof(int));
        M2[i] = (int *)malloc(newSize * sizeof(int));
        M3[i] = (int *)malloc(newSize * sizeof(int));
        M4[i] = (int *)malloc(newSize * sizeof(int));
        M5[i] = (int *)malloc(newSize * sizeof(int));
        M6[i] = (int *)malloc(newSize * sizeof(int));
        M7[i] = (int *)malloc(newSize * sizeof(int));
        tempA[i] = (int *)malloc(newSize * sizeof(int));
        tempB[i] = (int *)malloc(newSize * sizeof(int));
    }

    // Dividing the matrices into sub-matrices
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + newSize];
            A21[i][j] = A[i + newSize][j];
            A22[i][j] = A[i + newSize][j + newSize];
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + newSize];
            B21[i][j] = B[i + newSize][j];
            B22[i][j] = B[i + newSize][j + newSize];
        }
    }

    // Calculating M1 to M7
    add(A11, A22, tempA, newSize); // tempA = A11 + A22
    add(B11, B22, tempB, newSize); // tempB = B11 + B22
    strassen(tempA, tempB, M1, newSize); // M1 = (A11 + A22) * (B11 + B22)

    add(A21, A22, tempA, newSize); // tempA = A21 + A22
    strassen(tempA, B11, M2, newSize); // M2 = (A21 + A22) * B11

    subtract(B12, B22, tempB, newSize); // tempB = B12 - B22
    strassen(A11, tempB, M3, newSize); // M3 = A11 * (B12 - B22)

    subtract(B21, B11, tempB, newSize); // tempB = B21 - B11
    strassen(A22, tempB, M4, newSize); // M4 = A22 * (B21 - B11)

    add(A11, A12, tempA, newSize); // tempA = A11 + A12
    strassen(tempA, B22, M5, newSize); // M5 = (A11 + A12) * B22

    subtract(A21, A11, tempA, newSize); // tempA = A21 - A11
    add(B11, B12, tempB, newSize); // tempB = B11 + B12
    strassen(tempA, tempB, M6, newSize); // M6 = (A21 - A11) * (B11 + B12)

    subtract(A12, A22, tempA, newSize); // tempA = A12 - A22
    add(B21, B22, tempB, newSize); // tempB = B21 + B22
    strassen(tempA, tempB, M7, newSize); // M7 = (A12 - A22) * (B21 + B22)

    // Calculating C11, C12, C21, C22
    add(M1, M4, tempA, newSize); // tempA = M1 + M4
    subtract(tempA, M5, tempB, newSize); // tempB = M1 + M4 - M5
    add(tempB, M7, C11, newSize); // C11 = M1 + M4 - M5 + M7

    add(M3, M5, C12, newSize); // C12 = M3 + M5

    add(M2, M4, C21, newSize); // C21 = M2 + M4

    add(M1, M3, tempA, newSize); // tempA = M1 + M3
    subtract(tempA, M2, tempB, newSize); // tempB = M1 + M3 - M2
    add(tempB, M6, C22, newSize); // C22 = M1 + M3 - M2 + M6

    // Combining the results into one matrix
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            C[i][j] = C11[i][j];
            C[i][j + newSize] = C12[i][j];
            C[i + newSize][j] = C21[i][j];
            C[i + newSize][j + newSize] = C22[i][j];
        }
    }

    // Free allocated memory
    for (int i = 0; i < newSize; i++) {
        free(A11[i]);
        free(A12[i]);
        free(A21[i]);
        free(A22[i]);
        free(B11[i]);
        free(B12[i]);
        free(B21[i]);
        free(B22[i]);
        free(C11[i]);
        free(C12[i]);
        free(C21[i]);
        free(C22[i]);
        free(M1[i]);
        free(M2[i]);
        free(M3[i]);
        free(M4[i]);
        free(M5[i]);
        free(M6[i]);
        free(M7[i]);
        free(tempA[i]);
        free(tempB[i]);
    }
    free(A11); free(A12); free(A21); free(A22);
    free(B11); free(B12); free(B21); free(B22);
    free(C11); free(C12); free(C21); free(C22);
    free(M1); free(M2); free(M3); free(M4); free(M5); free(M6); free(M7);
    free(tempA); free(tempB);
}

// Helper function to allocate a 2D matrix
int **allocateMatrix(int size) {
    int **matrix = (int **)malloc(size * sizeof(int *));
    for (int i = 0; i < size; i++) {
        matrix[i] = (int *)malloc(size * sizeof(int));
    }
    return matrix;
}

// Helper function to free a 2D matrix
void freeMatrix(int **matrix, int size) {
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Helper function to print a matrix
void printMatrix(int **matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Main function to test Strassen's Algorithm
int main() {
    int size = 4;
    int **A = allocateMatrix(size);
    int **B = allocateMatrix(size);
    int **C = allocateMatrix(size);

    // Initialize matrices A and B
    int counter = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A[i][j] = counter++;
            B[i][j] = counter++;
        }
    }

    // Print matrices A and B
    printf("Matrix A:\n");
    printMatrix(A, size);
    printf("\nMatrix B:\n");
    printMatrix(B, size);

    // Perform Strassen's multiplication
    strassen(A, B, C, size);

    // Print the result matrix C
    printf("\nMatrix C (Result):\n");
    printMatrix(C, size);

    // Free allocated memory
    freeMatrix(A, size);
    freeMatrix(B, size);
    freeMatrix(C, size);

    return 0;
}
