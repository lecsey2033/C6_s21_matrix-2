#include "s21_matrix.h"

// main functions
int s21_create_matrix(int rows, int columns, matrix_t *result) {
    int status = 0;

    result->rows = rows;
    result->columns = columns;

    if (rows < 1 || columns < 1) {
        status = 1;
        result->matrix = NULL;
    } else {
        result->matrix = (double **)calloc(rows, sizeof(double *));
        for (int i = 0; i < rows; i++) {
            result->matrix[i] = (double *)calloc(columns, sizeof(double));
        }
    }
    return status;
}

void s21_remove_matrix(matrix_t *A) {
    if (A->rows && A->columns) {
        for (int i = 0; i < A->rows; i++) {
        free(A->matrix[i]);
        }
        free(A->matrix);
        A->rows = 0;
        A->columns = 0;
    }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int status = SUCCESS;
    if (incorrect_matrix(A) || incorrect_matrix(B)) {
        status = FAILURE;
    } else if (A->rows != B->rows || A->columns != B->columns) {
        status = FAILURE;
    } else {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= 1e-7) {
                    status = FAILURE;
                    break;
                }
            }
        }
    }
    return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int status;

    status = sum_or_sub(A, B, result, 1);
    return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int status;

    status = sum_or_sub(A, B, result, -1);
    return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int status = 0;

    if (incorrect_matrix(A)) {
        status = 1;
    } else {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = number * A->matrix[i][j];
            }
        }
    }
    return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int status = 0;

    if (incorrect_matrix(A) || incorrect_matrix(B)) {
        status = 1;
    } else if (A->columns != B->rows) {
        status = 2;
    } else {
        s21_create_matrix(A->rows, B->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < B->columns; j++) {
                for (int k = 0; k < A->columns; k++) {
                    result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
                }
            }
        }
    }
    return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
    int status = 0;

    if (incorrect_matrix(A)) {
        status = 1;
    } else {
        s21_create_matrix(A->columns, A->rows, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[j][i] = A->matrix[i][j];
            }
        }
    }
    return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
    int status = 0;

    if (incorrect_matrix(A)) {
        status = 1;
    } else if (A->columns != A->rows || A->columns == 1) {
        status = 2;
    } else {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = find_of_minor(A, i, j) * pow((-1), (i + j));
            }
        }
    }
    return status;
}

int s21_determinant(matrix_t *A, double *result) {
    int status = 0;

    if (incorrect_matrix(A)) {
        status = 1;
    } else if (A->columns != A->rows) {
        status = 2;
    } else {
        if (A->rows == 2) {
            *result = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
        } else if (A->rows == 1) {
            *result = A->matrix[0][0];
        } else {
            for (int n = 0; n < A->columns; n++) {
                *result += A->matrix[0][n] * find_of_minor(A, 0, n) * pow((-1), n);
            }
        }
    }
    return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int status = 0;
    double det = 0.0;
    matrix_t tmp_matrix, trans;

    s21_determinant(A, &det);
    if (incorrect_matrix(A)) {
        status = 1;
    } else if (A->columns != A->rows || det == 0.0) {
        status = 2;
    } else {
        if (A->rows == 1 && A->matrix[0][0] != 0) {
            result->matrix[0][0] = 1 / A->matrix[0][0];
        } else {
            s21_calc_complements(A, &tmp_matrix);
            s21_transpose(&tmp_matrix, &trans);
            s21_mult_number(&trans, 1.0 / det, result);
            s21_remove_matrix(&tmp_matrix);
            s21_remove_matrix(&trans);
        }
    }
    return status;
}

// helper functions
int incorrect_matrix(matrix_t *m) {
    int status = 0;
    if (m == NULL || m->matrix == NULL || m->rows < 1 || m->columns < 1) {
        status = 1;
    }
    return status;
}

int sum_or_sub(matrix_t *A, matrix_t *B, matrix_t *result, int sum_or_sub) {
    int status = 0;
    if (incorrect_matrix(A) || incorrect_matrix(B)) {
        status = 1;
    } else if (A->rows != B->rows || A->columns != B->columns) {
        status = 2;
    } else {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = A->matrix[i][j] + sum_or_sub * B->matrix[i][j];
            }
        }
    }
    return status;
}

double find_of_minor(matrix_t *A, int i, int j) {
    double minor = 0;
    int m, n, i_copy, j_copy;
    matrix_t result;
    s21_create_matrix(A->rows - 1, A->columns - 1, &result);
    for (i_copy = 0, m = 0; m < result.rows; m++, i_copy++) {
        for (j_copy = 0, n = 0; n < result.columns; n++, j_copy++) {
            if (i_copy == i) {
                    i_copy++;
            }
            if (j_copy == j) {
                    j_copy++;
            }
            result.matrix[m][n] = A->matrix[i_copy][j_copy];
        }
    }
    if (m == 2 && n == 2) {
        minor = result.matrix[0][0] * result.matrix[1][1] - result.matrix[0][1] * result.matrix[1][0];
    } else if (result.rows == 1) {
        minor = result.matrix[0][0];
    } else {
        for (n = 0; n < result.columns; n++) {
            minor += result.matrix[0][n] * find_of_minor(&result, 0, n) * pow((-1), n);
        }
    }
    s21_remove_matrix(&result);
    return minor;
}
