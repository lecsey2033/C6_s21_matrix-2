#include <check.h>
#include <string.h>
#include "s21_matrix.h"

void s21_gen_matrix(matrix_t *A) {
    double num = 1;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++, num++) {
            A->matrix[i][j] = num;
        }
    }
}

void s21_anti_gen_matrix(matrix_t *A) {
    double num = 32;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++, num--) {
            A->matrix[i][j] = num;
        }
    }
}

START_TEST(test_s21_eq_matrix) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t temp3;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &temp2);
    s21_create_matrix(3, 3, &temp3);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    ck_assert_int_eq(s21_eq_matrix(&temp, &temp2), SUCCESS);
    ck_assert_int_eq(s21_eq_matrix(&temp, &temp3), FAILURE);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
    s21_remove_matrix(&temp3);
}
END_TEST

START_TEST(test_s21_eq_matrix_2) {
    matrix_t temp;
    matrix_t temp2;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    temp2.matrix[1][1] = 0;
    temp2.matrix[1][2] = 77;
    ck_assert_int_eq(s21_eq_matrix(&temp, &temp2), FAILURE);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
}
END_TEST

START_TEST(test_s21_eq_matrix_3) {
    matrix_t temp;
    matrix_t temp2;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(2, 3, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    ck_assert_int_eq(s21_eq_matrix(&temp, &temp2), FAILURE);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
}
END_TEST

START_TEST(test_s21_eq_matrix_4) {
    matrix_t temp;
    matrix_t temp2;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 21, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    ck_assert_int_eq(s21_eq_matrix(&temp, &temp2), FAILURE);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
}
END_TEST

START_TEST(test_s21_eq_matrix_5) {
    matrix_t temp;
    matrix_t temp2;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    temp.matrix[0][1] = 4.67773;
    temp2.matrix[0][1] = 4.67774;
    ck_assert_int_eq(s21_eq_matrix(&temp, &temp2), FAILURE);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
}
END_TEST

START_TEST(test_s21_sum_matrix) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    matrix_t result_defolt;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &temp2);
    s21_create_matrix(3, 3, &result_defolt);
    s21_gen_matrix(&temp);
    s21_anti_gen_matrix(&temp2);
    for (int i = 0; (i < temp.rows); i++) {
        for (int j = 0; j < temp.columns; j++) {
            result_defolt.matrix[i][j] = temp.matrix[i][j] + temp2.matrix[i][j];
        }
    }
    ck_assert_int_eq(s21_sum_matrix(&temp, &temp2, &result), 0);
    ck_assert_int_eq(s21_eq_matrix(&result_defolt, &result), SUCCESS);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
    s21_remove_matrix(&result);
    s21_remove_matrix(&result_defolt);
}
END_TEST

START_TEST(test_s21_sum_matrix_2) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 0, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    ck_assert_int_eq(s21_sum_matrix(&temp, &temp2, &result), 1);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_sub_matrix) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    matrix_t result_defolt;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &temp2);
    s21_create_matrix(3, 3, &result_defolt);
    s21_gen_matrix(&temp);
    s21_anti_gen_matrix(&temp2);
    for (int i = 0; (i < temp.rows); i++) {
        for (int j = 0; j < temp.columns; j++) {
            result_defolt.matrix[i][j] = temp.matrix[i][j] - temp2.matrix[i][j];
        }
    }
    ck_assert_int_eq(s21_sub_matrix(&temp, &temp2, &result), 0);
    ck_assert_int_eq(s21_eq_matrix(&result_defolt, &result), SUCCESS);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
    s21_remove_matrix(&result);
    s21_remove_matrix(&result_defolt);
}
END_TEST

START_TEST(test_s21_sub_matrix_2) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(0, 3, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    ck_assert_int_eq(s21_sub_matrix(&temp, &temp2, &result), 1);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_mult_number) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    matrix_t result_defolt;
    double num = 777.777;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &temp2);
    s21_create_matrix(3, 3, &result_defolt);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    for (int i = 0; (i < temp.rows); i++) {
        for (int j = 0; j < temp.columns; j++) {
            result_defolt.matrix[i][j] = temp.matrix[i][j] * 777.777;
        }
    }
    ck_assert_int_eq(s21_mult_number(&temp, num, &result), 0);
    ck_assert_int_eq(s21_eq_matrix(&result_defolt, &result), SUCCESS);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
    s21_remove_matrix(&result);
    s21_remove_matrix(&result_defolt);
}
END_TEST

START_TEST(test_s21_mult_matrix_2) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    matrix_t result_defolt;
    s21_create_matrix(3, 2, &temp);
    s21_create_matrix(2, 3, &temp2);
    s21_create_matrix(3, 3, &result_defolt);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    result_defolt.matrix[0][0] = 9;
    result_defolt.matrix[0][1] = 12;
    result_defolt.matrix[0][2] = 15;
    result_defolt.matrix[1][0] = 19;
    result_defolt.matrix[1][1] = 26;
    result_defolt.matrix[1][2] = 33;
    result_defolt.matrix[2][0] = 29;
    result_defolt.matrix[2][1] = 40;
    result_defolt.matrix[2][2] = 51;
    ck_assert_int_eq(s21_mult_matrix(&temp, &temp2, &result), 0);
    ck_assert_int_eq(s21_eq_matrix(&result_defolt, &result), SUCCESS);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
    s21_remove_matrix(&result);
    s21_remove_matrix(&result_defolt);
}
END_TEST

START_TEST(test_s21_mult_matrix_3) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(0, 0, &temp2);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    ck_assert_int_eq(s21_mult_matrix(&temp, &temp2, &result), 1);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_mult_matrix_4) {
    matrix_t temp;
    matrix_t temp2;
    matrix_t result;
    matrix_t result_defolt;
    s21_create_matrix(2, 3, &temp);
    s21_create_matrix(3, 2, &temp2);
    s21_create_matrix(2, 2, &result_defolt);
    s21_gen_matrix(&temp);
    s21_gen_matrix(&temp2);
    result_defolt.matrix[0][0] = 22;
    result_defolt.matrix[0][1] = 28;
    result_defolt.matrix[1][0] = 49;
    result_defolt.matrix[1][1] = 64;
    ck_assert_int_eq(s21_mult_matrix(&temp, &temp2, &result), 0);
    ck_assert_int_eq(s21_eq_matrix(&result_defolt, &result), SUCCESS);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&temp2);
    s21_remove_matrix(&result);
    s21_remove_matrix(&result_defolt);
}
END_TEST

START_TEST(test_s21_transpose) {
    matrix_t temp;
    matrix_t result;
    matrix_t result_defolt;
    s21_create_matrix(3, 2, &temp);
    s21_create_matrix(2, 3, &result_defolt);
    s21_gen_matrix(&temp);
    result_defolt.matrix[0][0] = 1;
    result_defolt.matrix[0][1] = 3;
    result_defolt.matrix[0][2] = 5;
    result_defolt.matrix[1][0] = 2;
    result_defolt.matrix[1][1] = 4;
    result_defolt.matrix[1][2] = 6;
    ck_assert_int_eq(s21_transpose(&temp, &result), 0);
    ck_assert_int_eq(s21_eq_matrix(&result_defolt, &result), SUCCESS);
    s21_remove_matrix(&temp);
    s21_remove_matrix(&result);
    s21_remove_matrix(&result_defolt);
}
END_TEST

START_TEST(test_s21_determinant) {
    matrix_t temp;
    double result;
    int result_defolt;
    s21_create_matrix(3, 2, &temp);
    s21_gen_matrix(&temp);
    result_defolt = s21_determinant(&temp, &result);
    ck_assert_int_eq(result_defolt, 2);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_determinant_2) {
    matrix_t temp;
    double result;
    int result_defolt;
    s21_create_matrix(3, 3, &temp);
    s21_gen_matrix(&temp);
    result_defolt = s21_determinant(&temp, &result);
    ck_assert_int_eq(result_defolt, 0);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_determinant_3) {
    matrix_t temp;
    double result;
    int result_defolt;
    s21_create_matrix(3, 3, &temp);
    temp.matrix[0][0] = 0;
    temp.matrix[0][1] = 9;
    temp.matrix[0][2] = 5;

    temp.matrix[1][0] = 4;
    temp.matrix[1][1] = 3;
    temp.matrix[1][2] = -5;

    temp.matrix[2][0] = -1;
    temp.matrix[2][1] = 6;
    temp.matrix[2][2] = -4;
    result_defolt = s21_determinant(&temp, &result);
    ck_assert_int_eq(result_defolt, 0);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_determinant_4) {
    matrix_t temp;
    double result;
    int result_defolt;
    s21_create_matrix(4, 4, &temp);
    temp.matrix[0][0] = 1;
    temp.matrix[0][1] = 0;
    temp.matrix[0][2] = 2;
    temp.matrix[0][3] = -6;

    temp.matrix[1][0] = 9;
    temp.matrix[1][1] = 5;
    temp.matrix[1][2] = 4;
    temp.matrix[1][3] = 5;

    temp.matrix[2][0] = 7;
    temp.matrix[2][1] = 8;
    temp.matrix[2][2] = -4;
    temp.matrix[2][3] = 1;

    temp.matrix[3][0] = 2;
    temp.matrix[3][1] = 0;
    temp.matrix[3][2] = -1;
    temp.matrix[3][3] = 7;
    result_defolt = s21_determinant(&temp, &result);
    ck_assert_int_eq(result_defolt, 0);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_determinant_5) {
    matrix_t temp;
    double result;
    int result_defolt;
    s21_create_matrix(3, 3, &temp);
    temp.matrix[0][0] = 1;
    temp.matrix[1][1] = 1;
    temp.matrix[2][2] = 1;

    result_defolt = s21_determinant(&temp, &result);
    ck_assert_int_eq(result_defolt, 0);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_calc_complements) {
    matrix_t test, result;
    s21_create_matrix(3, 2, &test);
    int result_status;
    s21_gen_matrix(&test);
    result_status = s21_calc_complements(&test, &result);
    ck_assert_int_eq(result_status, 2);
    s21_remove_matrix(&test);
}
END_TEST

START_TEST(test_s21_calc_complements_2) {
    matrix_t test, result;
    s21_create_matrix(3, 3, &test);
    s21_gen_matrix(&test);
    int result_status = s21_calc_complements(&test, &result);
    ck_assert_int_eq(result_status, 0);
    s21_remove_matrix(&test);
    s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_calc_complements_3) {
    matrix_t test, test_2, result;
    s21_create_matrix(3, 3, &test);
    test.matrix[0][0] = 1;
    test.matrix[0][1] = 2;
    test.matrix[0][2] = 3;
    test.matrix[1][0] = 0;
    test.matrix[1][1] = 4;
    test.matrix[1][2] = 2;
    test.matrix[2][0] = 5;
    test.matrix[2][1] = 2;
    test.matrix[2][2] = 1;
    s21_calc_complements(&test, &result);
    s21_create_matrix(3, 3, &test_2);
    test_2.matrix[0][0] = 0;
    test_2.matrix[0][1] = 10;
    test_2.matrix[0][2] = -20;
    test_2.matrix[1][0] = 4;
    test_2.matrix[1][1] = -14;
    test_2.matrix[1][2] = 8;
    test_2.matrix[2][0] = -8;
    test_2.matrix[2][1] = -2;
    test_2.matrix[2][2] = 4;
    for (int i = 0; i < result.rows; i++) {
        for (int j = 0; j < result.columns; j++) {
            ck_assert_uint_eq(result.matrix[i][j], test_2.matrix[i][j]);
        }
    }
    s21_remove_matrix(&test);
    s21_remove_matrix(&test_2);
    s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_inverse_matrix) {
    matrix_t temp, result_defolt, result_defolt_2;
    s21_create_matrix(3, 3, &temp);
    s21_create_matrix(3, 3, &result_defolt);
    temp.matrix[0][0] = 2;
    temp.matrix[0][1] = 5;
    temp.matrix[0][2] = 7;
    temp.matrix[1][0] = 6;
    temp.matrix[1][1] = 3;
    temp.matrix[1][2] = 4;
    temp.matrix[2][0] = 5;
    temp.matrix[2][1] = -2;
    temp.matrix[2][2] = -3;
    result_defolt.matrix[0][0] = 1;
    result_defolt.matrix[0][1] = -1;
    result_defolt.matrix[0][2] = 1;
    result_defolt.matrix[1][0] = -38;
    result_defolt.matrix[1][1] = 41;
    result_defolt.matrix[1][2] = -34;
    result_defolt.matrix[2][0] = 27;
    result_defolt.matrix[2][1] = -29;
    result_defolt.matrix[2][2] = 24;

    s21_inverse_matrix(&temp, &result_defolt_2);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_uint_eq(result_defolt.matrix[i][j],
            result_defolt_2.matrix[i][j]);
        }
    }
    s21_remove_matrix(&temp);
    s21_remove_matrix(&result_defolt);
    s21_remove_matrix(&result_defolt_2);
}
END_TEST

START_TEST(test_s21_inverse_matrix_2) {
    matrix_t temp, result;
    s21_create_matrix(3, 3, &temp);
    s21_gen_matrix(&temp);
    int result_status = s21_inverse_matrix(&temp, &result);
    ck_assert_int_eq(result_status, 2);
    s21_remove_matrix(&temp);
}
END_TEST

START_TEST(test_s21_inverse_matrix_3) {
    matrix_t temp, result;
    s21_create_matrix(0, 0, &temp);
    s21_gen_matrix(&temp);
    s21_inverse_matrix(&temp, &result);
    int result_status = s21_inverse_matrix(&temp, &result);
    ck_assert_int_eq(result_status, 1);
}
END_TEST

START_TEST(test_s21_inverse_matrix_4) {
    matrix_t temp, result;
    s21_create_matrix(3, 3, &temp);
    temp.matrix[0][0] = 1;
    temp.matrix[1][1] = 1;
    temp.matrix[2][2] = 1;

    s21_inverse_matrix(&temp, &result);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_uint_eq(temp.matrix[i][j],
            result.matrix[i][j]);
        }
    }
    s21_remove_matrix(&temp);
    s21_remove_matrix(&result);
}
END_TEST

int main() {
    Suite *s1 = suite_create("s21_matrix: ");
    TCase *tc1_1 = tcase_create("s21_matrix: ");
    SRunner *sr = srunner_create(s1);
    int result;

    suite_add_tcase(s1, tc1_1);
    tcase_add_test(tc1_1, test_s21_eq_matrix);
    tcase_add_test(tc1_1, test_s21_eq_matrix_2);
    tcase_add_test(tc1_1, test_s21_eq_matrix_3);
    tcase_add_test(tc1_1, test_s21_eq_matrix_4);
    tcase_add_test(tc1_1, test_s21_eq_matrix_5);

    tcase_add_test(tc1_1, test_s21_sum_matrix);
    tcase_add_test(tc1_1, test_s21_sum_matrix_2);

    tcase_add_test(tc1_1, test_s21_sub_matrix);
    tcase_add_test(tc1_1, test_s21_sub_matrix_2);

    tcase_add_test(tc1_1, test_s21_mult_number);

    tcase_add_test(tc1_1, test_s21_mult_matrix_2);
    tcase_add_test(tc1_1, test_s21_mult_matrix_3);
    tcase_add_test(tc1_1, test_s21_mult_matrix_4);

    tcase_add_test(tc1_1, test_s21_transpose);

    tcase_add_test(tc1_1, test_s21_determinant);
    tcase_add_test(tc1_1, test_s21_determinant_2);
    tcase_add_test(tc1_1, test_s21_determinant_3);
    tcase_add_test(tc1_1, test_s21_determinant_4);
    tcase_add_test(tc1_1, test_s21_determinant_5);

    tcase_add_test(tc1_1, test_s21_calc_complements);
    tcase_add_test(tc1_1, test_s21_calc_complements_2);
    tcase_add_test(tc1_1, test_s21_calc_complements_3);

    tcase_add_test(tc1_1, test_s21_inverse_matrix);
    tcase_add_test(tc1_1, test_s21_inverse_matrix_2);
    tcase_add_test(tc1_1, test_s21_inverse_matrix_3);
    tcase_add_test(tc1_1, test_s21_inverse_matrix_4);

    srunner_run_all(sr, CK_ENV);
    result = srunner_ntests_failed(sr);
    srunner_free(sr);
    return result == 0 ? 0 : 1;
}
