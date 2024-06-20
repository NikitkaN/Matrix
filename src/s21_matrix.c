#include "s21_matrix.h"

/**
 *  Создание матриц
 *  0 - ОК
 *  1 - Ошибка, некорректная матрица
 */
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error_code = OK;

  result->rows = rows;
  result->columns = columns;

  // Если создание матрицы физически возможно, то создаем
  if (rows > 0 && columns > 0) {
    result->matrix = calloc(result->rows, sizeof(double *));
    if (result->matrix != NULL) {
      for (int i = 0; i < result->rows; i++) {
        result->matrix[i] = calloc(result->columns, sizeof(double));
        if (result->matrix[i] == NULL) {
          error_code = INCORRECT_MATRIX;
          break;
        }
      }
    } else {
      error_code = INCORRECT_MATRIX;
    }
  } else {
    result->matrix = NULL;
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

// Очистка матриц
void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->rows = 0;
    A->columns = 0;
    A = NULL;
  }
}

/**
 * Сравнение матриц
 * 1 - SUCCESS
 * 0 - FAILUER
 */
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error_code = SUCCESS;

  // Проверка на равенство строк и столбцов матриц А и В
  if (A != NULL && A->matrix != NULL && B != NULL && B->matrix != NULL &&
      A->rows == B->rows && A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        // Нахождение разности значений матриц А и В
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
          error_code = FAILURE;
          break;
        }
      }
      if (!error_code) {
        break;
      }
    }
  } else {
    error_code = FAILURE;
  }
  return error_code;
}

/**
 * Сложение матриц
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 * 2 - Ошибка вычисления
 */
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = OK;

  // Проверка на корректность матриц
  if (s21_check_matrix(A) == OK && s21_check_matrix(B) == OK) {
    // Проверка на равенство строк и столбцов матриц А и В
    if (A->rows != B->rows || A->columns != B->columns) {
      error_code = CALC_ERROR;
    } else {
      error_code = s21_create_matrix(A->rows, A->columns, result);

      if (error_code == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
          }
        }
      }
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Сложение матриц
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 * 2 - Ошибка вычисления
 */
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = OK;

  // Проверка на корректность матриц
  if (A != NULL && B != NULL && s21_check_matrix(A) == OK &&
      s21_check_matrix(B) == OK) {
    // Проверка на равенство строк и столбцов матриц А и В
    if (A->rows != B->rows || A->columns != B->columns) {
      error_code = CALC_ERROR;
    } else {
      error_code = s21_create_matrix(A->rows, A->columns, result);

      if (error_code == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
          }
        }
      }
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Умножение матрицы на число
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 */
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error_code = OK;

  // Проверка на корректность матрицы
  if (A != NULL && s21_check_matrix(A) == OK) {
    error_code = s21_create_matrix(A->rows, A->columns, result);
    if (error_code == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }

  return error_code;
}

/**
 * Умножение матриц
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 * 2 - Ошибка вычисления
 */
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = OK;

  // Проверка на корректность матриц
  if (!s21_check_matrix(A) && !s21_check_matrix(B)) {
    // Проверка на возможность вычисления
    if (A->columns == B->rows) {
      error_code = s21_create_matrix(A->rows, B->columns, result);
      if (error_code == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < B->columns; j++) {
            for (int k = 0; k < B->rows; k++) {
              result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
            }
          }
        }
      }
    } else {
      error_code = CALC_ERROR;
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Транспонирование матрицы
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 */
int s21_transpose(matrix_t *A, matrix_t *result) {
  int error_code = OK;

  // Проверка на корректность матрицы
  if (!s21_check_matrix(A)) {
    error_code = s21_create_matrix(A->columns, A->rows, result);

    if (error_code == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[j][i] = A->matrix[i][j];
        }
      }
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Минор матрицы и матрица алгебраических дополнений
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 * 2 - Ошибка вычисления
 */
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error_code = OK;
  // Проверка на корректность матриц
  if (!s21_check_matrix(A)) {
    // Проверка на возможность вычисления
    if (A->rows == A->columns) {
      error_code = s21_create_matrix(A->rows, A->columns, result);

      if (error_code == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            double determinant = 0.0;
            matrix_t minor;
            error_code = s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
            // Нахождение алгебраического дополнения
            if (error_code == OK) {
              s21_fill_minor(A, &minor, i, j);
              s21_determinant(&minor, &determinant);

              result->matrix[i][j] = pow(-1, i + j) * determinant;

              s21_remove_matrix(&minor);
            } else {
              error_code = INCORRECT_MATRIX;
            }
          }
        }
      } else {
        error_code = INCORRECT_MATRIX;
      }
    } else {
      error_code = CALC_ERROR;
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Определитель матрицы
 * 0 - ОК
 * 1 - Ошибка, некорректная матрица
 * 2 - Ошибка вычисления
 */
int s21_determinant(matrix_t *A, double *result) {
  int error_code = OK;

  // Проверка на корректность матрицы
  if (!s21_check_matrix(A)) {
    // Проверка на возможность вычисления
    if (A->rows == A->columns) {
      int count_of_iter = A->rows;
      double buffer = 0.0;

      // Если матрица состоит из 1 элемента, то она является определителем
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2 || A->rows == 3) {
        // Матрице размером 2х2 необходима одна итерация main_diagonal -
        // side_diagonal
        if (A->rows == 2) {
          count_of_iter--;
        }

        // Нахождение определителя методом "треугольника"
        for (int i = 0; i < count_of_iter; i++) {
          double main_diagonal = 1;
          double side_diagonal = 1;

          for (int j = 0; j < A->rows; j++) {
            main_diagonal *= A->matrix[j % A->rows][(j + i) % A->rows];
          }

          for (int j = 0; j < A->rows; j++) {
            side_diagonal *=
                A->matrix[j % A->rows][(A->rows - 1 - j + i) % A->rows];
          }

          buffer += main_diagonal - side_diagonal;
        }
        *result = buffer;
      } else {
        for (int i = 0; i < A->rows; i++) {
          double determinant = 0;

          // Рекурсивный вызов функции до тех пор, пока матрица не будет
          // размером 3х3
          matrix_t matrix_buffer;
          s21_create_matrix(A->rows - 1, A->columns - 1, &matrix_buffer);
          s21_fill_minor(A, &matrix_buffer, i, 0);
          s21_determinant(&matrix_buffer, &determinant);

          determinant *= A->matrix[i][0];

          if (i % 2 == 1) {
            *result -= determinant;
          } else {
            *result += determinant;
          }

          s21_remove_matrix(&matrix_buffer);
        }
      }
    } else {
      error_code = CALC_ERROR;
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Обратная матрица
 * 0 - OK
 * 1 - INCORRECT_MATRIX
 * 2 - Ошибка вычисления
 */
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error_code = OK;

  // Проверка на корректность матрицы
  if (!s21_check_matrix(A)) {
    // Проверка на возможность вычисления
    if (A->rows == A->columns) {
      double determinant = 0;
      s21_determinant(A, &determinant);

      if (fabs(determinant) > 1e-7) {
        matrix_t buffer;
        matrix_t buffer2;

        s21_calc_complements(A, &buffer);
        s21_transpose(&buffer, &buffer2);

        s21_mult_number(&buffer2, 1.0 / determinant, result);

        s21_remove_matrix(&buffer);
        s21_remove_matrix(&buffer2);
      } else {
        error_code = CALC_ERROR;
      }
    } else {
      error_code = CALC_ERROR;
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

/**
 * Проверка матрицы на корректность
 * 0 - OK
 * 1 - INCORRECT_MATRIX
 */
int s21_check_matrix(matrix_t *A) {
  int error_code = OK;

  if (A == NULL || A->matrix == NULL || A->rows < 1 || A->columns < 1) {
    error_code = INCORRECT_MATRIX;
  }

  return error_code;
}

// Нахождение минора матрицы
void s21_fill_minor(matrix_t *A, matrix_t *minor, int i, int j) {
  int minor_row_index = 0;
  int minor_column_index = 0;
  for (int m = 0; m < A->rows; m++) {
    for (int n = 0; n < A->columns; n++) {
      if (m != i && n != j) {
        minor->matrix[minor_row_index][minor_column_index] = A->matrix[m][n];

        if (minor_column_index + 1 != minor->columns) {
          minor_column_index++;
        } else {
          minor_row_index++;
          minor_column_index = 0;
        }
      }
    }
  }
}
