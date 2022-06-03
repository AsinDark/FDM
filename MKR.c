#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const int time_end = 15;
const double height = 6;
const double width = 8;
const double right_border = 3;
const double top_border = 5;
const double h = 0.25;
const double temp_init = 0;
const double alpha = 1;
int x_nodes = 0;
int y_nodes = 0;
int nodes_amount = 0;
double *gauss(double **matr, const double *y, int n) {
    double d, s;
    double *x = (double *) calloc(sizeof(double), n);
    double *b = (double *) calloc(sizeof(double), n);
    double **a = (double **) calloc(sizeof(double *), n);
    for (int i = 0; i < nodes_amount; ++i) {
        a[i] = (double *) calloc(sizeof(double), nodes_amount);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = matr[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        b[i] = y[i];
    }
    for (int k = 1; k < n; k++) { // прямой ход
        for (int j = k + 1; j < n; j++) {
            d = a[j][k] / a[k][k];
            for (int i = k; i < n; i++) {
                a[j][i] = a[j][i] - d * a[k][i];
            }
            b[j] = b[j] - d * b[k];
        }
    }
    for (int k = n - 1; k >= 1; k--) { // обратный ход
        d = 0;
        for (int j = k + 1; j < n; j++) {
            s = a[k][j] * x[j];
            d = d + s;
        }
        x[k] = (b[k] - d) / a[k][k];
    }
    free(b);
    for (int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    return x;
}
int print_matrix(double **A) {
    FILE *output = fopen("output_matr.dat", "w");
    if (!output) {
        fprintf(stderr, "output open error\n");
        return -1;
    }
    for (int i = 0; i < nodes_amount; i++) {
        for (int j = 0; j < nodes_amount; j++)
            fprintf(output, "%7.3lf ", A[i][j]);
        fprintf(output, "\n");
    }
    fclose(output);
    return 0;
}
double *solution(double *y, double **A) {
    double *x = (double *) calloc(sizeof(double), nodes_amount);
    for (int k = 0; k < time_end; ++k) {
        for (int i = 0; i < y_nodes; ++i) {
            for (int j = 0; j < x_nodes; ++j) {
// Учет генерации тепла (Q в векторе правых частей на каждой итерации)13
                if (i == y_nodes / 2 && j == x_nodes / 2) {
                    y[i * x_nodes + j] = x[i * x_nodes + j];
                }
                if (j > 0 && i > 0 && i < y_nodes - 1 && i > (j - top_border / h) && j < (x_nodes - 1) && i != y_nodes / 2 && j != x_nodes / 2) {
                    y[i * x_nodes + j] = x[i * x_nodes + j];
                }
            }
        }
        free(x);
        x = gauss(A, y, nodes_amount);
    }
    return x;
}
void initialise(double **vector, double ***matrix) {
    for (int i = 0; i < y_nodes; i++) {
        for (int j = 0; j < x_nodes; j++) {
            if (i == 0 && j <= top_border / h) { // Верхняя граница - ГУ 1-го рода
                (*vector)[i * x_nodes + j] = 50; // Само значение в векторе температур
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = 1; // Указание этого уравнения в СЛАУ
            } else if (j == 0) { // Левая граница - ГУ 3-го рода
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = -1 / h - 1;
                (*matrix)[i * x_nodes + j][i * x_nodes + (j + 1)] = 1 / h;
                (*vector)[i * x_nodes + j] = temp_init; // начальное условие нестационарной задачи
            } else if (i == y_nodes - 1) { // Нижняя граница - ГУ 1-го рода
                (*vector)[i * x_nodes + j] = 100;
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = 1;
            } else if (j >= top_border / h && i == (j - top_border / h)) { // Правая скошенная граница - ГУ 3-го рода
                double cos_45 = sqrt(2) / 2;
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = - 1 / (sqrt(2) * h)
                                                              - 1;
                (*matrix)[i * x_nodes + j][(i + 1) * x_nodes + j - 1] = 1 /
                                                                        (sqrt(2) * h);
                (*vector)[i * x_nodes + j] = temp_init;
            } else if (i * h <= (height - right_border) && j >= top_border / h &&
                       i <= (j - top_border / h)) { // Пустая область
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = 1;
                (*vector)[i * x_nodes + j] = 0;
            } else if (j == x_nodes - 1) { // Правая граница - ГУ 1-го рода
                (*vector)[i * x_nodes + j] = 50;
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = 1;
            } else if (j < x_nodes - 1) { // Остальные (внутренние) узлы
                (*matrix)[i * x_nodes + j][i * x_nodes + j] = 1 + 2 / (h * h) + 2 / (h * h);
                (*matrix)[i * x_nodes + j][i * x_nodes + j + 1] = -1 / (h * h);
                (*matrix)[i * x_nodes + j][i * x_nodes + j - 1] = -1 / (h * h);
                (*matrix)[i * x_nodes + j][(i + 1) * x_nodes + j] = -1 / (h * h);
                (*matrix)[i * x_nodes + j][(i - 1) * x_nodes + j] = -1 / (h * h);
            }
        }
    }
}
int output(double *x, int x_n, int y_n) {
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) {
        fprintf(stderr, "gnuplot open error\n");
        return -1;
    }
// Вывод
    FILE *output = fopen("output.dat", "w");
    if (!output) {
        fprintf(stderr, "output open error\n");
        return -1;
    }
    for (int i = 0; i < y_n; i++) {
        for (int j = 0; j < x_n; j++)
            fprintf(output, "%7.3lf ", x[i * x_n + j]);
        fprintf(output, "\n");
    }
    fprintf(gp, "set pm3d map\n"
                "set pm3d interpolate 2,2\n"
                "set cbrange [0:100]\n"
                "set palette defined (0 \"#000000\", 1 \"#001aff\", 3 \"#00f2ff\", 5 \"#00ffaa\", 7 \"#f2ff00\", 10 \"#ff2200\")\n;"
                "set autoscale fix\n"
                "set cbtics 10\n"
                 "splot '-' matrix\n");
    for (int i = y_n - 1; i >= 0; i--) {
        for (int j = 0; j < x_n; j++) {
            fprintf(gp, "%lf ", x[i * x_n + j]);
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
    fclose(gp);
    fclose(output);
    return 0;
}
int main() {
    y_nodes = (int) (height / h + 1);
    x_nodes = (int) (width / h + 1);
    nodes_amount = (int) x_nodes * y_nodes;
    double *temperatures = (double *) calloc(sizeof(double), nodes_amount);
    double **matrix = (double **) malloc(sizeof(double *) * nodes_amount);
    for (int i = 0; i < nodes_amount; ++i) {
        matrix[i] = (double *) calloc(sizeof(double), nodes_amount);
    }
    initialise(&temperatures, &matrix);
    double *x = solution(temperatures, matrix);
    if (output(x, x_nodes, y_nodes) < 0) {
        fprintf(stderr, "output err");
        return -1;
    }
    free(temperatures);
    for (int i = 0; i < nodes_amount; ++i) {
        free(matrix[i]);
    }
    free(matrix);
    free(x);
    return 0;
}
