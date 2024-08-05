#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const float PI = 3.1415926535;

float** allocate(int N) {
    int sz = N + 1;
    float** grid = malloc(sz * sizeof(float*));

    for (int i = 0; i < sz; ++i) 
        grid[i] = malloc(sz * sizeof(float));

    for (int i = 0; i < sz; ++i) 
        for (int j = 0; j < sz; ++j) 
            grid[i][j] = 0.0;

    return grid;
}

void deallocate(float** grid, int N) {
    for (int i = 0; i <= N; ++i) free(grid[i]);
    free(grid);
}

void print_grid(float** grid, int N) {
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            printf("u(%d,%d) = %f\n", i, j, grid[i][j]);
        }
    }
}

void print_latex_array(float** grid, int N) {
    printf("\\begin{bmatrix}\n");
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N-1; ++j) {
            printf("%f & ", i, j, grid[i][j]);
        }
        printf("%f \\\\\n", grid[i][N]);
    }
    printf("\\end{bmatrix}\n");
}

void print_dat_array(float** grid, int N) {
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N-1; ++j) {
            printf("%f\t", grid[i][j]);
        }
        printf("%f\n", grid[i][N]);
    }
}

void print_csv_array(float** grid, int N) {
    for (int i = 0; i <= N; ++i) {
        printf("%f", grid[i][0]);
        for (int j = 1; j <= N; ++j) {
            printf(", %f", grid[i][j]);
        }
        printf("\n");
    }
}

void print_python_array(float** grid, int N) {
    printf("grid = np.array([[%f", grid[0][0]);
    for (int j = 1; j <= N; ++j) {
        printf(", %f", grid[0][j]);
    }

    for (int i = 1; i <= N; ++i) {
        printf("], [%f", grid[i][0]);
        for (int j = 1; j <= N; ++j) {
            printf(", %f", grid[i][j]);
        }
    }

    printf("]])\n");
}

void copy_first_into_second(float** first, float** second, int N) {
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            second[i][j] = first[i][j];
        }
    }
}


void apply_boundary(float** grid, int N) {
    for (int i = 0; i <= N; ++i) {
        float x = (float)i / (float)(N);
        grid[i][0] = sin(PI * x) * sin(PI * x);
        grid[i][N] = 0;
    }
}

void gauss_seidel_relaxation(float** grid, int N) {
    for (int n = 0; n < 100; ++n) {
        for (int i = 1; i <= N-1; ++i) {
            for (int j = 1; j <= N-1; ++j) {
                grid[i][j] = 0.25 * (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1]);
            }
        }
    }
}

void jacobi_relaxation(float** grid, int N) {
    float** aux = allocate(N);
    copy_first_into_second(grid, aux, N);

    for (int k = 0; k < 100; ++k) {
        for (int i = 1; i <= N-1; ++i) {
            for (int j = 1; j <= N-1; ++j) {
                grid[i][j] = 0.25 * (aux[i-1][j] + aux[i+1][j] + aux[i][j-1] + aux[i][j+1]);
            }
        }

        copy_first_into_second(grid, aux, N);
    }

    deallocate(aux, N);
}

void analytic_solution(float** grid, int N) {
    for (int n = 1; n < 1000; ++n) {
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                float x = (float)i / (float)(N);
                float y = (float)j / (float)(N);
                float fn = (float)n;

                float arg1 = (2*fn - 1) * PI * (1-y);
                float arg2 = (2*fn - 1) * PI;

                float ratio = (1 - exp(-2*arg1)) / (1 - exp(-2*arg2));
                float factor = exp(arg1 - arg2) * ratio;
                // float factor = sinh(arg1) / sinh(arg2);

                float den = (4*fn*fn - 1) * (2*fn-3);
                float value = -8/PI / den * factor * sin((2*fn-1) * PI * x);
                grid[i][j] += value;
            }
        }
    }
}

int main() {
    int N = 3;
    float** grid = allocate(N);

    apply_boundary(grid, N);

    gauss_seidel_relaxation(grid, N);
    // jacobi_relaxation(grid, N);

    // print_grid(grid, N);
    // print_latex_array(grid, N);
    // print_python_array(grid, N);
    print_csv_array(grid, N);
    // print_dat_array(grid, N);

    deallocate(grid, N);
    return 0;
}
