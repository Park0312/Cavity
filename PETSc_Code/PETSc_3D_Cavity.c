#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <petscsys.h>
#include <petscksp.h>

#define N 100           // 격자 크기
#define Re 1000.0        // 레이놀즈 수
#define dt 0.01         // 시간 간격
#define T 100.0          // 총 시뮬레이션 시간

double u[N+1][N+2][N+2], v[N+2][N+1][N+2], w[N+2][N+2][N+1], p[N+2][N+2][N+2];
double u_star[N+1][N+2][N+2], v_star[N+2][N+1][N+2], w_star[N+2][N+2][N+1];
double dx = 1.0 / N;
double dy = 1.0 / N;
double dz = 1.0 / N;

////////////////// [First Function]  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

void initialize() {
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N+1; j++) {
            for (int k = 0; k <= N+1; k++) {
                u[i][j][k] = 0.0;
                u_star[i][j][k] = 0.0;
            }
            
        }
    }
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N+1; k++) {
                v[i][j][k] = 0.0;
                v_star[i][j][k] = 0.0;
            }
            
        }
    }
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N+1; j++) {
            for (int k = 0; k <= N; k++) {
                w[i][j][k] = 0.0;
                w_star[i][j][k] = 0.0;
            }
            
        }
    }
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N+1; j++) {
            for (int k = 0; k <= N+1; k++) {
                p[i][j][k] = 0.0;
            }
        }
    }
    // 상단 벽 경계 조건 설정
    for (int i = 0; i <= N; i++) {
        for (int k = 0; k <= N+1; k++) {
            u[i][N+1][k] = 1.0;  // 상단 벽에서 u = 1
        }
        
    }
}

////////////////// [Second Function]  ///////////////////////////////////////////////////////////////////////////////////////////////////////////


double compute_convection_u(int i, int j, int k) {
    double du2dx = ((u[i][j][k] + u[i+1][j][k]) * (u[i][j][k] + u[i+1][j][k]) - (u[i-1][j][k] + u[i][j][k]) * (u[i-1][j][k] + u[i][j][k])) / (4.0 * dx);
    double duvdy = ((v[i][j][k] + v[i+1][j][k]) * (u[i][j][k] + u[i][j+1][k]) - (v[i][j-1][k] + v[i+1][j-1][k]) * (u[i][j-1][k] + u[i][j][k])) / (4.0 * dy);
    double duwdz = ((w[i][j][k] + w[i+1][j][k]) * (u[i][j][k] + u[i][j][k+1]) - (w[i][j][k-1] + w[i+1][j][k-1]) * (u[i][j][k-1] + u[i][j][k])) / (4.0 * dz);
    return du2dx + duvdy + duwdz;
}

double compute_convection_v(int i, int j, int k) {
    double duvdx = ((u[i][j][k] + u[i][j+1][k]) * (v[i][j][k] + v[i+1][j][k]) - (u[i-1][j][k] + u[i-1][j+1][k]) * (v[i-1][j][k] + v[i][j][k])) / (4.0 * dx);
    double dv2dy = ((v[i][j][k] + v[i][j+1][k]) * (v[i][j][k] + v[i][j+1][k]) - (v[i][j-1][k] + v[i][j][k]) * (v[i][j-1][k] + v[i][j][k])) / (4.0 * dy);
    double dvwdz = ((w[i][j][k] + w[i][j+1][k]) * (v[i][j][k] + v[i][j][k+1]) - (w[i][j][k-1] + w[i][j+1][k-1]) * (v[i][j][k-1] + v[i][j][k])) / (4.0 * dz);
    return duvdx + dv2dy + dvwdz;
}

double compute_convection_w(int i, int j, int k) {
    double duwdx = ((u[i][j][k] + u[i][j][k+1]) * (w[i][j][k] + w[i+1][j][k]) - (u[i-1][j][k] + u[i-1][j][k+1]) * (w[i-1][j][k] + w[i][j][k])) / (4.0 * dx);
    double duvdy = ((v[i][j][k] + v[i][j][k+1]) * (w[i][j][k] + w[i][j+1][k]) - (v[i][j-1][k] + v[i][j-1][k+1]) * (w[i][j-1][k] + w[i][j][k])) / (4.0 * dy);
    double dw2dz = ((w[i][j][k] + w[i][j][k+1]) * (w[i][j][k] + w[i][j][k+1]) - (w[i][j][k-1] + w[i][j][k]) * (w[i][j][k-1] + w[i][j][k])) / (4.0 * dz);
    return duwdx + duvdy + dw2dz;
}

double compute_diffusion_u(int i, int j, int k) {
    double d2udx2 = (u[i+1][j][k] - 2.0 * u[i][j][k] + u[i-1][j][k]) / (dx * dx);
    double d2udy2 = (u[i][j+1][k] - 2.0 * u[i][j][k] + u[i][j-1][k]) / (dy * dy);
    double d2udz2 = (u[i][j][k+1] - 2.0 * u[i][j][k] + u[i][j][k-1]) / (dz * dz);
    return d2udx2 + d2udy2 + d2udz2;
}

double compute_diffusion_v(int i, int j, int k) {
    double d2vdx2 = (v[i+1][j][k] - 2.0 * v[i][j][k] + v[i-1][j][k]) / (dx * dx);
    double d2vdy2 = (v[i][j+1][k] - 2.0 * v[i][j][k] + v[i][j-1][k]) / (dy * dy);
    double d2vdz2 = (v[i][j][k+1] - 2.0 * v[i][j][k] + v[i][j][k-1]) / (dz * dz);
    return d2vdx2 + d2vdy2 + d2vdz2;
}

double compute_diffusion_w(int i, int j, int k) {
    double d2wdx2 = (w[i+1][j][k] - 2.0 * w[i][j][k] + w[i-1][j][k]) / (dx * dx);
    double d2wdy2 = (w[i][j+1][k] - 2.0 * w[i][j][k] + w[i][j-1][k]) / (dy * dy);
    double d2wdz2 = (w[i][j][k+1] - 2.0 * w[i][j][k] + w[i][j][k-1]) / (dz * dz);
    return d2wdx2 + d2wdy2 + d2wdz2;
}

void compute_intermediate_velocity() {
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++) {
                double convection = compute_convection_u(i, j, k);
                double diffusion = compute_diffusion_u(i, j, k);
                u_star[i][j][k] = u[i][j][k] + dt * (-convection + (1.0 / Re) * diffusion);
            }
            
        }
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j < N; j++) {
            for (int k = 1; k <= N; k++) {
                double convection = compute_convection_v(i, j, k);
                double diffusion = compute_diffusion_v(i, j, k);
                v_star[i][j][k] = v[i][j][k] + dt * (-convection + (1.0 / Re) * diffusion);
            }
            
        }
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k < N; k++) {
                double convection = compute_convection_w(i, j, k);
                double diffusion = compute_diffusion_w(i, j, k);
                w_star[i][j][k] = w[i][j][k] + dt * (-convection + (1.0 / Re) * diffusion);
            }
            
        }
    }
}

void solve_pressure_poisson() {
    Mat A;          // 행렬 A (라플라시안)
    Vec b, x;       // 벡터 b (rhs)와 x (압력)
    KSP ksp;        // Krylov Subspace Solver
    PetscErrorCode ierr;

    // 행렬 및 벡터 크기 정의
    PetscInt n = (N + 2) * (N + 2) * (N + 2); // 전체 격자 크기

    // 행렬과 벡터 생성
    ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetFromOptions(A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetUp(A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, n); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecSetFromOptions(b); CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = VecDuplicate(b, &x); CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // 행렬 A 설정 (라플라시안)
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                PetscInt row = i * (N + 2) * (N + 2) + j * (N + 2) + k;
                MatSetValue(A, row, row, -6.0, INSERT_VALUES); // 중앙

                // 이웃 설정
                if (k < N) MatSetValue(A, row, row + 1, 1.0, INSERT_VALUES);   // +x
                if (k > 1) MatSetValue(A, row, row - 1, 1.0, INSERT_VALUES);   // -x
                if (j < N) MatSetValue(A, row, row + (N + 2), 1.0, INSERT_VALUES); // +y
                if (j > 1) MatSetValue(A, row, row - (N + 2), 1.0, INSERT_VALUES); // -y
                if (i < N) MatSetValue(A, row, row + (N + 2) * (N + 2), 1.0, INSERT_VALUES); // +z
                if (i > 1) MatSetValue(A, row, row - (N + 2) * (N + 2), 1.0, INSERT_VALUES); // -z
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // 벡터 b 설정 (rhs)
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                PetscInt idx = i * (N + 2) * (N + 2) + j * (N + 2) + k;
                double rhs = ((u_star[i][j][k] - u_star[i-1][j][k]) / dx +
                              (v_star[i][j][k] - v_star[i][j-1][k]) / dy +
                              (w_star[i][j][k] - w_star[i][j][k-1]) / dz) / dt;

                VecSetValue(b, idx, rhs, INSERT_VALUES);
            }
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // KSP Solver 설정
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // Poisson 방정식 해결
    ierr = KSPSolve(ksp, b, x); CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // 압력 벡터를 p 배열로 복사
    PetscScalar *x_array;
    VecGetArray(x, &x_array);
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                PetscInt idx = i * (N + 2) * (N + 2) + j * (N + 2) + k;
                p[i][j][k] = x_array[idx];
            }
        }
    }
    VecRestoreArray(x, &x_array);

    // 메모리 정리
    KSPDestroy(&ksp);
    VecDestroy(&b);
    VecDestroy(&x);
    MatDestroy(&A);
}


void update_velocity() {
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++){
                u[i][j][k] = u_star[i][j][k] - dt * (p[i+1][j][k] - p[i][j][k]) / dx;
            }
        }
    }
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++){
                v[i][j][k] = v_star[i][j][k] - dt * (p[i][j+1][k] - p[i][j][k]) / dy;
            }
        }
    }
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++){
                w[i][j][k] = w_star[i][j][k] - dt * (p[i][j][k+1] - p[i][j][k]) / dz;
            }
        }
    }
}

////////////////// [Third Function]  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
void save_results() {
    FILE *file = fopen("3D_results_1000seconds.csv", "w");
    fprintf(file, "x,y,z,u,v,w\n"); // CSV 헤더 추가
    int step = 5; // 데이터 저장 간격
    for (int i = 0; i <= N; i += step) {
        for (int j = 0; j <= N; j += step) {
            for (int k = 0; k <= N; k += step) {
                fprintf(file, "%f,%f,%f,%f,%f,%f\n", i * dx, j * dy, k * dz, u[i][j][k], v[i][j][k], w[i][j][k]);
            }
        }
    }
    printf("File Updated!\n");
    fclose(file);
}

void simulate() {
    double time = 0.0;
    while (time < T) {
        compute_intermediate_velocity();
        solve_pressure_poisson();
        update_velocity();
        
        
        
        
        

        // 경계 조건 적용
        for (int i = 0; i <= N; i++) {
            for (int k = 0; k <= N+1; k++) {
                u[i][N][k] = 1.0;  // 상단 벽에서 u = 1
                u[i][0][k] = 0.0;   // 하단 벽에서 u = 0
            }
        }

        for (int j = 0; j <= N+1; j++) {
            for (int k = 0; k <= N+1; k++) {
                u[0][j][k] = 0.0;    // 좌측 벽에서 u = 0
                u[N][j][k] = 0.0;    // 우측 벽에서 u = 0
            }
        }

        for (int i = 0; i <= N+1; i++) {
            for (int j = 0; j <= N+1; j++) {
                w[i][j][0] = 0.0;    // 전면 벽에서 w = 0
                w[i][j][N] = 0.0;    // 후면 벽에서 w = 0
            }
        }

        for (int i = 0; i <= N+1; i++) {
            for (int k = 0; k <= N+1; k++) {
                v[i][0][k] = 0.0;    // 아래 벽에서 v = 0
                v[i][N][k] = 0.0;    // 위 벽에서 v = 0
            }
        }
        for (int i = 0; i <= N+1; i++) {
            for (int j = 0; j <= N+1; j++) {
                u[i][j][0] = 0.0;           // 전면 벽에서 속도 0
                u[i][j][N+1] = 0.0;         // 후면 벽에서 속도 0
                v[i][j][0] = 0.0;           // 전면 벽에서 속도 0
                v[i][j][N+1] = 0.0;         // 후면 벽에서 속도 0
                w[i][j][0] = 0.0;           // 전면 벽에서 속도 0
                w[i][j][N] = 0.0;           // 후면 벽에서 속도 0
            }
        }
        // 경계 설정
        for (int i = 0; i <= N+1; i++) {
            for (int j = 0; j <= N+1; j++) {
                p[i][0][j] = p[i][1][j];  // 전면 경계
                p[i][N+1][j] = p[i][N][j]; // 후면 경계
            }
        }

        time += dt;

        printf("Simulation time : %f\n", time);

        if (time/1 == 0) {
            save_results();
        }
        
    }
}



int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, "#d Lid Driven Cavity Simulation with PETSc");
    
    initialize();
    simulate();
    save_results();

    PetscFinalize();
    return 0;
}
