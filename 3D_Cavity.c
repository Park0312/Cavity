#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 100           // 격자 크기
#define Re 1000.0        // 레이놀즈 수
#define dt 0.001         // 시간 간격
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
    int max_iter = 10; // 최대 반복 횟수
    double tol = 1e-8;  // 잔차 허용 오차
    for (int iter = 0; iter < max_iter; iter++) {
        double residual = 0.0;

        // 압력 갱신
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < N; j++) {
                for (int k = 1; k < N; k++) {
                    double rhs = ((u_star[i][j][k] - u_star[i-1][j][k]) / dx +
                                  (v_star[i][j][k] - v_star[i][j-1][k]) / dy +
                                  (w_star[i][j][k] - w_star[i][j][k-1]) / dz) / dt;

                    double p_new = (1.0 / 6.0) * (p[i+1][j][k] + p[i-1][j][k] +
                                                  p[i][j+1][k] + p[i][j-1][k] +
                                                  p[i][j][k+1] + p[i][j][k-1] -
                                                  rhs * dx * dy * dz);

                    residual += fabs(p_new - p[i][j][k]);
                    p[i][j][k] = p_new;
                }
            }
        }

        // 경계 조건 적용
        for (int i = 0; i <= N+1; i++) {
            for (int j = 0; j <= N+1; j++) {
                p[i][j][0] = p[i][j][1];           // 전면 경계
                p[i][j][N+1] = p[i][j][N];         // 후면 경계
            }
        }
        for (int j = 0; j <= N+1; j++) {
            for (int k = 0; k <= N+1; k++) {
                p[0][j][k] = p[1][j][k];           // 좌측 경계
                p[N+1][j][k] = p[N][j][k];         // 우측 경계
            }
        }
        for (int i = 0; i <= N+1; i++) {
            for (int k = 0; k <= N+1; k++) {
                p[i][0][k] = p[i][1][k];           // 아래 경계
                p[i][N+1][k] = p[i][N][k];         // 위 경계
            }
        }

        // 잔차 확인
        if (residual < tol) {
            break;
        }
    }
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
    FILE *file = fopen("3D_results.csv", "w");
    fprintf(file, "x,y,z,u,v,w\n"); // CSV 헤더 추가
    int step = 5; // 데이터 저장 간격
    for (int i = 0; i <= N; i += step) {
        for (int j = 0; j <= N; j += step) {
            for (int k = 0; k <= N; k += step) {
                fprintf(file, "%f,%f,%f,%f,%f,%f\n", i * dx, j * dy, k * dz, u[i][j][k], v[i][j][k], w[i][j][k]);
            }
        }
    }
    printf("File Updated!");
    fclose(file);
}

void simulate() {
    double time = 0.0;
    while (time < T) {
        compute_intermediate_velocity();
        solve_pressure_poisson();
        update_velocity();
        
        
        printf("Simulation time : %f\n", time);

        if (time/5 == 0) {
            save_results();
        }
        
        

        // 경계 조건 적용
        for (int i = 0; i <= N; i++) {
            for (int k = 0; k <= N+1; k++) {
                u[i][N+1][k] = 1.0;  // 상단 벽에서 u = 1
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
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= N; j++) {
                p[i][0][j] = p[i][1][j];  // 전면 경계
                p[i][N+1][j] = p[i][N][j]; // 후면 경계
            }
        }

        time += dt;
        
    }
}



int main() {
    initialize();
    simulate();
    save_results();
    return 0;
}
