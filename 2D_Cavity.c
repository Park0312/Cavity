#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 120            // 격자 크기
#define Re 1000.0        // 레이놀즈 수
#define dt 0.001         // 시간 간격
#define T 100.0           // 총 시뮬레이션 시간

double u[N+1][N+2], v[N+2][N+1], p[N+2][N+2];
double u_star[N+1][N+2], v_star[N+2][N+1];
double dx = 1.0 / N;
double dy = 1.0 / N;

void initialize();
double compute_convection_u(int i, int j);
double compute_convection_v(int i, int j);
double compute_diffusion_u(int i, int j);
double compute_diffusion_v(int i, int j);
void compute_intermediate_velocity();
void solve_pressure_poisson();
void update_velocity();
void simulate();
void save_results();

int main() {
    initialize();
    simulate();
    save_results();
    return 0;
}


////////////////// [First Function]  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

void initialize() {
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N+1; j++) {
            u[i][j] = 0.0;
            u_star[i][j] = 0.0;
        }
    }
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N; j++) {
            v[i][j] = 0.0;
            v_star[i][j] = 0.0;
        }
    }
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N+1; j++) {
            p[i][j] = 0.0;
        }
    }
    // 상단 벽 경계 조건 설정
    for (int i = 0; i <= N; i++) {
        u[i][N+1] = 1.0;  // 상단 벽에서 u = 1
    }
}

////////////////// [Second Function]  ///////////////////////////////////////////////////////////////////////////////////////////////////////////


double compute_convection_u(int i, int j) {
    double du2dx = ((u[i][j] + u[i+1][j]) * (u[i][j] + u[i+1][j]) - (u[i-1][j] + u[i][j]) * (u[i-1][j] + u[i][j])) / (4.0 * dx);
    double duvdy = ((v[i][j] + v[i+1][j]) * (u[i][j] + u[i][j+1]) - (v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] + u[i][j])) / (4.0 * dy);
    return du2dx + duvdy;
}

double compute_convection_v(int i, int j) {
    double duvdx = ((u[i][j] + u[i][j+1]) * (v[i][j] + v[i+1][j]) - (u[i-1][j] + u[i-1][j+1]) * (v[i-1][j] + v[i][j])) / (4.0 * dx);
    double dv2dy = ((v[i][j] + v[i][j+1]) * (v[i][j] + v[i][j+1]) - (v[i][j-1] + v[i][j]) * (v[i][j-1] + v[i][j])) / (4.0 * dy);
    return duvdx + dv2dy;
}

double compute_diffusion_u(int i, int j) {
    double d2udx2 = (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / (dx * dx);
    double d2udy2 = (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / (dy * dy);
    return d2udx2 + d2udy2;
}

double compute_diffusion_v(int i, int j) {
    double d2vdx2 = (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j]) / (dx * dx);
    double d2vdy2 = (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]) / (dy * dy);
    return d2vdx2 + d2vdy2;
}

void compute_intermediate_velocity() {
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            double convection = compute_convection_u(i, j);
            double diffusion = compute_diffusion_u(i, j);
            u_star[i][j] = u[i][j] + dt * (-convection + (1.0 / Re) * diffusion);
        }
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j < N; j++) {
            double convection = compute_convection_v(i, j);
            double diffusion = compute_diffusion_v(i, j);
            v_star[i][j] = v[i][j] + dt * (-convection + (1.0 / Re) * diffusion);
        }
    }
}

void solve_pressure_poisson() {
    int max_iter = 1000;
    double tol = 5e-3;
    for (int iter = 0; iter < max_iter; iter++) {
        double residual = 0.0;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                double rhs = ((u_star[i][j] - u_star[i-1][j]) / dx + (v_star[i][j] - v_star[i][j-1]) / dy) / dt;
                double p_new = 0.25 * (p[i+1][j] + p[i-1][j] + p[i][j+1] + p[i][j-1] - rhs * dx * dy); // because dx == dy
                residual += fabs(p_new - p[i][j]);
                p[i][j] = p_new;
            }
        }
        if (residual < tol) break;
    }
}

void update_velocity() {
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            u[i][j] = u_star[i][j] - dt * (p[i+1][j] - p[i][j]) / dx;
        }
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j < N; j++) {
            v[i][j] = v_star[i][j] - dt * (p[i][j+1] - p[i][j]) / dy;
        }
    }
}

void simulate() {
    double time = 0.0;
    while (time < T) {
        compute_intermediate_velocity();
        solve_pressure_poisson();
        update_velocity();

        // 경계 조건 적용
        for (int i = 0; i <= N; i++) {
            u[i][N] = 1.0;
            u[i][1] = 0.0;
            v[i][N] = 0.0;
            v[i][1] = 0.0;

            p[i][0] = p[i][1];
            p[i][N+1] = p[i][N];
        }
        for (int j = 0; j <= N; j++) {
            v[1][j] = 0.0;
            v[N][j] = 0.0;
            u[1][j] = 0.0;
            u[N][j] = 0.0;

            p[0][j] = p[1][j];
            p[N+1][j] = p[N][j];
        }

        time += dt;
    }
}
////////////////// [Third Function]  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
void save_results() {
    FILE *file = fopen("2D_results.csv", "w");
    fprintf(file, "x,y,u,v\n"); // CSV 헤더 추가
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            fprintf(file, "%f,%f,%f,%f\n", i * dx, j * dy, u[i][j], v[i][j]);
        }
    }
    fclose(file);
}


