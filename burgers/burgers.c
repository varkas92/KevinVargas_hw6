#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *reserva_ones(int n_points);

int main(){
  
  int nx = 41;
  int ny = 41;
  int nt = 500;
  float dx = 2.0 / (nx - 1);
  float dy = 2.0 / (ny - 1);
  float sigma = 0.0009;
  float nu = 0.01;
  float dt = sigma * dx * dy / nu;

  float *u;
  float *v;

  float *un;
  float *vn;

  int i;
  int j;
  int k;
  int l;
  int t;

  u = reserva_ones(ny * nx);
  v = reserva_ones(ny * nx);

  un = reserva_ones(ny * nx);
  vn = reserva_ones(ny * nx);

  for(i = 0; i < 10; i++){

    for(j = 0; j < 10; j++){

      u[(j + 10) * ny + 10 + i] = 2.0;

      v[(j + 10) * nx + 10 + i] = 2.0;

      un[(j + 10) * ny + 10 + i] = 2.0;

      vn[(j + 10) * nx + 10 + i] = 2.0;

    }

  }

  for(i = 0; i < ny; i++){

    for(j = 0; j < nx; j++){

        printf("%f ", u[i * ny + j]);
    
    }

    for(j = 0; j < nx; j++){

        printf("%f ", v[i * ny + j]);

    }

    printf("\n");

}

  for(t = 1; t < (nt + 1); t++){

    for(k = 1; k < (ny - 1); k++){

      for(l = 1; l < (nx - 1); l++){

	un[k * ny + l] =  u[k * ny + l];

	vn[k * ny + l] =  v[k * ny + l];

      }

      for(l = 1; l < (nx - 1); l++){

	u[k * ny + l] = un[k * ny + l] - dt / dx * un[k * ny + l] * (un[k * ny + l] - un[k * ny + (l - 1)]) - dt / dy * vn[k * ny + l] * (un[k * ny + l] - un[(k - 1) * ny + l]) + nu * dt / pow(dx, 2) * (un[k * ny + (l + 1)] - 2 * un[k * ny + l] + un[k * ny + (l - 1)]) + nu * dt / pow(dy, 2) * (un[(k + 1) * ny + l] - 2 * un[k * ny + l] + un[(k - 1) * ny + l]);

	v[k * ny + l] = vn[k * ny + l] - dt / dx * un[k * ny + l] * (vn[k * ny + l] - vn[k * ny + (l - 1)]) - dt / dy * vn[k * ny + l] * (vn[k * ny + l] - vn[(k - 1) * ny + l]) + nu * dt / pow(dx, 2) * (vn[k * ny + (l + 1)] - 2 * vn[k * ny + l] + vn[k * ny + (l - 1)]) + nu * dt / pow(dy, 2) * (vn[(k + 1) * ny + l] - 2 * vn[k * ny + l] + vn[(k - 1) * ny + l]);

      }

    }

    for(i = 0; i < ny; i++){

      for(j = 0; j < nx; j++){

        printf("%f ", u[i * ny + j]);
    
      }

      for(j = 0; j < nx; j++){

        printf("%f ", v[i * ny + j]);

      }

      printf("\n");

    }

  }

  return 0;

}

float *reserva_ones(int n_points){

  float *array;

  int i;

  array = malloc(n_points * sizeof(float));

  for(i = 0; i < n_points; i++){

    array[i] = 1.0;

  }

  return array;

}
