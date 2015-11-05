#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define USAGE "./a.out valor"

double *reserva_ones(int n_points);

double pi = 3.14159265359;

int main(int argc, char **argv){

  double K_0;
  double alpha;

  if(argc != 3){

    printf("USAGE: %s\n", USAGE);

    exit(1);

  }

  K_0 = atof(argv[1]);

  alpha = atof(argv[2]);
  
  double R_T = 6378.1; //Km
  double B_0 = 3.0 * pow(10.0, -5); //Teslas
  double m = 938.272013; //MeV/c^2
  double q = 1.0; //carga e
  double c = 299792.458; //Km/s

  double dt = 0.01; //s
  int nt = (int)(100.0 / dt);
  int i;

  double x_0 = 2.0 * R_T;
  double y_0 = 0.0;
  double z_0 = 0.0;

  double a = (double) alpha * pi / 180.0;
    
  double *t;
  double *x;
  double *y;
  double *z;

  t = reserva_ones(nt + 1);
  x = reserva_ones(nt + 1);
  y = reserva_ones(nt + 1);
  z = reserva_ones(nt + 1);
  
  t[0] = 0.0;
  x[0] = x_0;
  y[0] = y_0;
  z[0] = z_0;
    
  double gamma_0 = K_0 / m + 1.0;
    
  double v_0 = c * sqrt(1 - pow(gamma_0, -2)); //Km/s

  double vx = 0.0;
  double vy = v_0 * sin(a);
  double vz = v_0 * cos(a);

  double vx_anterior;
  double vy_anterior;
  double vz_anterior;

  double r;
  double cte;

  double Bx;
  double By;
  double Bz;

  double delta_x;
  double delta_y;
  double delta_z;
    
  for(i = 1; i < (nt + 1); i++){
        
    t[i] = t[i - 1] + dt;
        
    x[i] = x[i - 1] + vx * dt;
        
    y[i] = y[i - 1] + vy * dt;
        
    z[i] = z[i - 1] + vz * dt;
        
    r = sqrt(pow(x[i - 1], 2) + pow(y[i - 1], 2) + pow(z[i - 1], 2));
        
    cte = -B_0 * pow(R_T, 3) / pow(r, 5);
        
    Bx = cte * 3.0 * x[i - 1] * z[i - 1];
        
    By = cte * 3.0 * y[i - 1] * z[i - 1];
        
    Bz = cte * (2.0 * pow(z[i - 1], 2) - pow(x[i - 1], 2) - pow(y[i - 1], 2));
        
    delta_x = q * (vy * Bz - vz * By);
        
    delta_y = q * (-vx * Bz + vz * Bx);
        
    delta_z = q * (vx * By - vy * Bx);
        
    vx_anterior = vx;
        
    vy_anterior = vy;
        
    vz_anterior = vz;
        
    vx = delta_x * dt / (gamma_0 * m) + vx_anterior;
        
    vy = delta_y * dt / (gamma_0 * m) + vy_anterior;
        
    vz = delta_z * dt / (gamma_0 * m) + vz_anterior;

  }

  for(i = 0; i < nt; i++){

    printf("%f %f %f %f\n", t[i], x[i], y[i], z[i]);

  }

}

double *reserva_ones(int n_points){

  double *array;

  int i;

  array = malloc(n_points * sizeof(double));

  for(i = 0; i < n_points; i++){

    array[i] = 1.0;

  }

  return array;

}
