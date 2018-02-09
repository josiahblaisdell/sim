#include <stdio.h>
#define IX(i,j) ((i)+(7)*(j))
void add_source(int, float*, float*, float);
void diffuse( int, int, float*, float*, float, float);
void advect(int, int, float*, float*, float*, float*, float);

//Adds the source to the density
void add_source(int N, float* x, float* s, float dt){
    int i, size = 49;
    for(i=0; i<size; i++) {
        x[i] += dt*s[i];
    }
}

void diffuse(int N, float* x, float* x0, float diff, float dt){
    int i,j,k;
    float a = dt*diff*N*N;
    for(k=0; k<20; k++){
        for(i=1; i<=N; i++){
            for(j=1; j<=N; j++){
                x[IX(i,j)] += a*(x[IX(i-1,j)] + x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)])/(1+4*a);

            }
        }
    }
    set_bnd(N,b,x);
}

void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt){
    int i,j,i0,j0,i1,j1;
    float x,y,s0,t0,s1,t1,dt0;
    dt0 = dt*N;
    for(i=1; i<=N; i++){
        for(j=1;j<=N;j++){
            x = i-dt0*u[IX(i,j)];
            y = j-dt0*v[IX(i,j)];
            if(x<.5) x = .5;
            if(x>N+.5) x = N+.5;
            i0 = (int)x;
            i1 = i0+1;
            j0 = (int)y;
            j1 = j0+1;
            s1 = x-i0;
            s0 = 1-s1;
            t1 = y-j0;
            t0 = 1-t1;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)]) 
                            + s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i0,j1)]);
        }
    }
    set_bnd(N,b,d);
}