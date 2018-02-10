#include <stdio.h>
#include <E:\projects\sim\sim.h>
#define size 49
int main(){

    //what are these statics? 
    //shouldn't "static" just limit the scope of the variable?

    // In this context (in main), "static" will have no meaningful effect.
    // In the function "f" below, the static variable x will stay around after the function call, so that it increments with each call.
    // For local variables like these, static means that the variable doesn't live in the function's "stack" 
    // (temporary memory that is created for a function call and thrown away when the call completes)
    // but lives outside in memory allocated to the program itself.  
    // However, its scope is still limited to the local context.  So you couldn't access these variable outside main().
    // Since main() is not a function you're going to call anywhere, static doesn't have any practical effect here.

    // Also, I'm surprised C lets you get away without specifying the data type of these variables.  Does this compile?

    static float u[size], u_prev[size], v[size], v_prev[size];
    static float dens[size], dens_prev[size];
    float dt = .05;
    float diff = 0.1;
    float visc = 0;
    float force = 1;
    float source = 10;
    memset(u_prev,0,49);
    memset(u,0,49);
    memset(v_prev,0,49);
    memset(v,0,49);
    memset(dens_prev,0,49);
    memset(dens,0,49);

    vel_step(5,u,v,u_prev,v_prev,visc,dt);
    for(int i=0;i<=49;i++){
        printf("%f, ",u[i]);
    }
    int x = 1;
    return x;
}

void f() {
    static int x = 0; // the initialization will only happen once.
    ++x; // x will continuously increment with each call to f().
}