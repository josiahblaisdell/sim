#include <stdio.h>
#include <sim.h>




int main(){
    //N+2 = 7 => N=5
    static size = 49;
    //what are these statics? 
    //shouldn't "static" just limit the scope of the variable?
    static u[size], u_prev[size], v_prev[size];
    static dens[size];


    return 0;
}