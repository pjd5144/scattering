#include <math.h>
#include <complex.h>
#include <stdio.h>

int main() {
    double t_film = 800;
    double k0 = 2*M_PI/1.54;
    int num_slice = 150.0;
    double slice_thick = t_film/num_slice;
    double layers[num_slice];
    int i;
    double out;
    // initialize layers array values
    for (i = 0; i < num_slice+1; i++)
    {
        layers[i] = -i*slice_thick;
        printf("%f\n",layers[i]);
    }
    
}
