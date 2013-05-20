#include <stdio.h>
#include <stdlib.h> // dynamic allocation
#include <math.h>

// Floating point precision
//typedef float Float;
typedef double Float;

// 3D vector
typedef struct {
    Float x;
    Float y;
    Float z;
} Float3;



//// SIMULATION PARAMETERS

// Number of dimensions
const int n = 3;

// Grid dims
const unsigned int nx = 3;
const unsigned int ny = 4;
const unsigned int nz = 2;

// Grid cell width
const Float dx = 2.0;

// Number of flow vectors in each cell
const int m = 19;

// Time step length
//const Float dt = 1.0;
const Float dt = 1.0e-4;

// Simulation end time
//const Float t_end = 2.0e-4;
const Float t_end = 1.0;

// Fluid dynamic viscosity
const Float nu = 8.9e-4;
//const Float nu = 8.9e-0;


// Gravitational acceleration
const Float3 g = {0.0, 0.0, -10.0};
//const Float3 g = {0.0, 0.0, 0.0};

// Lattice speed of sound for D2Q9 and D3Q19 (1.0/sqrt(3.0))
const Float c2_s = 1.0/1.7320508075688772;


//// FUNCTION DEFINITIONS

Float3 MAKE_FLOAT3(Float x, Float y, Float z)
{
    Float3 v;
    v.x = x; v.y = y; v.z = z;
    return v;
}

// Dot product of two Float3 vectors
Float dot(Float3 a, Float3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// Relaxation parameter, linked to fluid dynamic viscosity
Float tau(void) {

    //Float tau = (6.0*nu*dt/(dx*dx) + 1.0)/2.0;

    Float tau = nu/c2_s + 0.5;
    if (tau < 0.5) {
        fprintf(stderr, "Error: For positive viscosity: tau >= 0.5, ");
        fprintf(stderr, "but tau = %f.\n Increase the dynamic viscosity ", tau);
        fprintf(stderr, "(nu).\n");
        exit(1);
    } else
        return tau;
    //Float c2_s = 1.0/sqrtf(3.0);  // for D2Q9 and D3Q19
    //return nu/c2_s + 1.0/2.0;
}

// Get i-th value from cell x,y,z
unsigned int idx(
        unsigned int x,
        unsigned int y,
        unsigned int z)
{
    return x + nx*y + nx*ny*z;
}

// Get i-th value from cell x,y,z
unsigned int idxi(
        unsigned int x,
        unsigned int y,
        unsigned int z,
        unsigned int i)
{
    return x + ((y + z*ny)*nx) + (nx*ny*nz*i);
}

// Get i-th weight
Float w(unsigned int i)
{
    if (n == 3 && m == 19) {
        if (i == 0)
            return 1.0/3.0;
        else if (i > 0 && i < 7)
            return 1.0/18.0;
        else
            return 1.0/36.0;
    } else {
        fprintf(stderr, "Error in w: m = %d != 19", m);
        fprintf(stderr, ", n = %d != 3\n", n);
        exit(1);
    }
}

void set_e_values(Float3 *e)
{
    if (n == 3 && m == 19) {
        e[0] = MAKE_FLOAT3( 0.0, 0.0, 0.0); // zero vel.
        e[1] = MAKE_FLOAT3( 1.0, 0.0, 0.0); // face +x
        e[2] = MAKE_FLOAT3(-1.0, 0.0, 0.0); // face -x
        e[3] = MAKE_FLOAT3( 0.0, 1.0, 0.0); // face +y
        e[4] = MAKE_FLOAT3( 0.0,-1.0, 0.0); // face -y
        e[5] = MAKE_FLOAT3( 0.0, 0.0, 1.0); // face +z
        e[6] = MAKE_FLOAT3( 0.0, 0.0,-1.0); // face -z
        e[7] = MAKE_FLOAT3( 1.0, 1.0, 0.0); // edge +x,+y
        e[8] = MAKE_FLOAT3(-1.0,-1.0, 0.0); // edge -x,-y
        e[9] = MAKE_FLOAT3(-1.0, 1.0, 0.0); // edge -x,+y
        e[10] = MAKE_FLOAT3( 1.0,-1.0, 0.0); // edge +x,-y
        e[11] = MAKE_FLOAT3( 1.0, 0.0, 1.0); // edge +x,+z
        e[12] = MAKE_FLOAT3(-1.0, 0.0,-1.0); // edge -x,-z
        e[13] = MAKE_FLOAT3( 0.0, 1.0, 1.0); // edge +y,+z
        e[14] = MAKE_FLOAT3( 0.0,-1.0,-1.0); // edge -y,-z
        e[15] = MAKE_FLOAT3(-1.0, 0.0, 1.0); // edge -x,+z
        e[16] = MAKE_FLOAT3( 1.0, 0.0,-1.0); // edge +x,-z
        e[17] = MAKE_FLOAT3( 0.0,-1.0, 1.0); // edge -y,+z
        e[18] = MAKE_FLOAT3( 0.0, 1.0,-1.0); // edge +y,-z
    } else {
        fprintf(stderr, "Error in set_e_values: m = %d != 19", m);
        fprintf(stderr, ", n = %d != 3\n", n);
        exit(1);
    }
}

void init_fluid(Float* f, Float* rho, Float3* v)
{
    unsigned int x, y, z, i;
    const Float rho_init = 1.0;

    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                v[idx(x,y,z)].x = 0.0;
                v[idx(x,y,z)].y = 0.0;
                v[idx(x,y,z)].z = 0.0;
                rho[idx(x,y,z)] = rho_init;
                for (i=0; i<m; i++)
                    f[idxi(x,y,z,i)] = w(i) * rho_init;
            }
        }
    }
}

// Equilibrium distribution along flow vector e,
// Obtained from the local Maxwell-Boltzmann SPDF
// He and Luo, 1997
Float feq(
        Float rho,
        Float w,
        Float3 e,
        Float3 u)
{
    // Propagation speed on the lattice, squared
    Float c2 = dx/dt;
    return rho*w * (1.0 + 3.0*dot(e,u)/c2
            + 9.0/2.0*dot(e,u)*dot(e,u)/(c2*c2)
            - 3.0/2.0*dot(u,u)/c2);
}

// Bhatnagar-Gross-Kroop approximation collision operator,
// Bhatnagar et al., 1954
Float bgk(
        Float f,
        Float tau,
        Float rho,
        Float w,
        Float3 e,
        Float3 u)
{
    // Without gravitational drag
    //return f - (f - feq(rho, w, e, u))/tau;

    // With gravitational drag
    Float f_ext;    // External force along e
    Float m_f = dx*dx*dx*rho;   // Fluid mass
    Float3 f_g = {m_f*g.x, m_f*g.y, m_f*g.z}; // Gravitational force
    f_ext = dot(f_g, e);    // Drag force along e
    return f - (f - feq(rho, w, e, u))/tau
        + (2.0*tau - 1)/(2.0*tau)*3.0/w*f_ext;
}

// Cell fluid density
Float find_rho(
        Float* f,
        unsigned int x,
        unsigned int y,
        unsigned int z)
{
    int i;
    Float rho = 0.0;
    for (i=0; i<m; i++)
        rho += f[idxi(x,y,z,i)];
    return rho;
}

// Cell fluid velocity
Float3 find_v(
        Float* f,
        Float rho,
        Float3* e,
        unsigned int x,
        unsigned int y,
        unsigned int z)
{
    Float3 v = {0.0, 0.0, 0.0};
    Float f_i;
    unsigned int i;
    for (i=0; i<m; i++) {
        f_i = f[idxi(x,y,z,i)];
        v.x += f_i*e[i].x/rho;
        v.y += f_i*e[i].y/rho;
        v.z += f_i*e[i].z/rho;
    }
    return v;
}

// Lattice-Boltzmann collision step.
// Fluid distributions are modified towards the cell equilibrium.
// Values are read from f, and written to f_new.
void collide(
        Float* f,
        Float* rho,
        Float3* v,
        Float3* e)
{
    unsigned int x, y, z, i;
    Float rho_new;
    Float3 v_new;

    // Parallelize this with OpenMP
    // For each cell
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {

                // Calculate macroscopic parameters
                rho_new = find_rho(f, x, y, z);
                v_new = find_v(f, rho_new, e, x, y, z);

                // Find new f values by fluid particle collision
                for (i=0; i<m; i++) {
                    f[idxi(x,y,z,i)] =
                        bgk(f[idxi(x,y,z,i)], tau(), rho_new,
                                w(i), e[i], v_new);
                }

                // Store macroscopic parameters
                rho[idx(x,y,z)] = rho_new;
                v[idx(x,y,z)] = v_new;

            }
        }
    }
}

// Lattice-Boltzmann streaming step.
// Propagate fluid flows to cell neighbors.
// Boundary condition: Bounce back
void stream(Float* f, Float* f_new)
{

    unsigned int x, y, z;

    // For each cell
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                
                // Face 0
                f_new[idxi(x,y,z,0)] = fmax(0.0, f[idxi(x, y, z, 0)]);

                // Face 1 (+x): Bounce back
                if (x < nx-1)
                    f_new[idxi(x+1,  y,  z,  1)] = fmax(0.0, f[idxi(x, y, z, 1)]);
                else
                    f_new[idxi(  x,  y,  z,  2)] = fmax(0.0, f[idxi(x, y, z, 1)]);

                // Face 2 (-x): Bounce back
                if (x > 0)
                    f_new[idxi(x-1,  y,  z,  2)] = fmax(0.0, f[idxi(x, y, z, 2)]);
                else
                    f_new[idxi(  x,  y,  z,  1)] = fmax(0.0, f[idxi(x, y, z, 2)]);

                // Face 3 (+y): Bounce back
                if (y < ny-1)
                    f_new[idxi(  x,y+1,  z,  3)] = fmax(0.0, f[idxi(x, y, z, 3)]);
                else
                    f_new[idxi(  x,  y,  z,  4)] = fmax(0.0, f[idxi(x, y, z, 3)]);

                // Face 4 (-y): Bounce back
                if (y > 0)
                    f_new[idxi(  x,y-1,  z,  4)] = fmax(0.0, f[idxi(x, y, z, 4)]);
                else
                    f_new[idxi(  x,  y,  z,  3)] = fmax(0.0, f[idxi(x, y, z, 4)]);

                // Face 5 (+z): Bounce back
                if (z < nz-1)
                    f_new[idxi(  x,  y,z+1,  5)] = fmax(0.0, f[idxi(x, y, z, 5)]);
                else
                    f_new[idxi(  x,  y,  z,  6)] = fmax(0.0, f[idxi(x, y, z, 5)]);

                // Face 6 (-z): Bounce back
                if (z > 0)
                    f_new[idxi(  x,  y,z-1,  6)] = fmax(0.0, f[idxi(x, y, z, 6)]);
                else
                    f_new[idxi(  x,  y,  z,  5)] = fmax(0.0, f[idxi(x, y, z, 6)]);

                
                // Edge 7 (+x,+y): Bounce back
                if (x < nx-1 && y < ny-1)
                    f_new[idxi(x+1,y+1,  z,  7)] = fmax(0.0, f[idxi(x, y, z, 7)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z,  9)] = fmax(0.0, f[idxi(x, y, z, 7)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z, 10)] = fmax(0.0, f[idxi(x, y, z, 7)]);
                else
                    f_new[idxi(  x,  y,  z,  8)] = fmax(0.0, f[idxi(x, y, z, 7)]);

                // Edge 8 (-x,-y): Bounce back
                if (x > 0 && y > 0)
                    f_new[idxi(x-1,y-1,  z,  8)] = fmax(0.0, f[idxi(x, y, z, 8)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z,  9)] = fmax(0.0, f[idxi(x, y, z, 8)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z, 10)] = fmax(0.0, f[idxi(x, y, z, 8)]);
                else
                    f_new[idxi(  x,  y,  z,  7)] = fmax(0.0, f[idxi(x, y, z, 8)]);

                // Edge 9 (-x,+y): Bounce back
                if (x > 0 && y < ny-1)
                    f_new[idxi(x-1,y+1,  z,  9)] = fmax(0.0, f[idxi(x, y, z, 9)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z,  8)] = fmax(0.0, f[idxi(x, y, z, 9)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z,  7)] = fmax(0.0, f[idxi(x, y, z, 9)]);
                else
                    f_new[idxi(  x,  y,  z, 10)] = fmax(0.0, f[idxi(x, y, z, 9)]);

                // Edge 10 (+x,-y): Bounce back
                if (x < nx-1 && y > 0)
                    f_new[idxi(x+1,y-1,  z, 10)] = fmax(0.0, f[idxi(x, y, z, 10)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z,  8)] = fmax(0.0, f[idxi(x, y, z, 10)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z,  7)] = fmax(0.0, f[idxi(x, y, z, 10)]);
                else
                    f_new[idxi(  x,  y,  z,  9)] = fmax(0.0, f[idxi(x, y, z, 10)]);

                // Edge 11 (+x,+z): Bounce back
                if (x < nx-1 && z < nz-1)
                    f_new[idxi(x+1,  y,z+1, 11)] = fmax(0.0, f[idxi(x, y, z, 11)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z, 16)] = fmax(0.0, f[idxi(x, y, z, 11)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 15)] = fmax(0.0, f[idxi(x, y, z, 11)]);
                else
                    f_new[idxi(  x,  y,  z, 12)] = fmax(0.0, f[idxi(x, y, z, 11)]);

                // Edge 12 (-x,-z): Bounce back
                if (x > 0 && z > 0)
                    f_new[idxi(x-1,  y,z-1, 12)] = fmax(0.0, f[idxi(x, y, z, 12)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z, 15)] = fmax(0.0, f[idxi(x, y, z, 12)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 16)] = fmax(0.0, f[idxi(x, y, z, 12)]);
                else
                    f_new[idxi(  x,  y,  z, 11)] = fmax(0.0, f[idxi(x, y, z, 12)]);

                // Edge 13 (+y,+z): Bounce back
                if (y < ny-1 && z < nz-1)
                    f_new[idxi(  x,y+1,z+1, 13)] = fmax(0.0, f[idxi(x, y, z, 13)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z, 18)] = fmax(0.0, f[idxi(x, y, z, 13)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 17)] = fmax(0.0, f[idxi(x, y, z, 13)]);
                else
                    f_new[idxi(  x,  y,  z, 14)] = fmax(0.0, f[idxi(x, y, z, 13)]);

                // Edge 14 (-y,-z): Bounce back
                if (y > 0 && z > 0)
                    f_new[idxi(  x,y-1,z-1, 14)] = fmax(0.0, f[idxi(x, y, z, 14)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z, 17)] = fmax(0.0, f[idxi(x, y, z, 14)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 18)] = fmax(0.0, f[idxi(x, y, z, 14)]);
                else
                    f_new[idxi(  x,  y,  z, 13)] = fmax(0.0, f[idxi(x, y, z, 14)]);

                // Edge 15 (-x,+z): Bounce back
                if (x > 0 && z < nz-1)
                    f_new[idxi(x-1,  y,z+1, 15)] = fmax(0.0, f[idxi(x, y, z, 15)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z, 12)] = fmax(0.0, f[idxi(x, y, z, 15)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 11)] = fmax(0.0, f[idxi(x, y, z, 15)]);
                else
                    f_new[idxi(  x,  y,  z, 16)] = fmax(0.0, f[idxi(x, y, z, 15)]);

                // Edge 16 (+x,-z)
                if (x < nx-1 && z > 0)
                    f_new[idxi(x+1,  y,z-1, 16)] = fmax(0.0, f[idxi(x, y, z, 16)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z, 11)] = fmax(0.0, f[idxi(x, y, z, 16)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 12)] = fmax(0.0, f[idxi(x, y, z, 16)]);
                else
                    f_new[idxi(  x,  y,  z, 15)] = fmax(0.0, f[idxi(x, y, z, 16)]);

                // Edge 17 (-y,+z)
                if (y > 0 && z < nz-1)
                    f_new[idxi(  x,y-1,z+1, 17)] = fmax(0.0, f[idxi(x, y, z, 17)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z, 14)] = fmax(0.0, f[idxi(x, y, z, 17)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 13)] = fmax(0.0, f[idxi(x, y, z, 17)]);
                else
                    f_new[idxi(  x,  y,  z, 18)] = fmax(0.0, f[idxi(x, y, z, 17)]);

                // Edge 18 (+y,-z)
                if (y < ny-1 && z > 0)
                    f_new[idxi(  x,y+1,z-1, 18)] = fmax(0.0, f[idxi(x, y, z, 18)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z, 13)] = fmax(0.0, f[idxi(x, y, z, 18)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 14)] = fmax(0.0, f[idxi(x, y, z, 18)]);
                else
                    f_new[idxi(  x,  y,  z, 17)] = fmax(0.0, f[idxi(x, y, z, 18)]);
        
            }
        }
    }
}

// Swap Float pointers
void swapFloats(Float* a, Float* b)
{
    Float* tmp = a;
    a = b;
    b = tmp;
}

// Print density values to file stream (stdout, stderr, other file)
void print_rho(Float* rho, FILE* stream)
{
    unsigned int x, y, z;
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                fprintf(stream, "%f\t", rho[idx(x,y,z)]);
            }
            fprintf(stream, "\n");
        }
        fprintf(stream, "\n");
    }
}

// Print velocity values to file stream (stdout, stderr, other file)
void print_v(Float3* v, FILE* stream)
{
    unsigned int x, y, z;
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                fprintf(stream, "%.1ex%.1ex%.1e\t",
                        v[idx(x,y,z)].x,
                        v[idx(x,y,z)].y,
                        v[idx(x,y,z)].z);
            }
            fprintf(stream, "\n");
        }
        fprintf(stream, "\n");
    }
}


int main(int argc, char** argv)
{
    printf("### Lattice-Boltzman D%dQ%d test ###\n", n, m);


    // Print parameter vals
    unsigned int ncells = nx*ny*nz;
    printf("Grid dims: nx = %d, ny = %d, nz = %d: %d cells\n",
            nx, ny, nz, ncells);
    printf("tau = %f\n", tau());

    // Set cell flow vector values
    Float3 e[m]; set_e_values(e);

    // Particle distributions
    Float* f = malloc(ncells*m*sizeof(Float));
    Float* f_new = malloc(ncells*m*sizeof(Float));

    // Cell densities
    Float* rho = malloc(ncells*sizeof(Float));

    // Cell flow velocities
    Float3* v = malloc(ncells*sizeof(Float3));

    init_fluid(f, rho, v);


    double t;
    for (t = 0.0; t < t_end; t += dt) {
        collide(f, rho, v, e);
        stream(f, f_new);
        swapFloats(f, f_new);
    }

    fprintf(stdout, "rho\n");
    print_rho(rho, stdout);

    fprintf(stdout, "v\n");
    print_v(v, stdout);

    free(f);
    free(f_new);
    free(rho);
    free(v);

    return 0;
}
