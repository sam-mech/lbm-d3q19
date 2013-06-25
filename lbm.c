#include <stdio.h>
#include <stdlib.h> // dynamic allocation
#include <math.h>

// Courant kriterie for tidsskridt
// No slip top og bund
// Periodiske sider

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
//const unsigned int nx = 3;
//const unsigned int ny = 6;
//const unsigned int nz = 3;
const unsigned int nx = 37;
const unsigned int ny = 37;
const unsigned int nz = 37;

// Grid cell width
const Float dx = 1.0;

// Number of flow vectors in each cell
const int m = 19;

// Time step length
//const double dt = 1.0;
const double dt = 1.0e-3;
//const double dt = 0.01;

// Simulation end time
//const Float t_end = 1.5e-4;
const double t_end = 2.0;
//const double t_end = 1.0;
//const Float t_end = 10.1;

const double t_file = 0.01;

// Fluid dynamic viscosity
const Float nu = 8.9e-4;

// Gravitational acceleration
//const Float3 g = {0.0, 0.0, -10.0};
const Float3 g = {0.0, 0.0, 0.0};

// Initial cell fluid density (dimensionless)
const Float rho0 = 1.0;

// Inital cell fluid velocity (dimensionless)
const Float3 u0 = {0.0, 0.0, 0.0};

// Courant criteria limit
const Float C_max = 1.0;


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

// Viscosity parameter
Float tau(void) {
    return (6.0*nu*dt/(dx*dx) + 1.0)/2.0;
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
        exit(EXIT_FAILURE);
    }
}

void set_e_values(Float3 *e)
{
    if (n == 3 && m == 19) {
        e[0]  = MAKE_FLOAT3( 0.0, 0.0, 0.0); // zero vel.
        e[1]  = MAKE_FLOAT3( 1.0, 0.0, 0.0); // face +x
        e[2]  = MAKE_FLOAT3(-1.0, 0.0, 0.0); // face -x
        e[3]  = MAKE_FLOAT3( 0.0, 1.0, 0.0); // face +y
        e[4]  = MAKE_FLOAT3( 0.0,-1.0, 0.0); // face -y
        e[5]  = MAKE_FLOAT3( 0.0, 0.0, 1.0); // face +z
        e[6]  = MAKE_FLOAT3( 0.0, 0.0,-1.0); // face -z
        e[7]  = MAKE_FLOAT3( 1.0, 1.0, 0.0); // edge +x,+y
        e[8]  = MAKE_FLOAT3(-1.0,-1.0, 0.0); // edge -x,-y
        e[9]  = MAKE_FLOAT3(-1.0, 1.0, 0.0); // edge -x,+y
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
        exit(EXIT_FAILURE);
    }
}

// Equilibrium distribution along flow vector e
Float feq(
        Float rho,
        Float w,
        Float3 e,
        Float3 u)
{
    Float c2 = dx/dt;
    return rho*w * (1.0 + 3.0/c2*dot(e,u)
            + 9.0/(2.0*c2*c2)*dot(e,u)*dot(e,u)
            - 3.0/(2.0*c2)*dot(u,u)*dot(u,u));
}

// Initialize cell densities, velocities, and flow vectors
void init_rho_v(Float* rho, Float3* u)
{
    unsigned int x, y, z;
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {

                // Set velocity to u0
                u[idx(x,y,z)] = MAKE_FLOAT3(u0.x, u0.y, u0.z);

                // Set density to rho0
                rho[idx(x,y,z)] = rho0;
            }
        }
    }
}

void init_f(Float* f, Float* f_new, Float* rho, Float3* u, Float3* e)
{
    unsigned int x, y, z, i;
    Float f_val;

    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                for (i=0; i<m; i++) {

                    // Set fluid flow vectors to v0
                    f_val = feq(rho[idx(x,y,z)], w(i), e[i], u[idx(x,y,z)]);
                    f[idxi(x,y,z,i)] = f_val;
                    f_new[idxi(x,y,z,i)] = f_val;
                }
            }
        }
    }
}

// Bhatnagar-Gross-Kroop approximation collision operator
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
Float3 find_u(
        Float* f,
        Float rho,
        Float3* e,
        unsigned int x,
        unsigned int y,
        unsigned int z)
{
    Float3 u = {0.0, 0.0, 0.0};
    Float f_i;
    unsigned int i;
    for (i=0; i<m; i++) {
        f_i = f[idxi(x,y,z,i)];
        u.x += f_i*e[i].x/rho;
        u.y += f_i*e[i].y/rho;
        u.z += f_i*e[i].z/rho;
    }

    // Check the Courant-Frederichs-Lewy condition
    if ((u.x*dt/dx + u.y*dt/dx + u.z*dt/dx) > C_max) {
        fprintf(stderr, "Error, the Courant-Friderichs-Lewy condition is not ");
        fprintf(stderr, "satisfied.\nTry one or more of the following:\n");
        fprintf(stderr, "- Decrease the timestep (dt)\n");
        fprintf(stderr, "- Increase the cell size (dx)\n");
        fprintf(stderr, "- Decrease the fluid viscosity (nu)\n");
        fprintf(stderr, "- Decrease the fluid density (rho)\n");
        exit(EXIT_FAILURE);
    }

    return u;
}

// Lattice-Boltzmann collision step.
// Fluid distributions are modified towards the cell equilibrium.
// Values are read from f, and written to rho and u.
void collide(
        Float* f,
        Float* rho,
        Float3* u,
        Float3* e)
{
    unsigned int x, y, z, i;
    Float rho_new;
    Float3 u_new;

    // Parallelize this with OpenMP
    // For each cell
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {

                // Calculate macroscopic parameters
                rho_new = find_rho(f, x, y, z);
                u_new = find_u(f, rho_new, e, x, y, z);

                // Store macroscopic parameters
                rho[idx(x,y,z)] = rho_new;
                u[idx(x,y,z)] = u_new;

                // Find new f values by fluid particle collision
                for (i=0; i<m; i++) {
                    f[idxi(x,y,z,i)] =
                        bgk(f[idxi(x,y,z,i)], tau(), rho_new,
                                w(i), e[i], u_new);
                }
            }
        }
    }
}

// Lattice-Boltzmann streaming step.
// Propagate fluid flows to cell neighbors.
// Boundary condition: Bounce back
void stream(Float* f, Float* f_new)
{
    // For each cell
    unsigned int x, y, z;
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                
                // Face 0
                f_new[idxi(x,y,z,0)] = fmax(0.0, f[idxi(x, y, z, 0)]);

                // Face 1 (+x): Bounce back
                if (x < nx-1)
                    f_new[idxi(x+1,  y,  z,  1)]
                        = fmax(0.0, f[idxi(x, y, z, 1)]);
                else
                    f_new[idxi(  x,  y,  z,  2)]
                        = fmax(0.0, f[idxi(x, y, z, 1)]);

                // Face 2 (-x): Bounce back
                if (x > 0)
                    f_new[idxi(x-1,  y,  z,  2)]
                        = fmax(0.0, f[idxi(x, y, z, 2)]);
                else
                    f_new[idxi(  x,  y,  z,  1)]
                        = fmax(0.0, f[idxi(x, y, z, 2)]);

                // Face 3 (+y): Bounce back
                if (y < ny-1)
                    f_new[idxi(  x,y+1,  z,  3)]
                        = fmax(0.0, f[idxi(x, y, z, 3)]);
                else
                    f_new[idxi(  x,  y,  z,  4)]
                        = fmax(0.0, f[idxi(x, y, z, 3)]);

                // Face 4 (-y): Bounce back
                if (y > 0)
                    f_new[idxi(  x,y-1,  z,  4)]
                        = fmax(0.0, f[idxi(x, y, z, 4)]);
                else
                    f_new[idxi(  x,  y,  z,  3)]
                        = fmax(0.0, f[idxi(x, y, z, 4)]);

                // Face 5 (+z): Bounce back
                if (z < nz-1)
                    f_new[idxi(  x,  y,z+1,  5)]
                        = fmax(0.0, f[idxi(x, y, z, 5)]);
                else
                    f_new[idxi(  x,  y,  z,  6)]
                        = fmax(0.0, f[idxi(x, y, z, 5)]);

                // Face 6 (-z): Bounce back
                if (z > 0)
                    f_new[idxi(  x,  y,z-1,  6)]
                        = fmax(0.0, f[idxi(x, y, z, 6)]);
                else
                    f_new[idxi(  x,  y,  z,  5)]
                        = fmax(0.0, f[idxi(x, y, z, 6)]);

                
                // Edge 7 (+x,+y): Bounce back
                if (x < nx-1 && y < ny-1)
                    f_new[idxi(x+1,y+1,  z,  7)]
                        = fmax(0.0, f[idxi(x, y, z, 7)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z,  9)]
                        = fmax(0.0, f[idxi(x, y, z, 7)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z, 10)]
                        = fmax(0.0, f[idxi(x, y, z, 7)]);
                else
                    f_new[idxi(  x,  y,  z,  8)]
                        = fmax(0.0, f[idxi(x, y, z, 7)]);

                // Edge 8 (-x,-y): Bounce back
                if (x > 0 && y > 0)
                    f_new[idxi(x-1,y-1,  z,  8)]
                        = fmax(0.0, f[idxi(x, y, z, 8)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z,  9)]
                        = fmax(0.0, f[idxi(x, y, z, 8)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z, 10)]
                        = fmax(0.0, f[idxi(x, y, z, 8)]);
                else
                    f_new[idxi(  x,  y,  z,  7)]
                        = fmax(0.0, f[idxi(x, y, z, 8)]);

                // Edge 9 (-x,+y): Bounce back
                if (x > 0 && y < ny-1)
                    f_new[idxi(x-1,y+1,  z,  9)]
                        = fmax(0.0, f[idxi(x, y, z, 9)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z,  8)]
                        = fmax(0.0, f[idxi(x, y, z, 9)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z,  7)]
                        = fmax(0.0, f[idxi(x, y, z, 9)]);
                else
                    f_new[idxi(  x,  y,  z, 10)]
                        = fmax(0.0, f[idxi(x, y, z, 9)]);

                // Edge 10 (+x,-y): Bounce back
                if (x < nx-1 && y > 0)
                    f_new[idxi(x+1,y-1,  z, 10)]
                        = fmax(0.0, f[idxi(x, y, z, 10)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z,  8)]
                        = fmax(0.0, f[idxi(x, y, z, 10)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z,  7)]
                        = fmax(0.0, f[idxi(x, y, z, 10)]);
                else
                    f_new[idxi(  x,  y,  z,  9)]
                        = fmax(0.0, f[idxi(x, y, z, 10)]);

                // Edge 11 (+x,+z): Bounce back
                if (x < nx-1 && z < nz-1)
                    f_new[idxi(x+1,  y,z+1, 11)]
                        = fmax(0.0, f[idxi(x, y, z, 11)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z, 16)]
                        = fmax(0.0, f[idxi(x, y, z, 11)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 15)]
                        = fmax(0.0, f[idxi(x, y, z, 11)]);
                else
                    f_new[idxi(  x,  y,  z, 12)]
                        = fmax(0.0, f[idxi(x, y, z, 11)]);

                // Edge 12 (-x,-z): Bounce back
                if (x > 0 && z > 0)
                    f_new[idxi(x-1,  y,z-1, 12)]
                        = fmax(0.0, f[idxi(x, y, z, 12)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z, 15)]
                        = fmax(0.0, f[idxi(x, y, z, 12)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 16)]
                        = fmax(0.0, f[idxi(x, y, z, 12)]);
                else
                    f_new[idxi(  x,  y,  z, 11)]
                        = fmax(0.0, f[idxi(x, y, z, 12)]);

                // Edge 13 (+y,+z): Bounce back
                if (y < ny-1 && z < nz-1)
                    f_new[idxi(  x,y+1,z+1, 13)]
                        = fmax(0.0, f[idxi(x, y, z, 13)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z, 18)]
                        = fmax(0.0, f[idxi(x, y, z, 13)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 17)]
                        = fmax(0.0, f[idxi(x, y, z, 13)]);
                else
                    f_new[idxi(  x,  y,  z, 14)]
                        = fmax(0.0, f[idxi(x, y, z, 13)]);

                // Edge 14 (-y,-z): Bounce back
                if (y > 0 && z > 0)
                    f_new[idxi(  x,y-1,z-1, 14)]
                        = fmax(0.0, f[idxi(x, y, z, 14)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z, 17)]
                        = fmax(0.0, f[idxi(x, y, z, 14)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 18)]
                        = fmax(0.0, f[idxi(x, y, z, 14)]);
                else
                    f_new[idxi(  x,  y,  z, 13)]
                        = fmax(0.0, f[idxi(x, y, z, 14)]);

                // Edge 15 (-x,+z): Bounce back
                if (x > 0 && z < nz-1)
                    f_new[idxi(x-1,  y,z+1, 15)]
                        = fmax(0.0, f[idxi(x, y, z, 15)]);
                else if (x > 0)
                    f_new[idxi(x-1,  y,  z, 12)]
                        = fmax(0.0, f[idxi(x, y, z, 15)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 11)]
                        = fmax(0.0, f[idxi(x, y, z, 15)]);
                else
                    f_new[idxi(  x,  y,  z, 16)]
                        = fmax(0.0, f[idxi(x, y, z, 15)]);

                // Edge 16 (+x,-z)
                if (x < nx-1 && z > 0)
                    f_new[idxi(x+1,  y,z-1, 16)]
                        = fmax(0.0, f[idxi(x, y, z, 16)]);
                else if (x < nx-1)
                    f_new[idxi(x+1,  y,  z, 11)]
                        = fmax(0.0, f[idxi(x, y, z, 16)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 12)]
                        = fmax(0.0, f[idxi(x, y, z, 16)]);
                else
                    f_new[idxi(  x,  y,  z, 15)]
                        = fmax(0.0, f[idxi(x, y, z, 16)]);

                // Edge 17 (-y,+z)
                if (y > 0 && z < nz-1)
                    f_new[idxi(  x,y-1,z+1, 17)]
                        = fmax(0.0, f[idxi(x, y, z, 17)]);
                else if (y > 0)
                    f_new[idxi(  x,y-1,  z, 14)]
                        = fmax(0.0, f[idxi(x, y, z, 17)]);
                else if (z < nz-1)
                    f_new[idxi(  x,  y,z+1, 13)]
                        = fmax(0.0, f[idxi(x, y, z, 17)]);
                else
                    f_new[idxi(  x,  y,  z, 18)]
                        = fmax(0.0, f[idxi(x, y, z, 17)]);

                // Edge 18 (+y,-z)
                if (y < ny-1 && z > 0)
                    f_new[idxi(  x,y+1,z-1, 18)]
                        = fmax(0.0, f[idxi(x, y, z, 18)]);
                else if (y < ny-1)
                    f_new[idxi(  x,y+1,  z, 13)]
                        = fmax(0.0, f[idxi(x, y, z, 18)]);
                else if (z > 0)
                    f_new[idxi(  x,  y,z-1, 14)]
                        = fmax(0.0, f[idxi(x, y, z, 18)]);
                else
                    f_new[idxi(  x,  y,  z, 17)]
                        = fmax(0.0, f[idxi(x, y, z, 18)]);
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
void print_rho(FILE* stream, Float* rho)
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

// Print velocity values from y-plane to file stream
void print_rho_yplane(FILE* stream, Float* rho, unsigned int y)
{
    unsigned int x, z;
    for (z=0; z<nz; z++) {
        for (x=0; x<nx; x++) {
            fprintf(stream, "%f\t", rho[idx(x,y,z)]);
        }
        fprintf(stream, "\n");
    }
}


// Print velocity values to file stream (stdout, stderr, other file)
void print_u(FILE* stream, Float3* u)
{
    unsigned int x, y, z;
    for (z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
            for (x=0; x<nx; x++) {
                fprintf(stream, "%.1ex%.1ex%.1e\t",
                        u[idx(x,y,z)].x,
                        u[idx(x,y,z)].y,
                        u[idx(x,y,z)].z);
            }
            fprintf(stream, "\n");
        }
        fprintf(stream, "\n");
    }
}

// Print velocity values from y-plane to file stream
void print_u_yplane(FILE* stream, Float3* u, unsigned int y)
{
    unsigned int x, z;
    for (z=0; z<nz; z++) {
        for (x=0; x<nx; x++) {
            fprintf(stream, "%.1ex%.1ex%.1e\t",
                    u[idx(x,y,z)].x,
                    u[idx(x,y,z)].y,
                    u[idx(x,y,z)].z);
        }
        fprintf(stream, "\n");
    }
}


int main(int argc, char** argv)
{
    printf("### Lattice-Boltzman D%dQ%d test ###\n", n, m);

    FILE* frho;
    char filename[40];

    // Print parameter vals
    //printf("Grid dims: nx = %d, ny = %d, nz = %d: %d cells\n",
            //nx, ny, nz, ncells);

    // Set cell flow vector values
    Float3 e[m]; set_e_values(e);

    // Particle distributions
    unsigned int ncells = nx*ny*nz;
    Float* f = malloc(ncells*m*sizeof(Float));
    Float* f_new = malloc(ncells*m*sizeof(Float));

    // Cell densities
    Float* rho = malloc(ncells*sizeof(Float));

    // Cell flow velocities
    Float3* u = malloc(ncells*sizeof(Float3));

    // Set densities, velocities and flow vectors
    init_rho_v(rho, u);
    rho[idx(nx/2,ny/2,nz/2)] *= 1.0001;
    init_f(f, f_new, rho, u, e);

    // Temporal loop
    double t = 0.0;
    double t_file_elapsed = 0.0;

    // Save initial state
    sprintf(filename, "out/rho_y%d_t%.2f.txt", ny/2, t);
    if ((frho = fopen(filename, "w"))) {
        print_rho_yplane(frho, rho, ny/2);
        fclose(frho);
    } else {
        fprintf(stderr, "Error: Could not open output file ");
        fprintf(stderr, filename);
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
    }

    // Temporal loop
    for (t = 0.0; t < t_end; t += dt, t_file_elapsed += dt) {

        // Report time to stdout
        printf("\rt = %.1fs./%.1fs., %.1f%% done", t, t_end, t/t_end*100.0);

        // LBM collision and streaming
        collide(f, rho, u, e);
        stream(f, f_new);

        // Swap f and f_new
        Float* tmp = f;
        f = f_new;
        f_new = tmp;

        // Print x-z plane to file
        if (t_file_elapsed >= t_file) {
            sprintf(filename, "out/rho_y%d_t%.2f.txt", ny/2, t);
            if ((frho = fopen(filename, "w"))) {
                print_rho_yplane(frho, rho, ny/2);
                fclose(frho);
            } else {
                fprintf(stderr, "Error: Could not open output file ");
                fprintf(stderr, filename);
                fprintf(stderr, "\n");
                exit(EXIT_FAILURE);
            }
            t_file_elapsed = 0.0;
        }
    }
    printf("\n");

    // Report values to stdout
    //fprintf(stdout, "rho\n");
    //print_rho(stdout, rho);
    //fprintf(stdout, "u\n");
    //print_u(stdout, u);

    // Clear memory
    free(f);
    free(f_new);
    free(rho);
    free(u);

    return EXIT_SUCCESS;
}
