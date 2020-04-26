#include <iostream>
#include <cmath>
#include <SFML/Audio.hpp>
#include <SFML/Graphics.hpp>


using namespace std;

int size = 64;
int iter = 4;


#define IX(x, y, z) ((x) + (y) * size + (z) * size * size)


void set_bnd(bool b, float *x);
void lin_solve(int b, float *x, float *x0, float a, float c);
void diffuse (int b, float *x, float *x0, float diff, float dt);
void project(float *velocX, float *velocY, float *velocZ, float *p, float *div);
void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float *velocZ, float dt);

class Fluid{
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    public :float *density;
    
    float *Vx;
    float *Vy;
    float *Vz;

    float *Vx0;
    float *Vy0;
    float *Vz0;
    
    public : Fluid(int size0 ,float diffusion0,float viscosity0 ,float dt0){
        int N =size0;
        size = size0;
        diff = diffusion0;
        visc = viscosity0 ;
        dt = dt0;
        
        s = (float*)calloc(N * N * N, sizeof(float));
        density = (float*)calloc(N * N * N, sizeof(float));
        
        Vx = (float*)calloc(N * N * N, sizeof(float));
        Vy = (float*)calloc(N * N * N, sizeof(float));
        Vz = (float*)calloc(N * N * N, sizeof(float));
        
        Vx0 = (float*)calloc(N * N * N, sizeof(float));
        Vy0 = (float*)calloc(N * N * N, sizeof(float));
        Vz0 = (float*)calloc(N * N * N, sizeof(float));
        
    }
    public :~Fluid(){
        delete []s;
        delete []density;
        
        delete []Vx;
        delete []Vy;
        delete []Vz;
        
        delete []Vx0;
        delete []Vy0;
        delete []Vz0;
        
    }
    
    
    public : void addDensity(int x, int y, int z, float amount)
    {
        density[IX(x, y, z)] += amount;
    }
    public : void addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ)
    {
        int index = IX(x, y, z);
        
        Vx[index] += amountX;
        Vy[index] += amountY;
        Vz[index] += amountZ;
    }
    public : void Step(){
        
        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);
        diffuse(3, Vz0, Vz, visc, dt);
        
        project(Vx0, Vy0, Vz0, Vx, Vy);
        
        advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt);
        advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt);
        
        project(Vx, Vy, Vz, Vx0, Vy0);
        
        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, Vz,dt);
    }
    
};

// b would be 0 for x dimension and 1 for y

void set_bnd(int b, float *x){
    
    int N = size;
    
    for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0  )] = b == 3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
            x[IX(i, j, N-1)] = b == 3 ? -x[IX(i, j, N-2)] : x[IX(i, j, N-2)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  , k)] = b == 2 ? -x[IX(i, 1  , k)] : x[IX(i, 1  , k)];
            x[IX(i, N-1, k)] = b == 2 ? -x[IX(i, N-2, k)] : x[IX(i, N-2, k)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j, k)] = b == 1 ? -x[IX(1  , j, k)] : x[IX(1  , j, k)];
            x[IX(N-1, j, k)] = b == 1 ? -x[IX(N-2, j, k)] : x[IX(N-2, j, k)];
        }
    }
    
    x[IX(0, 0, 0)]       = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
    x[IX(0, N-1, 0)]     = 0.33f * (x[IX(1, N-1, 0)]
                                  + x[IX(0, N-2, 0)]
                                  + x[IX(0, N-1, 1)]);
    x[IX(0, 0, N-1)]     = 0.33f * (x[IX(1, 0, N-1)]
                                  + x[IX(0, 1, N-1)]
                                  + x[IX(0, 0, N)]);
    x[IX(0, N-1, N-1)]   = 0.33f * (x[IX(1, N-1, N-1)]
                                  + x[IX(0, N-2, N-1)]
                                  + x[IX(0, N-1, N-2)]);
    x[IX(N-1, 0, 0)]     = 0.33f * (x[IX(N-2, 0, 0)]
                                  + x[IX(N-1, 1, 0)]
                                  + x[IX(N-1, 0, 1)]);
    x[IX(N-1, N-1, 0)]   = 0.33f * (x[IX(N-2, N-1, 0)]
                                  + x[IX(N-1, N-2, 0)]
                                  + x[IX(N-1, N-1, 1)]);
    x[IX(N-1, 0, N-1)]   = 0.33f * (x[IX(N-2, 0, N-1)]
                                  + x[IX(N-1, 1, N-1)]
                                  + x[IX(N-1, 0, N-2)]);
    x[IX(N-1, N-1, N-1)] = 0.33f * (x[IX(N-2, N-1, N-1)]
                                  + x[IX(N-1, N-2, N-1)]
                                  + x[IX(N-1, N-1, N-2)]);
    
    
    
}

void lin_solve(int b, float *x, float *x0, float a, float c)
{
    float cRecip = 1.0 / c;
    int N = size;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m)] =
                        (x0[IX(i, j, m)]
                            + a*(    x[IX(i+1, j  , m  )]
                                    +x[IX(i-1, j  , m  )]
                                    +x[IX(i  , j+1, m  )]
                                    +x[IX(i  , j-1, m  )]
                                    +x[IX(i  , j  , m+1)]
                                    +x[IX(i  , j  , m-1)]
                           )) * cRecip;
                }
            }
        }
        set_bnd(b, x);
    }
}


void diffuse (int b, float *x, float *x0, float diff, float dt)
{
    int N = size; 
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

void project(float *velocX, float *velocY, float *velocZ, float *p, float *div)
{
    int  N = size;
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k)] = -0.5f*(
                         velocX[IX(i+1, j  , k  )]
                        -velocX[IX(i-1, j  , k  )]
                        +velocY[IX(i  , j+1, k  )]
                        -velocY[IX(i  , j-1, k  )]
                        +velocZ[IX(i  , j  , k+1)]
                        -velocZ[IX(i  , j  , k-1)]
                    )/N;
                p[IX(i, j, k)] = 0;
            }
        }
    }
    
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);
    
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k)] -= 0.5f * (  p[IX(i+1, j, k)]
                                                -p[IX(i-1, j, k)]) * N;
                velocY[IX(i, j, k)] -= 0.5f * (  p[IX(i, j+1, k)]
                                                -p[IX(i, j-1, k)]) * N;
                velocZ[IX(i, j, k)] -= 0.5f * (  p[IX(i, j, k+1)]
                                                -p[IX(i, j, k-1)]) * N;
            }
        }
    }
    
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
    set_bnd(3, velocZ);
}

void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float *velocZ, float dt)
{
    int N=size;
    
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;
    
    float Nfloat = N;
    float ifloat, jfloat, kfloat;
    int i, j, k;
    
    for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j, k)];
                tmp2 = dty * velocY[IX(i, j, k)];
                tmp3 = dtz * velocZ[IX(i, j, k)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                z    = kfloat - tmp3;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floorf(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floorf(y);
                j1 = j0 + 1.0f; 
                if(z < 0.5f) z = 0.5f;
                if(z > Nfloat + 0.5f) z = Nfloat + 0.5f;
                k0 = floorf(z);
                k1 = k0 + 1.0f;
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                u1 = z - k0;
                u0 = 1.0f - u1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;
                int k0i = k0;
                int k1i = k1;
                
                d[IX(i, j, k)] = 
                
                    s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i)]
                                +u1 * d0[IX(i0i, j0i, k1i)])
                        +( t1 * (u0 * d0[IX(i0i, j1i, k0i)]
                                +u1 * d0[IX(i0i, j1i, k1i)])))
                   +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i)]
                                +u1 * d0[IX(i1i, j0i, k1i)])
                        +( t1 * (u0 * d0[IX(i1i, j1i, k0i)]
                                +u1 * d0[IX(i1i, j1i, k1i)])));
            }
        }
    }
    set_bnd(b, d);
}

float float_rand( float min, float max )
{
    float scale = rand() / (float) RAND_MAX; /* [0, 1.0] */
    return min + scale * ( max - min );      /* [min, max] */
}

int int_rand(int lower, int upper) 
{ 
 
        return (rand() % (upper - lower + 1)) + lower; 
        
     
}

float interpolate(float a,float b,float x){
    float ft = x * 3.14159265359;
    float f = ( 1- cos(ft))*0.5;

    return a*( 1-f ) + b * f ;
}

class Noise{
    public:float * Array;
    
    public:Noise(int No,float min,float max){
        
        Array = (float*)calloc(No,sizeof(float));
        float prevRandomNo = float_rand(min,max);
        float randomNo = float_rand(min,max);
        
        Array[0]=prevRandomNo;
        int fraction =1;
        
        for(int i=1;i<No;i++){

            if(i%10 == 0){

                prevRandomNo=randomNo;
                randomNo=float_rand(min,max);
                Array[i]=prevRandomNo;
                int fraction = 1;
            }
            Array[i]=interpolate(prevRandomNo,randomNo, i/(No*1.0) );
            fraction+=1;
        }
    }
    
    public:~Noise(){
        delete []Array;
    }
};



int main() {
    int N = size;
    int scale = 4;
        const int H = N*scale;
        const int W = N*scale;

        sf::RenderTexture RT;
        RT.create(W, H);

        sf::RenderWindow window(sf::VideoMode(W, H, 60), "Test");

        //Creates the sf::Vertex matrix
        sf::Vertex **pixMat = new sf::Vertex*[H];
        
        
        
        for (int i=0;i<H;++i)
                pixMat[i] = new sf::Vertex[W];
        
        sf::Color color(0,0,0);

        
        

        //Loads the texture and creates the sf::Sprite we want to draw
        sf::Texture tex;
        
        //initialize Fluid
        
        Fluid cube(N,0.0001,0.00001,0.1);
        Noise noise(50,-5,5);
        
        while (window.isOpen()) {
                sf::Event event;
                while (window.pollEvent(event)) {
                        if (event.type == sf::Event::Closed)
                                window.close();
                }
                int index = 10;
                for(int i=-2;i<=2;i++)
                    for(int j=-2;j<=2;j++)
                        for(int k=-2;k<=2;k++)
                        {
                            cube.addDensity(N/2+i,N/2+j,N/2+k,int_rand(20,200));
                        
                            cube.addVelocity(N/2+i,N/2+j,N/2+k,float_rand(-5,5),float_rand(-5,5),0);
                        }
                
                cube.Step();
                for (int i=0;i<H;++i) {
                    for (int j=0;j<W;++j)
                    {
                        color.r=cube.density[IX(i/scale,j/scale,N/2)];
                        color.g=0;
                        color.b=0;
                        pixMat[i][j] = sf::Vertex(sf::Vector2f(j+.5f, i+.5f), color);
                    }
                }

                //Draws the pixel matrix
                for (int i=0;i<H;++i)
                        RT.draw(pixMat[i], W, sf::Points);


                //And finally draws the RenderTexture to the RenderWindow
                RT.display();
                window.draw(sf::Sprite(RT.getTexture()));
                window.display();
        }
}



