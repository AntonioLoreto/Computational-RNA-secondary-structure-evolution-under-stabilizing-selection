#ifndef random_hpp
#define random_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

gsl_rng* r;

class aleatorio{
    public:
    void sid(int seed);
    bool moneda();
    double uno();
    int tres();
    int length(int a);
    double random_tresdeci();
    
};

void aleatorio::sid(int seed){
    gsl_rng_env_setup();
     r = gsl_rng_alloc (gsl_rng_ranlxs2); //Esto crea una instancia (r) del tipo de generador ranlsx2
    gsl_rng_set(r, seed);
}


bool aleatorio::moneda(){
    double d;
    d = gsl_rng_uniform(r);
    bool moneda = true;
    if(d<0.49){
        moneda = false;
    }
    return moneda;
}

//checar
double aleatorio::uno(){
    return gsl_rng_uniform(r); //returns random int between 0 and 0.99
}

int aleatorio::tres(){
    return gsl_rng_uniform_int(r, 4);
}

int aleatorio::length(int a){
    return gsl_rng_uniform_int(r, a);//[0,n-1]
}

double aleatorio::random_tresdeci(){
    double d;

    d = gsl_rng_uniform_int(r,1000); //numero entre 0 - 999
    d = d/1000.0;
    return d;
}
#endif /* random_hpp */
