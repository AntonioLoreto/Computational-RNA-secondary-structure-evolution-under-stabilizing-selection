#ifndef basics_hpp
#define basics_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_eigen.h>


using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

class basic{
public:
        int contador_tam(std::string a, int inicio, int finall);
        void high(vector<double>& l, double* ma, int* lugar);
        void columna(vector<vector <double> >& a, vector<double>& col, int m);
        void print_matrix(vector<vector <float> >& a, int n, int m);
        void print_vector(vector <float>& v);
        void gauss(vector<vector <double> > a, vector<double> b , vector<double>& col, vector<double>& s, int n); //agregar & hace que aunque la funcion sea void se pueda modificar lo que va despues de &	
        void ceros(int a[], int n);
        void ceros(float a[], int n);
        float SD(std::vector<float> v);

private: 
};

float basic::SD(std::vector<float> v){
        float   prom            = 0;
        float suma              = 0;
        float cuantos = 0;

        for(int i = 0; i < v.size(); ++i){
                prom += v[i];
                ++cuantos;
        }
        prom = prom/cuantos;

        for(int i = 0; i < v.size(); ++i){
                suma += pow((v[i] - prom), 2);
        }
        suma = suma/cuantos;
        suma = sqrt(suma);

        return suma;
}

void basic::ceros(int a[], int n){
        for(int i = 0; i < n; ++i){
                a[i] = 0; 
        }
}
void basic::ceros(float a[], int n){
        for(int i = 0; i < n; ++i){
                a[i] = 0; 
        }
}

void basic::gauss(vector<vector <double> > a, vector<double> b , vector<double>& col, vector<double>& s, int n){
        double  max;
        double  mult;
        double  x;
        int     lugar;
        int     iwal;
        int     contador;
        for(int j = 0; j < n - 1; j++){
                col.clear();
                columna(a, col, j); 
                max = 0;
                lugar = 0;
                high(col, &max, &lugar);
                a[j].swap(a[lugar+j]); //cambiar hileras de a
                iter_swap(b.begin() + j, b.begin() + (lugar+j)); //cambiar elementos de b
                for(int i = j+1; i < n; i++){
                        mult = a[i][j]/a[j][j];
                        for(int k = j; k < n; k++){
                                a[i][k] = a[i][k] - mult*a[j][k];
                        }   
                        b[i] = b[i] - mult*b[j];
                }   
        }
        //cout<<a[n-1][n-1]<<endl;
        if(fabs(a[n-1][n-1]) < 0.000001){
                cout<<"SISTEMA INCONSISTENTE O CON INFINITAS SOLUCIONES"<<endl;
                return;
        }  

        x = 0;
        iwal = 0;
        //Back substitution
        for(int i = n; i --> 0;){ //hace un cuenta de mayor a menor (i = n es en realidad i = n-1)
                x = b[i];
                contador = 0;
                for(int j = n; j --> 0;){
                        if(i == (n-1)){
                                iwal = j;
                                break;
                        }
                        if(j == i){
                                iwal = j;
                                continue;
                        }
                        if(j < i){
                                break;
                        }
                        x -= a[i][j]*s[contador];
                        contador = contador + 1;
                }
                x = x*(1/a[iwal][iwal]);
                s.push_back(x);
        }

        reverse(s.begin(), s.end());
}

void basic::high(vector<double>& l, double* ma, int* lugar){
        
        int lug = 0;
        double m = l[0];
        for(int i = 0; i < l.size(); i++){
                if(fabs(l[i]) > fabs(m)){
                        m = l[i];
                        lug = i;
                }
        }
        *ma = m;
        *lugar = lug;
}

void basic::columna(vector<vector <double> >& a, vector<double>& col, int m){
        for(int i = m; i < a.size(); i++){
                col.push_back(a[i][m]);
        }
}

void basic::print_matrix(vector<vector <float> >& a, int n, int m){
        for(int i = 0; i < n; i++){
                for(int j = 0; j < m; j++){
                        cout<<a[i][j]<<"|";
                }
                cout<<endl;
        }
}

void basic::print_vector(vector <float>& v){
        cout<<"{";
        for(int i = 0; i < v.size(); i++){
                cout<<v[i]<<", ";
        }
        cout<<"}"<<endl;
}

int basic::contador_tam(std::string a, int inicio, int finall){
    int contador = 0;
    for(int j = inicio; j < finall ; ++j){
        if(a[j] != ' '){
            ++contador;
        }
        else{
            break;
        }
    }
    return contador;
}
#endif
