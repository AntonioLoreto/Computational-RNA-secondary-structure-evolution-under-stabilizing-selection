#ifndef rna_hpp
#define rna_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <utility>
#include <numeric>
#include "/Users/neo/Desktop/experimentos/Simul_RNA_chidas/reactor/evo_red_neutral/v7/redneu_v7/babel_v6/random.hpp"

extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/dist_vars.h>
#include <ViennaRNA/RNAstruct.h>
#include <ViennaRNA/stringdist.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/string_utils.h>
#include <ViennaRNA/subopt.h>
#include <ViennaRNA/treedist.h>
#include <ViennaRNA/plot_structure.h>
}

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

/*Esta seccion es para la callback function de la funcion rna_subopt_cb*/

struct losfen{
  string estrus;
  float eners;
  float prob;
};

bool probDescendente(const losfen& x, const losfen& y) { return x.prob > y.prob; }
std::vector<losfen> otipos;

int conta = 0; 

void subopt_callback(const char* structure, float energy, void* data){
    if(structure){
        otipos.push_back(losfen());
        otipos[otipos.size()-1].estrus = structure;
        otipos[otipos.size()-1].eners = energy; 
        //conta++;
    }
}

struct muestra{
    int         freq;
    int         dist;
    float       energy;
    float       probabily;
    std::string estrus;
};

struct encontrar_estrus{
    std::string estrus;
    encontrar_estrus(std::string estrus) : estrus(estrus) {}
    bool operator() ( const muestra& m ) const{
        return m.estrus == estrus;
    }
};

class rna{
public:
    aleatorio est;
    rna();
    rna (aleatorio& jaso);
    int find_ss(std::vector<muestra > m, std::string objetivo);
    int mutar_seq_r(std::string seq, std::string &b); 
    int contador_tam(std::string a, int inicio, int finall);
    void start_r (aleatorio& jaso);
    void foto_seqs(std::vector<string> seqs, std::vector< std::pair<vector<string>, vector<int> > >& v);
    void quitar_spacios(std::string &a);
    void neutral_walk_muest(std::string &a, int l, const char* structure, int numero_pasos, int &neutrales); //solo aceptas mutaciones nuetrales pero cuentas todas las mutaciones
    void generador_mutaciones (string a, int tam, string *aa, int *e);
    float probab_reperplastic_mfe(std::string seq, std::string obje, int del);
    double distprom(std::vector<string> seqs, std::vector<int> &dis );
    double funcion_seleccion(double dis, int l);
    double funcion_seleccion_simple (int dis_ham, int tam);
    double adecuacion(std::string secuencia, std::string ss_objetivo, float delta);
    double robustez_mutacional(std::string &secuencia, std::string objetivo);
    double funcion_seleccion_simple_tree(double dis_tree);
    double funcion_seleccion_tres (int dis_ham, int tam);
    std::string una_mutacion (std::string &a);
    std::string generador_mutaciones1 (string a, int tam);
    std::string convertString (char* a, int size);
    std::string caminata_neutral(std::string a, int lengt, const char* structure, int numero_pasos); //solo acepta mutaciones neutrales y solo cuenta mutaciones neutrales
    std::string neutral_walk_repertorio(std::string b, std::string objetivo_B, std::string objetivo_C,  int length, int numero_pasos, int BoC ); /*B = 0, C = 1*/
    //string poblaciones(string txt, string obj_A, string obj_B);
    //int find_tam(string a);
    std::string random_rna_seq (int length);
    std::string inverse_fold (std::string strus_objetivo);
    std::string caminata_restriccion_repertorio(std::string b, std::string objetivo_B, std::string objetivo_C, std::vector<std::vector<int> > &h );
    std::string caminata_restriccion_repertorio_alea(std::string b, std::string objetivo_B, std::string objetivo_C, bool siono);//siono = 1 entonces B en repertorio, lo contrario significa lo contrario 
    std::string mutar_seq (std::string &a);
private:
};

rna::rna()
{
}

rna::rna(aleatorio& jaso)//creates instance of class rna and assigns r to jaso
{
    start_r(jaso);
}

void rna::start_r(aleatorio& jaso)//assigns r to jaso
{
    est = jaso;
}

int rna::find_ss(std::vector<muestra > m, std::string objetivo){
    int index = m.size();
    for(int i = 0; i < m.size(); ++i){
        if(m[i].estrus == objetivo){
            index = i;
        }
    }
    return index;
}

float rna::probab_reperplastic_mfe(std::string seq, std::string obje, int del){
    int                         num_subopt;  
    int                         inservible  = 0;
    int                         index;
    int                         delta       = del;
    float                       p;
    double                      particion; 
    double                      kt;
    vrna_fold_compound_t        *foco; 
    std::vector<muestra>        m;

    kt      = (25 + 273.15) * 1.98717 / 1000;
    foco    = vrna_fold_compound(seq.c_str(), NULL, VRNA_OPTION_DEFAULT);
    vrna_subopt_cb(foco, delta, &subopt_callback, (void *)&inservible);
    vrna_fold_compound_free(foco);
    num_subopt = otipos.size();
    particion  = 0;
    for(int i = 0; i < num_subopt; ++i){
        particion += exp((-1.0*otipos[i].eners)/(kt)); //El repertorio plastico que se saca con subopt_cb tambien tiene a la estrus y eners de la mfe ss
    }   

    for(int i = 0; i < num_subopt; ++i){
        m.push_back(muestra());
        m[i].estrus     = otipos[i].estrus;
        m[i].energy     = otipos[i].eners;
        m[i].freq       = 1;
        m[i].probabily  = exp((-1.0*otipos[i].eners)/(kt))/particion;
    }
    

    index = find_ss(m, obje);
    if(index == m.size()){
        p = 0; 
    }else{
        p = m[index].probabily;
    }
    
    otipos.clear();
    m.clear();

    return p;
}

void rna::foto_seqs(std::vector<string> seqs, std::vector< std::pair<vector<string>, vector<int> > >& v){

    std::vector<int> num;
    std::pair<vector<string>, vector<int> > p;
    for(int i = 0; i < seqs.size(); ++i){
            num.push_back(1);
            for(int j = i + 1; j < seqs.size();){
                if(seqs[i] == seqs[j]){
                    num[i] += 1;
                    seqs.erase(seqs.begin() + j);
                }
                else ++j; 
            }
        }
        p.first     = seqs;
        p.second    = num;
        v.push_back(p);
}

double rna::distprom(std::vector<string> seqs, std::vector<int> &dis){
    int d           = 0;
    double s        = 0;
    double sumdist  = 0;
    
    for(int i = 0; i < (seqs.size() - 1); ++i){
        for(int j = i + 1; j < seqs.size(); ++j){
            d = vrna_hamming_distance(seqs[i].c_str(), seqs[j].c_str());
            dis.push_back(d);
        }
    }
    for(int i = 1; i < seqs.size(); ++i){
        s += (seqs.size() - i);
    }
    sumdist = std::accumulate(dis.begin(), dis.end(), 0);
    return (sumdist/s) ;
}

double rna::funcion_seleccion(double dis, int l){
    return 1/(0.01+(dis/l));
}

double rna::funcion_seleccion_simple_tree(double dis_tree){
    return 1/(0.01 + dis_tree);//0.99\aproxf<=100;
}

double rna::adecuacion(std::string secuencia, std::string ss_objetivo, float delta){
    int  l; 
    Tree *T1;
    Tree *T2;
    char *st1;
    char *st2;
    char *strus;
    double particion;
    double adecuacion = 0;
    double dis_tree;
    double kt = (25 + 273.15) * 1.98717 / 1000.0;  //kT in kcal/mol

    l       = secuencia.length();
    strus   = strcpy(new char[ss_objetivo.length() + 1], ss_objetivo.c_str()); //string a char*
    st1     = expand_Full(strus);
    T1      = make_tree(st1);
    free(st1);

    vrna_fold_compound_t *foco = vrna_fold_compound(secuencia.c_str(), NULL, VRNA_OPTION_DEFAULT);
    vrna_subopt_cb(foco, delta, &subopt_callback, (void *)&otipos);
    vrna_fold_compound_free(foco);

    particion = 0;
    for(int i = 0; i < otipos.size(); ++i){
        particion += exp((-1.0*otipos[i].eners)/(kt));
    }

    for(int i = 0; i < otipos.size(); ++i){
        st2      = expand_Full(otipos[i].estrus.c_str());
        T2       = make_tree(st2);
        dis_tree = tree_edit_distance(T1, T2); 
        dis_tree = dis_tree/(2*l);
        adecuacion += funcion_seleccion_simple_tree(dis_tree)*(exp((-1.0*otipos[i].eners)/kt)/particion);
        free_tree(T2);
        free(st2);
    }
    otipos.clear();
    return adecuacion;
}

int rna::contador_tam(std::string a, int inicio, int finall){
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


int rna::mutar_seq_r(std::string seq, std::string &b){
    double d, frecuencia_mutacion = 0.001;
    int y = 0, r = 0;
    const char nucleotidos[4] = {'A', 'C', 'G', 'U'};
    std::string copia = seq;

    for(int j = 0; j < copia.size(); ++j){
        //Generar numero aleatorio entre 0 y 0.999
        d = est.random_tresdeci();
       if(d < frecuencia_mutacion){
            r = 1;
           //Genera numeros entre 0 y 3
           y =  est.tres();
 
           while(copia[j]==nucleotidos[y]){
               y = est.tres();
           }
       copia[j] = nucleotidos[y];
       }
   }
   b = copia;
   return r;
}

std::string rna::una_mutacion(std::string& a){
    int y = 0, pos;
    const char nucleotidos[4] = {'A', 'C', 'G', 'U'};
    std::string b = a; 
    pos = est.length(a.size());
    y =  est.tres();
    while(a[pos]==nucleotidos[y]){
       y = est.tres();
    }
    b[pos] = nucleotidos[y];
   return b;
}

std::string rna::mutar_seq(std::string& a){
    int      y = 0;
    double   d;
    double   frecuencia_mutacion = 0.001;
    const char nucleotidos[4] = {'A', 'C', 'G', 'U'};
    for(int j = 0; j < a.size(); ++j){
        d = est.random_tresdeci(); //Generar numero aleatorio entre 0 y 0.999
       if(d < frecuencia_mutacion){
           y =  est.tres(); //Genera numeros entre 0 y 3
           while(a[j] == nucleotidos[y]){
               y = est.tres();
           }
       a[j] = nucleotidos[y];
       }
   }
   return a;
}

std::string rna::generador_mutaciones1(string a, int tam){
    double d, frecuencia_mutacion = 0.001;
    int y = 0;
    const char nucleotidos[4] = {'A', 'C', 'G', 'U'};

    for(int j=0; j<tam;++j){
        //Generar numero aleatorio entre 0 y 0.999
        d = est.random_tresdeci();
       if(d < frecuencia_mutacion){
           //Genera numeros entre 0 y 3
           y =  est.tres();
 
           while(a[j]==nucleotidos[y]){
               y =  est.tres();
           }
       a[j] = nucleotidos[y];
       }
   }
return a;
}

void rna::generador_mutaciones(std::string a, int tam, std::string *aa, int *e ){
    int y       = 0;
    int mien    = 0;
    double d;
    double frecuencia_mutacion = 0.001;
    const char nucleotidos[4] = {'A', 'C', 'G', 'U'};

    for(int j = 0; j < tam; ++j){
        //Generar numero aleatorio entre 0 y 0.999
        d = est.random_tresdeci();
        
       if(d < frecuencia_mutacion){
           mien = 1;
           //Genera numeros entre 0 y 3
           y =  est.tres();
 
           while(a[j]==nucleotidos[y]){
               y =  est.tres();
           }
           a[j] = nucleotidos[y];
       }
        *e = mien;
   }
        *aa = a;
}

double rna::funcion_seleccion_simple(int dis_ham, int tam){
    double f = 0;
    f = 1/(0.01+(dis_ham/double(tam)));
    return f;
}

double rna::funcion_seleccion_tres(int dis_ham, int tam){
    double f = 0;
    f = 1/(0.03+(dis_ham/double(tam)));
    return f;
}
string rna::convertString(char* a, int size){
    string s = ""; //Tal vez robe memoria
    for (int i = 0; i < size; i++) {
        s = s + a[i];
    }
    return s;
}

std::string rna::caminata_neutral(std::string a, int lengt, const char* structure, int numero_pasos){ 
    int dis, y, d;
    double mfe1;
    char nucleotidos[4] = {'A', 'C', 'G', 'U'}, j;
    char* structure0;
    std::string p = a;
    structure0 = new char[lengt+1];
    for(int i=0; i<numero_pasos;++i){
        dis = 1;
        while(dis!=0){
            d = est.length(lengt); //[0,n-1]
            y = est.tres();
            while(p[d]==nucleotidos[y]){
                y = est.tres();
            }
            j = p[d];
            p[d] = nucleotidos[y];
            mfe1 = vrna_fold(p.c_str(), structure0);
            dis = vrna_hamming_distance(structure, structure0);
            if(dis!=0){
                p[d] = j;
            }
        }
    }
    delete[] structure0;
return p;
}

void rna::neutral_walk_muest(std::string &a, int l, const char* structure, int numero_pasos, int &neutrales){
    int dis, y, d, no_neutral = 0; 
    double mfe1;
    char nucleotidos[4] = {'A', 'C', 'G', 'U'}, j;
    char* structure0;
    structure0 = new char[l+1];
    for(int i=0; i < numero_pasos;++i){
        d = est.length(l);
        y = est.tres();
        while(a[d]==nucleotidos[y]){ //esto sesga la caminata
            y = est.tres();
        }
        j = a[d];
        a[d] = nucleotidos[y];
        mfe1 = vrna_fold(a.c_str(), structure0);
        dis = vrna_hamming_distance(structure, structure0);
        /*if(dis==0){
            ++c;
            cout<<c<<endl;
        }*/
        if(dis!=0){
            a[d] = j;
            ++no_neutral;
        }   
    }
    neutrales = numero_pasos - no_neutral;
    delete[] structure0;
}

std::string rna::neutral_walk_repertorio (std::string b, std::string objetivo_B, std::string objetivo_C,  int length, int numero_pasos, int BoC ){ //BoC = 0 buscas B en repertorio, BoC = 1 buscas seq sin B ni C en repertorio
    int dis, y, d, ya_B, ya_C, dis_B, dis_C, k = 0;
    double delta = 300;
    double mfe1;
    char nucleotidos[4] = {'A', 'C', 'G', 'U'}, j;
    char* structure0;
    structure0 = new char[length+1];
    std::string a = b;
    mfe1 = vrna_fold(a.c_str(), structure0);
    std::string nativa(structure0);
    int inservible = 0; 
    vrna_fold_compound_t *foco; 
    for(int i = 0; i<numero_pasos;++i){
        dis = 1;
        while(dis!=0){
            d = est.length(length); //[0,n-1]
            y = est.tres();
            while(a[d]==nucleotidos[y]){
                y = est.tres();
            }
            j = a[d];
            a[d] = nucleotidos[y];
            mfe1 = vrna_fold(a.c_str(), structure0);
            dis = vrna_hamming_distance(nativa.c_str(), structure0);
            if(dis!=0){
                a[d] = j;
            }
        }
        //conta = 0;

        foco = vrna_fold_compound(a.c_str(), NULL, VRNA_OPTION_DEFAULT);
        vrna_subopt_cb(foco, delta, &subopt_callback, (void *)&inservible);
        vrna_fold_compound_free(foco);
        
        if(BoC == 0){
            ya_B = 10;
            ya_C = 10;
            for(int i = 0; i < otipos.size(); ++i){
                dis_B = vrna_hamming_distance(otipos[i].estrus.c_str(), objetivo_B.c_str());
                dis_C = vrna_hamming_distance(otipos[i].estrus.c_str(), objetivo_C.c_str());
                if(dis_B == 0){
                    ya_B = 1;
                }
                if(dis_C ==0){
                    ya_C = 1;
                }
            }
            if(ya_B != 1 && ya_C != 0){
                a[d] = nucleotidos[y];
                --i;
                k++;
            }
        }
        else{
            for(int i = 0; i < otipos.size(); ++i){
                dis_B = vrna_hamming_distance(otipos[i].estrus.c_str(), objetivo_B.c_str());
                dis_C = vrna_hamming_distance(otipos[i].estrus.c_str(), objetivo_C.c_str());
                
                if(dis_B == 0){
                    ya_B = 1;
                }
                if(dis_C ==0){
                    ya_C = 1;
                }
            }
            if(ya_B != 0 && ya_C != 0){
                a[d] = nucleotidos[y];
                --i;
                k++;
            }
        }
       otipos.clear();
       //conta = 0;
   }
    delete[] structure0;
    return a;
}
std::string rna::random_rna_seq (int length){
    const char symbols[4] = {'A','U','G','C'};
    string seq;
    for(int i=0;i<length;++i){
        int a = est.tres();
        seq.push_back(symbols[a]);
    }
    return seq;
}

std::string rna::inverse_fold (std::string strus_objetivo){
    int y = 0, u = 0, distancia = 1, distancia2 = 1, contador = 0, fracaso = 0, a;
    char b;
    const char nucleotidos[4] = {'A', 'C', 'G', 'U'};
    int tam = strus_objetivo.size();
    std::string seq;
    for(int i=0;i<tam;++i){
        a = est.tres();
        seq.push_back(nucleotidos[a]);
    }
    char *structure1 = new char[tam+1];
    double mfe;
    mfe = vrna_fold(seq.c_str(), structure1);
    int max = pow(tam,4);

    while (distancia2 != 0){
        mfe = vrna_fold(seq.c_str(), structure1);
        distancia2 = vrna_hamming_distance(strus_objetivo.c_str(), structure1);
        
        u = est.length(tam);
        y = est.tres();
        
        while(seq[u]==nucleotidos[y]){
            y =  est.tres();
        }
        
        b = seq[u];
        seq[u] = nucleotidos[y];
        
        mfe = vrna_fold(seq.c_str(), structure1);
        
        distancia = vrna_hamming_distance(strus_objetivo.c_str(), structure1);
        
        if(distancia2 < distancia){
            seq[u] = b;
        }
        
        //lo siguiente no es del todo cierto, tam^4 intentos no indican necesariamente que no existe una secuencia, pero me vale madres y aqui se para todo alv
        
        if(contador == max){
            if(distancia < distancia2){
                fracaso = distancia;
            }
            else{
                fracaso = distancia2;
            }
            cout<<"No hay una secuencia con esa estructura mi pana, la mas cerca es esta "<<seq<<" y esta a "<<fracaso<<" posiciones"<<endl;
            break;
        }
        contador++;
    }
    delete []structure1;

    return seq;
}

std::string rna::caminata_restriccion_repertorio(std::string b, std::string objetivo_B, std::string objetivo_C, std::vector< std::vector<int> > &h ){
    char nucleotidos[4] = {'A', 'C', 'G', 'U'}, j;
    int tam = b.size(), dis, dis_B, dis_C, ya_B, ya_C, terminaste = 0, delta = 300, size;
    double mfe;
    char* structure0 = new char[tam+1];
    char* structure1 = new char[tam+1];
    std::string a;
    mfe = vrna_fold(b.c_str(), structure1);
    bool tarea;
    vrna_fold_compound_t *foco; 
    int inservible = 0; 
    for(int i = 0; i < tam; ++i){
        a = b;
        for(int k = 0; k < 4; ++k){
            tarea = 0;
            a = b;
            dis = 0;
            if(h.size() != 0 ){
                for(int m = 0; m < h.size(); ++m){
                    if( i == h[m][0] && k == h[m][1]){
                        tarea = 1;
                        break;
                    }
                }
            }
            if(tarea){
                break; //no se si es mejor un break o un continue
            }
            if (a[i] == nucleotidos[k]){
                continue;
            }
            j = a[i];
            a[i] = nucleotidos[k];
            mfe = vrna_fold(a.c_str(), structure0);
            dis = vrna_hamming_distance(structure0, structure1);
            if(dis!=0){
                a[i] = j;
                continue;
            }
            //conta = 0;
            foco = vrna_fold_compound(a.c_str(), NULL, VRNA_OPTION_DEFAULT);
            vrna_subopt_cb(foco, delta, &subopt_callback, (void *)&inservible);
            vrna_fold_compound_free(foco);
            ya_B = 0;
            ya_C = 0;
            
            for(int l = 0; l < otipos.size(); ++l){
                dis_B = vrna_hamming_distance(otipos[l].estrus.c_str(), objetivo_B.c_str());
                dis_C = vrna_hamming_distance(otipos[l].estrus.c_str(), objetivo_C.c_str());
            
                if(dis_B == 0){
                    ya_B = 1;
                }
                if(dis_C ==0){
                    ya_C = 1;
                }
            }
            if(ya_B == 1 && ya_C == 0){
                terminaste = 1;
                h.push_back(std::vector<int> ());
                size = h.size();
                h[(size - 1)].push_back(i);
                h[(size - 1)].push_back(k);
                break;
            }
            otipos.clear();
        }
        if(terminaste == 1){
            break;
        }
    }
    delete []structure0;
    delete []structure1;
    
    return a;
}

std::string rna::caminata_restriccion_repertorio_alea(std::string b, std::string objetivo_B, std::string objetivo_C, bool siono){
    char nucleotidos[4] = {'A', 'C', 'G', 'U'}, j;
    int tam = b.size(), dis, dis_B, dis_C, ya_B, ya_C, terminaste = 0, delta = 300, d, y, cuantos=0;
    double mfe;
    char* structure0 = new char[tam+1];
    char* structure1 = new char[tam+1];
    std::string a;
    a = b;
    mfe = vrna_fold(b.c_str(), structure1);
    int inservible = 0; 
    vrna_fold_compound_t *foco;
    while(terminaste != 1){
        cout<<"Intentos fallidos "<<cuantos<<endl;
        //terminaste = 0;
        dis = 1;
        while(dis != 0){
            d = est.length(tam); //[0,n-1]
            y = est.tres();
            while(a[d] == nucleotidos[y]){
                y = est.tres();
            }
            j = a[d];
            a[d] = nucleotidos[y];
            mfe = vrna_fold(a.c_str(), structure0);
            dis = vrna_hamming_distance(structure0, structure1);
            if(dis != 0){
                a[d] = j;
            }
        }
        
        //conta = 0;
        foco  = vrna_fold_compound(a.c_str(), NULL, VRNA_OPTION_DEFAULT);
        vrna_subopt_cb(foco, delta, &subopt_callback, (void *)&inservible);
        vrna_fold_compound_free(foco);
        ya_B = 0;
        ya_C = 0;
        for(int l = 0; l < otipos.size(); ++l){
            dis_B = vrna_hamming_distance(otipos[l].estrus.c_str(), objetivo_B.c_str());
            dis_C = vrna_hamming_distance(otipos[l].estrus.c_str(), objetivo_C.c_str());
            
            if(dis_B == 0){
                ya_B = 1;
                break;
            }
            if(dis_C ==0){
                ya_C = 1;
            }
        }
        if(siono == 1){
            if(ya_B == 1 /*&& ya_C ==*/ ){
                terminaste = 1;
            }else{
                a[d] = j;
                ++cuantos;
            }
        }else{
            if(ya_B == 0 /*&& ya_C == 0*/){
                terminaste = 1;
            }else{
                a[d] = j;
                ++cuantos;
            }
        }
        otipos.clear();
    }
    delete []structure0;
    delete []structure1;
    return a;
}

void rna::quitar_spacios(std::string &a){
    a.erase(remove(a.begin(), a.end(), ' '), a.end());
}
double rna::robustez_mutacional(std::string &secuencia, std::string objetivo){
    char    nucleotidos[4] = {'A', 'C', 'G', 'U'};
    char    bolsita;
    char    *estrus; 
    int     tam;
    int     dis;
    int     i;
    int     j; 
    float   mfe;
    double  t;
    double  robuz; 
    double  neutrales = 0; 
    
    tam         = secuencia.size();
    t           = tam*3;
    estrus      = new char[tam + 1];

    for(i = 0; i < tam; ++i){
        bolsita = secuencia[i];
        for(j = 0; j < 4; ++j){
            secuencia[i] = bolsita;
            if(secuencia[i] == nucleotidos[j]){
                continue;
            }
            secuencia[i]    = nucleotidos[j];
            mfe             = vrna_fold(secuencia.c_str(), estrus);
            dis             = vrna_hamming_distance(estrus, objetivo.c_str());
            if(dis == 0){
                neutrales += 1;
            }
        }
        secuencia[i] = bolsita;
    }
    robuz = neutrales/t;
    delete[] estrus; 

    return robuz;
}

#endif
