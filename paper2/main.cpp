/**
 * Authors:
 *          Heredia Parillo, JuanPablo Andrew
 *          Rendon Zuniga, Luis Francisco
 * 
 * File:    
 *          main.cpp (sssp)
 * Purpose:
 *          Find the shortest path from a source node to the rest.
 *          Implementation of "Scalable Single Source Shortest Path Algorithms 
 *          for Massively Parallel Systems".
 * Compile:
 *          g++ -Wall -fopenmp -I. -o sssp main.cpp -lm (makefile)
 * Run:
 *          ./sssp <cant_threads> <cant_nodes> <source_node> <delta>
 * Notes:
 *          -
 * 
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>
#include <utility>
#include <math.h>

#define MAX_WEIGHT 10
#define MAX_DEGREE 10
#define MAX_DISTANCE 1000
#define INFINITO 10000

typedef unsigned int uint;

int cant_threads;
int delta;

std::vector< std::vector< std::pair<int,int> > > graph;         // Matriz dispersa, cada nodo contiene sus aristas (destiny,weight)
std::vector< int > distances;
std::vector< int > node_predecessor;
std::vector< std::vector< int > > buckets;                      // Contenedor de nodos (indices al grafo)

void thread_work (int first_index, int last_index);
bool relax (int u, int v);

int main (int argc, char** argv){
    if (argc != 5){
        printf ("Incorrect parameters \n");
        return 0;
    }
    cant_threads        = atoi (argv[1]);
    uint cant_nodes     = atoi (argv[2]);
    uint cant_edges     = 0;
    uint node_source    = atoi (argv[3]);
    delta               = atoi (argv[4]);

    srand (0);

 /////////////////////////////////////////////////////////////////////////////
 // Inicializar
 /////////////////////////////////////////////////////////////////////////////

    std::vector< std::pair<int,int> >  temp;
    for (int i = 0; i < cant_nodes; ++i){
        graph.push_back (temp);
 // ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^  Con Poco hicimos una manera para rezize vector, si recuerdas lo cambias
        distances.push_back(INFINITO);
        node_predecessor.push_back(-1);
    }

    for (int i = 0; i < cant_nodes; ++i){
        int weight;
        int node_aux;
        int d = 1 + (rand() % MAX_DEGREE);
        
        for (int j = 0; j < d; ++j){
            weight = 1 + (rand () % MAX_WEIGHT);
            node_aux = rand () % cant_nodes;
        
            while(node_aux == i) node_aux = rand () % cant_nodes;
        
            graph[i].push_back (std::make_pair (node_aux, weight) );
            graph[node_aux].push_back (std::make_pair (i, weight) );
            cant_edges++;
        }
    }
 // ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ Puede mejorar, o variar para generar RMAT o para recivirlo externamente

    uint cant_buckets = (MAX_DISTANCE / delta) +1;

    std::vector< int >  temp2;
    for (int i = 0; i < cant_buckets; ++i)
        buckets.push_back(temp2);

    buckets[0].push_back(node_source);
    distances[node_source] = 0;

    for (int i = 0; i < cant_nodes; ++i)
        buckets[cant_buckets].push_back (i);
    
 /////////////////////////////////////////////////////////////////////////////
 // Separación del trabajo
 /////////////////////////////////////////////////////////////////////////////

    int thread_cant_nodes = ceil ((double)cant_nodes / (double)cant_threads);    

    #pragma omp parallel//num_threads(cant_threads)
    for(int i = 0; i < cant_nodes; i += thread_cant_nodes){
        #pragma omp task
        thread_work (i,i + (thread_cant_nodes-1) );
    }
    //#pragma omp taskwait
    #pragma omp barrier

 /////////////////////////////////////////////////////////////////////////////
 // Finalize
 /////////////////////////////////////////////////////////////////////////////

    for(size_t i = 0; i < cant_nodes; ++i)
        printf ("%d : (%d, %d)", i, node_predecessor[i], distances[i]);

    return 0;
}

/////No esta usando los indices, se me olvido, revisalo
void thread_work (
        int first_index         /*indice de inicio del bloque de nodos*/, 
        int last_index          /*indice de final del bloque de nodos*/){

    #pragma omp parallel num_threads(cant_threads)

    for (int k = 0; k != (MAX_DISTANCE/delta); ){         /* Epoca*/
        /*process*/
        while ( !buckets[k].empty() ){                 /*Fase*/
            std::vector<int> aux_bucket;
            
            #pragma omp for
            for (int u = 0; u < buckets[k].size (); ++u){ //// [ solo nodos entre first y last]
                for (int v = 0; v < graph[u].size (); ++v){
                    /*relax*/
                    if (relax (buckets[k][u], graph[u][v].first)){
                        aux_bucket.push_back (graph[u][v].first);
                    }
                }
            }
            //Intersección
            std::vector<int> temp_bucket(buckets[k]);
            buckets[k].clear ();
            for (int i = 0; i < temp_bucket.size (); ++i){
                for (int j = 0; j < aux_bucket.size (); ++j){
                    if (temp_bucket[i] == aux_bucket[j])
                        buckets[k].push_back (temp_bucket[i]);
                }
            }
        }
    }
}

/// Revisa si debe tener una parte critica
bool relax (int u, int v){
    uint old_bucket = floor ((double)distances[v] / (double)delta);
    int temp = distances[v];
    distances[v] = std::min (distances[v], (distances[u] + graph[u][v].second));
    uint new_bucket = floor ((double)distances[v] / (double)delta);

    if (new_bucket < old_bucket){
        buckets[new_bucket].push_back (v);
        for (int i = 0; i < buckets[old_bucket].size (); ++i){
            if (buckets[old_bucket][i] == v){
                buckets[old_bucket].erase (buckets[old_bucket].begin () + i);
                break;
            }
        }
    }
    return (temp != distances[v]);
}