/*
 * Authors: Heredia Parillo, JuanPablo Andrew   &
 *          Rendon Zuniga, Luis Francisco
 * 
 * File:     main.c (sssp)
 * Purpose:  Implement "Scalable Single Source Shortest Path Algorithms 
 *           for Massively Parallel Systems"
 * Input:
 *      Number of nodes, delta parameter & source node for paths
 * Output:
 *      for now, the execution time
 *
 * Compile:  mpicc -Wall -o sssp main.c dyrray.c -I (makefile)
 * Run:
 *    mpiexec -n <p> sssp <nodes> <root> <delta> 
 *       - p: the number of processes
 *       - nodes: number of nodes
 *       - root: source for paths
 *       - delta: delta
 *
 */

#include <stdio.h>
#include "dyrray.h"
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define MAX_WEIGHT 10
#define MAX_DEGREE 10
#define INFI 1000

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void Get_args(int argc,char* argv[],int* global_cant_nodes,int* global_root_node,int* local_cant_nodes,int delta,int my_rank,int p,MPI_Comm comm);

int main(int argc, char** argv){
    int my_rank, comm_sz;
    double start, finish, loc_elapsed, elapsed;
    MPI_Comm comm;

    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    int global_cant_nodes;      // Cantidad total de los nodos del grafo
    int local_cant_nodes;       // Cantidad de nodos para cada proceso
    int global_root;            // Nodo Raiz
    int delta;                  // Parametro delta
    int *local_nodes_u;         // Arreglo que contendrá el identificador (indice) del nodo source (u)
    int *local_nodes_v;         // Arreglo que contendrá el identificador (indice) del nodo destino (v)
    int *local_weights;         // Arreglo que contendrá el peso o distancia entre el nodo u y v

    Get_args(argc, argv, &global_cant_nodes, &global_root, &local_cant_nodes, &delta, my_rank, comm_sz, comm);
    
    /**Arreglos que contendran todas las aristas*/
    int *temp_global_nodes_u, *temp_global_nodes_v, *temp_global_weights;
    int *last_global_nodes_u, *last_global_nodes_v, *last_global_weights;

    int truly_cant_elements;        // Cantidad total de aristas
    int local_cant_edges;           // Cantidad de aristas asignadas por proceso

    if(my_rank == 0){

        temp_global_nodes_u = new_array(last_global_nodes_u);
        temp_global_nodes_v = new_array(last_global_nodes_v);
        temp_global_weights = new_array(last_global_weights);
        int k,v,w;
        /**Creación del Grafo**/
        for(int i=0;i<global_cant_nodes;++i){
            k = 1 + (random() % MAX_DEGREE);
            for (int j = 0; j < k; ++j){
                v = random() % global_cant_nodes;
                while(v == i)
                    v = random() % global_cant_nodes;
                w = 1+(random()%MAX_WEIGHT);

                push_element(temp_global_nodes_u, i, last_global_nodes_u);
                push_element(temp_global_nodes_v, v, last_global_nodes_v);
                push_element(temp_global_weights, w, last_global_weights);
            }
        }
        /*****************************/

        truly_cant_elements = last_global_weights  - temp_global_nodes_u;
        local_cant_edges = truly_cant_elements / comm_sz;
        
        MPI_Bcast(local_cant_edges,1,MPI_INT,0,comm);
        MPI_Scatter(temp_global_nodes_u,local_cant_edges,MPI_Init,
                    local_nodes_u,local_cant_edges,MPI_INT,0,comm);
        MPI_Scatter(temp_global_nodes_v,local_cant_edges,MPI_Init,
                    local_nodes_v,local_cant_edges,MPI_INT,0,comm);
        MPI_Scatter(temp_global_weights,local_cant_edges,MPI_Init,
                    local_weights,local_cant_edges,MPI_INT,0,comm);

        free(temp_global_nodes_u);  free(temp_global_nodes_v);  free(temp_global_weights);

    }else{
        MPI_Scatter(temp_global_nodes_u,local_cant_edges,MPI_Init,
                    local_nodes_u,local_cant_edges,MPI_INT,0,comm);
        MPI_Scatter(temp_global_nodes_v,local_cant_edges,MPI_Init,
                    local_nodes_v,local_cant_edges,MPI_INT,0,comm);
        MPI_Scatter(temp_global_weights,local_cant_edges,MPI_Init,
                    local_weights,local_cant_edges,MPI_INT,0,comm);
    }
    MPI_Barrier(comm);
    /* Todos los procesos ya tienen la cantidad de nodos/aristas desigandas*/

    int *distances = malloc(global_cant_nodes * sizeof(int));       // Arreglo de las distancias a cada nodo 
    int *predecessors = malloc(global_cant_nodes * sizeof(int));    // Arreglo de los predecesores de cada nodo
    int **Bucks , **last_of_buck; // Arreglo de Buckets
    
    Bucks = malloc((INFI/delta + 1)*sizeof(int*));
    last_of_buck = malloc((INFI/delta + 1)*sizeof(int*));
    
    for(int i=0;i<(INFI/delta + 1);++i) 
        Bucks[i] = new_array(last_of_buck[i]);      // Todos los buckets son inicializados como nulos

    distances[global_root] = 0;
    predecessors[global_root] = global_root;

    push_element(Bucks[0],global_root,last_of_buck[0]);     // El bucket 0 contendrá al nodo raiz
    
    for (int i = 0; i < global_cant_nodes; ++i){
        if(i == global_root) continue;          // El resto de nodos que no son raiz, estan en el bucket infinito
        push_element(Bucks[INFI/delta],i,last_of_buck[INFI/delta]); 
    }
    
    MPI_Barrier(comm);
    start = MPI_Wtime();
    for(int k=0;k!=(INFI/delta);){            //////////////Epocas
        /*process bucket k*/
        int * Aux = Bucks[k], *last_aux = last_of_buck[k];      // Bucket Auxiliar
        int *last_Relax, *Relaxeds = new_array(last_Relax);     // Bucket auxiliar de relajaciones
        
        while(!array_is_empty(Aux,last_aux)){           ////////////////// Fases
            
            for(int i=0;i<(last_aux-Aux);++i){ // i en los nodos del Bucket Auxiliar
                
                for(int j=0;j<local_cant_edges;j++){ 
                                                            // Se revisan todas las aristas del nodo del Bucket Auxiliar
                    if(local_nodes_u[j] == *(Aux+i)){
                        /*relax*/
                        int old_buk = floor((double)distances[local_nodes_v[j]]/(double)delta);
                        int dis = distances[local_nodes_v[j]];
                        
                        distances[local_nodes_v[j]] = MIN(distances[local_nodes_v[j]],distances[local_nodes_u[j]]+local_weights[j]);
                        
                        int new_buk = floor((double)distances[local_nodes_v[j]]/(double)delta);
                        
                        if(new_buk < old_buk){
                            push_element(Bucks[new_buk], local_nodes_v[j],last_of_buck[new_buk]);
                            delete_element(Bucks[old_buk], local_nodes_v[j],last_of_buck[old_buk]);
                            /*
                            * Se comunica el estado de los buckets
                            */
                        }
                        /*end relax*/

                        // Si la distancia al nodod v fue actualizada se añade el nodo al bucket de relajaciones
                        if(dis != distances[local_nodes_v[j]]) 
                            push_element(Relaxeds,local_nodes_v[j],last_Relax);
                    }
                }
            }
            // El bucket auxiliar será la intersección entre el Bucket_k original y el bucket de relajaciones
            Intersec(Aux, last_aux, Bucks[k],last_of_buck[k],Relaxeds,last_Relax);
        }
        /*end process*/
       
        k++;
        // Bucket_k sera el siguiente bucket no nulo
        while(! array_is_empty(Bucks[k],last_of_buck[k]) && k < (INFI/delta))
            k++;
    }
    finish = MPI_Wtime();
    loc_elapsed = finish-start;
    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    free(local_nodes_u);    free(local_nodes_v);    free(local_weights);
    
    for(int i=0; i<(INFI/delta);++i)
        free(Bucks[i]);
    free(Bucks);
    
    if (my_rank == 0)
      printf("Tiempo transcurrido de ejecucion = %e\n", elapsed);

    MPI_Finalize();
    return 0;
}


void Get_args(
        int argc,               /*argc*/
        char* argv[],           /*argv*/
        int* global_cant_nodes, /*cantidad total de nodos*/
        int* global_root_node,  /*nodo raiz*/
        int* local_cant_nodes,  /*cantidad de nodos para cada proceso*/
        int* delta,             /*delta*/
        int my_rank,            /*rank identificador*/
        int p,                  /*cantidad de procesos (comm_sz)*/
        MPI_Comm comm           /*comm*/ ){

    if(my_rank ==0){
        if(argc != 4){
//            fflush(stderr);
            *global_cant_nodes = -1;
        }else{
            *global_cant_nodes = atoi(argv[1]);
            *global_root_node = atoi(argv[2]);
            *delta = atoi(argv[3]);
        }
    }
    MPI_Bcast(global_cant_nodes,1,MPI_INT,0,comm);
    MPI_Bcast(global_root_node,1,MPI_INT,0,comm);
    MPI_Bcast(delta,1,MPI_INT,0,comm);
    if(*global_cant_nodes <= 0){
        MPI_Finalize();
        exit(-1);
    }
    *local_cant_nodes = *global_cant_nodes / p;
}  /* Get_args */
