#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <unordered_set>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
//#define VERBOSE true
#define FRONTIER_TOP_BFS_THRESHOLD 700000
#define TOP_BFS true
#define DOWN_BFS false
#define MAX_FRONTIER_NODES_FOUND 1000

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            if (distances[outgoing] == NOT_VISITED_MARKER) {
                distances[outgoing] = distances[node] + 1;
                int index = new_frontier->count++;
                new_frontier->vertices[index] = outgoing;
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;
    int** frontier_subarrays = new int*[graph->num_nodes];
    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++){
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    int* frontier_subarray_counts = new int [graph->num_nodes];
    int* frontier_partial_sums = new int [graph->num_nodes];
    
    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
        
#pragma omp parallel for schedule(guided)
        for(int i = 0; i < frontier->count; i++){
            int node = frontier->vertices[i];
            int start_edge = graph->outgoing_starts[node];
        int end_edge = (node == graph->num_nodes - 1)
                           ? graph->num_edges
                           : graph->outgoing_starts[node + 1];
        int local_frontier_count = 0;
        int local_frontier[end_edge - start_edge];;
        
        //frontier_subarrays[i] =  local_frontier;
        
            for(int i = start_edge; i < end_edge; i++){
                int outgoing_edge = graph->outgoing_edges[i];
                

                if( sol->distances[outgoing_edge] == NOT_VISITED_MARKER){
                    
                    //if(__sync_bool_compare_and_swap(& sol->distances[outgoing_edge], NOT_VISITED_MARKER, sol->distances[node] + 1)){
                    sol->distances[outgoing_edge] = sol->distances[node] + 1;
                        local_frontier[local_frontier_count++] = outgoing_edge;
               
                /*
                while(!__sync_bool_compare_and_swap(&(new_frontier->count), index, index+1)){
                    index = new_frontier->count;
                }*/
                        
                        
                //new_frontier->vertices[index] = outgoing_edge;
                        //}
                        
                        //}
                }
            }
            int start_i = __sync_fetch_and_add(&new_frontier->count, local_frontier_count);
            int end_i = start_i + local_frontier_count;
            memcpy(new_frontier->vertices + start_i, local_frontier, local_frontier_count * sizeof(int));
            //free(local_frontier);
            
        }
        //frontier_subarray_counts[i] = local_frontier_count;
            
    
        
    /*int sum = 0;
#pragma parallel for reduction(+:sum)
        for(int i =0; i < frontier->count; i++){
            frontier_partial_sums[i] = sum;
#pragma omp scan exclusive(sum)
            sum += frontier_subarray_counts[i];
        }
        new_frontier->count = frontier_partial_sums[frontier->count -1] + frontier_subarray_counts[frontier->count-1];
        
#pragma omp parallel for schedule(guided)
        for(int i = 0; i < frontier->count; i++){
            int start_i = frontier_partial_sums[i];
            int end_i = i == frontier->count - 1? new_frontier->count : frontier_partial_sums[i+1];
            

            for(int j = start_i, ii = 0; j < end_i; j++, ii++){
                new_frontier->vertices[j] = frontier_subarrays[i][ii];
            }
            delete [] frontier_subarrays[i];
        }
    */        
        vertex_set*  temp = new_frontier;
            new_frontier = frontier;
            frontier = temp;
            
        vertex_set_clear(new_frontier);

        //top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

    }
}


void bfs_bottom_up(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.
    int frontier_count  = 1;
    bool* frontier = new bool [graph->num_nodes];
    bool* new_frontier = new bool[graph->num_nodes];
    
    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++){
        sol->distances[i] = NOT_VISITED_MARKER;
        frontier[i] = false;
        new_frontier[i] = false;
    }

    
    sol->distances[ROOT_NODE_ID] = 0;
    frontier[ROOT_NODE_ID] = true;
    int curr_distance = 0;
    

    while(frontier_count != 0){
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
#pragma omp parallel for schedule(guided)
        for(int i = 0; i < graph->num_nodes; i++){
            int node = i;
            int incoming_start = graph->incoming_starts[node];
        int incoming_end = (node == graph->num_nodes - 1)
                           ? graph->num_edges
                           : graph->incoming_starts[node + 1];
        bool ancestor_found = false;
        if(sol->distances[node] == NOT_VISITED_MARKER){
            for(int  ancestor_i = incoming_start; ancestor_i < incoming_end; ancestor_i++){
                int ancestor = graph->incoming_edges[ancestor_i];
                if(frontier[ancestor] == true){
                    sol->distances[node] = sol->distances[ancestor] + 1;
                    ancestor_found = true;
                    break; 
                }
            }
        }
        
         new_frontier[node] = ancestor_found;
            
        }
        bool* temp = frontier;
        frontier = new_frontier;
        new_frontier = temp;
        frontier_count = 0;
        
#pragma omp parallel for reduction(+:frontier_count)
        for(int i = 0; i < graph->num_nodes; i++){
            frontier_count += frontier[i];
        }
#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier_count, end_time - start_time);
#endif

    }
    delete [] frontier;
    delete [] new_frontier;
    
}

int bfs_top_down_step(Graph graph, solution* sol, vertex_set* frontier, vertex_set* new_frontier){
#pragma omp parallel for schedule(guided)
        for(int i = 0; i < frontier->count; i++){
            int node = frontier->vertices[i];
            int start_edge = graph->outgoing_starts[node];
        int end_edge = (node == graph->num_nodes - 1)
                           ? graph->num_edges
                           : graph->outgoing_starts[node + 1];


            for(int i = start_edge; i < end_edge; i++){
                int outgoing_edge = graph->outgoing_edges[i];
                

                if( sol->distances[outgoing_edge] == NOT_VISITED_MARKER){
                    if ( __sync_bool_compare_and_swap( & (sol->distances[outgoing_edge]), NOT_VISITED_MARKER, (sol->distances[node]) + 1)){
                        int index = 0;
                        
                        #pragma omp atomic capture
                        
                        {
                            index = new_frontier->count;
                            new_frontier->count++;
                        }
                        
                        new_frontier->vertices[index] = outgoing_edge;
                        
                    }
                }
            
     }
        }
        return new_frontier->count;
}

int bfs_bottom_up_step(Graph graph, solution* sol, bool* frontier, bool* new_frontier){
    
#pragma omp parallel for schedule(guided)
        for(int i = 0; i < graph->num_nodes; i++){
            int node = i;
            int incoming_start = graph->incoming_starts[node];
        int incoming_end = (node == graph->num_nodes - 1)
                           ? graph->num_edges
                           : graph->incoming_starts[node + 1];
        bool ancestor_found = false;
        if(sol->distances[node] == NOT_VISITED_MARKER){
            for(int  ancestor_i = incoming_start; ancestor_i < incoming_end; ancestor_i++){
                int ancestor = graph->incoming_edges[ancestor_i];
                if(frontier[ancestor] == true){
                    sol->distances[node] = sol->distances[ancestor] + 1;
                    ancestor_found = true;
                    break; 
                }
            }
        }
        
         new_frontier[node] = ancestor_found;
            
        }
        bool* temp = frontier;
        frontier = new_frontier;
        new_frontier = temp;
        int frontier_count = 0;
        
#pragma omp parallel for reduction(+:frontier_count)
        for(int i = 0; i < graph->num_nodes; i++){
            frontier_count += frontier[i];
        }
        return frontier_count;
}
void transformFrontierToDownBFS(Graph graph, vertex_set* frontier, bool* frontier_bool){
    #pragma omp parallel for
    for(int i = 0; i < graph->num_nodes; i++){
        frontier_bool[i] = false;
    }

    #pragma omp parallel for
    for(int i = 0; i < frontier->count; i++){
        frontier_bool[frontier->vertices[i]] = true;
    }
}
void transformFrontierToTopBFS(Graph graph, bool* frontier_bool, vertex_set* frontier, int frontier_count){
    frontier->count = frontier_count;
    int index = 0;
    #pragma omp parallel for
    for(int i = 0; i < graph->num_nodes; i++){
        if (frontier_bool[i] == true){
            int v_index = 0;
            #pragma omp atomic capture
            {
                v_index = index;
                index++;
            }
            frontier->vertices[v_index] = i;
        }
    }
    
}
void bfs_hybrid(Graph graph, solution* sol)
{
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;
    bool* frontier_bool = new bool [graph->num_nodes];
    bool* new_frontier_bool = new bool [graph->num_nodes];

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;
    int frontier_count = frontier->count;
    bool prev_step = TOP_BFS;
    while(frontier_count != 0){
        
        if(frontier_count > FRONTIER_TOP_BFS_THRESHOLD){
            if(prev_step == TOP_BFS){
                transformFrontierToDownBFS(graph, frontier, frontier_bool);
            }
            frontier_count = bfs_bottom_up_step(graph, sol, frontier_bool, new_frontier_bool);
            prev_step = DOWN_BFS;
            bool* temp = frontier_bool;
            frontier_bool = new_frontier_bool;
            new_frontier_bool = temp;
        }
        else{
            if(prev_step == DOWN_BFS){
                transformFrontierToTopBFS(graph, frontier_bool, frontier, frontier_count);
            }
            frontier_count = bfs_top_down_step(graph, sol, frontier, new_frontier);
            vertex_set* temp = frontier;
            frontier = new_frontier;
            new_frontier = temp;
            prev_step = TOP_BFS;
            vertex_set_clear(new_frontier);
        }
        
    }
    delete [] frontier_bool;
    delete [] new_frontier_bool;
    
}
