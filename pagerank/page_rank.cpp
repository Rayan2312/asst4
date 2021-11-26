#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"


// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double* solution, double damping, double convergence)
{


  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
  double* prev_sol = new double [numNodes];
  
  #pragma omp parallel for
  for (int i = 0; i < numNodes; ++i) {
    solution[i] = equal_prob;
    prev_sol[i] = equal_prob;
  }

  #pragma omp barrier
  int numIters = 0;
  
  bool converged = false;
  while( !converged ){
      
      double outgoing_sum = 0;
#pragma omp parallel for reduction(+: outgoing_sum)
      for(int i = 0; i < numNodes; i++){
          if(outgoing_size(g, i) == 0) outgoing_sum += damping*prev_sol[i]/numNodes;
          else outgoing_sum += 0;
      }
      #pragma omp barrier
      
      #pragma omp parallel for
      for(int i = 0; i < numNodes; i++){
          double  sum = 0;
          const Vertex* Start = incoming_begin(g, i);
          const Vertex* End = incoming_end(g, i);
          for(const Vertex* v = Start; v != End; v++){
              sum +=  (prev_sol[*v]/ outgoing_size(g, *v) );
          }
          solution[i] = (sum*damping) + (1.0 - damping)/numNodes;
          sum = 0;
          
          /* #pragma omp parallel for reduction(+:sum)
          for(int j = 0; j < numNodes; j++){
              if (outgoing_size(g, j) == 0) sum = (damping*prev_sol[j]/numNodes);
          } 
          */
          solution[i] += outgoing_sum;
          
      }
      #pragma omp barrier
      
      double* temp = prev_sol;
      prev_sol = solution;
      solution = temp;
      
      
      double global_diff = 0;
      
#pragma omp parallel for reduction(+: global_diff)
              for(int i = 0; i < numNodes; i++){
                  global_diff += abs(solution[i] - prev_sol[i]);
              }
          converged = global_diff <= convergence;
          numIters++;
          
          if(converged && numIters % 2 == 0){
              #pragma omp parallel for
              for(int i = 0; i < numNodes; i++){
                  solution[i] = prev_sol[i];
              }
              delete [] prev_sol;
          }
          else if(converged)  delete [] solution;
      
  }
  
  
  /*
     CS149 students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes v in graph with no outgoing edges
                          { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */
}
