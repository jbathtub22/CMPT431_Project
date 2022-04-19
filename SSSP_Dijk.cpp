#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <vector>

void display_result(int total_n, std::vector<int> local_dist, std::vector<int> local_pred) {

    std::cout << "Source vertex : 0" << std::endl;  
    std::cout << "Vertex id,  Shortest distance to source,  Predecessor vertex" << std::endl;
    for (uintV i = 0; i < total_n; i++)
    {
        if (local_dist[i] == INT_MAX)
        {
            std::cout << i << "No path to source vertex"<< std::endl;
        } else {
            std::cout << i << ", " << local_dist[i] << ", " << local_pred[i] << std::endl;
        }            
    }
    
}

intV find_local_min(int local_n, std::vector<int>& local_dist, std::vector<bool>& marked) {
    intV min_v = -1; 
    int local_min = INT_MAX;
    for (uintV i = 0; i < local_n; i++) {
        if (marked[i] == false && local_dist[i] < local_min) {
            local_min = local_dist[i];
            min_v = i;
        }
    }
    return min_v;
}

void SSSP(Graph &g, uintV source_vertex) 
{
    //----- Initialization -----

    int total_n = g.n_;
    std::vector<bool> marked(total_n, false);
    std::vector<int> local_dist(total_n, INT_MAX);
    std::vector<int> local_pred(total_n, source_vertex);

    marked[source_vertex] = true;
    local_dist[source_vertex] = 0;

    uintE out_degree = g.vertices_[source_vertex].getOutDegree();
    for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[source_vertex].getOutNeighbor(i);
        local_dist[v] = 1;
    }
    
    /*
    std::cout << g.vertices_[source_vertex].getOutDegree() << g.vertices_[source_vertex].getInDegree() << std::endl;

    for (int m=0; m < local_dist.size(); m++) {
        std::cout << local_dist[m] << " ";
    }
    std::cout << std::endl;
    */

 
    intV local_min_v; // local vertex that is closest to the source vertex
    //----- End of Initialization -----------------

    //----- Dijkstra Algorithm -----
    for (int i = 0; i < total_n-1; i++) {
        local_min_v = find_local_min(total_n, local_dist, marked);
        if (local_min_v == -1) {
            break;
        } 

        marked[local_min_v] = true;

        // update local dist for each process     
        //std::cout << "current local_min vetex: " << local_min_v << std::endl; 
        //std::cout << "neighbors" << std::endl;
        uintE out_degree = g.vertices_[local_min_v].getOutDegree();
        int update_dist;
        for (uintE j = 0; j < out_degree; j++) {
            uintV k = g.vertices_[local_min_v].getOutNeighbor(j);
            //std::cout << k << std::endl;
            if (marked[k] == false) {
                update_dist = local_dist[local_min_v] + 1;
                //std::cout << "new_dist: " << update_dist << std::endl;
                if (update_dist < local_dist[k]) {
                    local_dist[k] = update_dist;
                    local_pred[k] = local_min_v;
                }                 
            }
        }                          
    }
    display_result(total_n, local_dist, local_pred);
}

int main (int argc, char *argv[]){
    
    cxxopts::Options options("SSSP_Dijk", "Calculate SSSP using Dijk");
    options.add_options(
        "",
        {
            {"source", "Source vertex",
            cxxopts::value<uintV>()->default_value("1")},
            {"inputFile", "Input graph file path",
            cxxopts::value<std::string>()->default_value(
            "/scratch/input_graphs/roadNet-CA")}
        });

    auto cl_options = options.parse(argc, argv);
    uintV source_vertex = cl_options["source"].as<uintV>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";


    SSSP(g, source_vertex);    

    return 0;
}