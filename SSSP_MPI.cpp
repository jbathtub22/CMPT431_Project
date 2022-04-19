#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

#define SOURCE "0"
#define Inf 999999

void display_result(uintV source_vertex, int total_n, int* global_dist, int* global_pred) {

    std::cout << "Source vertex : " << source_vertex << std::endl;  
    std::cout << "Vertex_id,  min_distance,  Predecessor" << std::endl;
    for (uintV i = 0; i < total_n; i++)
    {
        if (global_dist[i] == Inf)
        {
            std::cout << i << "  No path to source vertex"<< std::endl;
        } else {
            std::cout << i << ",      " << global_dist[i] << ",      " << global_pred[i] << std::endl;
        }            
    }
    
}

// find the closest vertex to source in each process
int find_local_min(int local_n, int* local_dist, bool* marked) {
    int min_v = -1; 
    int local_min = INT_MAX;
    for (int i = 0; i < local_n; i++) {
        if (marked[i] == false && local_dist[i] < local_min) {
            local_min = local_dist[i];
            min_v = i;
        }
    }
    return min_v;
}  


void SSSP(Graph &g, int source_vertex, std::string y_or_n) 
{
    //----- Initialization -----

    timer t1;
    double time;

	int world_size;
	int world_rank;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int local_n; 
    int total_n = g.n_;

    // Equally distribute vertices among processes
    int* count_arr = new int[world_size];
    for(int i = 0; i < world_size-1; i++) {
        count_arr[i] = total_n/world_size;
    }
    count_arr[world_size-1] = total_n - (world_size-1)*(total_n/world_size);
    local_n = count_arr[world_rank];

    int* displc_arr = new int[world_size];
    for (int j = 0; j < world_size; j++) {
        displc_arr[j] = j * (total_n/world_size);
    }
    
    bool* marked = new bool[local_n];
    int* local_dist = new int[local_n];
    int* local_pred = new int[local_n];
    for (int i = 0; i < local_n; i++) {
        marked[i] = false;
        local_dist[i] = Inf;
        local_pred[i] = source_vertex;
    }
    
    // Initialization for source vertex
    // case of the last process
    if (source_vertex >= (world_size-1)*(total_n/world_size)) {
        int i = source_vertex - world_rank*(total_n/world_size);
        marked[i] = true;
        local_dist[i] = 0;
    } else {
    // case of the other processes
        if (world_rank == source_vertex/(total_n/world_size)) {
            int i = source_vertex - world_rank*(total_n/world_size);
            marked[i] = true;
            local_dist[i] = 0;
        }
    }

    // Initialization for neighbors of source vertex
    uintE out_degree = g.vertices_[source_vertex].getOutDegree();
    for (uintE j = 0; j < out_degree; j++) {
        uintV glb_v = g.vertices_[source_vertex].getOutNeighbor(j);

        if (glb_v >= (world_size-1)*(total_n/world_size)) {
            int local_v = glb_v - world_rank*(total_n/world_size);
            local_dist[local_v] = 1;
        } else {
            if (world_rank == glb_v/(total_n/world_size)) {
                int local_v = glb_v - world_rank*(total_n/world_size);
                local_dist[local_v] = 1;
            }
        }
    }

    int local_min_v; // local vertex that is closest to the source vertex
    int global_min_v; // gloval vertex that is cloest to the source vertex
    int global_min_dist;
    int* local_min_pair = new int[2];
    int* global_min_pair = new int[2];
    //----- End of Initialization -----------------


    t1.start();
    //----- Dijkstra Algorithm -----
    for (int step = 0; step < total_n-1; step++) {
        local_min_v = find_local_min(local_n, local_dist, marked); 
        //std::cout << world_rank << "  step:  " << step << "local_min " << local_min_v << std::endl;
        if (local_min_v != -1) {
            local_min_pair[0] = local_dist[local_min_v];
            local_min_pair[1] = local_min_v + world_rank * (total_n/world_size);
        } else {
            local_min_pair[0] = Inf;
            local_min_pair[1] = -1;           
        }
        
        // find the global minumum distance and boradcast it to all the processes
        MPI_Allreduce(local_min_pair, global_min_pair, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        global_min_dist = global_min_pair[0];
        global_min_v = global_min_pair[1];
             
        if (global_min_v == -1) {
            break;
        }

        // update marked
        if (global_min_v >= (world_size-1)*(total_n/world_size)) {
            int local_i = global_min_v - world_rank*(total_n/world_size);
            marked[local_i] = true;
        } else {
            if (world_rank == global_min_v/(total_n/world_size)) {
                int local_i = global_min_v - world_rank*(total_n/world_size);
                marked[local_i] = true;
            }
        }

        // update local dist for each process       
        uintE out_degree = g.vertices_[global_min_v].getOutDegree();
        int update_dist;
        for (uintE j = 0; j < out_degree; j++) {
            int glb_k = g.vertices_[global_min_v].getOutNeighbor(j);
            //std::cout << k << std::endl;
            if (glb_k >= (world_size-1)*total_n/world_size) {
                int local_k = glb_k - world_rank*(total_n/world_size);
                if (marked[local_k] == false) {
                    int update_dist = global_min_dist + 1;
                    if (update_dist < local_dist[local_k]) {
                        local_dist[local_k] = update_dist;
                        local_pred[local_k] = global_min_v;
                    }
                }
            } else {
                if (world_rank == glb_k/(total_n/world_size)) {
                    int local_k = glb_k - world_rank*(total_n/world_size);
                    if (marked[local_k] == false) {
                        int update_dist = global_min_dist + 1;
                        if (update_dist < local_dist[local_k]) {
                            local_dist[local_k] = update_dist;
                            local_pred[local_k] = global_min_v;
                        }
                    }
                }
            }
        }
    }

    int* global_dist;
    int* global_pred;

    if (world_rank == 0) {
        global_dist = new int[total_n];
        global_pred = new int[total_n];
    }
    // root process gather all the information
    MPI_Gatherv(local_dist, local_n, MPI_INT, global_dist, count_arr, displc_arr, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(local_pred, local_n, MPI_INT, global_pred, count_arr, displc_arr, MPI_INT, 0, MPI_COMM_WORLD);
    
    time = t1.stop();
    
    //----- Print result -----
    if (world_rank == 0) {
        if (y_or_n == "yes") {
            display_result(source_vertex, total_n, global_dist, global_pred);
        }

        delete[] global_dist;
        delete[] global_pred;
    }
    if (world_rank == 0) {
        std::cout << "running time: " << time << std::endl;
    }

    delete[] marked;
    delete[] local_dist;
    delete[] local_pred;
    delete[] local_min_pair;
    delete[] global_min_pair;
    delete[] count_arr;
    delete[] displc_arr;

    MPI_Finalize();
}

int main (int argc, char *argv[]){
    
    cxxopts::Options options("SSSP_MPI", "Calculate SSSP using MPI");
    options.add_options(
        "",
        {
            {"sourceVertex", "Source vertex",
            cxxopts::value<int>()->default_value(SOURCE)},
            {"inputFile", "Input graph file path",
            cxxopts::value<std::string>()->default_value(
            "/scratch/input_graphs/roadNet-CA")},
            {"y_or_n", "displayOutput",
            cxxopts::value<std::string>()->default_value("no")}
        });

    auto cl_options = options.parse(argc, argv);
    int source_vertex = cl_options["sourceVertex"].as<uintV>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    std::string y_or_n = cl_options["y_or_n"].as<std::string>();

    Graph g;
    g.readGraphFromBinary<int>(input_file_path);
    SSSP(g, source_vertex, y_or_n);    

    return 0;
}