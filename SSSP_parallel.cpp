#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <queue>
#include <thread>
#include <atomic> 
#include "common/non_blocking_queue.h"

typedef struct thread_obj{
  NonBlockingQueue<uintV> *g_queue_1;
  NonBlockingQueue<uintV> *g_queue_2;
  Graph *sub_gr;
  CustomBarrier *barrier;
  uint *dist_array;
  uintV *prev_vertex_array;
  uint pid;  
} thread_obj; 


void dijkstra_task(thread_obj *thread, std::atomic<bool> *check_if_empty){
    NonBlockingQueue<uintV> *g_queue_1 = thread->g_queue_1;
    NonBlockingQueue<uintV> *g_queue_2 = thread->g_queue_2;
    Graph *graph = thread->sub_gr; 
    CustomBarrier *barrier = thread->barrier;
    uint pid = thread->pid;

    uintV n = graph->n_;
    uint *dist_array = thread->dist_array;
    uintV *prev_vertex_array = thread->prev_vertex_array;
    uintV current_vertex=0;

    while(check_if_empty->load()==false){
    while (g_queue_1->dequeue(&current_vertex))
    {

        uintE out_degree = graph->vertices_[current_vertex].getOutDegree();
        for (uintE i = 0; i < out_degree; i++)
        {
            uintV adjacent_vertex = graph->vertices_[current_vertex].getOutNeighbor(i);
            if (dist_array[adjacent_vertex] == INT_MAX)
                {
                    dist_array[adjacent_vertex] = dist_array[current_vertex] + 1;
                    prev_vertex_array[adjacent_vertex] = current_vertex;
                    g_queue_2->enqueue(adjacent_vertex);

                
                }

        }
    }

    //Swapping two queues after finishing one layer
    barrier->wait(); 
    NonBlockingQueue<uintV> *g_queue_temp = g_queue_1; 
    g_queue_1 = g_queue_2;
    g_queue_2 = g_queue_temp; 

    //If first process, check if queue is empty.
    if(pid==0){
        if(g_queue_1->dequeue(&current_vertex)){
            uintE out_degree = graph->vertices_[current_vertex].getOutDegree();
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV adjacent_vertex = graph->vertices_[current_vertex].getOutNeighbor(i);
                if (dist_array[adjacent_vertex] == INT_MAX)
                {
                    dist_array[adjacent_vertex] = dist_array[current_vertex] + 1;
                    prev_vertex_array[adjacent_vertex] = current_vertex;
                    g_queue_2->enqueue(adjacent_vertex);
                
                }

            }
        }
        //If queue is empty, notify other processes
        else{
             check_if_empty->store(true);
        }
    }
    barrier->wait(); 
    }

}


void parallel_dijkstra(Graph *graph, uintV nThreads, uint *dist_array, uintV *prev_vertex_array, uintV source_vertex, double *time_taken){

    //Creating two queues. One queue for dequeuing from current layer, other queue for enqueueing in the next layer 
    NonBlockingQueue<uintV> *g_queue_1 = new NonBlockingQueue<uintV>();
    NonBlockingQueue<uintV> *g_queue_2 = new NonBlockingQueue<uintV>();

    uint number_vertex = graph->n_;
    std::thread threads[nThreads];

    //Initializing thread function parameters object 
    thread_obj threads_para[nThreads]; 

    //Barrier
    CustomBarrier barrier(nThreads);

    //Initializing queues
    g_queue_1->initQueue(number_vertex);
    g_queue_2->initQueue(number_vertex);

    timer t1;
    g_queue_1->enqueue(source_vertex);
    dist_array[source_vertex] = 0;
    prev_vertex_array[source_vertex] = source_vertex;
    std::atomic<bool> check_if_empty;
    t1.start();
    check_if_empty.store(false); 

    
    for(uint j = 0; j < nThreads; j++){
      threads_para[j].g_queue_1 = g_queue_1; 
      threads_para[j].g_queue_2 = g_queue_2;
      threads_para[j].sub_gr = graph;
      threads_para[j].barrier = &barrier;
      threads_para[j].pid = j;
      threads_para[j].dist_array = dist_array; 
      threads_para[j].prev_vertex_array = prev_vertex_array;
      threads[j] = std::thread(dijkstra_task, &threads_para[j], &check_if_empty);
      

  }

    for(int i = 0; i < nThreads; i++){
        threads[i].join();
    }

    *time_taken = t1.stop();
}


void display_result(uint *dist_array, uintV *prev_vertex_array, uintV number_vertex, double time_taken){
    std::cout << "Vertex id, Shortest distance to source, Previous vertex on the path" << std::endl;
    for (uintV i = 0; i < number_vertex; i++)
    {
        if (dist_array[i] == INT_MAX)
        {
            std::cout << i << " No path, No previous vertex"<< std::endl;
        }
        else
        {
            std::cout << i << ", " << dist_array[i] << ", " << prev_vertex_array[i] << std::endl;
        }
    }
    std::cout <<"\n"; 
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "SSSP_serial",
      "Calculate SSSP using serial execution");
  options.add_options(
      "",
      {
          {"nThreads", "Number of threads",
           cxxopts::value<uintV>()->default_value(DEFAULT_MAX_ITER)},
          {"sourceVertex", "Maximum number of iterations",
           cxxopts::value<uintV>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
          {"displayOutput", "Indicating printing output or not",
           cxxopts::value<std::string>()->default_value("no")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uintV>();
  uintV source_vertex = cl_options["sourceVertex"].as<uintV>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  std::string y_or_n = cl_options["displayOutput"].as<std::string>();


  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  std::cout << "Number of threads: " << n_threads << "\n"; 
  std::cout << "Graph being tested: " << input_file_path << "\n";
  uintV number_vertex = g.n_;
  double time_taken = 0.0;
  uint *dist_array = new uint[number_vertex];
  uintV *prev_vertex_array = new uintV[number_vertex];
  for (uintV i = 0; i < number_vertex; i++)
  {
      dist_array[i] = INT_MAX;
  }

  
  parallel_dijkstra(&g, n_threads, dist_array, prev_vertex_array, source_vertex, &time_taken);

  if(y_or_n == "yes"){
  display_result(dist_array, prev_vertex_array, number_vertex, time_taken);
  }

  std::cout <<"Total time taken: " << time_taken << std::endl;
  std::cout <<"\n"; 
}
