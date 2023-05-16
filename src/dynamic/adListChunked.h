#ifndef ADLISTCHUNKED_H_
#define ADLISTCHUNKED_H_

#include <algorithm>
#include <chrono>
#include <queue>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>
#include "adListPerChunk.h"

using namespace std;

bool compare_and_swap(bool &x, const bool &old_val, const bool &new_val);

template <typename T> class Neighborhood;
template <typename T> class neighborhood_iter;

// T can be either node or nodeweight
template <typename T>
class adListChunked: public dataStruc {
    friend class neighborhood_iter<adListChunked<T>>;
    friend class Neighborhood<adListChunked<T>>;
    private:
        class partition {
        	friend adListChunked<T>;
            friend class neighborhood_iter<adListChunked<T>>;
            friend class Neighborhood<adListChunked<T>>;

          private:
            int64_t label;
            int64_t num_directed_partitions;

            adListPerChunk<T>* partAdList;
        	EdgeQueue q;
        	std::mutex q_mutex;
            void insert(EdgeList el);

          public:
        	partition(int64_t label, int64_t _num_partitions, bool w, bool d, int64_t _num_nodes);
	        ~partition();
	        inline void enqueue(Edge const &e);
        };
      int64_t num_nodes_initialize; // The max amount of nodes we would initialize.
      // static const int64_t num_partitions = 8;      
      int64_t num_partitions;
      vector<unique_ptr<partition>> in, out;
      static void dequeue_loop(partition *pt, volatile bool& done);
      inline int32_t pt_hash(NodeID const &n) const;
      inline int32_t hash_within_chunk(NodeID const &n) const;
      
    public:  
      adListChunked(bool w, bool d, int64_t _num_nodes, int64_t _num_parts);
      //adListChunked(bool w, bool d, int64_t _num_nodes);
      ~adListChunked();   
      void update(const EdgeList& el) override;
      void print() override;
      int64_t in_degree(NodeID n) override;
      int64_t out_degree(NodeID n) override;      
};

// dfa----------------------------------Partition----------------------------------------
template <typename T>
adListChunked<T>::partition::partition(int64_t label, int64_t _num_partitions, bool w, bool d, int64_t _num_nodes): label(label), num_directed_partitions(_num_partitions) {
    int64_t num_nodes_per_partition = ceil(((double) _num_nodes) / ((double) _num_partitions));
    partAdList = new adListPerChunk<T>(w, d, num_nodes_per_partition);
}

template <typename T>
adListChunked<T>::partition::~partition() {
    delete partAdList;
}

template <typename T>
void adListChunked<T>::partition::insert(EdgeList el) {
    for (auto& e:el) 
        e.source = (int) e.source/num_directed_partitions;
    partAdList -> update(el);
}

template <typename T>
void adListChunked<T>::partition::enqueue(Edge const &e) {
    lock_guard<std::mutex> guard(q_mutex);
    q.push(e);
}

// // ---------------------------------adListChunked--------------------------------------
template <typename T>
adListChunked<T>::adListChunked(bool w, bool d, int64_t _num_nodes, int64_t _num_parts)
: dataStruc(w, d){

#ifdef _OPENMP
		if(_num_parts > 0){
			omp_set_num_threads(_num_parts);
		}
#endif

    num_nodes_initialize = _num_nodes;
    if (_num_parts % 2)
        num_partitions = _num_parts - 1;
    else
        num_partitions = _num_parts;
    cout << "Num parts: " << _num_parts << endl;
    // initialize 1) property 2) affected 3) vertices vectors 4) markers
    property.resize(num_nodes_initialize, -1);
    affected.resize(num_nodes_initialize); affected.fill(false);

    for (int i = 0; i < num_partitions / 2; i++) {
        if (directed) {
	        in.push_back(
                unique_ptr<partition>(new partition(i, num_partitions / 2, w, d, _num_nodes))
            );
            
	        out.push_back(
                unique_ptr<partition>(new partition(i, num_partitions / 2, w, d, _num_nodes))
            );
        }
        else {
            out.push_back(
                unique_ptr<partition>(new partition(2*i, num_partitions, w, d, _num_nodes))
            );
        
	        out.push_back(
                unique_ptr<partition>(new partition(2*i + 1, num_partitions, w, d, _num_nodes))
            );
        }
    }
} 

//template <typename T>
//adListChunked<T>::adListChunked(bool w, bool d, int64_t _num_nodes)
//: dataStruc(w, d){
//    num_nodes_initialize = _num_nodes;
//    num_partitions = 16;
//
//    // initialize 1) property 2) affected 3) vertices vectors 4) markers
//    property.resize(num_nodes_initialize, -1);
//    affected.resize(num_nodes_initialize); affected.fill(false);
//
//    for (int i = 0; i < num_partitions / 2; i++) {
//        if (directed) {
//	        in.push_back(
//                unique_ptr<partition>(new partition(i, num_partitions / 2, w, d, _num_nodes))
//            );
//
//	        out.push_back(
//                unique_ptr<partition>(new partition(i, num_partitions / 2, w, d, _num_nodes))
//            );
//        }
//        else {
//            out.push_back(
//                unique_ptr<partition>(new partition(2*i, num_partitions, w, d, _num_nodes))
//            );
//
//	        out.push_back(
//                unique_ptr<partition>(new partition(2*i + 1, num_partitions, w, d, _num_nodes))
//            );
//        }
//    }
//}

template <typename T>
adListChunked<T>::~adListChunked(){

}

template <typename T>
int32_t adListChunked<T>::pt_hash(NodeID const &n) const {
    if (directed)
        return n % (num_partitions/2);
    else
        return n % num_partitions;
}

template <typename T>
int32_t adListChunked<T>::hash_within_chunk(NodeID const &n) const {
    if (directed)
        return  (int) n/(num_partitions / 2);
    else
        return  (int) n/(num_partitions);
}

template <typename T>
void adListChunked<T>::dequeue_loop(partition *partPtr, volatile bool& done){
	LIKWID_MARKER_START("udp");
    EdgeList el;
    partPtr->q_mutex.lock();    
    while (!done || !partPtr->q.empty()) {
        if (partPtr->q.empty()) {
            // no edges, keep waiting
            partPtr->q_mutex.unlock();
            this_thread::sleep_for(chrono::milliseconds(1));
        }
        else {
            // read out all the edges
            while (!partPtr->q.empty()) {
                el.push_back(partPtr->q.front());
                partPtr->q.pop();
            }
            // unlock            
            partPtr->q_mutex.unlock();
            partPtr->insert(el);
            el.clear();
        }
        partPtr->q_mutex.lock();
    }
    partPtr->q_mutex.unlock();
    LIKWID_MARKER_STOP("udp");
}

template <typename T>
void adListChunked<T>::update(const EdgeList& el) {
    bool done = false;
    vector<unique_ptr<thread>> dqs;
    
    //#######............. thread pinning.............#########
    int64_t count = 0;
    if (directed) {
        //cout << "IN" << endl;
        for(auto& ptr:in){        
            //cout << "Count: " << count << " Cpu: " << count+2 << endl;
            //unique_ptr<thread> t= new thread t(dequeue_loop, ptr.get(), ref(done));
            thread* t = new thread(dequeue_loop, ptr.get(), ref(done));
            dqs.push_back(unique_ptr<thread>(t));

            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(count, &cpuset);

            int rc = pthread_setaffinity_np(t->native_handle(), sizeof(cpu_set_t), &cpuset);
        
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }       
            count++;                  
        }

        //assert(count == 22);
        //count = 26;
        
        //cout << "OUT" << endl;
        for (auto& ptr:out){
            //cout << "Count: " << count << " Cpu: " << count+2 << endl;
            thread* t = new thread(dequeue_loop, ptr.get(), ref(done));
            dqs.push_back(unique_ptr<thread>(t));

            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(count, &cpuset);

            int rc = pthread_setaffinity_np(t->native_handle(), sizeof(cpu_set_t), &cpuset);
        
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }        
            count++;             
        } 
    } 
    else {
        //cout << "NO threads for IN." << endl
             //<< "OUT" << endl;
        for (auto& ptr:out){
            //cout << "Count: " << count << " Cpu: " << count + 2 << endl;
            thread* t = new thread(dequeue_loop, ptr.get(), ref(done));
            dqs.push_back(unique_ptr<thread>(t));

            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(count, &cpuset);

            int rc = pthread_setaffinity_np(t->native_handle(), sizeof(cpu_set_t), &cpuset);
        
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }        
            count++;
            //if(count == 22) count = 26;
        } 
    }
    //#######............. thread pinning.............#########
            
    /*for (auto& ptr:in)
        dqs.push_back(
            unique_ptr<thread>(new thread(dequeue_loop, ptr.get(), ref(done)))
        );
    //for (unsigned int k = 0; k < num_partitions; k++)
    for (auto& ptr:out) 
        dqs.push_back(
            unique_ptr<thread>(new thread(dequeue_loop, ptr.get(), ref(done)))
        );*/

    int o_ix, i_ix;
    LIKWID_MARKER_START("q");
    for(unsigned int i=0; i<el.size(); i++){
        Edge e_reverse = el[i].reverse();
        
        // update affected
        num_edges += 2;
        if (!el[i].sourceExists)
            num_nodes++;
        if (!el[i].destExists)
            num_nodes++;
        affected[el[i].source] = true;
        affected[el[i].destination] = true;

	    o_ix = pt_hash(el[i].source);
	    i_ix = pt_hash(e_reverse.source);

        if (!directed) {
            out[o_ix] -> enqueue(el[i]);
            out[i_ix] -> enqueue(e_reverse);
        }
        else {
            out[o_ix] -> enqueue(el[i]);
            in[i_ix] -> enqueue(e_reverse);
        }
    }
    LIKWID_MARKER_STOP("q");

    done = true;
    for(auto& dq: dqs) {
	dq->join();
    }     
}

template <typename T>
int64_t adListChunked<T>::in_degree(NodeID n) {
    int32_t part_ix = pt_hash(n);
    int64_t sub_idx = hash_within_chunk(n); 
    if(directed) 
        return in[part_ix]->partAdList->degree(sub_idx);
    else 
        return out[part_ix]->partAdList->degree(sub_idx);
}

template <typename T>
int64_t adListChunked<T>::out_degree(NodeID n) {   
    int part_ix = pt_hash(n);
    int64_t sub_idx = hash_within_chunk(n);
    return out[part_ix]->partAdList->degree(sub_idx);
}


template <typename T>
void adListChunked<T>::print(void) {
    u64 totMem = 0;
    totMem += in.capacity()*sizeof(unique_ptr<partition>);
    totMem += out.capacity()*sizeof(unique_ptr<partition>);
    for (const auto &p : in) {
        totMem += p->q.size() * sizeof(Edge);
        totMem += sizeof(adListPerChunk<T>);
        totMem += p->partAdList->neighbors.capacity()*sizeof(std::vector<T>);
        totMem += sizeof(partition);
        for (const auto &nei : p->partAdList->neighbors) {
            totMem += sizeof(T)*nei.capacity();
        }
    }

    for (const auto &p : out) {
        totMem += p->q.size() * sizeof(Edge);
        totMem += p->partAdList->neighbors.capacity()*sizeof(std::vector<T>);
        totMem += sizeof(adListPerChunk<T>);
        totMem += sizeof(partition);
        for (const auto &nei : p->partAdList->neighbors) {
            totMem += sizeof(T)*nei.capacity();
        }
    }

    std::cout << "Total Memory adListChunked: " << totMem << std::endl;

//	for(const auto& p : in){
//		insTot += p->partAdList->insTot;
//		insSucc += p->partAdList->insSucc;
//		delTot += p->partAdList->delTot;
//		delSucc += p->partAdList->delSucc;
//	}
//	for(const auto& p : out){
//		insTot += p->partAdList->insTot;
//		insSucc += p->partAdList->insSucc;
//		delTot += p->partAdList->delTot;
//		delSucc += p->partAdList->delSucc;
//	}
//
//	std::cout << "Inserts--------------------" << std::endl;
//	std::cout << "    Total: " << insTot << std::endl;
//	std::cout << "    Succ : " << insSucc << std::endl;
//	std::cout << "    Fail : " << insTot - insSucc << std::endl;
//	std::cout << std::endl;
//
//	std::cout << "Deletes--------------------" << std::endl;
//	std::cout << "    Total: " << delTot << std::endl;
//	std::cout << "    Succ : " << delSucc << std::endl;
//	std::cout << "    Fail : " << delTot - delSucc << std::endl;
//	std::cout << std::endl;
//
//	std::cout << "Final number of edges: " << insSucc - delSucc << std::endl;
}
#endif  // ADLISTCHUNKED_H_
