#ifndef DARHH_H
#define DARHH_H

#include <algorithm>
#include <chrono>
#include <queue>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

#include "abstract_data_struc.h"

//#include "darhh_ld.h"
//#include "darhh_hd.h"
#include "rhh.h"
#include "common.h"
#include "omp.h"

/* Data Structure: Degree-Aware Hashing */

template<typename U> class neighborhood;
template<typename U> class neighborhood_iter;

template<typename T>
class darhh: public dataStruc {
	friend neighborhood<darhh<T>> ;
	friend neighborhood_iter<darhh<T>> ;

	using super = dataStruc;
	using HDHashMap = rhh<NodeID, Weight>;
	using LDHashMap = rhh<EdgeID, Weight>;
	//using HDHashMap = unordered_map<NodeID, Weight>;
	//using LDHashMap = unordered_map<EdgeID, Weight>;

	static const u64 LD_THRESHOLD = 2;

	typedef struct{
		i64 nodeCnt = 0;
		i64 edgeCnt = 0;
		LDHashMap inLDHash;
		LDHashMap outLDHash;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	void transferL2H(NodeID src, LDHashMap* __restrict ldMap, HDHashMap* __restrict hdMap){
		u32 i = src % ldMap->arr.capacity();
		const u32 cap = ldMap->arr.capacity();
		while(true){
			if(ldMap->arr[i].empty()){
				break;
			}
			else if(ldMap->arr[i].key.first == src){
				NodeID dst = ldMap->arr[i].key.second;
				Weight val = ldMap->arr[i].val;
				ldMap->arr[i].mark_deleted();
				ldMap->sz--;
				hdMap->insert(dst, val);
			}
			i++;
			if(i == cap){
				i = 0;
			}
		}
	}

	void insertOutEdge(u64 thId, const Edge& e){
		u64& degree = vArray[e.source].outDegree;
		u32 prevSize, newSize;
		if(degree < LD_THRESHOLD){	//use low degree map
			LDHashMap* __restrict hashMap = &thInfo[thId].outLDHash;
			prevSize = hashMap->size();
			hashMap->insert({e.source,e.destination}, e.weight);
			//hashMap->at({e.source,e.destination}) = e.weight;
			newSize = hashMap->size();
		}
		else{ //use high degree map
			if(!vArray[e.source].outMap){
				vArray[e.source].outMap = new HDHashMap;
				//move from low to high degree map
				transferL2H(e.source, &thInfo[thId].outLDHash, vArray[e.source].outMap);
			}
			HDHashMap* hashMap = vArray[e.source].outMap;
			prevSize = hashMap->size();
			hashMap->insert(e.destination, e.weight);
			//hashMap->at(e.destination) = e.weight;
			newSize = hashMap->size();
		}
		if(prevSize != newSize){
			degree++;
			thInfo[thId].edgeCnt++;
		}
	}

	void insertInEdge(u64 thId, const Edge& e){
		u64& degree = vArray[e.destination].inDegree;
		u32 prevSize, newSize;
		if(degree < LD_THRESHOLD){	//use low degree map
			LDHashMap* __restrict hashMap = &thInfo[thId].inLDHash;
			prevSize = hashMap->size();
			hashMap->insert({e.destination,e.source}, e.weight);
			//hashMap->at({e.source,e.destination}) = e.weight;
			newSize = hashMap->size();
		}
		else{ //use high degree map
			if(!vArray[e.destination].inMap){
				vArray[e.destination].inMap = new HDHashMap;
				transferL2H(e.destination, &thInfo[thId].inLDHash, vArray[e.destination].inMap);
			}
			HDHashMap* hashMap = vArray[e.destination].inMap;
			prevSize = hashMap->size();
			hashMap->insert(e.source, e.weight);
			//hashMap->at(e.source) = e.weight;
			newSize = hashMap->size();
		}
		if(prevSize != newSize){
			degree++;
		}
	}

	void deleteOutEdge(u64 thId, const Edge& e){
		u32 prevSize, newSize;
		if(!vArray[e.source].outMap){
			//low degree vertex
			LDHashMap* __restrict hashMap = &thInfo[thId].outLDHash;
			prevSize = hashMap->size();
			hashMap->erase({e.source,e.destination});
			newSize = hashMap->size();
		}
		else{
			//high degree vertex
			HDHashMap* hashMap = vArray[e.source].outMap;
			prevSize = hashMap->size();
			hashMap->erase(e.destination);
			newSize = hashMap->size();
		}

		if(prevSize != newSize){
			 vArray[e.source].outDegree--;
			 thInfo[thId].edgeCnt--;
		}
	}

	void deleteInEdge(u64 thId, const Edge& e){
		u32 prevSize, newSize;
		if(!vArray[e.destination].inMap){
			//low degree vertex
			LDHashMap* __restrict hashMap = &thInfo[thId].inLDHash;
			prevSize = hashMap->size();
			hashMap->erase({e.destination,e.source});
			newSize = hashMap->size();
		}
		else{
			//high degree vertex
			HDHashMap* hashMap = vArray[e.destination].inMap;
			prevSize = hashMap->size();
			hashMap->erase(e.source);
			newSize = hashMap->size();
		}

		if(prevSize != newSize){
			 vArray[e.destination].inDegree--;
		}
	}

	typedef struct {
		u64 inDegree = 0;
		u64 outDegree = 0;
		HDHashMap* __restrict inMap = nullptr;
		HDHashMap* __restrict outMap = nullptr;
	} Vertex;

	//static void dequeue_loop(partition *pt, volatile bool &done);
	//inline int32_t pt_hash(NodeID const &n) const;

	const int64_t init_num_nodes;
	const u64 num_threads;
	//std::vector<std::unique_ptr<partition>> in, out;
	Vertex* vArray;

public:
	darhh(bool w, bool d, int64_t init_nn, int64_t nt);
	void update(EdgeList const &el) override;
	int64_t out_degree(NodeID n) override;
	int64_t in_degree(NodeID n) override;
	void print() override;
	std::string to_string() const;
};

//template<typename T>
//darhh<T>::partition::partition(darhh *parent) :
//		parent(parent) {
//	ld = new ld_rhh<T>();
//	hd = new hd_rhh<T>();
//}
//
//template<typename T>
//darhh<T>::partition::~partition() {
//	delete ld;
//	delete hd;
//}
//
//template<typename T>
//void darhh<T>::partition::transfer_low_to_high(NodeID const &n) {
//	for (auto it = ld->begin(n), end = ld->end(n); it != end; ++it) {
//		EdgeID id(n, it.cursor->getNodeID());
//		hd->insert_elem(id, it.cursor->getWeight());
//		ld->arr[it.pos].mark_deleted();
//	}
//}
//
//template<typename T>
//void darhh<T>::partition::insertInEdge(NodeID src, T neigh) {
//	u64 degree = vArray[src]
//
//
//	if (!e.sourceExists) {
//		ld->insert_elem(e);
//		//std::lock_guard<std::mutex> guard(parent->num_nodes_mutex);
//		//++parent->super::num_nodes;
//	} else {
//		int deg = ld->get_degree(e.source);
//		if (deg > 0 && deg < ld_threshold) {
//			ld->insert_elem(e);
//		} else if (deg == ld_threshold) {
//			transfer_low_to_high(e.source);
//			hd->insert_elem(e);
//		} else {
//			if (hd->get_degree(e.source) > 0)
//				hd->insert_elem(e);
//			else
//				ld->insert_elem(e);
//		}
//	}
//}

//template<typename T>
//void darhh<T>::partition::erase(Edge const &e) {
//	//__sync_fetch_and_add(&parent->delTot, 1);
//	ld->delete_elem(EdgeID(e.source, e.destination));
//	//hd->delete_elem(EdgeID(e.source, e.destination));
//}
//
//template<typename T>
//void darhh<T>::partition::enqueue(Edge const &e) {
//	std::lock_guard<std::mutex> guard(q_mutex);
//	q.push(e);
//}

template<typename T>
darhh<T>::darhh(bool w, bool d, int64_t init_nn, int64_t nt) :
		super(w, d), init_num_nodes(init_nn), num_threads(nt){
	super::property.resize(init_num_nodes, -1);
	super::affected.resize(init_num_nodes);
	super::affected.fill(false);
	std::cout << "Vertex class size: " << sizeof(Vertex) << endl;
	std::cout << "ThreadInfo class size: " << sizeof(ThreadInfo) << endl;
	vArray = (Vertex*)aligned_alloc(64, sizeof(Vertex) * init_num_nodes);
	memset(vArray, 0, sizeof(Vertex) * init_num_nodes);
	for(int i = 0; i < num_threads; i++){
		thInfo[i].edgeCnt = 0;
		thInfo[i].nodeCnt = 0;
	}
	omp_set_num_threads(num_threads);
}

template<typename T>
void darhh<T>::update(EdgeList const &el) {
	const u64 batchSize = el.size();


	#pragma omp parallel
	for(u64 i = 0; i < batchSize; i++){
		const i64 src = el[i].source;
		const i64 dst = el[i].destination;
		const i64 actualTh = omp_get_thread_num();

		i64 targetTh = (src / 64) % num_threads;
		if(targetTh == actualTh){
			if(!el[i].sourceExists){
				thInfo[actualTh].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[actualTh].nodeCnt++;
			}
			if(!affected[src]){
				affected[src] = true;
			}

			if(!el[i].isDelete){
				//insert out edge
				insertOutEdge(targetTh, el[i]);
			}
			else{
				//delete out edge
				deleteOutEdge(targetTh, el[i]);
			}
		}

		targetTh = (dst / 64) % num_threads;
		if(targetTh == actualTh){
			if(!affected[dst]){
				affected[dst] = true;
			}

			if(!el[i].isDelete){
				//insert in edge
				insertInEdge(targetTh, el[i]);
			}
			else{
				//delete in edge
				deleteInEdge(targetTh, el[i]);
			}
		}
	}

	for(int i = 0; i < num_threads; i++){
		num_nodes += thInfo[i].nodeCnt;
		num_edges += thInfo[i].edgeCnt;
		thInfo[i].edgeCnt = 0;
		thInfo[i].nodeCnt = 0;
	}
}

template<typename T>
int64_t darhh<T>::in_degree(NodeID n) {
	return vArray[n].inDegree;
}

template<typename T>
int64_t darhh<T>::out_degree(NodeID n) {
	return vArray[n].outDegree;
}

template<typename T>
std::string darhh<T>::to_string() const {
	std::ostringstream os;
	os << *this;
	return os.str();
}

template<typename T>
void darhh<T>::print() {
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

#endif
