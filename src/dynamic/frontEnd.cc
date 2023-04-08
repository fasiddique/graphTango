#include <unistd.h>
#include <fstream>
#include <cstring>
#include <mutex>
#include <thread>
#include <pthread.h>
#include <sstream>

#include "topAlg.h"
#include "topDataStruc.h"
#include "fileReader.h"
#include "parser.h"
#include "../common/timer.h"

using namespace std;
/* Main thread that launches everything else */

ofstream algF("Alg.csv");

int main(int argc, char *argv[]) {
	cmd_args opts = parse(argc, argv);
	ifstream file(opts.filename);
	if (!file.is_open()) {
		cout << "Couldn't open file " << opts.filename << endl;
		exit(-1);
	}

	int batch_id = 0;
	NodeID lastAssignedNodeID = -1;
	MapTable VMAP;
	EdgeList el;
	el.reserve(opts.batch_size);

	Timer t;
	dataStruc *ds = createDataStruc(opts.type, opts.weighted, opts.directed,
			opts.num_nodes, opts.num_threads);
	Algorithm alg(opts.algorithm, ds, opts.type);

	ofstream updF("Update.csv");

	while (!file.eof()) {
		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted,
				VMAP, lastAssignedNodeID);

		t.Start();
		ds->update(el);
		t.Stop();

		updF << t.Seconds() << endl;
		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes
				<< ", Edges " << ds->num_edges << endl;

		batch_id++;
	}
	file.close();

#ifndef ENABLE_PROFILING
	alg.performAlg();
#endif

	if (opts.enDeleteEdges) {
		ofstream updD("delUpdate.csv");
		stringstream ss;
		ss << opts.filename << ".del";
		file.open(ss.str());
		if (!file.is_open()) {
			cout << "Couldn't open file for delete" << ss.str() << endl;
			exit(-1);
		}

		batch_id = 0;
		el.clear();
		while (!file.eof()) {
			readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted,
					VMAP, lastAssignedNodeID);

			t.Start();
			ds->update(el);
			t.Stop();

			updD << t.Seconds() << endl;
			cout << "Deleted Batch " << batch_id << ": Nodes " << ds->num_nodes
					<< ", Edges " << ds->num_edges << endl;

			batch_id++;
		}
		file.close();
		updD.close();
	}
	ds->print();
	if (ds) {
		delete ds;
	}

#ifdef CALC_EDGE_TOUCHED
	cout << "EDGES TOUCHED: " << g_edge_touched << endl;
#endif
}

