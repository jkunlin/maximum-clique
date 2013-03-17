/*
   Copyright 2007-2012 Janez Konc 

   If you use this program, please cite: 
   Janez Konc and Dusanka Janezic. An improved branch and bound algorithm for the 
   maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.

   More information at: http://www.sicmm.org/~konc

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include "mcqd.h"

using namespace std;

void read_file(string name, bool ** &conn, int &size) {
	ifstream f(name.c_str());
	char buffer[80];
	assert(f.is_open());
	set<int> v;
	multimap<int, int> e;

	while(!f.eof()) {
		f.getline (buffer, 70);
		if(buffer[0] == 'e') {
			int vi, vj;
			sscanf(buffer, "%*c %d %d", &vi, &vj);
			v.insert(vi);
			v.insert(vj);
			e.insert(make_pair(vi, vj));
		}
	}
	size = v.size();

	conn = new bool *[size];
	for(size_t i = 0; i < size; i++) {
		conn[i] = new bool[size];
		memset(conn[i], 0, size * sizeof(bool));
	}
	for(multimap<int, int>::iterator it = e.begin(); it != e.end(); it++) {
		conn[it->first - 1][it->second -1] = true;
		conn[it->second - 1][it->first - 1] = true;
	}

	cout << "|E| = "<< e.size() << " |V| = " << v.size() << " p = " << (double) e.size() / (v.size() * (v.size() - 1) /2) <<endl;
	f.close();
}

int main(int argc, char *argv[]) {
	assert(argc == 2);
	cout << "args = " << argv[1] << endl;
	bool **conn;
	int size;
	read_file(argv[1], conn, size);

	clock_t start1, start2;
	int *qmax;
	int qsize;

	cout << "---------- Example 1: run max clique with improved coloring ----------------"<<endl;
	start1 = time(NULL);
	start2 = clock();
	Maxclique m(conn, size);
	m.mcq(qmax, qsize);  // run max clique with improved coloring
	cout << "Size = " << qsize << endl;
	cout << "Number of steps = " << m.steps() << endl;
	cout << "Time = " << difftime(time(NULL), start1) << endl;
	cout << "Time (precise) = " << ((double) (clock() - start2)) / CLOCKS_PER_SEC << endl << endl;
	delete [] qmax;

	cout << "---------- Example 2: run max clique with improved coloring and dynamic sorting of vertices ----------------"<<endl;
	start1 = time(NULL);
	start2 = clock();
	Maxclique md(conn, size, 0.025);  //(3rd parameter is optional - default is 0.025 - this heuristics parameter enables you to use dynamic resorting of vertices (time expensive)
	// on the part of the search tree that is close to the root - in this case, approximately 2.5% of the search tree -
	// you can probably find a more optimal value for your graphs
	md.mcqdyn(qmax, qsize);  // run max clique with improved coloring and dynamic sorting of vertices 
	cout << "Size = " << qsize << endl;
	cout << "Number of steps = " << md.steps() << endl;
	cout << "Time = " << difftime(time(NULL), start1) << endl;
	cout << "Time (precise) = " << ((double) (clock() - start2)) / CLOCKS_PER_SEC << endl << endl;

	delete [] qmax;
	for (int i=0;i<size;i++)
		delete [] conn[i];
	delete [] conn;
	return 0;
}
