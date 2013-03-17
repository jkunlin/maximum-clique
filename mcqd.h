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

#ifndef MCQD_H
#define MCQD_H

#include <iostream>
#include <algorithm>
#include <string>
#include <assert.h>

class Maxclique {
	const bool* const* e;
	int pk, level;
	const float Tlimit;
	class Vertices {
		class Vertex {
			int i, d;
			public:
			void set_i(const int ii)  { i = ii; }
			int get_i() const { return i; }
			void set_degree(int dd) { d = dd; }
			int get_degree() const { return d; }
		};
		Vertex *v;
		int sz;
		static bool desc_degree(const Vertex vi, const Vertex vj) { return (vi.get_degree() > vj.get_degree()); }
		public:
#ifndef NDEBUG
		void dbg_v(const std::string msg="") const {
			std::cout << msg << " Vertices: [";
			for (int i=0; i < sz; i++) 
				std::cout << "(" << v[i].get_i() << "," << v[i].get_degree() << ") ";
			std::cout << "]" << std::endl;
		}
#endif
		Vertices(int size) : sz(0) { v = new Vertex[size]; }
		~Vertices () {}
		void dispose() { if (v) delete [] v; }
		void sort() { std::sort(v, v+sz, desc_degree); }
		void init_colors();

//		void MCR_sort(Maxclique &m);
//		void MCR_init_colors();

		void set_degrees(Maxclique&);
		int size() const { return sz; }
		void push(const int ii) { v[sz++].set_i(ii); };
		void pop() { sz--; };
		void clean() { sz = 0; }
		Vertex& at(const int ii) const { return v[ii]; };
		Vertex& end() const { return v[sz - 1]; };
	};
	class ColorClass {
		int *i;
		int sz;
		public:
#ifndef NDEBUG
		void dbg_i(const std::string msg="") const {
			std::cout << msg << " Class: [";
			for (int ii=0; ii < sz; ii++) 
				std::cout << i[ii] << " ";
			std::cout << "]" << std::endl;
		}
#endif
		ColorClass() : sz(0), i(0) {}
		ColorClass(const int sz) : sz(0), i(0) { i = new int[sz]; }
		~ColorClass() { if (i) delete [] i;
		}
		void init(const int sz) { i = new int[sz]; rewind(); }
		void push(const int ii) { i[sz++] = ii; };
		void pop() { sz--; };
		void rewind() { sz = 0; };
		int size() const { return sz; }
		int& at(const int ii) const { return i[ii]; }
		ColorClass& operator=(const ColorClass& dh) {
			for (int j = 0; j < dh.sz; j++) i[j] = dh.i[j];
			sz = dh.sz;
			return *this;
		}
	};
	Vertices V;
	ColorClass *C, QMAX, Q;
	class StepCount {
		int i1, i2;
		public:
		StepCount() : i1(0), i2(0) {}
		void set_i1(const int ii)  { i1 = ii; }
		int get_i1() const { return i1; }
		void set_i2(const int ii)  { i2 = ii; }
		int get_i2() const { return i2; }
		void inc_i1()  { i1++; }
	};
	StepCount *S;
	bool connection(const int i, const int j) const { return e[i][j]; }
	bool conflict(const int, const ColorClass&);
	void cut(const Vertices&, Vertices&);
	void color_sort(Vertices&);
	void expand(Vertices &);
	void expand_dyn(Vertices &);
	void _mcq(int*&, int&, bool);
	void degree_sort(Vertices &R) { R.set_degrees(*this); R.sort(); }
	void MCR_sort(Vertices &);
	public:
#ifndef NDEBUG
	void dbg_C() const {
		for (int i=0; i < V.size(); i++) {
			std::cout << "C["<< i << "] : ";
			C[i].dbg_i();
		}
	}
	void dbg_conn() const {
		for (int i=0; i < V.size(); i++) {
			for (int j=0; j < V.size(); j++) {
				std::cout <<e[i][j];
			}
			std::cout<< std::endl;
		}
	}
#endif
	Maxclique(const bool* const*, const int, const float=0.025);
	int steps() const { return pk; }
	void mcq(int* &maxclique, int &sz) { _mcq(maxclique, sz, false); }
	void mcqdyn(int* &maxclique, int &sz) { _mcq(maxclique, sz, true); }
	~Maxclique() {
		if (C) delete [] C;
		if (S) delete [] S;
		V.dispose();
	};
};

Maxclique::Maxclique (const bool* const* conn, const int sz, const float tt) : pk(0), level(1), Tlimit(tt), V(sz), Q(sz), QMAX(sz) {
	assert(conn!=0 && sz>0);
	for (int i=0; i < sz; i++) V.push(i);
	e = conn;
	C = new ColorClass[sz];
	for (int i=0; i < sz; i++) C[i].init(sz);
	S = new StepCount[sz];
}

void Maxclique::_mcq(int* &maxclique, int &sz, bool dyn) { 
//	V.sort();
//	V.init_colors();
	MCR_sort(V);
	color_sort(V);

	if (dyn) {
		/*	
			for (int i=0; i < V.size() + 1; i++) {
			S[i].set_i1(0);
			S[i].set_i2(0);
			}
		 */
		expand_dyn(V);
	}
	else
		expand(V);
	maxclique = new int[QMAX.size()]; 
	for (int i=0; i<QMAX.size(); i++) { 
		maxclique[i] = QMAX.at(i);
	}
	sz = QMAX.size();
}
/*********init is no a good idea******************
void Maxclique::Vertices::init_colors() { 
	const int max_degree = v[0].get_degree();
	for (int i = 0; i < max_degree; i++)
		v[i].set_degree(i + 1);
	for (int i = max_degree; i < sz; i++)
		v[i].set_degree(max_degree + 1);
}
 */
/*
void Maxclique::Vertices::MCR_sort(Maxclique &m) {
	set_degrees(m);
	int min_degree = sz;
	int min_index = sz - 1;
	Vertices R(sz);
	for(int i=0; i < sz; i++) {
		//find the vertex with min degree
		R.clean();
		min_degree = v[0].get_degree();
		min_index = 0;
		R.push(v[0].get_i());

		for(int j=1; j < sz - i; j++) {
			if(v[j].get_degree() == min_degree) {
				R.push(v[j].get_i());
			}
			else if(v[j].get_degree() < min_degree) {
				min_degree = v[j].get_degree();
				min_index = j;
				R.clean();
				R.push(v[j].get_i());
			}
		}
		//regular subgraph
		if(R.size() == sz - i) 
			break;
		std::swap(v[min_index], v[sz - i - 1]);

		//update degrees
		for(int j = 0; j < sz - i - 1; j++) {
			if(m.connection(v[sz - i - 1].get_i(), v[j].get_i()))
				v[j].set_degree(v[j].get_degree() - 1);
		}
	}
}
*/
void Maxclique::MCR_sort(Vertices &V) {
	V.set_degrees(*this);
	int sz = V.size();
	int min_degree = sz;
	int min_index = sz - 1;
	Vertices R(sz);
	for(int i=0; i < sz; i++) {
		//find the vertex with min degree
		R.clean();
		min_degree = V.at(0).get_degree();
		min_index = 0;
		R.push(V.at(0).get_i());

		for(int j=1; j < sz - i; j++) {
			if(V.at(j).get_degree() == min_degree) {
				R.push(V.at(j).get_i());
			}
			else if(V.at(j).get_degree() < min_degree) {
				min_degree = V.at(j).get_degree();
				min_index = j;
				R.clean();
				R.push(V.at(j).get_i());
			}
		}
		//regular subgraph
		if(R.size() == sz - i) 
			break;
		std::swap(V.at(min_index), V.at(sz - i - 1));

		//update degrees
		for(int j = 0; j < sz - i - 1; j++) {
			if(connection(V.at(sz - i - 1).get_i(), V.at(j).get_i()))
				V.at(j).set_degree(V.at(j).get_degree() - 1);
		}
	}
}
//use inc_degree maybe quicker
void Maxclique::Vertices::set_degrees(Maxclique &m) { 
	int d;
	for (int i=0; i < sz; i++) {
		d = 0;
		for (int j=0; j < sz; j++)
			if (m.connection(v[i].get_i(), v[j].get_i())) d++;
		v[i].set_degree(d);
	}
}

//judge if vertex pi is in conflict with ColorClass A
bool Maxclique::conflict(const int pi, const ColorClass &A) {
	for (int i = 0; i < A.size(); i++)
		if (connection(pi, A.at(i)))
			return true;
	return false;
}

void Maxclique::cut(const Vertices &A, Vertices &B) {
	for (int i = 0; i < A.size() - 1; i++) {
		if (connection(A.end().get_i(), A.at(i).get_i()))
			B.push(A.at(i).get_i());
	}
}

void Maxclique::color_sort(Vertices &R) {
	int j = 0;
	int maxno = -1;
	int min_k = QMAX.size() - Q.size() - 1;
	C[0].rewind();
	int k = 0;
	for (int i=0; i < R.size(); i++) {
		int pi = R.at(i).get_i();
		k = 0;
		while (conflict(pi, C[k]))
			k++;
		if (k > maxno) {
			maxno = k;
			C[maxno + 1].rewind();
		}
		C[k].push(pi);
		if (k < min_k) {
			R.at(j++).set_i(pi);
		}
	}
	if (j > 0) R.at(j-1).set_degree(0);
	if (min_k <= 0) min_k = 0;
	for (k = min_k; k <= maxno; k++)
		for (int i = 0; i < C[k].size(); i++) {
			R.at(j).set_i(C[k].at(i));
			R.at(j++).set_degree(k);
		}
}

void Maxclique::expand(Vertices &R) {
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Rp(R.size());
			cut(R, Rp);
			if (Rp.size()) {
				color_sort(Rp);
				pk++;
				expand(Rp);
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
			Rp.dispose();
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	}
}

void Maxclique::expand_dyn(Vertices &R) {
	S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
	S[level].set_i2(S[level - 1].get_i1());
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Rp(R.size());
			cut(R, Rp);
			if (Rp.size()) {
				if ((float)S[level].get_i1()/++pk < Tlimit) {
					degree_sort(Rp);
				}
				color_sort(Rp);
				S[level].inc_i1();
				level++;
				expand_dyn(Rp);
				level--;
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
			Rp.dispose();
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	}
}

#endif
