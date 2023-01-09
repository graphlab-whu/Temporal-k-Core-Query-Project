#include<iostream>
#include<fstream>
#include<cstring>
#include<string>
#include<vector>
#include<cstdio>
#include<unordered_map>
#include<unordered_set>
#include<set>
#include<map>
#include<algorithm>
#include<chrono>
#include<array>
using namespace std;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;
using std::chrono::system_clock;
using std::chrono::high_resolution_clock;

template<typename T>
void hash_combine(size_t& seed, T const& v)
{
	seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash
{
	template<typename T1, typename T2>
	size_t operator()(std::pair<T1, T2>const& p) const
	{
		size_t seed = 0;
		hash_combine(seed, p.first);
		hash_combine(seed, p.second);
		return seed;
	}
};

const int VMAX = 2610000;
const int EMAX = 64000000;
const int TMAX = 239800000;

int vern = 0;
int arcn = 0;
int tmax = 0;

struct arc {
	int src, dst, t;
}arcs[EMAX];

map<int, unordered_set<int>> tarcs;

const int QMAX = 100;
struct Q {
	int l, r, k;
}q[QMAX];
int qcnt = 0;

void loadgraph(const char* name)
{
	ifstream fin(name, ios::in);
#ifdef _DEBUG
	if (fin.is_open() == false) { printf("open graph %s fail\n", name); exit(1); }
#endif
	vector<int> v;
	int tmin = 0x7fffffff;

	string l;
	while (getline(fin, l))
	{
		int uvt[3] = { 0 };
		int p = -1, np = -1;
		for (int i = 0; i < 3; ++i)
		{
			p = np + 1, np = l.find(' ', np + 1);
			if (np == -1) np = l.size();
			uvt[i] = stoi(l.substr(p, np - p));
		}
		v.push_back(uvt[0]), v.push_back(uvt[1]);
		tmin = min(tmin, uvt[2]);
		arcs[arcn++] = { uvt[0], uvt[1], uvt[2] };
	}

	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	vern = v.size();

	auto get = [&](int k) {
		return (lower_bound(v.begin(), v.end(), k) - v.begin()) + 1;
	};

	for (int i = 0; i < arcn; ++i)
	{
		arcs[i].t -= tmin;
		arcs[i].src = get(arcs[i].src), arcs[i].dst = get(arcs[i].dst);
		if (arcs[i].src > arcs[i].dst) swap(arcs[i].src, arcs[i].dst);
		tmax = max(tmax, arcs[i].t);
	}

}

void loadtest(const char* name)
{
	ifstream fin(name, ios::in);
#ifdef _DEBUG
	if (fin.is_open() == false) { printf("open test %s fail\n", name); exit(1); }
#endif

	string l;
	while (getline(fin, l))
	{
		int lrk[3] = { 0 };
		int p = -1, np = -1;
		for (int i = 0; i < 3; ++i)
		{
			p = np + 1, np = l.find('\t', np + 1);
			if (np == -1) np = l.size();
			lrk[i] = stoi(l.substr(p, np - p));
		}
		q[qcnt++] = { lrk[0], lrk[1], lrk[2] };
	}
}

void tdump(int l, int r)
{
	for (int i = 0; i < arcn; ++i)
		if(arcs[i].t >= l && arcs[i].t <= r)
	tarcs[arcs[i].t].insert(i);
}

int h[VMAX];
int nxt[EMAX], nbr[EMAX], val[EMAX];
int idx = 0;

void addarc(int src, int dst, int t)
{
	nbr[idx] = dst, val[idx] = t, nxt[idx] = h[src], h[src] = idx++;
}

void buildll(int l, int r)
{
	idx = 0;
	memset(h, -1, sizeof h);
	for (int i = 0; i < arcn; ++i)
	{
		int src = arcs[i].src;
		int dst = arcs[i].dst;
		int t = arcs[i].t;
		if (t >= l && t <= r)
		{
			addarc(src, dst, t);
			addarc(dst, src, t);
		}
	}
}

unordered_map<int, int> Mv;
unordered_map<pair<int, int>, int, pair_hash> Mc;
set<pair<int, int>> Hv;
int dmax = 0;

int *core = new int[VMAX];

bool cAdd(int src, int dst)
{
	if (Mc.count({ src, dst }) == 0)
	{
		Mc[{src, dst}] = 1;
		return true;
	}
	Mc[{src, dst}] ++;
	return false;
}

void vAdd(int v)
{
	if (Mv.count(v) == 0) Mv[v] = 0;
	Hv.erase({ Mv[v], v });
	Mv[v] ++;
	Hv.insert({ Mv[v], v });
}


void initMH(int l, int r)
{
	Mv.clear();
	Mc.clear();
	Hv.clear();
	for (int i = 0; i < arcn; ++i)
	{
		int src = arcs[i].src;
		int dst = arcs[i].dst;
		int t = arcs[i].t;
		if (t < l || t > r) continue;
		if (cAdd(src, dst))
		{
			vAdd(src);
			vAdd(dst);
		}
	}
}

void vUpd(int v)
{
	if (Mv.count(v) == 0)
return;
	int d = Mv[v];
	Hv.erase({ d, v });
	Hv.insert({ d - 1, v });
	Mv[v] --;
}

void coredecomp(int l, int r)
{
	memset(core, 0, sizeof(int)*VMAX);
	initMH(l, r);
	if (Hv.empty()) return;
	dmax = (--Hv.end())->first;
	for (int k = 2; k <= dmax + 1; ++k)
	{
		while (Hv.size() && (Hv.begin()->first < k))
		{
			auto nv = *(Hv.begin());
			Hv.erase(Hv.begin());
			int n = nv.first, v = nv.second;
			Mv.erase(v);
			core[v] = k - 1;
			unordered_set<int> nbrs;
			for (int i = h[v]; ~i; i = nxt[i])
			{
				int u = nbr[i];
				int t = val[i];
				if (t >= l && t <= r)
			nbrs.insert(u);
			}
			for (auto u : nbrs) vUpd(u);
		}
	}
}


const int NMAX = 142000;
const int KMAX = 16;
unordered_map<int, int>* CN = new unordered_map<int, int>[VMAX];
auto PHC = new array<array<map<int, int>, KMAX>, VMAX>;

void initCN(int l, int r)
{
	initMH(l, r);
	for (int i = 1; i <= vern; ++i) CN[i].clear();
	for (auto vn : Mv)
	{
		int v = vn.first;
		unordered_set<int> nbrs;
		for (int i = h[v]; ~i; i = nxt[i])
		{
			int u = nbr[i];
			int t = val[i];
			if (t >= l && t <= r)
				nbrs.insert(u);
		}
		for (auto u : nbrs)
		{
			int c = v < u ? Mc[{v, u}] : Mc[{u, v}];
			if(core[u] >= core[v])
		CN[v][u] = c;
		}
	}
}

void updCN(int ver, int nbr)
{
#ifdef _DEBUG
	if (CN[ver].count(nbr) == 0) printf("updCN nbr missing\n");
#endif
	CN[ver][nbr] --;
	if (CN[ver][nbr] == 0)
CN[ver].erase(nbr);
}

int localcore(int ver, int l, int r)
{
	static int cnt[NMAX];
	memset(cnt, 0, sizeof cnt);
	unordered_set<int> nbrs;
	for (int i = h[ver]; ~i; i = nxt[i])
	{
		int u = nbr[i];
		int t = val[i];
		if (t >= l && t <= r)
	nbrs.insert(u);
	}
	for (auto u : nbrs) cnt[min(core[ver], core[u])]++;
	int c = 0;
	for (int k = core[ver]; k >= 1; --k)
	{
		c += cnt[k];
		if (c >= k) return k;
	}
	return 0;
}

void localCN(int ver, int l, int r)
{
	CN[ver].clear();
	unordered_map<int, int> cn;
	for (int i = h[ver]; ~i; i = nxt[i])
	{
		int u = nbr[i];
		int t = val[i];
		if (t >= l && t <= r)
			if (cn.count(u) == 0) cn[u] = 1;
			else cn[u] ++;
	}
	for (auto un : cn)
	{
		int u = un.first, n = un.second;
		if (core[u] >= core[ver])
	CN[ver][u] = n;
	}
}

int qq[VMAX] = {};
void deledge(int t, int ts, int te)
{
	int hh, tt;
	hh = 0, tt = -1;
	for (int eid:tarcs[t])
	{
		int u = arcs[eid].src;
		int v = arcs[eid].dst;
		if (core[u] <= core[v])
		{
			updCN(u, v);
			if (CN[u].size() < core[u]) qq[++tt] = u;
		}
		if (core[v] <= core[u])
		{
			updCN(v, u);
			if (CN[v].size() < core[v]) qq[++tt] = v;
		}
	}
	while (hh <= tt)
	{
		int u = qq[hh++];
		int oldcore = core[u];
		core[u] = localcore(u, ts, te);
		localCN(u, ts, te);
		for (int k = oldcore; k >= core[u] + 1; --k)
			if ((*PHC)[u][k].empty() || (--(*PHC)[u][k].end())->second < te + 1)
				if(k < KMAX)
			(*PHC)[u][k].insert({ ts, te + 1 });
		unordered_set<int> nbrs;
		for (int i = h[u]; ~i; i = nxt[i])
		{
			int v = nbr[i];
			int t = val[i];
			if (t >= ts && t <= te)
		nbrs.insert(v);
		}
		for(auto v:nbrs)
			if (oldcore >= core[v] && core[u] < core[v])
			{
				CN[v].erase(u);
				if (CN[v].size() < core[v])
			qq[++tt] = v;
			}
	}

}

void updCCN(int t, int ts, int te)
{
	int hh, tt;
	hh = 0, tt = -1;
	for (int eid : tarcs[t])
	{
		int u = arcs[eid].src;
		int v = arcs[eid].dst;
		if (core[u] <= core[v])
		{
			updCN(u, v);
			if (CN[u].size() < core[u]) qq[++tt] = u;
		}
		if (core[v] <= core[u])
		{
			updCN(v, u);
			if (CN[v].size() < core[v]) qq[++tt] = v;
		}
	}
	while (hh <= tt)
	{
		int u = qq[hh++];
		int oldcore = core[u];
		core[u] = localcore(u, ts, te);
		localCN(u, ts, te);
		unordered_set<int> nbrs;
		for (int i = h[u]; ~i; i = nxt[i])
		{
			int v = nbr[i];
			int t = val[i];
			if (t >= ts && t <= te)
				nbrs.insert(v);
		}
		for (auto v : nbrs)
			if (oldcore >= core[v] && core[u] < core[v])
			{
				CN[v].erase(u);
				if (CN[v].size() < core[v])
					qq[++tt] = v;
			}
	}
}

int *_core = new int[VMAX];
unordered_map<int, int>* _CN = new unordered_map<int, int>[VMAX];

void bckCCN()
{
	memcpy(_core, core, sizeof(int)*VMAX);
	for (int i = 1; i <= vern; ++i)
_CN[i] = CN[i];
}

void rstCCN()
{
	memcpy(core, _core, sizeof(int)*VMAX);
	for (int i = 1; i <= vern; ++i)
CN[i] = _CN[i];
}


void buildPHC(int l, int r)
{
	coredecomp(l, r);
	initCN(l, r);
	for (int ts = l; ts <= r; ++ts)
	{
		bckCCN();
		for (int t = r; t >= ts; --t)
	deledge(t, ts, t - 1);
		rstCCN();
		updCCN(ts, ts + 1, r);
	}
}

void cleanPHC()
{
	for (int i = 1; i <= vern; ++i)
		for (int j = 1; j < KMAX; ++j)
	(*PHC)[i][j].clear();
}

void fixPHC(int l, int r)
{
	for(int v = 1; v <= vern; ++ v)
		for (int k = 1; k < KMAX; ++k)
			if ((*PHC)[v][k].empty())
		(*PHC)[v][k].insert({ l, r + 1 });
			else {
				auto p = *(--(*PHC)[v][k].end());
				int ts = p.first;
				int te = p.second;
				coredecomp(ts, te);
				initCN(ts, te);
#ifdef _DEBUG 
				if (core[v] < k) printf("fixPHC:suspicious PHC item\n");
#endif
				int t = ts;
				while (t <= te)
				{
					updCCN(t, t + 1, te);
					t++;
					if (core[v] < k)
				break;
				}
				(*PHC)[v][k].insert({ t, r + 1 });
			}
}

int getCT(int v, int k, int ts)
{
	return (--(*PHC)[v][k].upper_bound(ts))->second;
}


unordered_set<int> V;
unordered_set<int> E;
set<pair<int, int>> Hct;
set<pair<int, int>> He;

void process(int ts, int te);

void baseline(int ql, int qr, int qk)
{
	for (int ts = ql; ts <= qr; ++ts)
	{
		Hct.clear(), He.clear(), V.clear(), E.clear();
		for (int v = 1; v <= vern; ++ v)
		{
			int ct = getCT(v, qk, ts);
			Hct.insert({ ct, v });
		}
		for (int i = 0; i < arcn; ++ i)
			if(arcs[i].t >= ts && arcs[i].t <= qr)
		He.insert({ arcs[i].t, i });

		for (int te = ts; te <= qr; ++te)
		{
			while (Hct.size() && Hct.begin()->first <= te)
			{
				V.insert(Hct.begin()->second);
				Hct.erase(Hct.begin());
			}
			vector<pair<int,int>> tempE;
			while (He.size() && He.begin()->first <= te)
			{
				tempE.push_back(*(He.begin()));
				He.erase(He.begin());
			}
			for (auto p : tempE)
			{
				int t = p.first;
				int id = p.second;
				int src = arcs[id].src, dst = arcs[id].dst;
				if (V.count(src) && V.count(dst))
			E.insert(id);
				else
			He.insert(p);
			}
		}
	}

}



void process(int ts, int te)
{
	printf("[%d,%d]\n", ts, te);
	for (auto v : V) printf("%d ", v);
	printf("\n");
	for (auto e : E) printf("(%d,%d,%d)", arcs[e].src, arcs[e].dst, arcs[e].t);
	printf("\n");
}

void pinfo(int qid, const char* graphname, long long clapse)
{
	printf("Graph Name:%s\n", graphname);
	printf("Platform:Omen Laptop 15\n");
	printf("Query:%d\t%d\t%d\n", q[qid].l, q[qid].r, q[qid].k);
	printf("Clapse(nanoseconds):%lld\n", clapse);
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("./baseline [graph] [query]\n");
		return 1;
	}
	loadgraph(argv[1]);
	loadtest(argv[2]);

	for (int i = 0; i < qcnt; ++i)
	{
		int l = q[i].l, r = q[i].r, k = q[i].k;
		buildll(l, r);
		tdump(l, r);
		buildPHC(l, r);
		fixPHC(l, r);
		auto t0 = system_clock::now();
		baseline(l, r, k);
		auto t1 = system_clock::now();
		auto clapse_ns = duration_cast<nanoseconds>(t1 - t0);
		pinfo(i, argv[1], clapse_ns.count());
		cleanPHC();
	}

	system("pause");
	return 0;
}