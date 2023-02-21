#include "definition.h"

class DYGraph_Unified {
public:
    explicit DYGraph_Unified(const long &t_end, vector<Edge> &intervalEdges) : t_e(0) {
        curV = {};
        eccTS = {};
        valid_check = {};
        Q_res = {};
        V_id = {};
        ecc_loc = {};
        eccTS.emplace_back(vector<vector<Vertex>>{});   // t == 0
        this->insert_Edges(intervalEdges);
    }

    void update_Edges(vector<Edge> & intervalEdges){
        this->insert_Edges(intervalEdges);
    }

    DYGraph_Unified &operator=(const DYGraph_Unified &obj) = delete;

    // TD: interface for new edges insertion by
    // t_e++
    // and insert_Edges

private:
    ~DYGraph_Unified() = default; // allocate on heap

    long t_e;
    vector<vector<bool>> valid_check; // cmpnt check where eccTS valid if true
    vector<vector<vector<Vertex>>> eccTS; // index i, t_e == t_s+i // t_s in reverse ordering
    queue<Edge> Q_res;
    unordered_map<Vertex, set<Vertex>> curV;   // adjacency
    vector<Vertex> V_id;
    unordered_map<Vertex, Vertex> ecc_loc;
    map<pair<Vertex, Vertex>, long> adj_ti_last;   //latest timestamp between u,v that u<v


    [[nodiscard]]auto BFS_SS(Vertex &src, long time_offset, auto &v_dist_t) -> vector<long> {

        priority_queue<adj_repr, vector<adj_repr>, adj_time_cmp> PQ;
        set < Vertex > S;
        S.insert(src);
        PQ.push({src, 0, t_e});

        // src of v_dist_t on the time_offset
        v_dist_t[ecc_loc.at(src)].emplace_back(make_pair(1, 0));

        Vertex MD = 0;
        long MDTime = t_e;
        map<long, Vertex> t_MD;

        while (!PQ.empty()) {
            auto cur = PQ.top();
            PQ.pop();
            // curV -> adj meanwhile adj_ti_last -> time

            if (cur.time > MDTime) [[unlikely]] puts("MDT error");

            if (cur.time != MDTime) {
                t_MD[MDTime] = MD;
                MD = 0;
                MDTime = cur.time;
            }

            S.insert(cur.v);
            //
            // for adj (v,t): v in curV, time in adj_ti_last
            //
            Vertex adj_u = cur.v;
            auto &adj = curV.at(adj_u);
            for (auto adj_v: adj) {
                // if adj_u>adj_v swap adj_u,adj_v
                long time_late = min(cur.time, adj_ti_last.at(make_pair(min(adj_u, adj_v), max(adj_u, adj_v))));
                adj_repr x = {adj_v, cur.dist + 1, time_late};

                if (time_late <= time_offset) {
                    continue;
                }

                if (!S.count(adj_v)) {   // virgin
                    S.insert(adj_v);
                    // PQ and adj_v 都 emplace_back
                    PQ.push(x);
                    v_dist_t[ecc_loc.at(adj_v)].emplace_back(make_pair(time_late, cur.dist + 1));
                    MD = max(MD, cur.dist + 1);
                } else {
                    // while t equal if dist decrease, replace v_dist_t && emplace PQ
                    // else if t decrease, dist cannot decrease, continue
                    // both t and dist decrease, emplace v_dist_t and PQ
                    auto &visit_adj_v = v_dist_t[ecc_loc.at(adj_v)];
                    auto &visit_v_last = visit_adj_v.at(visit_adj_v.size() - 1);
                    if (visit_v_last.first == time_late) {
                        if (visit_v_last.second <= cur.dist + 1) [[likely]] {
                            continue;
                        } else {
                            // PQ emplace_back && replace v_dist_t
                            visit_v_last.second = cur.dist + 1;
                            PQ.push(x);
                            MD = max(cur.dist + 1, MD);  // reset MD
                        }
                    } else if (visit_v_last.second > cur.dist + 1) {  // t in PQ must be <= that if not == then <
                        // PQ and v_dist_t do emplace_back
                        PQ.push(x);
                        // v_dist_t[ecc_loc.at(adj_v)].emplace_back(make_pair(time_late,cur.dist+1));
                        visit_adj_v.emplace_back(make_pair(time_late, cur.dist + 1));
                        MD = max(MD, cur.dist + 1);
                    } else [[unlikely]] { puts("error in PQ"); }
                }
            }
        }

        vector<Vertex> u_dist(t_e + 1, 0); // DIST
        long index_offset = time_offset+1;
        long last_d = 0;
        for (const auto &[t, d]: t_MD) {
            for (long i = index_offset; i < t; ++i) {
                u_dist[i] = last_d;
            }
            index_offset = t;
            last_d = d;
        }
        for (long i = index_offset; i < u_dist.size(); ++i) {
            u_dist[i] = last_d;
        }
        return u_dist;
        // check cmpt
    }

    void update_DY(Vertex u, Vertex v, long time_offset, vector<vector<Vertex>> &curEcc) {
        auto nV = static_cast<Vertex>(V_id.size());

        // BFS info can be stored, if no such info, do searching
        vector<vector<pair<long, Vertex>>> verti_dist_t_u(nV); // each v_r-> vector (t, dist)
        vector<vector<pair<long, Vertex>>> verti_dist_t_v(nV); // each v_r-> vector (t, dist)
        auto tMD_u = BFS_SS(u, time_offset, verti_dist_t_u);
        auto tMD_v = BFS_SS(v, time_offset, verti_dist_t_v);  // valid after offset

        /** combine eccTS with tMD before time offset **/
        auto ecc_u_last = curEcc.at(ecc_loc.at(u));
        auto ecc_v_last = curEcc.at(ecc_loc.at(v));

        for (int i = 0; i < time_offset; ++i) {
            tMD_u[i]= ecc_u_last[i];
            tMD_v[i]= ecc_v_last[i];
        }

        // bfs results
        unordered_map<Vertex, long> dirty_a;
        unordered_map<Vertex, long> dirty_b;
        unordered_map<Vertex, long> dirty_ab;
        Vertex j = 0;   // pair index

        for (int v_r = 0; v_r < nV; ++v_r) {

            vector<Vertex> dist2u_(t_e + 1, 0), dist2v_(t_e + 1, 0);
            Vertex x = V_id.at(v_r);

            // inherent
            // pair转换成vector，然后进行+1比较 存留比较后来的结果作为u，v的新vector of ecc
            for (long i = t_e; i > 0; --i) {
                if (i <= time_offset) break;
                if (i >= verti_dist_t_u[v_r][j].first) {
                    dist2u_[i] = verti_dist_t_u[v_r][j].second;
                } else {
                    j++;
                    if (i < verti_dist_t_u[v_r][j].first) [[unlikely]] {
                        printf("error dist2u_");
                    }
                }
                if (i >= verti_dist_t_v[v_r][j].first) {
                    dist2v_[i] = verti_dist_t_v[v_r][j].second;
                } else {
                    j++;
                    dist2v_[i] = verti_dist_t_v[v_r][j].second;
                    if (i < verti_dist_t_u[v_r][j].first) [[unlikely]] {
                        printf("error dist2v_");
                    }
                }

                // a' = min(b+1,a)
                // b' = min(a+1,b)
                // if a' != a  C_a
                auto dist2v = min(dist2v_[i], dist2u_[i] + 1);
                auto dist2u = min(dist2u_[i], dist2v_[i] + 1);
                // TD: 不同的t设置ca cb 或者整体c_a或c_b 所有边的dirty一起bfs
                // 逐个增加时候 详尽ab集合划分，即 a b ab 以及 各点不同的t_s offset时间的
                if (dist2v != dist2v_[i]) {
                    // TD: 找到最后一个不同的时间 然后 BFS
                    // 不考虑 u,t_i that t_i >= time_offset
                    // time offset
                    dirty_a[x] = i;
                    if (dirty_b.count(x)) {
                        dirty_ab[x] = i;
                    }
                }
                if (dist2u != dist2u_[i]) {
                    dirty_b[x] = i;
                    if (dirty_a.count(x)) {
                        dirty_ab[x] = i;
                    }
                }
                // update dist after insert
                dist2v_[i] = dist2v;
                dist2u_[i] = dist2u;
                // v possible in different set C-a or C-b vary t
            }
        }


        printf("dirty_a: %lu ", dirty_a.size());
        printf("dirty_b: %lu ", dirty_b.size());
        printf("dirty_ab: %lu ", dirty_ab.size());
        // bfs a+ab to reduce b-ab
        // bfs b

        curV[u].insert(v);
        curV[v].insert(u);
        adj_ti_last.at(make_pair(u, v)) = t_e;
    }

    void v_insert_DY(Edge &edge_cur) {

    }

    void insert_Edges(vector<Edge> &newE) {
        vector<vector<Vertex>> curEcc; // v_r -> t_x-> ecc
        t_e++;
        if (!curV.empty()) {
            curEcc = eccTS[eccTS.size() - 1];
            for (auto &e: newE) {
                Q_res.push(e);
            }
        } else [[unlikely]] {
            curEcc = {};
            //this->eccTS[this->eccTS.size()-1];

            for (auto &e: newE) {
                Q_res.push(e);
            }

            auto edge_0 = Q_res.front();
            Q_res.pop();

            curV[get<0>(edge_0)].insert(get<1>(edge_0));
            curV[get<1>(edge_0)].insert(get<0>(edge_0));

            ecc_loc[get<0>(edge_0)] = 0;
            ecc_loc[get<1>(edge_0)] = 1;

            V_id.emplace_back(get<0>(edge_0));
            V_id.emplace_back(get<1>(edge_0));
            auto x = get<0>(edge_0) < get<1>(edge_0) ?
                     make_pair(get<0>(edge_0), get<1>(edge_0)) : make_pair(get<1>(edge_0), get<0>(edge_0));
            adj_ti_last.at(x) = 1;

            curEcc.emplace_back(vector<Vertex>{0, 1}); // v0 -> ts = te = 1
            curEcc.emplace_back(vector<Vertex>{0, 1});
        }

        // TD
        // 如果 G_[1, t_e] 没见过 u,v 那么到任意起始时间都不能一个cc
        // 如果 G_[1, t_e] 只见过 u,v 的一个 那么 t_e +1 服从 t_e 的cc

        // adjacency
        // printf("insertE: adjacency preparing finished by pop 1 edge with rest %lu edges in Queue\n",Q_res.size());
        long s = 0;
        while (s != Q_res.size()) {
            s = (long) Q_res.size();
            // printf("insertE: cur Q_res size: %d\n",s);
            for (long i = 0; i < s; ++i) {
                auto edge_cur = Q_res.front();
                Q_res.pop();
                Vertex u = get<0>(edge_cur);
                Vertex v = get<1>(edge_cur);

                // printf("insertE: edge %lu %lu\n",u,v);
                // cout << curV.count(u)+curV.count(v) << endl;

                switch (curV.count(u) + curV.count(v)) {
                    case 2: {
                        if (u > v) swap(u, v);
                        long time_offset = 0;
                        if (curV.at(u).count(v)) {
                            time_offset = adj_ti_last.at(make_pair(u, v));
                        }
                        adj_ti_last.at(make_pair(u, v)) = time_offset;
                        update_DY(u, v, time_offset, curEcc); // the updating only happens within (time_offset,t_e]
                        break;
                    }
                    case 1: {
                        if (curV.count(v)) swap(u, v); // v is newly inserted
                        auto nv = static_cast<Vertex> (V_id.size());
                        ecc_loc[v] = nv;
                        V_id.emplace_back(v);
                        // curV[u].insert(v);
                        // curV[v].insert(u);
                        auto x = u < v ? make_pair(u, v) : make_pair(v, u);
                        adj_ti_last.at(x) = t_e;
                        //更改cur
                        v_insert_DY(edge_cur);
                        break;
                    }
                    case 0: {
                        Q_res.push(edge_cur);
                        break;
                    }
                }
            }
        }

        eccTS.emplace_back(curEcc);

        // append a valid_check array
//        if(!Q_res.empty()){
//            return;
//        }

    }


};