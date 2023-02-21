
#include "definition.h"


class DYGraph_Base{
public:
    long t_s,t_e_s;
    vector<bool> valid_check; // cmpnt check where eccTS valid if true
    vector<vector<pair<Vertex,Vertex>>> eccTS; // index i, t_e == t_s+i
    map<Vertex,set<Vertex>> curV;   // adjacency
    vector<Edge> curE;  // if u_1,v_1,t_1, u_2,v_2,t_2, when adjacency pair equals, t_1 != t_2
    queue<Edge> Q_res;
    explicit DYGraph_Base(const long & t_start,vector<Edge>& intervalEdges):t_s(t_start){
        curV={};curE={};eccTS={};valid_check={};t_e_s=0;
        this->init_E(intervalEdges);
    }
    explicit DYGraph_Base(const long & t_start,const long & t_end,vector<Edge>& intervalEdges):t_s(t_start),t_e_s(t_end){
        curV={};curE={};eccTS={};valid_check={};
        this->init_E(intervalEdges);
    }

    [[nodiscard]]consteval long get_curT_e() const{ // check current t_e once t_s could initially be an interval
        if(this->t_e_s){
            return t_s + (long)this->eccTS.size() + t_e_s;
        }
        return  t_s + (long)this->eccTS.size();
    }

    void add_E(vector<Edge>& newE){
        if ( curE.empty()||curV.empty()||eccTS.empty() ){
            printf("Need to initE firstly");
            return;
        }
        for(auto& e:newE){
            Q_res.push(e);
        }

        printf("current vertices:%lu \n",curV.size());
        vector<pair<Vertex,Vertex>> curEcc = this->eccTS[this->eccTS.size()-1];

        int s=0;
        while(s != Q_res.size()){
            s = (int)Q_res.size();
            printf("add_E: current s: %d \n",s);
            for (int i=0;i<s;++i){
                auto edge_cur = Q_res.front();
                Q_res.pop();
                Vertex u = get<0>(edge_cur);
                Vertex v = get<1>(edge_cur);
                printf("add_E: edge %lu %lu: ",u,v);
                if(curV.count(u)+curV.count(v)==2){

                    if (curV[u].count(v)&&curV[v].count(u)){
                        curE.emplace_back(edge_cur);
                        printf("existed\n");
                    } else{
                        curEcc = insertE(edge_cur,curEcc);
                        curE.emplace_back(edge_cur);
                        printf("2\n");
                    }
                } else if (curV.count(u)+curV.count(v)==0){
                    Q_res.push(edge_cur);
                    printf("0\n");
                } else if (curV.count(u)+curV.count(v)==1){
                    printf("1\n");
                    if (curV.count(u)){
                        auto& lastEcc = curEcc;
                        Vertex ecc_tmp=-1;
                        for(auto& p:lastEcc){
                            if(p.first==u){
                                ecc_tmp=p.second;
                            }
                        }
                        // ECC eccV(v,ecc_tmp+1);
                        curEcc.emplace_back(make_pair(v,ecc_tmp+1));
                        curE.emplace_back(edge_cur);
                        curV[u].insert(v);
                        curV[v].insert(u);
                        //BFS
                        bfs_MDist(curEcc,u,ecc_tmp);
                    } else{
                        auto& lastEcc = curEcc;
                        Vertex ecc_tmp=-1;
                        for(auto& p:lastEcc){
                            if(p.first==v){
                                ecc_tmp=p.second;
                            }
                        }
                        curEcc.emplace_back(make_pair(u,ecc_tmp+1));
                        curE.emplace_back(edge_cur);
                        curV[u].insert(v);
                        curV[v].insert(u);
                        //BFS
                        bfs_MDist(curEcc,v,ecc_tmp);
                    }
                } else{
                    printf("add edge_cur error");
                }
            }
        }
        this->eccTS.emplace_back(curEcc);
        if (Q_res.empty()){
            valid_check.emplace_back(true);
        } else{
            valid_check.emplace_back(false);
        }

    }

private:
    ~DYGraph_Base()=default;

    /**BFS algorithm below**/

    void bfs_MDist(vector<pair<Vertex,Vertex>>& curEcc, Vertex& src, Vertex& MDist){
        queue<Vertex> Q;
        set<Vertex> S;
        S.insert(src);
        Q.push(src);
        Vertex MD=0;
        Vertex dist = 0;
        while(!Q.empty()){
            ++dist;
            auto s = Q.size();
            for (int i = 0; i < s; ++i) {
                auto cur = Q.front();
                Q.pop();
                for(auto x:curV.at(cur)){
                    if (!S.count(x)){
                        S.insert(x);
                        Q.push(x);
                        // 如果 x的ecc 等于 现在的距离 则x的ecc加一
                        MD = dist;
                        for(auto& p:curEcc){
                            if(p.first==x){
                                if(p.second==dist){
                                    p.second++;
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    vector<pair<Vertex,Vertex>> bfs_curG(const Vertex& src){
        queue<Vertex> Q;
        set<Vertex> S;
        S.insert(src);
        Q.push(src);
        vector<pair<Vertex,Vertex>> src_dist;
        Vertex dist = 0;
        Vertex MD=0;
        src_dist.emplace_back(make_pair(src,MD));
        while(!Q.empty()){
            dist++;
            auto s = Q.size();
            for (int i = 0; i < s; ++i) {
                auto cur = Q.front();
                Q.pop();
                for(auto x:curV.at(cur)){
                    if (!S.count(x)){
                        S.insert(x);
                        Q.push(x);
                        src_dist.emplace_back(make_pair(x,dist));
                        MD = dist;
                    }
                }
            }
        }
        return src_dist;
    }

    static vector<pair<Vertex,Vertex>> bfs_lastG(const Vertex& src,const map<Vertex,set<Vertex>>& lastV){
        vector<pair<Vertex,Vertex>> src_dist;
        queue<Vertex> Q;
        set<Vertex> S;
        S.insert(src);
        Q.push(src);
        Vertex dist = 0;
        Vertex MD=0;
        src_dist.emplace_back(make_pair(src,MD));
        while(!Q.empty()){
            dist++;
            auto s = Q.size();
            for (int i = 0; i < s; ++i) {
                auto cur = Q.front();
                Q.pop();
                for(auto x:lastV.at(cur)){
                    if (!S.count(x)){
                        S.insert(x);
                        Q.push(x);
                        src_dist.emplace_back(make_pair(x,dist));
                        MD = dist;
                    }
                }
            }
        }
        return src_dist;
    }

    void BFS_LS(vector<ECC>& ecc_arr, const Vertex& src, const map<Vertex,Vertex>& V_id, const map<Vertex,Vertex>& ecc_loc){
        queue<Vertex> Q;
        set<Vertex> S;
        Q.push(src);
        S.insert(src);
        Vertex dist = 0;
        Vertex MD =0;
        while (!Q.empty()){
            dist++;
            auto s=Q.size();
            for (int i = 0; i < s; ++i) {
                auto cur = Q.front();
                Q.pop();
                for (auto x: curV.at(cur)) {
                    if (!S.count(x)) {
                        S.insert(x);
                        Q.push(x);
                        MD = dist;
                    }
                }
            }
        }
        ECC& src_e = ecc_arr.at(ecc_loc.at(src));
        src_e.lb = MD;
        src_e.ub = MD;
        src_e.eccv = MD;
    }

    Vertex BFS_pecc_a(const Vertex& src,const set<Vertex>& C_a,const map<Vertex,Vertex>& ecc_loc){
        queue<Vertex> Q;
        set<Vertex> S;
        Q.push(src);
        S.insert(src);
        Vertex dist = 0;
        Vertex MD=0;
        set<Vertex> pC_a = C_a;
        while (!Q.empty() || !pC_a.empty()){
            dist++;
            auto s= Q.size();
            for (int i = 0; i < s; ++i) {
                auto cur = Q.front();
                Q.pop();
                pC_a.erase(ecc_loc.at(cur));
                for (auto x: curV.at(cur)) {
                    if (!S.count(x)) {
                        S.insert(x);
                        Q.push(x);
                        MD=dist;
                    }
                }
            }
        }
        return MD;
    }

    static Vertex BFS_pecc_a_(const Vertex& src,const set<Vertex>& C_a,const map<Vertex,Vertex>& ecc_loc,
                              const map<Vertex,set<Vertex>>& lastV){
        queue<Vertex> Q;
        set<Vertex> S;
        Q.push(src);
        S.insert(src);
        Vertex dist = 0;
        Vertex MD=0;
        set<Vertex> pC_a = C_a;
        while (!Q.empty() || !pC_a.empty()){
            dist++;
            auto s = Q.size();
            for (int i = 0; i < s; ++i) {
                auto cur = Q.front();
                Q.pop();
                pC_a.erase(ecc_loc.at(cur));
                for (auto x: lastV.at(cur)) {
                    if (!S.count(x)) {
                        S.insert(x);
                        Q.push(x);
                        MD = dist;
                    }
                }
            }
        }
        return MD;
    }

    /**BFS algorithm ends**/

    vector<pair<Vertex,Vertex>> insertE(Edge& e,const vector<pair<Vertex,Vertex>>& curEcc){
        // DYECC
        //

        map<Vertex,Vertex> V_id;
        map<Vertex,Vertex> ecc_loc;
        Vertex V_repr = 0;
        for(auto& [k,_]:curV){
            V_id[V_repr]=k;
            ecc_loc[k]=V_repr;
            ++V_repr;
        }

        vector<ECC> ecc_arr (curV.size());
        for (Vertex i = 0; i < curV.size(); ++i) {
            ecc_arr[i].lb=0;
            ecc_arr[i].ub=POSINF;
            ecc_arr[i].eccv=-1;
        }

        Vertex u = get<0>(e);
        Vertex v = get<1>(e);
        const auto lastV = curV;
        vector<pair<Vertex,Vertex>> u_dist_ = bfs_lastG(u,lastV);
        vector<pair<Vertex,Vertex>> v_dist_ = bfs_lastG(v,lastV);

        curV.at(u).insert(v);
        curV.at(v).insert(u);
        vector<pair<Vertex,Vertex>> u_dist = bfs_curG(u);
        vector<pair<Vertex,Vertex>> v_dist = bfs_curG(v);
        // U_DIST storing V_id dist pairs from v

        set<Vertex> C_u; //vertex repr instead of id
        set<Vertex> C_v;

        vector<Vertex> u_d(curV.size(),-1);
        for(auto& [x,d]:u_dist){
            u_d[ecc_loc.at(x)]=d;
        }
        for(auto& [x_,d_]:u_dist_){
            if (d_!=u_d[ecc_loc.at(x_)]){
                C_u.insert(ecc_loc[x_]);
            }
        }

        vector<Vertex> v_d(curV.size(),-1);
        for(auto& [x,d]:v_dist){
            v_d[ecc_loc.at(x)]=d;
        }

        for(auto& [x_,d_]:v_dist_){
            if (d_!=v_d[ecc_loc.at(x_)]){
                C_v.insert(ecc_loc[x_]);
            }
        }

        auto ecc_last = curEcc;
        for(auto& [x_,ecc_]:ecc_last){
            if ( !C_u.count(ecc_loc.at(x_))&&!C_v.count(ecc_loc.at(x_)) ){
                ecc_arr.at(ecc_loc.at(x_)).eccv = ecc_;
                ecc_arr.at(ecc_loc.at(x_)).lb = ecc_;
                ecc_arr.at(ecc_loc.at(x_)).ub = ecc_;
            }
            else {
                ecc_arr.at(ecc_loc.at(x_)).eccv = -1;
                ecc_arr.at(ecc_loc.at(x_)).lb = 1;
                ecc_arr.at(ecc_loc.at(x_)).ub = ecc_;   // upperbound is old ecc;
            }
        }

        // determine which is C_a and BFS on C_a with local spread
        auto& C_a = C_v.size()<C_u.size()?C_v:C_u;
        auto& C_b = C_v.size()<C_u.size()?C_u:C_v;
        ecc_arr.at(ecc_loc.at(u)).eccv = u_dist[u_dist.size()-1].second;
        ecc_arr.at(ecc_loc.at(u)).lb = u_dist[u_dist.size()-1].second;
        ecc_arr.at(ecc_loc.at(u)).ub = u_dist[u_dist.size()-1].second;
        ecc_arr.at(ecc_loc.at(v)).eccv = v_dist[v_dist.size()-1].second;
        ecc_arr.at(ecc_loc.at(v)).lb = v_dist[v_dist.size()-1].second;
        ecc_arr.at(ecc_loc.at(v)).ub = v_dist[v_dist.size()-1].second;

        for (auto& x:C_a){
            BFS_LS(ecc_arr,V_id.at(x),V_id,ecc_loc);
            // possibly apply Local Spread here
        }

        // for each node v∈C_b do
        // Compute the partial eccentricity pecc(v|C_a) on G and pecc (v|C_a) on G' , respectively;

        set<Vertex> C_b_r = C_b;
        for (auto& x:C_b) {
            if (ecc_arr.at(x).eccv!=-1) continue;
            // BFS_pecc_a should be avoided by perform BFS_LS on C_a
            auto p_ecc_x = BFS_pecc_a(V_id.at(x), C_a,ecc_loc);
            // another pecc_a on last adjacency
            auto p_ecc_x_ = BFS_pecc_a_(V_id.at(x), C_a,ecc_loc,lastV);
            /*TD*/
            // using lemma to minimize C_b
            if (p_ecc_x_ < ecc_arr.at(x).ub){
                ecc_arr.at(x).lb = ecc_arr.at(x).ub;
                ecc_arr.at(x).eccv = ecc_arr.at(x).ub;
                C_b_r.erase(x);
                continue;
            }
            if (p_ecc_x_ == ecc_arr.at(x).ub == p_ecc_x){
                ecc_arr.at(x).lb = ecc_arr.at(x).ub;
                ecc_arr.at(x).eccv = ecc_arr.at(x).ub;
                C_b_r.erase(x);
                continue;
            }
            const auto neighbors = curV.at(V_id.at(x));
            for (auto& p:neighbors){
                auto p_loc = ecc_loc.at(p);
                if (ecc_arr.at(p_loc).eccv!=-1) continue;
                else if (ecc_arr.at(p_loc).eccv > ecc_arr.at(x).ub) {
                    ecc_arr.at(x).lb = ecc_arr.at(x).ub;
                    ecc_arr.at(x).eccv = ecc_arr.at(x).ub;
                    C_b_r.erase(x);
                    break;
                }
            }
        }

        // BFS on remaining C_b
        if (!C_b_r.empty()){
            for (auto& x:C_b){
                BFS_LS(ecc_arr,V_id.at(x),V_id,ecc_loc);
            }
        }

//        if (ecc_arr.size()!=curV.size()){
//            printf("insertE error in set C_a or C_b");
//        }
//        for (auto& x:ecc_arr){
//            if (x.eccv==-1){
//                printf("insertE error in set C_a or C_b");
//            }
//        }

        vector<pair<Vertex,Vertex>> new_ecc;
        for (int i = 0; i < ecc_arr.size(); ++i) {
            new_ecc.emplace_back(make_pair(V_id.at(i),ecc_arr.at(i).eccv));
        }
        /*TD*/
        // just put it into the new vector by pairs
        return new_ecc;
    }

    void init_E(vector<Edge>& newE){
        if ( !( curE.empty()&&curV.empty()&&eccTS.empty() ) ){
            printf("file to initE");
            return;
        }

        vector<pair<Vertex,Vertex>> curEcc = {};
        //this->eccTS[this->eccTS.size()-1];

        for(auto& e:newE){
            Q_res.push(e);
        }

        auto edge_0 = Q_res.front();
        Q_res.pop();
        curE.push_back(edge_0);

        curV[get<0>(edge_0)].insert(get<1>(edge_0));
        curV[get<1>(edge_0)].insert(get<0>(edge_0));
        curEcc.emplace_back(make_pair(get<0>(edge_0),1));
        curEcc.emplace_back(make_pair(get<1>(edge_0),1));
        // adjacency
        // printf("insertE: adjacency preparing finished by pop 1 edge with rest %lu edges in Queue\n",Q_res.size());
        int s=0;
        // cout<< boolalpha<< (s < Q_res.size()) <<endl;
        while(s != Q_res.size()){
            s = (int)Q_res.size();
            printf("insertE: cur Q_res size: %d\n",s);
            for (int i=0;i<s;++i){
                auto edge_cur = Q_res.front();
                Q_res.pop();
                Vertex u = get<0>(edge_cur);
                Vertex v = get<1>(edge_cur);
                // printf("insertE: edge %lu %lu\n",u,v);
                // cout << curV.count(u)+curV.count(v) << endl;
                if(curV.count(u)+curV.count(v)==2){
                    if (curV[u].count(v)&&curV[v].count(u)){
                        curE.emplace_back(edge_cur);
                        // printf("insertE: edge is redundant\n");
                    } else{
                        curEcc = insertE(edge_cur,curEcc);
                        curE.emplace_back(edge_cur);
                        // printf("insertE: edge is inserted\n");
                    }
                } else if (curV.count(u)+curV.count(v)==0){
                    Q_res.push(edge_cur);
                    // printf("insertE: edge is unmanageable\n");
                } else if (curV.count(u)+curV.count(v)==1){
                    if (curV.count(u)){
                        auto& lastEcc = curEcc;
                        Vertex ecc_tmp=-1;
                        for(auto& p:lastEcc){
                            if(p.first==u){
                                ecc_tmp=p.second;
                            }
                        }
                        // ECC eccV(v,ecc_tmp+1);
                        curEcc.emplace_back(make_pair(v,ecc_tmp+1));
                        curE.emplace_back(edge_cur);
                        curV[u].insert(v);
                        curV[v].insert(u);
                        //BFS
                        bfs_MDist(curEcc,u,ecc_tmp);

                        // printf("insertE: one novel vertex\n");
                    } else{
                        auto& lastEcc = curEcc;
                        Vertex ecc_tmp=-1;
                        for(auto& p:lastEcc){
                            if(p.first==v){
                                ecc_tmp=p.second;
                            }
                        }
                        curEcc.emplace_back(make_pair(u,ecc_tmp+1));
                        curE.emplace_back(edge_cur);
                        curV[u].insert(v);
                        curV[v].insert(u);
                        //BFS
                        bfs_MDist(curEcc,v,ecc_tmp);

                        // printf("insertE: one novel vertex\n");
                    }
                } else{
                    printf("add edge_cur error\n");
                }
            }
        }

        this->eccTS.emplace_back(curEcc);
        if (Q_res.empty()){
            valid_check.emplace_back(true);
        } else{
            valid_check.emplace_back(false);
        }

    }

};
