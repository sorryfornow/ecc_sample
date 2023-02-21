#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <queue>
#include <map>
#include "DYGraph_Base.hpp"
#include "DYGraph_Unified.hpp"

using namespace std;

/**declaration set in Unified file for code below**/

vector<Edge> edgesReceive(const int& paraNum, const char* dataPath){
    vector<Edge> edges;
    string curLine;
    int lineNum = 0;    // number of edges

    ifstream infile;
    infile.open(dataPath,ios::in);
    if (!infile.is_open()){
        throw fileOpenException(dataPath);
    }
    const char *pDelimiter = " ";
    char* pToken = nullptr;
    char* pSave = nullptr;
    while(getline(infile,curLine)){
        char * strs = new char[curLine.length() + 1];
        assert(strs);
        strcpy(strs, curLine.c_str());
        /* split the line into paraNum parts */
        pToken = strtok_r(strs,pDelimiter,&pSave);
        int pNum = 0;

        vector<Vertex> edge(paraNum);
        while(pToken && pNum<paraNum){
            edge[pNum] = stoi(pToken);
            pToken = strtok_r(nullptr, pDelimiter, &pSave);
            ++pNum;
        }
        Edge e(edge[0],edge[1],edge[2]);

        edges.emplace_back(e);
        ++lineNum;
        delete[] strs; // may cause memory leak somehow it interrupted
    }
    infile.close();
    return edges;
}

void baseline_algorithm(const vector<Edge>& edges){
    long t_max = get<2>(edges[edges.size()-1]); // last timestamp t_max
    vector<long> time_offset(t_max+1,0);   // offset array
    time_offset[t_max]=(long)edges.size();  // last timestamp
    // Edge tri;
    // tuple_element<2, decltype(tri)>::type time_stmp;
    long t_log = 1;
    long i_ = 0;
    for(auto& edge:edges){
        long timeslot = get<2>(edge);
        if (timeslot > t_log){
            time_offset[t_log]=i_;
            t_log = timeslot;
        }
        ++i_;
    }

    vector<Edge> edges_t;
    vector<vector<Edge>> time2edges(t_max+1,vector<Edge>({}));
    for (long i=1;i<t_max+1;++i){   // split as timestamp
        for (long j=time_offset[i-1]; j<time_offset[i];++j){
            time2edges[i].emplace_back(edges[j]);
        }
    }

    printf("start t_s to t_max:%lu\n",t_max);

    long t_s,t_e;
    vector<vector<vector<pair<Vertex,Vertex>>>> ecc_totoal; /** <t_s,t_e,v,ecc> **/
    for (t_s=1;t_s<t_max+1;++t_s){
        //DYGraph_Base curG(t_s,time2edges[t_s]);
        auto curG = new DYGraph_Base(t_s,time2edges[t_s]);
        printf("initalized at t_s: %ld\n",t_s);
        for (t_e=t_s+1;t_e<t_max+1;++t_e){
            curG->add_E(time2edges[t_e]);
            printf("done added edges at : %ld\n",t_e);
        }
        ecc_totoal.emplace_back(curG->eccTS);
        break;
        /** visualize result below **/
//        printf("\n");
//        int time_i=1;
//        for (auto& eccA:curG->eccTS) {
//            if (!(curG->valid_check[time_i-1])){
//                time_i++;
//                continue;
//            }
//            cout<<time_i++<<endl;
//            for (auto& [v,d]:eccA){
//                printf("%lu ",v);
//            }
//            printf("\n");
//            for (auto& [v,d]:eccA){
//                printf("%lu ",d);
//            }
//            printf("\n");
//        }
//        printf("ecc as following\n");
//        for (int src = 0; src < 21; ++src) {
//            queue<Vertex> Q;
//            set<Vertex> S;
//            Q.push(src);
//            S.insert(src);
//            Vertex dist = 0;
//            Vertex MD = 0;
//            // cout<<endl;
//            // cout<<"src:"<<src<<endl;
//            while (!Q.empty()){
//                dist++;
//                // cout<< dist <<": ";
//                auto s = Q.size();
//                // cout<<"("<<s<<")";
//                for (int i = 0; i < s; ++i) {
//                    auto cur = Q.front();
//                    Q.pop();
//                    for (auto x: curG->curV.at(cur)) {
//                        if (S.count(x)) {
//                            continue;
//                        }
//                        // cout<< x <<" ";
//                        S.insert(x);
//                        Q.push(x);
//                        MD=dist;
//                    }
//                }
//            }
//            cout<<MD<<" ";
//        }
//        break;
    }
}


void unified_algorithm(const vector<Edge>& edges) {
    long t_max = get<2>(edges[edges.size() - 1]); // last timestamp t_max
    vector<long> time_offset(t_max + 1, 0);   // offset array
    time_offset[t_max] = (long) edges.size();  // last timestamp
    // Edge tri;
    // tuple_element<2, decltype(tri)>::type time_stmp;
    long t_log = 1;
    long i_ = 0;
    for (auto &edge: edges) {
        long timeslot = get<2>(edge);
        if (timeslot > t_log) {
            time_offset[t_log] = i_;
            t_log = timeslot;
        }
        ++i_;
    }

    vector<Edge> edges_t;
    vector<vector<Edge>> time2edges(t_max + 1, vector<Edge>({}));
    for (long i = 1; i < t_max + 1; ++i) {   // split as timestamp
        for (long j = time_offset[i - 1]; j < time_offset[i]; ++j) {
            time2edges[i].emplace_back(edges[j]);
        }
    }

    printf("start t_s to t_max:%lu\n", t_max);


    vector<vector<vector<pair<Vertex, Vertex>>>> ecc_totoal; /** <t_s,t_e,v,ecc> **/

    auto* curG = new DYGraph_Unified(1,time2edges[1]);
    for (long i = 2; i < t_max + 1; ++i){
        // add new edges

    }

}
int main(int argc,char* argv[]) {
    /** data_path might be an array of paths respecting individual datasets **/


    // constexpr string_view dataPATH = R"(/Users/siqingzhang/Desktop/testPY/wiki-talk-temporal-fine-20-example.txt)";
    constexpr const char * dataPATH = R"(/Users/siqingzhang/Desktop/testPY/wiki-talk-temporal-fine-24818-example.txt)";
    //
    // constexpr const char * dataPath = argv[1];

    const int paraNums = 3;    // number of elements of one line in data doc : u, v, timestamp
    //const int paraNum = atoi(argv[2]);

    vector<Edge> temp_edges;
    try {
        temp_edges = edgesReceive(paraNums,dataPATH);
        puts("temporal edges received");
    } catch(fileOpenException& OE){
        cout<<"fail to open file: "<<OE.what()<<endl;
        //
    } catch(...){
        puts("unexpected error occurred during input");
    }

    /**
     * end input new edges from files
     * **/
    puts("baseline algorithm below");
    auto t = clock();
    baseline_algorithm(temp_edges);
    t = clock() - t;
    puts("baseline algorithm ends here");
    printf("total time=%0.3lf sec\n", (double)t * 1.0 / CLOCKS_PER_SEC);
    puts("~~~~~~");

    printf("unified algorithm below\n");
    auto t_ = clock();
    unified_algorithm(temp_edges);
    t_ = clock() - t_;
    printf("unified algorithm ends here\n");
    printf("total time=%0.3lf sec\n", (double)t_ * 1.0 / CLOCKS_PER_SEC);
    return 0;
}
