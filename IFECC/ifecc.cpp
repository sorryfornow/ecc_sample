//
// Created by Siqing Zhang on 4/12/2022.
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <queue>

using namespace std;

typedef long long Vertex;

typedef struct eccBound{
    Vertex u;
    Vertex lowerBound;
    Vertex upperBound;
}eccBound;

struct Label{
    Vertex v;
    int Dist;
};


void computeFFO(unordered_map<Vertex,set<Vertex>>& adjacency,vector<Label>& L, vector<Vertex>& LVisitOrder, Vertex src, Vertex nV){

    // perform BFS
    for(Vertex i=0; i<nV; ++i){
        L[i].v=i;
        L[i].Dist=-1;   // unvisited
    }
    queue<Label> Q;
    L[src].Dist=0;
    LVisitOrder.emplace_back(src);
    Q.push(L[src]);
    while (!Q.empty()){
        Label cur = Q.front();
        Q.pop();
        for (auto iter = adjacency[cur.v].cbegin(); iter != adjacency[cur.v].end();++iter)
        {
            if(L[*iter].Dist>-1){
                continue;
            }
            L[*iter].Dist=cur.Dist+1;
            Q.push(L[*iter]);
            LVisitOrder.emplace_back(L[*iter].v);
        }
    }
}


void Lemma3_1(Vertex v, Vertex z, vector<eccBound>& IFECC, vector<Label>& refNode){
    Vertex distVZ = refNode[v].Dist;
    // IFECC[z].lowerBound == IFECC[z].upperBound == ecc(z)
    IFECC[v].lowerBound = max(max(IFECC[v].lowerBound,distVZ),IFECC[z].lowerBound-distVZ);
    IFECC[v].upperBound = min(IFECC[v].upperBound,IFECC[z].lowerBound+distVZ);
}


void BFS4DIST(unordered_map<Vertex,set<Vertex>>& adjacency, vector<Vertex>& distance,Vertex src){
    distance[src] = 0;
    queue<int> Q;
    Q.push(src);
    while(!Q.empty()){
        int cur = Q.front();
        Q.pop();
        for (auto iter = adjacency[cur].cbegin(); iter != adjacency[cur].end();++iter)
        {
            if(distance[*iter]>-1){
                continue;
            }
            distance[*iter]=distance[cur]+1;
            Q.push(*iter);
        }
    }
    printf("BFS: %lld\n",src);
    for (auto&x:distance){
        printf("%lld ",x);
    }
    printf("\n");
}


vector<Vertex> computeIFECC(Vertex nV, vector<vector<Vertex>>& edges){
    // start to evaluate vector edges regardless timestamp
    printf("~~~starting IFECC computing~~~");
    unordered_map<Vertex,set<Vertex>> adjacency;
    for (auto& edge:edges){
        if (!adjacency[edge[0]].count(edge[1])){
            adjacency[edge[0]].emplace(edge[1]);
        }
        if (!adjacency[edge[1]].count(edge[0])){
            adjacency[edge[1]].emplace(edge[0]);
        }
    }


    vector<eccBound> IFECC(nV);
    for(Vertex i=0;i<nV;++i){
        IFECC[i].u=i;
        IFECC[i].lowerBound=0;
        IFECC[i].upperBound=nV;
    }


    vector<Vertex> degree(nV);
    int maxDeg =0;
    for (Vertex i = 0; i<nV; ++i){
        degree[i]= adjacency[i].size();
        if (degree[i]>maxDeg){
            maxDeg = degree[i];
        }
    }

    vector<Vertex> referenceNodeZ;
    for (Vertex i = 0; i<nV; ++i){
        if (degree[i]==maxDeg){
            referenceNodeZ.emplace_back(i);
        }
    }

    // show reference nodes with degree
//    for(int& z:referenceNodeZ){
//        cout<<z<<": "<<degree[z]<<endl;
//        for (auto iter = adjacency[z].cbegin(); iter != adjacency[z].end();++iter)
//        {
//            cout << *iter<<" ";
//        }
//        cout<<endl;
//    }

    // bfs from each reference node z
    unordered_map<Vertex,vector<Label>> distNodes;// distance index ordering
    vector<Vertex> assignment(nV,referenceNodeZ[0]);
    unordered_map<Vertex,vector<Label>> FFO;// distance reverse ordering

    for(Vertex& z:referenceNodeZ){
        vector<Label> L(nV);
        vector<Label> Lz(nV);
        vector<Vertex> LVisitOrder;
        computeFFO(adjacency,L,LVisitOrder,z,nV);

        IFECC[z].upperBound= L[LVisitOrder[nV-1]].Dist;
        IFECC[z].lowerBound= L[LVisitOrder[nV-1]].Dist;
        distNodes[z]=L;

        for(Vertex i=0;i<nV;++i){
            if (L[i].Dist <= distNodes[assignment[i]][i].Dist){
                assignment[i]=z;
            }
        }

        for(Vertex i=0;i<nV;++i){
            Lz[nV-i-1]=L[LVisitOrder[i]];
        }
        FFO[z]=Lz;

    }

    printf("V^z assignment");
    for(Vertex i=0;i<nV;++i){
        printf("%lld ",assignment[i]);
    }
    printf("\n");

    for(Vertex& z:referenceNodeZ) {
        printf("V^z: %lld\n",z);
        auto x = FFO[z];
        for (Vertex i = 0; i < nV; ++i) {
            printf("%lld ",x[i].v);
        }
        printf("\n");
        for (Vertex i = 0; i < nV; ++i) {
            printf("%d  ",x[i].Dist);
        }
        printf("\n");
    }

    unordered_map<Vertex,vector<Vertex>> VZ;
    for(Vertex i=0; i<nV; ++i){
        if (i==assignment[i]){
            continue;
        }
        VZ[assignment[i]].emplace_back(i);
        Lemma3_1(i,assignment[i],IFECC,distNodes[assignment[i]]);
    }

    printf("reference nodes probe have done.\n");

    vector<Vertex> checkConvergence(nV, -1);
    for(auto& z:referenceNodeZ){ // ecc of reference nodes
        checkConvergence[z]=IFECC[z].lowerBound;
    }
    unordered_map<Vertex,vector<Vertex>> BFSResult;

    for(auto& z:referenceNodeZ) {
        Vertex n_z = VZ[z].size();
        for (auto& V_i :FFO[z]){

            if(!BFSResult.count(V_i.v)){
                vector<Vertex> distV_i_(nV,-1);
                BFS4DIST(adjacency,distV_i_,V_i.v);
                BFSResult[V_i.v]=distV_i_;
            }
            vector<Vertex> distV_i = BFSResult[V_i.v];

            for (auto& v:VZ[z]){
                if (checkConvergence[v]>-1){
                    continue;
                }
                // Lemma3.1 && 3.3
                IFECC[v].lowerBound = max(IFECC[v].lowerBound,distV_i[v]);
                IFECC[v].upperBound = min(IFECC[v].upperBound,
                                          max(IFECC[v].lowerBound,distV_i[z]+distNodes[z][v].Dist));
                if(IFECC[v].lowerBound==IFECC[v].upperBound){
                    checkConvergence[v]=IFECC[v].lowerBound;
                    --n_z;
                }
            }
            if (!n_z){
                break;
            }
        }
    }
//    for(auto&x:IFECC){
//        cout<<x.u<<" "<<x.lowerBound<<" "<<x.upperBound<<endl;
//    }
    return checkConvergence;
}

int main(int argc,char* argv[]) {

    string dataPath;
    dataPath = R"(./wiki-talk-temporal-fine-20-example.txt)";
    // dataPath = argv[1];

    int paraNum;    // number of elements in one line: u, v, timestamp
    paraNum = 3;
    //paraNum = atoi(argv[2]);

    ifstream infile;
    infile.open(dataPath,ios::in);
    if (!infile.is_open()){
        cout << "fail to open file: " << dataPath << endl;
        return 1;
    }

//    vector<vector<int>> edges(lineCount,vector<int>(3));
    vector<vector<Vertex>> edges;
    string curLine;

    const char *pDelimiter = " ";
    char* pToken = nullptr;
    char* pSave = nullptr;

    int lineNum = 0;    // number of edges
//    while(lineNum<=lineCount){
//        getline(infile,curLine);
    set<Vertex> curVertices;

    while(getline(infile,curLine)){
        char * strs = new char[curLine.length() + 1];
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
        curVertices.insert(edge[0]);
        curVertices.insert(edge[1]);
        edges.emplace_back(edge);
        ++lineNum;
        delete[] strs;
    }
    infile.close();



    printf("initially having (%zu) vertices:\n",curVertices.size());
    for(auto x:curVertices){
        printf("%lld ",x);
    }
    printf("\n");

    Vertex nV = curVertices.size();

//    vector<int> IFEccentricities(nV, -1);
    vector<Vertex> IFEccentricities = computeIFECC(nV, edges);

    printf("eccentricities by IFECC:\n");
    for(auto& ecc:IFEccentricities){
        printf("%lld ",ecc);
    }
    printf("\n");

    /* Dynamic Updating*/


    return 0;
}
