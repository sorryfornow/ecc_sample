//
// Created by Siqing Zhang on 16/1/2023.
//
#pragma once

using namespace std;


// helper type for the visitor
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
// explicit deduction guide (not needed as of C++20)
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;


typedef long Vertex; // ecc dist is also Vertex type
// unsigned long is also applicable while -1 representing POSINF

template<typename T1, typename T2, typename T3>
using triplet = tuple<T1, T2, T3>;

typedef triplet<Vertex,Vertex,long> Edge;    //u, v, (long)timestamp

inline constexpr Vertex POSINF = 2'147'483'647;

struct fileOpenException : public exception
{
    const char * errorPath;
    explicit fileOpenException(const char * path_):errorPath(path_){}
    [[nodiscard]] const char * what () const noexcept override
    {
        return errorPath;
    }
};

struct ECC{
    Vertex eccv;
    Vertex ub,lb;
    //vector<Vertex> adj;
//    explicit ECC(const Vertex& x_);
//    explicit ECC(const Vertex& x_, const Vertex& e);
};
//ECC::ECC(const Vertex& x_):eccv(x_){ub=-1;lb=0;adj={};}
//ECC::ECC(const Vertex& x_, const Vertex& e):eccv(x_),ub(e),lb(e){adj={};}


struct adj_repr{
    Vertex v;
    Vertex dist;
    long time;
};

struct adj_time_cmp{
    inline bool operator()(adj_repr a, adj_repr b){
        return a.time==b.time ? a.dist>b.dist:a.time<b.time;
    }   // priority ordering: latest time, then shortest dist
};

/**definition ends**/