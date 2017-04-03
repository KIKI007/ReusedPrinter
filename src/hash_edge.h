//
// Created by 汪子琦 on 3/26/17.
//

#ifndef SUPPORTER_HASH_EDGE_H
#define SUPPORTER_HASH_EDGE_H

#include <vector>
#include <Eigen/Core>

//EdgeHash stores the edge (u, v) and its index
//usage:
//buildConnection(u, v, idx) : connect u -> v and the edge (u,v)'s index is idx
//find(u, v, idx): find the edge which starts at u, and output its v and index.
class HashEdge
{
public:

    HashEdge()
    {
        prime_n = 5003;
        source.resize(prime_n);
        target.resize(prime_n);
        edge.resize(prime_n);
    }

    ~HashEdge()
    {

    }

public:

    int compute_key(Eigen::RowVector2i u)
    {
        //compute hash key value
        //return (u[0] * 73856093) ^ (u[1] * 19349663) ^ (u[2] * 83492791) % prime_n;
        return (((u[0] * 73856093) ^ (u[1] * 19349663)) % prime_n + prime_n) % prime_n;
    }

    void add_info(int key, Eigen::RowVector2i u, Eigen::RowVector2i v, int idx)
    {
        source[key].push_back(u);
        target[key].push_back(v);
        edge[key].push_back(idx);
    }

    void build_connection(Eigen::RowVector2i u, Eigen::RowVector2i v, int idx)
    {
        int key = compute_key(u);
        if(source[key].empty())
            add_info(key, u, v, idx);
        else
        {
            bool duplicate = false;
            for(int kd = 0;kd < source[key].size(); kd++)
            {
                if(source[key][kd] == u)
                {
                    duplicate = true;
                    break;
                }
            }
            if(!duplicate)
                add_info(key, u, v, idx);
        }
    }

    bool find(Eigen::RowVector2i u, Eigen::RowVector2i &v, int &idx)
    {
        int key = compute_key(u);
        if(source[key].empty())
            return false;
        else
        {
            for(int kd = 0;kd < source[key].size(); kd++)
            {
                if(source[key][kd] == u)
                {
                    v = target[key][kd];
                    idx = edge[key][kd];
                    return true;
                }
            }
        }
        return false;
    }

public:

    int prime_n;
    std::vector<std::vector<Eigen::RowVector2i>> source;
    std::vector<std::vector<Eigen::RowVector2i>> target;
    std::vector<std::vector<int>> edge;
};


#endif //SUPPORTER_HASH_EDGE_H
