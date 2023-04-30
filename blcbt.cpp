#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <cmath>
#include <boost/functional/hash.hpp>
#include <cstdlib>  
#include <random>
#include <chrono>

int K;
double EPSILON;

using namespace std;

struct Vertex {
    int id;
    vector<int> neighbours;
};

struct Edge {
    int uid;
    int vid;
    int sign;
    int balanceSupport;
    int unbalanceSupport;
    unordered_set<pair<int, int>, boost::hash<pair<int, int>>> blooms; //all the the blooms containing the edge
    bool fixed;
    bool pruned() {
        return balanceSupport < (1 - EPSILON) * K || balanceSupport + unbalanceSupport < K;
    }
    bool isCandidate() {
        return unbalanceSupport * 1.0 / (balanceSupport + unbalanceSupport) > EPSILON;
    }
    int followers; // number of followers
};

struct Bloom {
    int uid; //The vertex with highest priority
    int wid; //The vertex a same side with u

    unordered_set<int> sWedges; 
    unordered_set<int> dWedges; 
};



struct Graph {
    long long balance_count;
    long long unbalance_count;
    vector<int> vID; //All IDs
    unordered_map<int, Vertex> vMap; //ID and the actual Vertex
    vector<pair<int, int>> eID; // All edge IDs
    unordered_map<pair<int, int>, Edge, boost::hash<pair<int, int>>> eMap; // Left ID, right ID, and the actual Edge
    unordered_set<pair<int, int>, boost::hash<pair<int, int>>> candSet; // set of candidate edges to be removed
    unordered_map<pair<int, int>, Bloom, boost::hash<pair<int, int>>> bIndex; // uid, wid, and the bloom
    vector<pair<int, int>> fixedEID;  // global variable for exact
    vector<Graph> partitions; // butterfly-connected subgraphs
    
    
    int currentBestBBSize;
    
    void addEdge(int uid, int vid, int sign) {
        if (eMap.find(make_pair(uid, vid)) != eMap.end()) {
            return;
        }
        // Left IDs are positive
        if (vMap.find(uid) == vMap.end()) {
            vID.push_back(uid);
            Vertex u;
            u.id = uid;
            vMap[uid] = u;
        } 
        // Right IDs are negative
        Vertex v;
        if (vMap.find(vid) == vMap.end()) { 
            vID.push_back(vid);
            Vertex v;
            v.id = vid;
            vMap[vid] = v;
        }
        vMap[uid].neighbours.push_back(vid);
        vMap[vid].neighbours.push_back(uid);
        eMap[make_pair(uid, vid)].uid = uid;
        eMap[make_pair(uid, vid)].vid = vid;
        eMap[make_pair(uid, vid)].sign = sign;
        eMap[make_pair(uid, vid)].balanceSupport = 0;
        eMap[make_pair(uid, vid)].unbalanceSupport = 0;
        eMap[make_pair(uid, vid)].fixed = false;
        eMap[make_pair(uid, vid)].followers = -1;

        //cout << eMap[make_pair(u.id, v.id)].uid << " " << eMap[make_pair(u.id, v.id)].vid  << " " << eMap[make_pair(u.id, v.id)].sign << endl;
    }

    /***
    * Operator for Vertex Priority of Two Nodes
    */
    bool vpComp(const Vertex u, const Vertex v) {
        if (u.neighbours.size() == v.neighbours.size()) {
            return u.id > v.id;
        }
        return u.neighbours.size() > v.neighbours.size();
    }

    /***
    * Reverse Operator for Vertex Priority of Two Nodes
    * Used for the purposes of sorting
    */

    bool vpRevComp(const Vertex u, const Vertex v) {
        if (u.neighbours.size() == v.neighbours.size()) {
            return u.id < v.id;
        }
        return u.neighbours.size() < v.neighbours.size();
    }


    

    void vpNeighbourSort(Vertex &v) {	
        sort(v.neighbours.begin(), v.neighbours.end(), [&](const int &u, const int &w) -> bool {
            return vpRevComp(vMap[u], vMap[w]);
        });
    }


    void initializeVp() {
        for (unordered_map<int, Vertex>::iterator i = vMap.begin(); i != vMap.end(); i++) {
            vpNeighbourSort(i->second);
        }
    }
    

    void graphStat() {
        int leftSize = 0, rightSize = 0, pEdgeSize = 0, nEdgeSize = 0;
        for (auto i: vMap) {
            if (i.second.neighbours.size() > 0) {
                if (i.first > 0) {
                    leftSize++;
                } else {
                    rightSize++;
                }
            }
        }

        for (auto i: eMap) {
            if (i.second.sign > 0) {
                pEdgeSize++;
            } else {
                nEdgeSize++;
            }
        }
        cout << "Left size: " << leftSize << endl;
        cout << "Right size: " << rightSize << endl; 
        cout << "Positive edge count: " << pEdgeSize << endl;
        cout << "Negative edge count: " << nEdgeSize << endl;
        cout << "Average degree: " << (pEdgeSize + nEdgeSize) * 2.0 / (leftSize + rightSize) << endl;
    }

    size_t indexSize() {
        size_t size = 0;
        for (auto i: eMap) {
            size += i.second.blooms.size() * sizeof(pair<pair<int, int>, pair<int, int>>);
        }
        size += bIndex.size() * sizeof(pair<pair<int, int>, Bloom>);
        for (auto i: bIndex) {
            size += (i.second.sWedges.size() + i.second.dWedges.size()) * sizeof(int);
        }
        return size;
    }

    int positiveEdgeCount() {
        int pEdgeSize = 0;
        for (auto i: eMap) {
            if (i.second.sign > 0) {
                pEdgeSize++;
            } 
        }
        return pEdgeSize;
    }

    int negativeEdgeCount() {
        int nEdgeSize = 0;
        for (auto i: eMap) {
            if (i.second.sign < 0) {
                nEdgeSize++;
            } 
        }
        return nEdgeSize;
    }

    double averageBalanceSpp() {
        return balance_count * 4.0 / eMap.size();
    }
    
    double averageUnbalanceSpp() {
        return unbalance_count * 4.0 / eMap.size();
    }

    void initializeGraph(string filename) {
        
        ifstream file(filename);
        string line;

        if (file.is_open()) {
            //cout << "file " << filename << " opened" << endl;
            while (getline(file, line)) {
                vector<int> nums;
                istringstream is(line);
                int n;
                while( is >> n ) {
                    nums.push_back( n );
                }
                
                addEdge(nums[0], -nums[1], nums[2]);
            }
            initializeVp();
        
        }
        file.close();
    }

    void printEdges() {
        for (auto i: eMap) {
            cout << i.second.uid << "\t" << i.second.vid << "\t" << i.second.sign << "\t" 
                << i.second.balanceSupport << "\t" << i.second.unbalanceSupport << "\t" << i.second.blooms.size() << endl;
        }
    }

    void countButterfliesBL() {
        balance_count = 0;
        unbalance_count = 0;
        unordered_map<pair<int, int>, vector<int>, boost::hash<pair<int, int>>> h;
        for (int u: vID) {
            for (int v: vMap[u].neighbours) {
                if (vpComp(vMap[u], vMap[v])) {
                    for (int w: vMap[v].neighbours) {
                        if (vpComp(vMap[u], vMap[w])) {
                            h[make_pair(u, w)].push_back(v);
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
            }

        }
        for (auto i : h) {
            int u = i.first.first;
            int w = i.first.second;
            for (auto it0 = i.second.begin(); it0 + 1 != i.second.end(); it0++) {
                for (auto it1 = it0 + 1; it1 != i.second.end(); it1++) {
                    if (u > 0) { // u and w are in left side
                        if (eMap[make_pair(u, *it0)].sign * eMap[make_pair(u, *it1)].sign * eMap[make_pair(w, *it0)].sign * eMap[make_pair(w, *it1)].sign == 1) {
                            balance_count++;
                        } else {
                            unbalance_count++;
                        }
                    } else { // u and w are in right side
                        if (eMap[make_pair(*it0, u)].sign * eMap[make_pair(*it1, u)].sign * eMap[make_pair(*it0, w)].sign * eMap[make_pair(*it1, w)].sign == 1) {
                            balance_count++;
                        } else {
                            unbalance_count++;
                        }
                    }
                }
            }
        }
    }

    void countButterflies() {
        balance_count = 0;
        unbalance_count = 0;
        // First stores s wedges, second stores d wedges
        unordered_map<pair<int, int>, pair<int, int>, boost::hash<pair<int, int>>> h;
        for (int u: vID) {
            if (u > 0) { // u and w are in left side
                for (int v: vMap[u].neighbours) {
                    if (vpComp(vMap[u], vMap[v])) {
                        for (int w: vMap[v].neighbours) {
                            if (vpComp(vMap[u], vMap[w])) {
                                if (eMap[make_pair(u, v)].sign ==  eMap[make_pair(w, v)].sign) {
                                    h[make_pair(u, w)].first++;
                                }  else {
                                    h[make_pair(u, w)].second++;
                                }
                                
                            } else {
                                break;
                            }
                        }
                    } else {
                        break;
                    }
                }
            } else { // u and w are in right side
                for (int v: vMap[u].neighbours) {
                    if (vpComp(vMap[u], vMap[v])) {
                        for (int w: vMap[v].neighbours) {
                            if (vpComp(vMap[u], vMap[w])) {
                                if (eMap[make_pair(v, u)].sign ==  eMap[make_pair(v, w)].sign) {
                                    h[make_pair(u, w)].first++;
                                }  else {
                                    h[make_pair(u, w)].second++;
                                }
                                
                            } else {
                                break;
                            }
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        for (auto i : h) {
            balance_count += (i.second.first * (i.second.first - 1)) / 2;
            balance_count += (i.second.second * (i.second.second - 1)) / 2;
            unbalance_count += i.second.first * i.second.second;
        }
    }


    void constructIndex() {
        balance_count = 0;
        unbalance_count = 0;
        currentBestBBSize = 0;
        for (int u: vID) {
            if (u > 0) { // u and w are in left side
                for (int v: vMap[u].neighbours) {
                    if (vpComp(vMap[u], vMap[v])) {
                        for (int w: vMap[v].neighbours) {
                            if (vpComp(vMap[u], vMap[w])) {
                                pair<int, int> uv = make_pair(u, v);
                                pair<int, int> vw = make_pair(w, v);
                                bIndex[make_pair(u, w)].uid = u;
                                bIndex[make_pair(u, w)].wid = w;
                                if (eMap[uv].sign ==  eMap[vw].sign) {
                                    bIndex[make_pair(u, w)].sWedges.insert(v);
                                }  else {
                                    bIndex[make_pair(u, w)].dWedges.insert(v);
                                }
                                eMap[uv].blooms.insert(make_pair(u, w));
                                eMap[vw].blooms.insert(make_pair(u, w));
                                
                            } else {
                                break;
                            }
                        }
                    } else {
                        break;
                    }
                }
            } else { // u and w are in right side
                for (int v: vMap[u].neighbours) {
                    if (vpComp(vMap[u], vMap[v])) {
                        for (int w: vMap[v].neighbours) {
                            if (vpComp(vMap[u], vMap[w])) {
                                pair<int, int> uv = make_pair(v, u);
                                pair<int, int> vw = make_pair(v, w);
                                bIndex[make_pair(u, w)].uid = u;
                                bIndex[make_pair(u, w)].wid = w;
                                if (eMap[uv].sign ==  eMap[vw].sign) {
                                    bIndex[make_pair(u, w)].sWedges.insert(v);
                                }  else {
                                    bIndex[make_pair(u, w)].dWedges.insert(v);
                                }
                                eMap[uv].blooms.insert(make_pair(u, w));
                                eMap[vw].blooms.insert(make_pair(u, w));
                                
                            } else {
                                break;
                            }
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        for (auto i :bIndex) {
            int u = i.second.uid;
            int w = i.second.wid;
            pair<int, int> uv, vw;
            if (u > 0) { // u and w are in left side
                for (auto v: i.second.sWedges) {
                    uv = make_pair(u, v);
                    vw = make_pair(w, v);
                    eMap[uv].balanceSupport += i.second.sWedges.size() - 1;
                    eMap[vw].balanceSupport += i.second.sWedges.size() - 1;
                    eMap[uv].unbalanceSupport += i.second.dWedges.size();
                    eMap[vw].unbalanceSupport += i.second.dWedges.size();
                }
                for (auto v: i.second.dWedges) {
                    uv = make_pair(u, v);
                    vw = make_pair(w, v);
                    eMap[uv].balanceSupport += i.second.dWedges.size() - 1;
                    eMap[vw].balanceSupport += i.second.dWedges.size() - 1;
                    eMap[uv].unbalanceSupport += i.second.sWedges.size();
                    eMap[vw].unbalanceSupport += i.second.sWedges.size();
                }
            } else { // u and w are in right side
                for (auto v: i.second.sWedges) {
                    uv = make_pair(v, u);
                    vw = make_pair(v, w);
                    eMap[uv].balanceSupport += i.second.sWedges.size() - 1;
                    eMap[vw].balanceSupport += i.second.sWedges.size() - 1;
                    eMap[uv].unbalanceSupport += i.second.dWedges.size();
                    eMap[vw].unbalanceSupport += i.second.dWedges.size();
                }
                for (auto v: i.second.dWedges) {
                    uv = make_pair(v, u);
                    vw = make_pair(v, w);
                    eMap[uv].balanceSupport += i.second.dWedges.size() - 1;
                    eMap[vw].balanceSupport += i.second.dWedges.size() - 1;
                    eMap[uv].unbalanceSupport += i.second.sWedges.size();
                    eMap[vw].unbalanceSupport += i.second.sWedges.size();
                }
            }
            
            balance_count += (i.second.sWedges.size() * (i.second.sWedges.size() - 1)) / 2;
            balance_count += (i.second.dWedges.size() * (i.second.dWedges.size() - 1)) / 2;
            unbalance_count += i.second.sWedges.size() * i.second.dWedges.size();
        }
    }

    void removeEdge(Edge & e) {
        for (auto i: e.blooms) {
            int u = i.first;
            int w = i.second;
            pair<int, int> twin, uv, vw;
            if (u > 0) { // u and w are in left side
                if (e.uid == u) {
                    twin = make_pair(w, e.vid);
                } else {
                    twin = make_pair(u, e.vid);
                }
                if (e.sign == eMap[twin].sign) {
                    bIndex[make_pair(u, w)].sWedges.erase(e.vid); 
                    //bIndex[make_pair(u, w)].sWedges.erase(remove(bIndex[make_pair(u, w)].sWedges.begin(), bIndex[make_pair(u, w)].sWedges.end(), e.vid), bIndex[make_pair(u, w)].sWedges.end());
                } else {
                    bIndex[make_pair(u, w)].dWedges.erase(e.vid); 
                    //bIndex[make_pair(u, w)].dWedges.erase(remove(bIndex[make_pair(u, w)].dWedges.begin(), bIndex[make_pair(u, w)].dWedges.end(), e.vid), bIndex[make_pair(u, w)].dWedges.end());
                }
                if (bIndex[make_pair(u, w)].sWedges.empty() && bIndex[make_pair(u, w)].dWedges.empty()){
                    bIndex.erase(make_pair(u, w));
                }

            } else { // u and w are in right side
                if (e.vid == u) {
                    twin = make_pair(e.uid, w);
                } else {
                    twin = make_pair(e.uid, u);
                }
                if (e.sign == eMap[twin].sign) {
                    bIndex[make_pair(u, w)].sWedges.erase(e.uid); 
                    //bIndex[make_pair(u, w)].sWedges.erase(remove(bIndex[make_pair(u, w)].sWedges.begin(), bIndex[make_pair(u, w)].sWedges.end(), e.uid), bIndex[make_pair(u, w)].sWedges.end());
                } else {
                    bIndex[make_pair(u, w)].dWedges.erase(e.uid); 
                    //bIndex[make_pair(u, w)].dWedges.erase(remove(bIndex[make_pair(u, w)].dWedges.begin(), bIndex[make_pair(u, w)].dWedges.end(), e.uid), bIndex[make_pair(u, w)].dWedges.end());
                }
                if (bIndex[make_pair(u, w)].sWedges.empty() && bIndex[make_pair(u, w)].dWedges.empty()){
                    bIndex.erase(make_pair(u, w));
                }
            }
            eMap[twin].blooms.erase(make_pair(u, w));
            for (auto v: bIndex[make_pair(u, w)].sWedges) {
                if (u > 0) { // u and w are in left side
                    uv = make_pair(u, v);
                    vw = make_pair(w, v);
                    
                } else { // u and w are in right side
                    uv = make_pair(v, u);
                    vw = make_pair(v, w);
                }
                if (e.sign == eMap[twin].sign) {
                    eMap[uv].balanceSupport -= 1;
                    eMap[vw].balanceSupport -= 1;
                    eMap[twin].balanceSupport -= 1;
                    balance_count -= 1;
                } else {
                    eMap[uv].unbalanceSupport -= 1;
                    eMap[vw].unbalanceSupport -= 1;
                    eMap[twin].unbalanceSupport -= 1;
                    unbalance_count -= 1;
                }
                if (eMap[uv].isCandidate()) {
                    candSet.insert(uv);
                }
                if (eMap[vw].isCandidate()) {
                    candSet.insert(vw);
                }
            }
            for (auto v: bIndex[make_pair(u, w)].dWedges) {
                if (u > 0) { // u and w are in left side
                    uv = make_pair(u, v);
                    vw = make_pair(w, v);
                } else { // u and w are in right side
                    uv = make_pair(v, u);
                    vw = make_pair(v, w);
                }
                if (e.sign == eMap[twin].sign) {
                    eMap[uv].unbalanceSupport -= 1;
                    eMap[vw].unbalanceSupport -= 1;
                    eMap[twin].unbalanceSupport -= 1;
                    unbalance_count -= 1;
                } else {
                    eMap[uv].balanceSupport -= 1;
                    eMap[vw].balanceSupport -= 1;
                    eMap[twin].balanceSupport -= 1;
                    balance_count -= 1;
                }
                if (eMap[uv].isCandidate()) {
                    candSet.insert(uv);
                }
                if (eMap[vw].isCandidate()) {
                    candSet.insert(vw);
                }
            }
            if (eMap[twin].isCandidate()) {
                candSet.insert(twin);
            }
        }
        eMap.erase(make_pair(e.uid, e.vid));
        candSet.erase(make_pair(e.uid, e.vid));
        vMap[e.uid].neighbours.erase(remove(vMap[e.uid].neighbours.begin(), vMap[e.uid].neighbours.end(), e.vid), vMap[e.uid].neighbours.end());
        vMap[e.vid].neighbours.erase(remove(vMap[e.vid].neighbours.begin(), vMap[e.vid].neighbours.end(), e.uid), vMap[e.vid].neighbours.end());
        //vMap[e.uid].neighbours.erase(e.vid);
        //vMap[e.vid].neighbours.erase(e.uid);
    }

    void recursiveRemoveEdge(queue<pair<int, int>> & eq) {
        while (eq.size() > 0) {
            if (eMap.find(eq.front()) != eMap.end()) {
                Edge & e = eMap[eq.front()];
                //cout << "(" << e.uid << ", " << e.vid << ", " << e.balanceSupport << ", " << e.unbalanceSupport << ") ";
                eq.pop();
                for (auto i: e.blooms) {
                    int u = i.first;
                    int w = i.second;
                    pair<int, int> twin, uv, vw;
                    if (u > 0) { // u and w are in left side
                        if (e.uid == u) {
                            twin = make_pair(w, e.vid);
                        } else {
                            twin = make_pair(u, e.vid);
                        }
                        if (e.sign == eMap[twin].sign) {
                            bIndex[make_pair(u, w)].sWedges.erase(e.vid);
                            //bIndex[make_pair(u, w)].sWedges.erase(remove(bIndex[make_pair(u, w)].sWedges.begin(), bIndex[make_pair(u, w)].sWedges.end(), e.vid), bIndex[make_pair(u, w)].sWedges.end());
                        } else {
                            bIndex[make_pair(u, w)].dWedges.erase(e.vid);
                            //bIndex[make_pair(u, w)].dWedges.erase(remove(bIndex[make_pair(u, w)].dWedges.begin(), bIndex[make_pair(u, w)].dWedges.end(), e.vid), bIndex[make_pair(u, w)].dWedges.end());
                        }
                        if (bIndex[make_pair(u, w)].sWedges.empty() && bIndex[make_pair(u, w)].dWedges.empty()){
                            bIndex.erase(make_pair(u, w));
                        }

                    } else { // u and w are in right side
                        if (e.vid == u) {
                            twin = make_pair(e.uid, w);
                        } else {
                            twin = make_pair(e.uid, u);
                        }
                        if (e.sign == eMap[twin].sign) {
                            bIndex[make_pair(u, w)].sWedges.erase(e.uid);
                            //bIndex[make_pair(u, w)].sWedges.erase(remove(bIndex[make_pair(u, w)].sWedges.begin(), bIndex[make_pair(u, w)].sWedges.end(), e.uid), bIndex[make_pair(u, w)].sWedges.end());
                        } else {
                            bIndex[make_pair(u, w)].dWedges.erase(e.uid);
                            //bIndex[make_pair(u, w)].dWedges.erase(remove(bIndex[make_pair(u, w)].dWedges.begin(), bIndex[make_pair(u, w)].dWedges.end(), e.uid), bIndex[make_pair(u, w)].dWedges.end());
                        }
                        if (bIndex[make_pair(u, w)].sWedges.empty() && bIndex[make_pair(u, w)].dWedges.empty()){
                            bIndex.erase(make_pair(u, w));
                        }
                    }
                    eMap[twin].blooms.erase(make_pair(u, w));
                    for (auto v: bIndex[make_pair(u, w)].sWedges) {
                        if (u > 0) { // u and w are in left side
                            uv = make_pair(u, v);
                            vw = make_pair(w, v);
                            
                        } else { // u and w are in right side
                            uv = make_pair(v, u);
                            vw = make_pair(v, w);
                        }
                        if (e.sign == eMap[twin].sign) {
                            eMap[uv].balanceSupport -= 1;
                            eMap[vw].balanceSupport -= 1;
                            eMap[twin].balanceSupport -= 1;
                            balance_count -= 1;
                        } else {
                            eMap[uv].unbalanceSupport -= 1;
                            eMap[vw].unbalanceSupport -= 1;
                            eMap[twin].unbalanceSupport -= 1;
                            unbalance_count -= 1;
                        }
                        if (eMap[uv].pruned()) {
                            eq.push(uv);
                        } else if (eMap[uv].isCandidate()) {
                            candSet.insert(uv);
                        } else {
                            candSet.erase(uv);
                        }
                        if (eMap[vw].pruned()) {
                            eq.push(vw);
                        } else if (eMap[vw].isCandidate()) {
                            candSet.insert(vw);
                        } else {
                            candSet.erase(vw);
                        }
                    }
                    for (auto v: bIndex[make_pair(u, w)].dWedges) {
                        if (u > 0) { // u and w are in left side
                            uv = make_pair(u, v);
                            vw = make_pair(w, v);
                        } else { // u and w are in right side
                            uv = make_pair(v, u);
                            vw = make_pair(v, w);
                        }
                        if (e.sign == eMap[twin].sign) {
                            eMap[uv].unbalanceSupport -= 1;
                            eMap[vw].unbalanceSupport -= 1;
                            eMap[twin].unbalanceSupport -= 1;
                            unbalance_count -= 1;
                        } else {
                            eMap[uv].balanceSupport -= 1;
                            eMap[vw].balanceSupport -= 1;
                            eMap[twin].balanceSupport -= 1;
                            balance_count -= 1;
                        }
                        if (eMap[uv].pruned()) {
                            eq.push(uv);
                        } else if (eMap[uv].isCandidate()) {
                            candSet.insert(uv);
                        } else {
                            candSet.erase(uv);
                        }
                        if (eMap[vw].pruned()) {
                            eq.push(vw);
                        } else if (eMap[vw].isCandidate()) {
                            candSet.insert(vw);
                        } else {
                            candSet.erase(vw);
                        }
                    }
                    if (eMap[twin].pruned()) {
                        eq.push(twin);
                    } else if (eMap[twin].isCandidate()) {
                        candSet.insert(twin);
                    } else {
                        candSet.erase(twin);
                    }
                }
            eMap.erase(make_pair(e.uid, e.vid));
            candSet.erase(make_pair(e.uid, e.vid));
            vMap[e.uid].neighbours.erase(remove(vMap[e.uid].neighbours.begin(), vMap[e.uid].neighbours.end(), e.vid), vMap[e.uid].neighbours.end());
            vMap[e.vid].neighbours.erase(remove(vMap[e.vid].neighbours.begin(), vMap[e.vid].neighbours.end(), e.uid), vMap[e.vid].neighbours.end());
            } else {
                eq.pop();
            }
        }
    }

    void countFollowers(Edge & e, int min = INT_MAX) {
        unordered_map<pair<int, int>, pair<int, int>, boost::hash<pair<int, int>>> originalSupport; // Original balance / unbalance support
        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> fSet; // set of followers
        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> vfSet; // set of visited followers
        queue<pair<int, int>> fq; // followers to be investigated
        fq.push(make_pair(e.uid, e.vid));
        fSet.insert(make_pair(e.uid, e.vid));
        while (fq.size() > 0) {
            Edge & f = eMap[fq.front()];
            vfSet.insert(fq.front());
            //cout << "(" << f.uid << ", " << f.vid << ", " << f.balanceSupport << ", " << f.unbalanceSupport << ") ";
            fq.pop();
            for (auto i: f.blooms) {
                int u = i.first;
                int w = i.second;
                pair<int, int> twin, uv, vw;
                if (u > 0) { // u and w are in left side
                    if (f.uid == u) {
                        twin = make_pair(w, f.vid);
                    } else {
                        twin = make_pair(u, f.vid);
                    }

                } else { // u and w are in right side
                    if (f.vid == u) {
                        twin = make_pair(f.uid, w);
                    } else {
                        twin = make_pair(f.uid, u);
                    }
                }
                if (vfSet.find(twin) != vfSet.end()) {
                    continue;
                }
                if (originalSupport.find(twin) == originalSupport.end()) {
                    originalSupport[twin] = make_pair(eMap[twin].balanceSupport, eMap[twin].unbalanceSupport);
                }
                for (auto v: bIndex[make_pair(u, w)].sWedges) {
                    if (u > 0) { // u and w are in left side
                        uv = make_pair(u, v);
                        vw = make_pair(w, v);
                        
                    } else { // u and w are in right side
                        uv = make_pair(v, u);
                        vw = make_pair(v, w);
                    }
                    if (vfSet.find(uv) != vfSet.end() || (vfSet.find(vw) != vfSet.end())) {
                        continue;
                    }
                    if (originalSupport.find(uv) == originalSupport.end()) {
                        originalSupport[uv] = make_pair(eMap[uv].balanceSupport, eMap[uv].unbalanceSupport);
                    }
                    if (originalSupport.find(vw) == originalSupport.end()) {
                        originalSupport[vw] = make_pair(eMap[vw].balanceSupport, eMap[vw].unbalanceSupport);
                    }
                    
                    if (f.sign == eMap[twin].sign) {
                        eMap[uv].balanceSupport -= 1;
                        eMap[vw].balanceSupport -= 1;
                        eMap[twin].balanceSupport -= 1;
                    } else {
                        eMap[uv].unbalanceSupport -= 1;
                        eMap[vw].unbalanceSupport -= 1;
                        eMap[twin].unbalanceSupport -= 1;
                    }
                    if (eMap[uv].pruned() && fSet.find(uv) == fSet.end()) {
                        fq.push(uv);
                        fSet.insert(uv);
                        
                        if (fSet.size() > min) {
                            e.followers = fSet.size();
                            // Recover Supports
                            for (auto i: originalSupport) {
                                eMap[i.first].balanceSupport = i.second.first;
                                eMap[i.first].unbalanceSupport = i.second.second;
                            }
                            return;
                        }
                        
                    }
                    if (eMap[vw].pruned() && fSet.find(vw) == fSet.end()) {
                        fq.push(vw);
                        fSet.insert(vw);
                        if (fSet.size() > min) {
                            e.followers = fSet.size();
                            // Recover Supports
                            for (auto i: originalSupport) {
                                eMap[i.first].balanceSupport = i.second.first;
                                eMap[i.first].unbalanceSupport = i.second.second;
                            }
                            return;
                        }
                    }
                }
                for (auto v: bIndex[make_pair(u, w)].dWedges) {
                    if (u > 0) { // u and w are in left side
                        uv = make_pair(u, v);
                        vw = make_pair(w, v);
                    } else { // u and w are in right side
                        uv = make_pair(v, u);
                        vw = make_pair(v, w);
                    }
                    if (vfSet.find(uv) != vfSet.end() || (vfSet.find(vw) != vfSet.end())) {
                        continue;
                    }
                    if (originalSupport.find(uv) == originalSupport.end()) {
                        originalSupport[uv] = make_pair(eMap[uv].balanceSupport, eMap[uv].unbalanceSupport);
                    }
                    if (originalSupport.find(vw) == originalSupport.end()) {
                        originalSupport[vw] = make_pair(eMap[vw].balanceSupport, eMap[vw].unbalanceSupport);
                    }
                    if (f.sign == eMap[twin].sign) {
                        eMap[uv].unbalanceSupport -= 1;
                        eMap[vw].unbalanceSupport -= 1;
                        eMap[twin].unbalanceSupport -= 1;
                    } else {
                        eMap[uv].balanceSupport -= 1;
                        eMap[vw].balanceSupport -= 1;
                        eMap[twin].balanceSupport -= 1;
                    }
                    if (eMap[uv].pruned() && fSet.find(uv) == fSet.end()) {
                        fq.push(uv);
                        fSet.insert(uv);
                        
                        if (fSet.size() > min) {
                            e.followers = fSet.size();
                            // Recover Supports
                            for (auto i: originalSupport) {
                                eMap[i.first].balanceSupport = i.second.first;
                                eMap[i.first].unbalanceSupport = i.second.second;
                            }
                            return;
                        }
                        
                    }
                    if (eMap[vw].pruned() && fSet.find(vw) == fSet.end()) {
                        fq.push(vw);
                        fSet.insert(vw);
                        
                        if (fSet.size() > min) {
                            e.followers = fSet.size();
                            // Recover Supports
                            for (auto i: originalSupport) {
                                eMap[i.first].balanceSupport = i.second.first;
                                eMap[i.first].unbalanceSupport = i.second.second;
                            }
                            return;
                        }
                        
                    }
                }
                if (eMap[twin].pruned() && fSet.find(twin) == fSet.end()) {
                    fq.push(twin);
                    fSet.insert(twin);
                    
                    if (fSet.size() > min) {
                        e.followers = fSet.size();
                        // Recover Supports
                        for (auto i: originalSupport) {
                            eMap[i.first].balanceSupport = i.second.first;
                            eMap[i.first].unbalanceSupport = i.second.second;
                        }
                        return;
                    }
                    
                }
            }
        }
        e.followers = fSet.size();
        // Recover Supports
        for (auto i: originalSupport) {
            eMap[i.first].balanceSupport = i.second.first;
            eMap[i.first].unbalanceSupport = i.second.second;
        }
    }



    void pruneEdges() {
        queue<pair<int, int>> to_delete;
        for (auto & e: eMap) {
            if (e.second.pruned()) {
                to_delete.push(make_pair(e.second.uid, e.second.vid));
            } else if (e.second.isCandidate()) {
                candSet.insert(e.first);
            }
        }
        recursiveRemoveEdge(to_delete); 
        /*
        auto i = to_delete.begin();
        while (i != to_delete.end()) {
            // eMap[ePair].balanceSupport < eMap[ePair].unbalanceSupport ||
            pair<int, int> ePair = make_pair(i->first, i->second);
            ++i;
            if (eMap.find(ePair) != eMap.end()) {
                queue<pair<int, int>> eq;
                eq.push(ePair);
                 
            }
            
        }
        */
    }

     void random() {
        while (!candSet.empty()) {
            //pruneEdges();
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(0, candSet.size() - 1);
            auto e = std::next(candSet.begin(), dis(gen));
            //cout << "(" << e->first << ", " << e->second << ") ";
            queue<pair<int, int>> eq;
            eq.push(*e);
            recursiveRemoveEdge(eq);
            //cout << candSet.size() << endl;
        } 
    }

    void supportsGreedy() {
        //pruneEdges();
        double max_ubp;
        int tuid, tvid;
        bool removing = true;
        while (removing) {
            max_ubp = 0;
            tuid = 0;
            tvid = 0;

            for (auto & p: candSet) {
                if (eMap[p].unbalanceSupport * 1.0 / (eMap[p].balanceSupport + eMap[p].unbalanceSupport) > max_ubp) {
                    max_ubp = eMap[p].unbalanceSupport * 1.0 / (eMap[p].balanceSupport + eMap[p].unbalanceSupport);
                    tuid = p.first;
                    tvid = p.second;
                }
            }
            pair<int, int> ePair = make_pair(tuid, tvid);
            /*
            for (auto bPair: eMap[ePair].blooms) {
                cout << bIndex[bPair].sWedges.size() << "\t" << bIndex[bPair].dWedges.size() << endl;
            }
            */
            if (tuid != 0 && tvid != 0 && max_ubp > EPSILON) { // max_ubp > EPSILON
                /*
                cout << tuid << " " << tvid << endl;
                cout << max_ubp << "\t" << eMap[ePair].balanceSupport << "\t" << eMap[ePair].unbalanceSupport <<  "\t" << eMap[ePair].blooms.size()<< endl;
                cout << "Remaining edge size: " << eMap.size() << endl;
                cout << "Balanced butterfly count: " << balance_count << endl;
                cout << "Unbalanced butterfly count: " << unbalance_count << endl;
                */
                // cout << "(" << ePair.first << ", " << ePair.second << ") ";
                queue<pair<int, int>> eq;
                eq.push(ePair);
                recursiveRemoveEdge(eq);
                //cout << candSet.size() << endl;
            } else {
                removing = false;
            }
        }
    }

    void followersGreedy() {
        //pruneEdges();
        bool removing = true;
        int min_f;
        while (removing) {
            pair<int, int> ePair = make_pair(0, 0);
            min_f = eMap.size(); // An upper bound of followers
            
            queue<pair<int, int>> eq;
            for (auto & p: candSet) {
                countFollowers(eMap[p], min_f);
                if (eMap[p].followers < min_f) {
                    min_f = eMap[p].followers;
                    ePair = p;
                }
                
                if (eMap[p].followers == 1) {
                     eq.push(p);
                }
                
                
            }

            if (ePair != make_pair(0, 0)) {
                //cout << "(" << ePair.first << ", " << ePair.second << ") ";
                if (min_f > 1) {
                    eq.push(ePair);
                }
                recursiveRemoveEdge(eq);
                int candCount = 0;
                //cout << candSet.size() << endl;
            } else {
                removing = false;
            }
        }
    }


    bool fixedEdgesDeleted() {
        for (auto p: fixedEID) {
            if (eMap.find(p) == eMap.end()) {
                return true;
            }
        }
        return false;
    }

    
    void exact() {
        if (partitions.size() > 1) {
            /*
            int partsum = 0;
            for (Graph & H: partitions){
                partsum += H.eMap.size();
            } 
            cout << "Sum of partition size: " << partsum << endl;
            */
            for (Graph & H: partitions){
                H.currentBestBBSize = 0;
                H.exact();
                currentBestBBSize += H.currentBestBBSize;
            }
            return;
        }

        if (candSet.empty()) {
            if (eMap.size() > currentBestBBSize) {
                currentBestBBSize = eMap.size();
            }
            return;
        } else if (candSet.size() == fixedEID.size()) { // No more edge to delete
            return;
        }
        /*
        for (auto ep: candSet) {
            cout << "(" << ep.first << ", " << ep.second << ")";
            if (eMap[ep].fixed) {
                cout << "fixed, ";
            } else {
                cout << ", ";
            }
        }
        cout << endl;
        */

        //cout << candSet.size() << endl;
        double max_ubp = 0;
        pair<int, int> ePair = make_pair(0, 0);
        for (auto & p: candSet) {
            if (eMap[p].unbalanceSupport * 1.0 / (eMap[p].balanceSupport + eMap[p].unbalanceSupport) > max_ubp && !eMap[p].fixed) {
                max_ubp = eMap[p].unbalanceSupport * 1.0 / (eMap[p].balanceSupport + eMap[p].unbalanceSupport);
                ePair = p;
            }
        }
        if (ePair == make_pair(0, 0)) {
            return;
        }
        Graph H = *this;
        queue<pair<int, int>> eq;
        eq.push(ePair);
        H.recursiveRemoveEdge(eq);
        if (H.eMap.size() < currentBestBBSize) {
            return;
        }
        if (!H.fixedEdgesDeleted()) {
            H.exact();
            currentBestBBSize = H.currentBestBBSize;
        }

        H = *this;
        H.eMap[ePair].fixed = true;
        H.fixedEID.push_back(ePair);
        H.exact();
        if (H.currentBestBBSize > currentBestBBSize) {
            currentBestBBSize = H.currentBestBBSize;
        }

    }

    void constructBCSubgraphs() {
        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> visitedEdges;
        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> visitedBlooms;
        queue<pair<int, int>> eq;
        for (auto & e: eMap) {
            if (visitedEdges.find(e.first) == visitedEdges.end()) {
                eq.push(e.first);
                visitedEdges.insert(e.first);
            }
            Graph H;
            H.balance_count = 0;
            H.unbalance_count = 0;
            while (eq.size() > 0) {
                pair<int, int> ePair = eq.front();
                eq.pop();
                H.addEdge(eMap[ePair].uid, eMap[ePair].vid, eMap[ePair].sign);
                H.eMap[ePair].balanceSupport = eMap[ePair].balanceSupport;
                H.eMap[ePair].unbalanceSupport = eMap[ePair].unbalanceSupport;
                H.eMap[ePair].blooms = eMap[ePair].blooms;
                H.balance_count += eMap[ePair].balanceSupport;
                H.unbalance_count += eMap[ePair].unbalanceSupport;
                if (eMap[ePair].isCandidate()) {
                    H.candSet.insert(ePair);
                }
                for (auto bid : eMap[ePair].blooms) {
                    if (visitedBlooms.find(bid) == visitedBlooms.end()) {
                        int u = bid.first;
                        int w = bid.second;
                        pair<int, int> uv, vw;
                        for (int v: bIndex[bid].sWedges) {
                            if (u > 0) { // u and w are in left side
                                uv = make_pair(u, v);
                                vw = make_pair(w, v);
                            } else { // u and w are in right side
                                uv = make_pair(v, u);
                                vw = make_pair(v, w);
                            }
                            if (visitedEdges.find(uv) == visitedEdges.end()) {
                                eq.push(uv);
                                visitedEdges.insert(uv);
                            }
                            if (visitedEdges.find(vw) == visitedEdges.end()) {
                                eq.push(vw);
                                visitedEdges.insert(vw);
                            }
                        }
                        for (int v: bIndex[bid].dWedges) {
                            if (u > 0) { // u and w are in left side
                                uv = make_pair(u, v);
                                vw = make_pair(w, v);
                            } else { // u and w are in right side
                                uv = make_pair(v, u);
                                vw = make_pair(v, w);
                            }
                            if (visitedEdges.find(uv) == visitedEdges.end()) {
                                eq.push(uv);
                                visitedEdges.insert(uv);
                            }
                            if (visitedEdges.find(vw) == visitedEdges.end()) {
                                eq.push(vw);
                                visitedEdges.insert(vw);
                            }
                        }
                        H.bIndex[bid] = bIndex[bid];
                        visitedBlooms.insert(bid);
                    }
                }
            }
            if (H.eMap.size() > 0) {
                H.balance_count /= 4;
                H.unbalance_count /= 4;
                partitions.push_back(H);  
            }
        }
    }

    void constructECSubgraphs() { // Edge-connected subgraphs
        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> visitedEdges;
        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> visitedBlooms;
        queue<pair<int, int>> eq;
        for (auto & e: eMap) {
            if (visitedEdges.find(e.first) == visitedEdges.end()) {
                eq.push(e.first);
                visitedEdges.insert(e.first);
            }
            Graph H;
            H.balance_count = 0;
            H.unbalance_count = 0;
            while (eq.size() > 0) {
                pair<int, int> ePair = eq.front();
                eq.pop();
                H.addEdge(eMap[ePair].uid, eMap[ePair].vid, eMap[ePair].sign);
                H.eMap[ePair].balanceSupport = eMap[ePair].balanceSupport;
                H.eMap[ePair].unbalanceSupport = eMap[ePair].unbalanceSupport;
                H.eMap[ePair].blooms = eMap[ePair].blooms;
                H.balance_count += eMap[ePair].balanceSupport;
                H.unbalance_count += eMap[ePair].unbalanceSupport;
                if (eMap[ePair].isCandidate()) {
                    H.candSet.insert(ePair);
                }
                for (auto bid : eMap[ePair].blooms) {
                    if (visitedBlooms.find(bid) == visitedBlooms.end()) {
                        H.bIndex[bid] = bIndex[bid];
                        visitedBlooms.insert(bid);
                    }
                }
                for (auto w: vMap[eMap[ePair].uid].neighbours) {
                    if (visitedEdges.find(make_pair(eMap[ePair].uid, w)) == visitedEdges.end()) {
                        eq.push(make_pair(eMap[ePair].uid, w));
                        visitedEdges.insert(make_pair(eMap[ePair].uid, w));
                    }
                }

                for (auto w: vMap[eMap[ePair].vid].neighbours) {
                    if (visitedEdges.find(make_pair(w, eMap[ePair].vid)) == visitedEdges.end()) {
                        eq.push(make_pair(w, eMap[ePair].vid));
                        visitedEdges.insert(make_pair(w, eMap[ePair].vid));
                    }
                }
            }
            if (H.eMap.size() > 0) {
                H.balance_count /= 4;
                H.unbalance_count /= 4;
                partitions.push_back(H);  
            }
        }
    }
};


int main(int argc, char *argv[]) {

    string algorithm;
    string usage = R"(Usage: ./program datafile algorithm k epsilon
Available algorithms are:
countBL: baseline balanced butterfly counting algorithm (does not require k or epsilon).
count: improved balanced butterfly counting algorithm (does not require k or epsilon).
supportsGreedy: greedy heuristics by balanced supports ratio.
followersGreedy: greedy heuristics by followers.
exact: exact algorithm.)";
    // Check if correct number of arguments provided
    if (argc == 3) {
        // Parse algorithm
        if (strcmp(argv[2], "countBL") == 0) {
            algorithm = "countBL";
        } else if (strcmp(argv[2], "count") == 0) {
            algorithm = "count";
        } else {
            cout << "Invalid algorithm selection or number of arguments." << endl;
            cout << usage;
            return 1;
        }

    } else if (argc == 5) {
        // Parse k
        try {
            K = stoi(argv[3]);
        }
        catch (invalid_argument) {
            cout << "Error: k must be an integer\n";
            return 1;
        }

        // Parse epsilon
        try {
            EPSILON = stod(argv[4]);
            if (EPSILON < 0 || EPSILON > 1) {
                throw invalid_argument("");
            }
        }
        catch (invalid_argument) {
            cout << "Error: epsilon must be a real value between 0 and 1\n";
            return 1;
        }

        // Parse algorithm
        if (strcmp(argv[2], "supportsGreedy") == 0) {
            algorithm = "supportsGreedy";
        } else if (strcmp(argv[2], "followersGreedy") == 0) {
            algorithm = "followersGreedy";
        } else if (strcmp(argv[2], "random") == 0) {
            algorithm = "random";
        } else if (strcmp(argv[2], "exact") == 0) {
            algorithm = "exact";
        } else {
            cout << "Invalid algorithm selection." << endl;
            cout << usage;
            return 1;
        }
    } else {
        cout << "Error: Invalid number of arguments\n";
        cout << usage;
        return 1;
    }

    // Parse datafile
    
    Graph G;
    G.initializeGraph(argv[1]);
    

    chrono::high_resolution_clock::time_point t0, t1;
    chrono::duration<float> fs;
    if (algorithm == "countBL") {
        G.graphStat();
        t0 = chrono::high_resolution_clock::now();
        G.countButterfliesBL();
        t1 = chrono::high_resolution_clock::now();
        fs = t1 - t0;
        cout << "Total runtime: " << fs.count() << "s\n";
        cout << "Balanced butterfly count: " << G.balance_count << endl;
        cout << "Unbalanced butterfly count: " << G.unbalance_count << endl;
        return 0;
    } else if (algorithm == "count") {
        G.graphStat();
        t0 = chrono::high_resolution_clock::now();
        G.countButterflies();
        t1 = chrono::high_resolution_clock::now();
        fs = t1 - t0;
        cout << "Total runtime: " << fs.count() << "s\n";
        cout << "Balanced butterfly count: " << G.balance_count << endl;
        cout << "Unbalanced butterfly count: " << G.unbalance_count << endl;
        return 0;
    }
    t0 = chrono::high_resolution_clock::now();
    G.constructIndex();
    G.graphStat();
    cout << "SBE index size: " << G.indexSize() * 1.0 / 1048576 << " MB" << endl;
    // G.printEdges();
    cout << "Balanced butterfly count: " << G.balance_count << endl;
    cout << "Unbalanced butterfly count: " << G.unbalance_count << endl;

    G.pruneEdges();

    if (algorithm == "supportsGreedy") { 
        G.supportsGreedy();
    }
    else if (algorithm == "followersGreedy") {
        G.followersGreedy();
    }
    else if (algorithm == "random") {
        G.random();
    }
    else if (algorithm == "exact") {
        G.constructBCSubgraphs();
        cout << G.partitions.size() << " Partitions." << endl;
        for (Graph & H : G.partitions) {
            cout << "Partition size: " << H.eMap.size() << "\t";
        }
        G.exact();
        cout << "Exact size: " << G.currentBestBBSize << endl;
    } 

    t1 = chrono::high_resolution_clock::now();
    fs = t1 - t0;
    cout << "Total runtime: " << fs.count() << "s\n";
    G.graphStat();
    cout << "Remaining edge size: " << G.eMap.size() << endl;
    cout << "Balanced butterfly count: " << G.balance_count << endl;
    cout << "Unbalanced butterfly count: " << G.unbalance_count << endl;
    return 0;
}

