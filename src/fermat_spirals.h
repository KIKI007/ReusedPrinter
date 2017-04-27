//
// Created by 汪子琦 on 4/6/17.
//

#ifndef SUPPORTER_FERMAT_SPIRALS_H
#define SUPPORTER_FERMAT_SPIRALS_H

#include <Eigen/Core>
#include "settings.h"
#include "clipper.hpp"

typedef struct tagMSTGraphNode
{
    tagMSTGraphNode()
    {
        key = -1;
        x = 0;
        t1 = t2 = 0;
        e1 = e2 = 0;
    }

    bool valid()
    {
        return key != -1;
    }

    double key;
    int x;
    double t1;
    double t2;
    int e1;
    int e2;
}MSTGraphNode;

typedef struct tagFermatEdge
{
    //edge is in clockwise
public:
    tagFermatEdge()
    {
        p0 = Eigen::Vector2d(0, 0);
        p1 = Eigen::Vector2d(0, 0);
        yh = 0;
    }

    tagFermatEdge(ClipperLib::IntPoint in_p0, ClipperLib::IntPoint in_p1, double h)
    {
        p0 = Eigen::Vector2d(settings.int2mm(in_p0.X), settings.int2mm(in_p0.Y));
        p1 = Eigen::Vector2d(settings.int2mm(in_p1.X), settings.int2mm(in_p1.Y));
        yh = h;
    }

public:

    bool on(Eigen::Vector2d p)
    {
        double cross = (p - p0)[0] * (p - p1)[1] - (p - p0)[1] * (p - p1)[0];
        if(std::abs(cross) < settings.ZERO_EPS && (p - p0).dot(p - p1) < settings.ZERO_EPS)
            return true;
        else
            return false;
    }

    double closest(Eigen::Vector2d p)
    {
        Eigen::Vector2d q0 = p0 - p,
                q1 = p1 - p;
        double t = (q0.dot(q0) -q0.dot(q1)) / (q0.dot(q0) + q1.dot(q1) - 2 * q0.dot(q1));
        if(t < 0)
            t = 0;
        else if (t > 1)
            t = 1;
        return t;
    }

    Eigen::Vector2d point(double t)
    {
        return p0 * (1 - t) + p1 * t;
    }

    Eigen::RowVector3d P0()
    {
        return Eigen::RowVector3d(p0.x(), yh, p0.y());
        //return Eigen::RowVector3d(p0.x(), 0, p0.y());
    }

    Eigen::RowVector3d P1() {
        return Eigen::RowVector3d(p1.x(), yh, p1.y());
        //return Eigen::RowVector3d(p1.x(), 0, p1.y());
    }
public:
    Eigen::Vector2d p0, p1;
    double yh;
    Settings settings;
}FermatEdge;

class Fermat_Level_Set
{

public:
    Fermat_Level_Set()
    {
        levelset.clear();
    }



    Fermat_Level_Set(ClipperLib::Path &outside, double pin_height);

    //return false means it is the center of all level sets
    // t1 means the point is e1.p0 * (1 - t) + e1.p1 * t
    bool gradient(int level, int e1, double t1, int &e2, double &t2);

    //return false means ds is too large for searching
    //only valid for small ds
    bool anticlockwise(int level, int e1, double t1, double ds ,int &e2, double &t2);

    bool clockwise(int level, int e1, double t1, double ds ,int &e2, double &t2);

    double closest(Fermat_Level_Set &poly, int &E1, double &T1, int &E2, double &T2);

    void reroute(int level,
                   int e0, double t0,
                   int e1, double t1,
                   std::list<FermatEdge> &upper, std::list<FermatEdge> &lower);

    void partition(int level,
                   int e0, double t0,
                   int e1, double t1,
                   std::list<FermatEdge> &upper, std::list<FermatEdge> &lower);

    void fermat_spiral(int e0,
                       double t0,
                       std::list<FermatEdge> &path);

    void fermat_spiral(int e0,
                       double t0,
                       std::vector<std::list<FermatEdge>> &path_cut,
                       std::vector<MSTGraphNode> &child);

public:

    void reverse(std::list< FermatEdge> &path);

    void connect(std::list<FermatEdge> &path1, std::list<FermatEdge> &path2);

    void cut_head(std::list<FermatEdge> &path,
             std::list<FermatEdge> &head,
             Eigen::Vector2d p);

public:

    bool empty(){return num_level() == 0;}

    int num_level_id(int id){
        assert(0 <= id && id <= levelset.size());
        return levelset[id].size();
    }

    int num_level(){return levelset.size();}

    void get_level(std::vector<FermatEdge> &polygon, int level)
    {
        assert(0 <= level && level <= levelset.size());
        polygon = levelset[level];
    }


    void get_point(Eigen::Vector2d &p, int level, int e, double t)
    {
        assert(0 <= level && level <= levelset.size());
        assert(0 <= e && e < levelset[level].size());
        assert(0 <= t && t <= 1);

        p =  levelset[level][e].point(t);
        return;
    }


private:

    Settings settings;
    std::vector< std::vector<FermatEdge>> levelset;
};

Fermat_Level_Set::Fermat_Level_Set(ClipperLib::Path &outside, double pin_height)
{
    double outside_area = 0;
    do
    {
        if(!ClipperLib::Orientation(outside))
            ClipperLib::ReversePath(outside);

        std::vector<FermatEdge> polygon;
        for(int jd = 0; jd < outside.size(); jd++)
        {
            ClipperLib::IntPoint p0 = outside[jd];
            ClipperLib::IntPoint p1 = outside[(jd + 1) % outside.size()];
            FermatEdge e(p0, p1, pin_height);
            polygon.push_back(e);
        }
        levelset.push_back(polygon);

        ClipperLib::ClipperOffset co;
        co.AddPath(outside, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
        ClipperLib::Paths solution;
        co.Execute(solution, -settings.mm2int(settings.extrusion_width * 2));
        ClipperLib::SimplifyPolygons(solution);

        outside.clear();
        outside_area = 0;
        for(int jd = 0; jd < solution.size(); jd++)
        {
            double sol_area =  std::abs(ClipperLib::Area(solution[jd]));
            if(outside_area <= sol_area)
            {
                outside_area = sol_area;
                outside = solution[jd];
            }
        }

        outside_area *= settings.UNIT * settings.UNIT;
    }while(outside_area > settings.support_center_area);
}

bool Fermat_Level_Set::gradient(int level, int e1, double t1, int &e2, double &t2) {
    assert(0 <= level && level < num_level());
    assert(0 <= e1 && e1 < levelset[level].size());
    assert(0 <= t1 && t1 <= 1);

    if (level + 1 == num_level())
        return false;

    Eigen::Vector2d p0 = levelset[level][e1].point(t1);
    double closest_dist = 0;
    bool first_closest_flag = true;
    for (int id = 0; id < levelset[level + 1].size(); id++)
    {
        double t = levelset[level + 1][id].closest(p0);
        if((levelset[level + 1][id].point(t) - p0).norm() < closest_dist ||
           first_closest_flag)
        {
            first_closest_flag = false;
            closest_dist = (levelset[level + 1][id].point(t) - p0).norm();
            e2 = id;
            t2 = t;
        }
    }

    return true;
}

void Fermat_Level_Set::fermat_spiral(int e0, double t0, std::vector<std::list<FermatEdge>> &path_cut, std::vector<MSTGraphNode> &child)
{
    std::list<FermatEdge> fermat_path;
    fermat_spiral(e0, t0, fermat_path);

    std::vector<Eigen::Vector2d> p0s, p1s;

    for(int id = 0;id < child.size(); id++)
    {
        Eigen::Vector2d p0, p1;
        get_point(p0, 0, child[id].e1, child[id].t1);
        int e2;double t2;
        clockwise(0, child[id].e1, child[id].t1, 2 * settings.fermat_cut_width, e2, t2);
        get_point(p1, 0, e2, t2);
        p0s.push_back(p0);
        p1s.push_back(p1);
    }

    std::list<FermatEdge> head;
    for(int id = 0; id < child.size(); id++)
    {
        cut_head(fermat_path, head, p0s[id]);
        path_cut.push_back(head);
        cut_head(fermat_path, head, p1s[id]);
    }
    path_cut.push_back(fermat_path);

    return;
}

void Fermat_Level_Set::cut_head(std::list<FermatEdge> &path,
                           std::list<FermatEdge> &head,
                           Eigen::Vector2d p)
{
    std::list<FermatEdge>::iterator it;
    head.clear();
    for(it = path.begin(); it != path.end();)
    {
        if(it->on(p))
        {
            //head
            FermatEdge fe = *it;
            fe.p1 = p;
            head.push_back(fe);

            //reset path
            fe.p0 = fe.p1;
            fe.p1 = it->p1;
            path.erase(it);
            path.push_front(fe);
            break;
        }
        else
        {
            head.push_back(*it);
            it = path.erase(it);
        }
    }
}


bool Fermat_Level_Set::anticlockwise(int level, int e1, double t1, double ds ,int &e2, double &t2) {
    assert(0 <= level && level < num_level());
    assert(0 <= e1 && e1 < levelset[level].size());
    assert(0 <= t1 && t1 <= 1);
    assert(0 <= ds);

    Eigen::Vector2d p0 = levelset[level][e1].point(t1);
    ds += (levelset[level][e1].p1 - p0).norm();

    int id = e1;
    do
    {
        double s = (levelset[level][id].p1 - levelset[level][id].p0).norm();
        if(s >= ds)
        {
            e2 = id;
            t2 = (s - ds) / s;
            break;
        }
        ds -= s;
        id = (id - 1 + levelset[level].size()) % levelset[level].size();
    }while(id != e1);

//    if(id == e1)
//        return false;

    return true;
}

bool Fermat_Level_Set::clockwise(int level, int e1, double t1, double ds ,int &e2, double &t2)
{
    assert(0 <= level && level < num_level());
    assert(0 <= e1 && e1 < levelset[level].size());
    assert(0 <= t1 && t1 <= 1);
    assert(0 <= ds);

    Eigen::Vector2d p0 = levelset[level][e1].point(t1);
    ds += (levelset[level][e1].p0 - p0).norm();

    int id = e1;
    t2 = t1;
    do
    {
        double s = (levelset[level][id].p1 - levelset[level][id].p0).norm();
        if(s >= ds)
        {
            e2 = id;
            t2 = ds / s;
            break;
        }
        ds -= s;
        id = (id + 1) % levelset[level].size();
    }while(id != e1);

    return true;
}

double Fermat_Level_Set::closest(Fermat_Level_Set &poly, int &E1, double &T1, int &E2, double &T2)
{
    std::vector<FermatEdge> outside1;
    poly.get_level(outside1, 0);

    double closest_dist = 0;
    bool first_closest_flag = true;

    for(int e1 = 0; e1 < levelset[0].size(); e1++)
    {
        for(double t1 = 0; t1 <= 1; t1 += 0.1)
        {
            Eigen::Vector2d p1 = levelset[0][e1].point(t1);
            for(int e2 = 0; e2 < outside1.size(); e2++)
            {
                double t2 = outside1[e2].closest(p1);
                Eigen::Vector2d p2 = outside1[e2].point(t2);
                if((p1 - p2).norm() < closest_dist
                   || ((p1 - p2).norm() == closest_dist && std::abs(t1 - 0.5) < std::abs(T1 - 0.5))
                   || first_closest_flag)
                {
                    first_closest_flag = false;
                    closest_dist = (p1 - p2).norm();

                    E1 = e1;
                    T1 = t1;
                    E2 = e2;
                    T2 = t2;
                }
            }
        }
    }

    return closest_dist;
}


void Fermat_Level_Set::reroute(int level,
             int e0, double t0,
             int e1, double t1,
             std::list<FermatEdge> &upper, std::list<FermatEdge> &lower)
{
    int e2, e3;
    double t2, t3;
    anticlockwise(level, e0, t0, settings.fermat_cut_width, e2, t2);
    anticlockwise(level, e1, t1, settings.fermat_cut_width, e3, t3);

    //upper
    FermatEdge fermatEdge;
    int size = levelset[level].size();
    if(e0 != e3 || t0 > t3)
    {
        fermatEdge = levelset[level][e0];
        fermatEdge.p0 = levelset[level][e0].point(t0);
        upper.push_back(fermatEdge);

        for(int id = (e0 + 1) % size; id != e3; id = (id + 1) % size)
        {
            upper.push_back(levelset[level][id]);
        }

        fermatEdge = levelset[level][e3];
        fermatEdge.p1 = levelset[level][e3].point(t3);
        upper.push_back(fermatEdge);
    }
    else
    {
        fermatEdge = levelset[level][e0];
        fermatEdge.p0 = levelset[level][e0].point(t0);
        fermatEdge.p1 = levelset[level][e0].point(t3);
        upper.push_back(fermatEdge);
    }

    //lower
    if(e1 != e2 || t1 > t2)
    {
        fermatEdge = levelset[level][e1];
        fermatEdge.p0 = levelset[level][e1].point(t1);
        lower.push_back(fermatEdge);

        for(int id = (e1 + 1) % size; id != e2; id = (id + 1) % size)
        {
            lower.push_back(levelset[level][id]);
        }

        fermatEdge = levelset[level][e2];
        fermatEdge.p1 = levelset[level][e2].point(t2);
        lower.push_back(fermatEdge);
    }
    else
    {
        fermatEdge = levelset[level][e1];
        fermatEdge.p0 = levelset[level][e1].point(t1);
        fermatEdge.p1 = levelset[level][e1].point(t2);
        lower.push_back(fermatEdge);
    }

    return;
}

void Fermat_Level_Set::partition(int level,
               int e0, double t0,
               int e1, double t1,
               std::list<FermatEdge> &upper, std::list<FermatEdge> &lower)
{

    //upper
    FermatEdge fermatEdge;
    int size = levelset[level].size();

    if(e0 != e1 || t0 > t1) {
        fermatEdge = levelset[level][e0];
        fermatEdge.p0 = levelset[level][e0].point(t0);
        upper.push_back(fermatEdge);

        for (int id = (e0 + 1) % size; id != e1; id = (id + 1) % size) {
            upper.push_back(levelset[level][id]);
        }

        fermatEdge = levelset[level][e1];
        fermatEdge.p1 = levelset[level][e1].point(t1);
        upper.push_back(fermatEdge);
    }
    else
    {
        fermatEdge = levelset[level][e0];
        fermatEdge.p0 = levelset[level][e0].point(t0);
        fermatEdge.p1 = levelset[level][e0].point(t1);
        upper.push_back(fermatEdge);
    }

    //lower
    if(e0 != e1 || t1 > t0)
    {
        fermatEdge = levelset[level][e1];
        fermatEdge.p0 = levelset[level][e1].point(t1);
        lower.push_back(fermatEdge);

        for (int id = (e1 + 1) % size; id != e0; id = (id + 1) % size) {
            lower.push_back(levelset[level][id]);
        }

        fermatEdge = levelset[level][e0];
        fermatEdge.p1 = levelset[level][e0].point(t0);
        lower.push_back(fermatEdge);

    }
    else
    {
        fermatEdge = levelset[level][e1];
        fermatEdge.p0 = levelset[level][e1].point(t1);
        fermatEdge.p1 = levelset[level][e1].point(t0);
        lower.push_back(fermatEdge);
    }
    return;
}

void Fermat_Level_Set::connect(std::list<FermatEdge> &path1, std::list<FermatEdge> &path2)
{
    if(!path1.empty())
    {
        FermatEdge fermatEdge = path1.back();
        fermatEdge.p0 = fermatEdge.p1;
        fermatEdge.p1 = path2.front().p0;
        path1.push_back(fermatEdge);
    }

    path1.insert(path1.end(), path2.begin(), path2.end());

    return;
}

void Fermat_Level_Set::reverse(std::list<FermatEdge> &path)
{
    path.reverse();
    std::list<FermatEdge>::iterator it;
    for(it = path.begin(); it != path.end(); it++)
    {
        Eigen::Vector2d tmp = it->p0;
        it->p0 = it->p1;
        it->p1 = tmp;
    }

    return;
}

void Fermat_Level_Set::fermat_spiral(int e0, double t0, std::list<FermatEdge> &path)
{
    std::vector< std::list<FermatEdge>> lowers, uppers;
    lowers.resize(num_level() - 1);
    uppers.resize(num_level() - 1);

    int level = 0 , e1 = 0;
    double t1 = 0;

    anticlockwise(0, e0, t0, 2 * settings.fermat_cut_width, e1, t1);

    if(num_level() <= 1)
    {
        std::list< FermatEdge> centers[2];
        partition(0, e0, t0, e1, t1, centers[0], centers[1]);
        path = centers[0];
        return;
    }

    while(level < num_level() - 1)
    {
        reroute(level, e0, t0, e1, t1, uppers[level], lowers[level]);
        double t2, t3;
        int e2, e3;
        gradient(level, e0, t0, e2, t2);
        anticlockwise(level + 1, e2, t2, 2 * settings.fermat_cut_width, e3, t3);
        e0 = e2;t0 = t2;
        e1 = e3;t1 = t3;
        level++;
    }


    std::list< FermatEdge> centers[2];
    partition(level, e0, t0, e1, t1, centers[0], centers[1]);

    //inward
    FermatEdge fermatEdge;
    for(level = 0; level < num_level() - 1; level ++)
    {
        if(level % 2 == 0) {
            connect(path, uppers[level]);
        }
        else{
            connect(path, lowers[level]);
        }

    }

    //center
    double dist_center[2] = {0, 0};
    std::list<FermatEdge>::iterator it;
    for(int id = 0; id < 2; id++)
    {
        for(it = centers[id].begin(); it != centers[id].end(); ++it)
        {
            dist_center[id] += (it->p0 - it->p1).norm();
        }
    }

    if(dist_center[0] > dist_center[1])
    {
        if(num_level() % 2 == 0) reverse(centers[0]);
        connect(path, centers[0]);
    }
    else
    {
        if(num_level() % 2 == 1) reverse(centers[1]);
        connect(path, centers[1]);
    }

    //outward
    for(level = num_level() - 2; level >= 0; level --)
    {
        if(level % 2 == 1)
        {
            reverse(uppers[level]);
            connect(path, uppers[level]);
        }
        else {
            reverse(lowers[level]);
            connect(path, lowers[level]);
        }
    }

    return;
}

//void Fermat_Level_Set::get_level(std::vector<FermatEdge> &polygon, int level, int e1, int t1, int e2, int t2)
//{
//    if(e1 == e2 && t1 < t2)
//    {

//        FermatEdge fe = levelset[level][e1];
//        fe.p0 = levelset[level][e1].point(t1);
//        fe.p1 = levelset[level][e1].point(t2);
//        polygon.push_back(fe);
//
//    } else{
//
//        FermatEdge fe = levelset[level][e1];
//        fe.p0 = levelset[level][e1].point(t1);
//        polygon.push_back(fe);
//
//        for(int e = e1 +1; e != e2; e++)
//        {
//            polygon.push_back(levelset[level][e]);
//        }
//
//        FermatEdge fe = levelset[level][e2];
//        fe.p1 = levelset[level][e2].point(t2);
//        polygon.push_back(fe);
//    }
//    return;
//}


#endif //SUPPORTER_FERMAT_SPIRALS_H
