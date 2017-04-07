//
// Created by 汪子琦 on 4/6/17.
//

#ifndef SUPPORTER_FERMAT_SPIRALS_H
#define SUPPORTER_FERMAT_SPIRALS_H

#include <Eigen/Core>
#include "settings.h"
#include "clipper.hpp"

typedef struct tagFermatEdge
{
    //edge is in clockwise
public:
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
        if(std::abs(cross) < settings.ZERO_EPS && (p - p0).dot(p - p1) < 0)
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
    }

    Eigen::RowVector3d P1() {
        return Eigen::RowVector3d(p1.x(), yh, p1.y());
    }
public:
    Eigen::Vector2d p0, p1;
    double yh;
    Settings settings;
}FermatEdge;

class Fermat_Level_Set
{

public:
    Fermat_Level_Set(ClipperLib::Path &outside, double pin_height)
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
            co.Execute(solution, -settings.mm2int(settings.support_width));
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
            std::cout << outside_area << std::endl;
        }while(outside_area > settings.support_center_area);
    }

    //return false means it is the center of all level sets
    // t1 means the point is e1.p0 * (1 - t) + e1.p1 * t
    bool gradient(int level, int e1, double t1, int &e2, double &t2) {
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

    //return false means ds is too large for searching
    //only valid for small ds
    bool anticlockwise(int level, int e1, double t1, double ds ,int &e2, double &t2) {
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
            id = (id - 1 + levelset[level].size()) % levelset[level].size();
        }while(id != e1);

        if(id == e1)
            return false;

        return true;
    }

    void closest(Fermat_Level_Set &poly, int &E1, double &T1, int &E2, double &T2)
    {
        std::vector<FermatEdge> outside1;
        poly.get_level(outside1, 0);

        for(int e1 = 0; e1 < levelset[0].size(); e1++)
        {
            double closest_dist = 0;
            bool first_closest_flag = true;
            for(double t1 = 0; t1 <= 1; t1 += 0.1)
            {
                Eigen::Vector2d p1 = levelset[0][e1].point(t1);
                for(int e2 = 0; e2 < outside1.size(); e2++)
                {
                    double t2 = outside1[e2].closest(p1);
                    Eigen::Vector2d p2 = outside1[e2].point(t2);
                    if((p1 - p2).norm() < closest_dist ||
                            first_closest_flag)
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

        return;
    }

    void partition(int level,
                   int e0, double t0,
                   int e1, double t1,
                   std::list<FermatEdge> &upper, std::list<FermatEdge> &lower)
    {
        int e2, e3;
        double t2, t3;
        anticlockwise(level, e0, t0, settings.fermat_cut_width, e2, t2);
        anticlockwise(level, e1, t1, settings.fermat_cut_width, e3, t3);

        //upper
        FermatEdge fermatEdge = levelset[level][e0];
        fermatEdge.p0 = levelset[level][e0].point(t0);
        upper.push_back(fermatEdge);

        int size = levelset[level].size();
        for(int id = e0 + 1; id != e3; id = (id + 1) % size)
        {
            upper.push_back(levelset[level][id]);
        }

        fermatEdge = levelset[level][e3];
        fermatEdge.p1 = levelset[level][e3].point(t3);
        upper.push_back(fermatEdge);

        //lower
        fermatEdge = levelset[level][e1];
        fermatEdge.p0 = levelset[level][e0].point(t0);
        upper.push_back(fermatEdge);

        for(int id = e1 + 1; id != e2; id = (id + 1) % size)
        {
            lower.push_back(levelset[level][id]);
        }

        fermatEdge = levelset[level][e2];
        fermatEdge.p1 = levelset[level][e2].point(t2);
        lower.push_back(fermatEdge);

        return;
    }

public:

    int num_level(){return levelset.size();}

    void get_level(std::vector<FermatEdge> &polygon, int level)
    {
        assert(0 <= level && level <= levelset.size());
        polygon = levelset[level];
    }

private:

    Settings settings;
    std::vector< std::vector<FermatEdge>> levelset;
};


#endif //SUPPORTER_FERMAT_SPIRALS_H
