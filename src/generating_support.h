//
// Created by 汪子琦 on 3/2/17.
//

#ifndef SUPPORTER_GENERATING_SUPPORT_H
#define SUPPORTER_GENERATING_SUPPORT_H

#include <limits>
#include <Eigen/Core>
#include <map>
#include <Eigen/StdVector>

#include "rendering_tree_support.h"
#include "detecting_overhangs.h"
#include "scene_organizer.h"
#include "settings.h"

class GeneratingSupport {
public:
    GeneratingSupport()
    {
        tmp_r = 0.5;

        platform_button_ = false;
    };

    ~GeneratingSupport()
    {

    };

public:

    void generate_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &sV, Eigen::MatrixXi &sF);


private:

    bool ConeMesh(Eigen::Vector3d p, Eigen::Vector3d &m, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    bool ConeCone(Eigen::Vector3d upper, Eigen::Vector3d lower, Eigen::Vector3d &point);

    Eigen::Vector3d ConeTriangle(Eigen::Vector3d cone, Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, int n, int k);

    Eigen::Vector3d PointTriangle(Eigen::Vector3d p, Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3);

private:

    inline double length(Eigen::Vector3d p){return p.norm();}

    inline double squared_length(Eigen::Vector3d p){return p.dot(p);}

    inline double dot(Eigen::Vector3d p, Eigen::Vector3d q){return p.dot(q);}

    inline Eigen::Vector3d cross(Eigen::Vector3d p, Eigen::Vector3d q){return p.cross(q);}

    void set_base( std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> &overhangs_);

    void set_platform(Eigen::MatrixXd &V);

    double get_platform_height(double x, double z);

private:

    void addBaseSupport(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iv);

    void addTrunk(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iu, unsigned long &iv);

    void addMeshSupport(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iv);

    void quick_sort(std::vector<Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> &vS, std::vector<unsigned long> &iS, int left, int right);

    void add_node(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iv);

private:

    double base_y_height_;

    double tmp_r;

public:

    Eigen::MatrixXd H_;

    bool platform_button_;
private:

    RenderingTreeSupport tree_;

    Settings settings_;

};

void GeneratingSupport::quick_sort(std::vector<Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> &vS,
                                    std::vector<unsigned long> &iS,
                                    int left,
                                    int right)
{
    int i = left, j = right;
    Eigen::Vector3d tmpV;
    unsigned long tmpI;
    double pivot = vS[(left + right) / 2](1);

    /* partition */
    while (i <= j) {
        while (vS[i](1) < pivot)
            i++;
        while (vS[j](1) > pivot)
            j--;
        if (i <= j) {
            //swap
            tmpV = vS[i];
            vS[i] = vS[j];
            vS[j] = tmpV;

            tmpI = iS[i];
            iS[i] = iS[j];
            iS[j] = tmpI;

            i++;
            j--;
        }
    }

    /* recursion */
    if (left < j)
        quick_sort(vS, iS, left, j);
    if (i < right)
        quick_sort(vS, iS, i, right);
}

void GeneratingSupport::set_platform(Eigen::MatrixXd &V)
{

    H_ = Eigen::MatrixXd::Zero(9, 11);
    Eigen::MatrixXi Flag = Eigen::MatrixXi::Zero(9, 11);
    for(size_t id = 0; id < V.rows() ;id ++)
    {
        int x_id = (int)(V(id, 0) / settings_.pad_size_);
        int z_id = (int)(V(id, 2) / settings_.pad_size_);
        if(H_(z_id,x_id) > V(id, 1) || Flag(z_id,x_id) == 0)
        {
            H_(z_id,x_id) = V(id, 1);
            Flag(z_id,x_id) = 1;
        }
    }
    return;
}

double GeneratingSupport::get_platform_height(double x, double z)
{
    int x_id = (int)(x / settings_.pad_size_);
    int z_id = (int)(z / settings_.pad_size_);
    return H_(z_id,x_id);
}

void GeneratingSupport::set_base( std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> &overhangs_)
{
    base_y_height_ = overhangs_[0](1);
    for(size_t kd = 0; kd < overhangs_.size(); kd++)
    {
        if(base_y_height_ > overhangs_[kd](1))
        {
            base_y_height_ = overhangs_[kd](1);
        }
    }

    std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> tmp_overhangs;
    for(size_t id = 0;id < overhangs_.size(); id++)
    {
        if(std::abs(overhangs_[id](1) - base_y_height_) > 0.1)
        {
            tmp_overhangs.push_back(overhangs_[id]);
        }
    }
    overhangs_ = tmp_overhangs;
    return;
}

void GeneratingSupport::generate_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &sV, Eigen::MatrixXi &sF)
{
    DetectingOverhangs overhangs_;
    std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> vS;
    overhangs_.sample(V, F, vS);
    std::vector<unsigned long>iS;

    set_base(vS);
    if(platform_button_)set_platform(V);
    iS.resize(vS.size(), ULONG_MAX);

    Eigen::Vector3d c0(0, 0, 0);

    bool isCC, isCM;
    std::cout << base_y_height_ << std::endl;
    while(!vS.empty())
    {
        //uniqueOverhangs(queue);
        quick_sort(vS, iS, 0, vS.size() - 1);
        std::cout << vS.size() << std::endl;
//      for(size_t kd = 0; kd < vS.size(); kd++)
//      {
//          std::cout << vS[kd].transpose() << ", " << iS[kd] << std::endl;
//      }

        isCC = false;
        isCM = false;

        Eigen::Vector3d c1(0, 0, 0);     //present highest apex point among these cones
        Eigen::Vector3d c2(0, 0, 0);     //the closest cone with cone c1
        Eigen::Vector3d c0(0, 0, 0);     //the intersection point of cone c1 and c0

        c1 = vS.back();
        double minDc1c0 = settings_.MAX_DOUBLE;
        size_t c1_index = vS.size() - 1, c2_index = 0;

        //Cone 2 Cone
        for(size_t id = 0; id < vS.size() - 1; id++)
        {
            Eigen::Vector3d tmpc0;
            //std::cout << "overhangs_[id]:\t" << overhang << std::endl;
            if(ConeCone(c1, vS[id], tmpc0))
            {
                isCC = true;
                if (length(c1 - tmpc0) < minDc1c0)
                {
                    c0 = tmpc0;
                    minDc1c0 = length(c1 - c0);
                    c2 = vS[id];
                    c2_index = id;
                }
            }
        }

//        std::cout <<acos((c1  - c0)(1) / (c1  - c0).norm())  / 3.1415926f * 180 <<
//                  "," << ( acos((c1  - c0)(1) / (c1  - c0).norm())  / 3.1415926f * 180 > 31 ) << std::endl;
//
//        std::cout << acos((c2  - c0)(1) / (c2  - c0).norm())  / 3.1415926f * 180 <<
//                  "," << ( acos((c2  - c0)(1) / (c2  - c0).norm())  / 3.1415926f * 180 > 31 ) << std::endl;

        //Cone 2 Mesh
        Eigen::Vector3d m0;      //the itersection point of mesh and cone c1
        isCM = ConeMesh(c1, m0, V, F);

        //Cone 2 Base
        Eigen::Vector3d b0;
        if(!platform_button_)
             b0 = Eigen::Vector3d(c1(0), base_y_height_, c1(2));
        else
             b0 = Eigen::Vector3d(c1(0), get_platform_height(c1(0),c1(2)) , c1(2));

//        std::cout << "c1:\t" << c1.transpose() << std::endl;
//        std::cout << "c2:\t" << c2.transpose() << std::endl;
//        std::cout << "m0:\t" << m0.transpose() << std::endl;
//        std::cout << "c0:\t" << c0.transpose() << std::endl;
//        std::cout << "b0:\t" << b0.transpose() << std::endl;

        if(isCM)
        {
            if(length(m0 - c1) > length(b0 - c1)) isCM = false;
            if(isCC && length(m0 - c1) > length(c0 - c1)) isCM = false;
        }

        if(isCC)
        {
            if(length(c0 - c1) > length(b0 - c1)) isCC = false;
            if(isCM && length(c0 - c1) > length(m0 - c1)) isCC = false;
        }


        unsigned long ic1, ic2;

        if(!isCC && !isCM)
        {
            ic1 = iS[c1_index];
            addBaseSupport(b0, c1, ic1);
            iS[c1_index] = ic1;
            vS.pop_back();
            iS.pop_back();
        }
        else if(isCM)
        {
            ic1 = iS[c1_index];
            addBaseSupport(m0, c1, ic1);
            iS[c1_index] = ic1;
            vS.pop_back();
            iS.pop_back();
        }
        else if(isCC)
        {
            unsigned long ic0 = ULONG_MAX - 1;
            ic1 = iS[c1_index];
            ic2 = iS[c2_index];
            addTrunk(c0, c1, ic0, ic1);
            addTrunk(c0, c2, ic0, ic2);
            iS[c1_index] = ic1;
            iS[c2_index] = ic2;

            vS.erase(vS.begin() + c2_index);
            iS.erase(iS.begin() + c2_index);
            vS.pop_back();
            iS.pop_back();
            vS.push_back(c0);
            iS.push_back(ic0);
        }
    }

    tree_.draw(sV, sF);
    return;
}

void GeneratingSupport::add_node(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iv)
{
    if(iv == ULONG_MAX)
    {
        unsigned iv0, iv1;
        iv0 = tree_.add_node(0.1 * u + 0.9 * v);
        iv1 = tree_.add_node(v);
        tree_.add_edge(iv0, iv1);
        iv = iv0;
    }
    else if(iv == ULONG_MAX - 1)
    {
        iv =  tree_.add_node(v);
    }
}

void GeneratingSupport::addBaseSupport(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iv)
{
    unsigned long iu = tree_.add_root(u);
    add_node(u, v, iv);
    tree_.add_edge(iu, iv);
    return;
}

void GeneratingSupport::addTrunk(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iu, unsigned long &iv)
{
    if(iu == ULONG_MAX || iu == ULONG_MAX - 1) {
        iu = tree_.add_node(u);
    }
    add_node(u, v, iv);
    tree_.add_edge(iu, iv);
    return;
}

void GeneratingSupport::addMeshSupport(Eigen::Vector3d &u, Eigen::Vector3d &v, unsigned long &iv)
{
    unsigned long iu = tree_.add_node(u);
    add_node(u, v, iv);
    tree_.add_edge(iu, iv);
    return;
}

bool  GeneratingSupport::ConeMesh(Eigen::Vector3d p, Eigen::Vector3d &m, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    bool flag = false;
    float minD = settings_.MAX_DOUBLE;
    m = Eigen::Vector3d(0,0,0);
    double Cangle = settings_.printing_max_angle;
    Eigen::Vector3d q[3];
    for(size_t fd = 0; fd < F.rows(); fd++)
    {
        q[0] = V.row(F(fd, 0)).transpose();
        q[1] = V.row(F(fd, 1)).transpose();
        q[2] = V.row(F(fd, 2)).transpose();
        if ( length(p - q[0]) < settings_.ZERO_EPS || length(p - q[1]) < settings_.ZERO_EPS || length(p - q[2]) < settings_.ZERO_EPS)
        {
            continue;                           //ignore the triangle containing cone vertex p
        }

        if ( q[0](1) > p(1) || q[1](1) > p(1) || q[2](1) > p(1))   //triangle is above point p
        {
            continue;
        }

        bool in_q[3] = {false, false ,false};
        for(size_t id = 0; id < 3; id++)
        {
            double cos_incl_angle = (p - q[id])(1) / length(q[id] - p);
            in_q[id] = (cos_incl_angle >= cos(Cangle));
        }
        if(in_q[0] && in_q[1] && in_q[2])
        {
            flag = true;
            Eigen::Vector3d mid = (q[0] + q[1] + q[2]) / 3;
            if(length(mid - p) < minD)
            {
                minD = length(mid - p);
                m = mid;
            }
        }
    }
    return flag;
}

bool GeneratingSupport::ConeCone(Eigen::Vector3d upper, Eigen::Vector3d lower, Eigen::Vector3d &point)
{
    Eigen::Vector3d p = upper; //p is equal to upper cone's apex point
    Eigen::Vector3d q = lower; //q is equal to lower cone's apex point
    double Cangle = settings_.printing_max_angle; //Cangle is equal to cone's apex angle

    double cos_incl_angle = (p - q)(1) / length(p - q);

    if (cos_incl_angle < cos(Cangle))
    {
        //q is out of the cone determined by p
        double R1 = (p(1) - base_y_height_) * tan(Cangle);
        double R2 = (q(1) - base_y_height_) * tan(Cangle);
        double d  = std::sqrt((p(0) - q(0)) * (p(0) - q(0))
                              +(p(2) - q(2)) * (p(2) - q(2)));

        if (R1 + R2 < d)             //the two cone have no intersection above the base line
        {
            point = Eigen::Vector3d(0, 0, 0);
            return false;
        }
        else if (R1 + R2 > d)       //the intersection of the two cones is above the base line
        {
            double h = 0.5 * (p(1) + q(1) - d / std::tan(Cangle) - 2 * base_y_height_);
            double y = base_y_height_ + h;

            double r1 = tan(Cangle) * (p(1) - y);
            double r2 = tan(Cangle) * (q(1) - y);

            double x = r2 * (p(0) - q(0)) / (r1 + r2) + q(0);
            double z = r2 * (p(2) - q(2)) / (r1 + r2) + q(2);
            point = Eigen::Vector3d(x, y, z);
//            std::cout << acos((p - point)(1)/length(p - point)) / PI_ * 180 << std::endl;
//            std::cout << acos((q - point)(1)/length(q - point)) / PI_ * 180 << std::endl;
            return true;
        }
        else                      //the intersection of the two cones is on the base line
        {
            double x = R2 * (p(0) - q(0))/(R1 + R2) + q(0);
            double z = R2 * (p(2) - q(2))/(R1 + R2) + q(2);
            double y = base_y_height_;
            point = Eigen::Vector3d(x, y, z);
            return true;
        }
    }
    else if (cos_incl_angle >= cos(Cangle))         //one cone contains the other, there is no intersection
    {
        point = Eigen::Vector3d(0, 0, 0);
        return false;
    }
    else                      //point q is on the surface of the cone determined by p (p is higher than q)
    {
        point = q;
        return true;
    }
}



//bool GeneratingSupport::ConeMesh(Eigen::Vector3d p, Eigen::Vector3d &m, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
//{
//
//    bool flag = false;           //flag is true if the cone is intersected with the mesh
//    float minD = MAX_FLOAT_;            //the shortest length from cone vertex p to intersection m, if there is a intersection
//    m = Eigen::Vector3d(0,0,0);
//    double Cangle = cone_apex_angle_; //Cangle is equal to cone's apex angle
//
//    Eigen::Vector3d q[3];
//
//    for(size_t fd = 0; fd < F.rows(); fd++)
//    {
//        q[0] = V.row(F(fd, 0)).transpose();
//        q[1] = V.row(F(fd, 1)).transpose();
//        q[2] = V.row(F(fd, 2)).transpose();
//
//        if ( length(p - q[0]) < settings_.ZERO_EPS || length(p - q[1]) < settings_.ZERO_EPS || length(p - q[2]) < settings_.ZERO_EPS)
//        {
//            continue;                           //ignore the triangle containing cone vertex p
//        }
//
//        if (!(q[0](1) < p(1)) || !(q[1](1) < p(1)) || !(q[2](1) < p(2)))   //triangle is above point p
//        {
//            continue;
//        }
//
//        float isq[3];            //judge q[i] is in cone or not, q[i] is in the cone if isq[i] is greater than 0
//        for (size_t id = 0; id < 3; id++)
//        {
//            isq[id] = std::sin(Cangle) * std::sin(Cangle) * ( q[id](2) - p(2)) * (q[id](2) - p(2))
//                      -std::cos(Cangle) * std::cos(Cangle) * ((q[id](0) - p(0)) * (q[id](0) - p(0)) +(q[id](1) - p(1)) * (q[id](1) - p(1)));
//        }
//
//        if ((isq[0] < 0 && isq[1] < 0 && isq[2] < 0))         //triangle q1q2q3 is out of the cone, there is no intersection
//        {
//            continue;
//        }
//        else if ((isq[0] * isq[1] * isq[2]) < settings_.ZERO_EPS)
//        {
//            for (size_t id = 0; id <3; id ++)
//            {
//                if ( std::abs(isq[id]) < settings_.ZERO_EPS && length(q[id] - p) < minD)   //q[i] is the intersection of cone to mesh ƒ⁄ or Õ‚£ø£ø£ø
//                {
//                    m = q[id];
//                    minD = length(q[id] - p);
//                    //cout<<fixed<<setprecision(6)<<minD<<"  1111  "<<m<<endl;
//                    flag = true;
//                }
//            }
//        }
//        else if ((isq[0] > 0  && isq[1] > 0 && isq[2] > 0))    //triangle q1q2q3 is in the cone, there exist intersections
//        {
//            if (length(p - PointTriangle(p, q[0], q[1], q[2])) < minD)
//            {
//                m = PointTriangle(p, q[0], q[1], q[2]);     //m is the nearest point in triangle q1q2q3 from point p
//                minD = length(p - m);
//                //cout<<fixed<<setprecision(6)<<minD<<"  2222  "<<m<<endl;
//                flag = true;
//            }
//        }
//        else
//        {
//            Eigen::Vector3d p1, p2, p3;
//            if (isq[0] * isq[1] * isq[2] > 0)      // 1 in and 2 out, sign(+, -, -)
//            {
//                size_t id;
//                for (id = 0; id <3; id ++)
//                {
//                    if (isq[id] > 0)
//                    {
//                        p1 = q[id];
//                        break;
//                    }
//                }
//                p2 = q[(id + 1) % 3];
//                p3 = q[(id + 2) % 3];
//            }
//            else if (isq[0] * isq[1] * isq[2] < 0)   // 2 in and 1 out, sign(+, +, -)
//            {
//                size_t id;
//                for (id = 0; id < 3; id ++)
//                {
//                    if (isq[id] < 0)
//                    {
//                        p1 = q[id];
//                        break;
//                    }
//                }
//                p2 = q[(id + 1) % 3];
//                p3 = q[(id + 2) % 3];
//            }
//
//            int n = 5;                  //sample distance, divide the edge of triangle n segments
//            Eigen::Vector3d m0;         //sample point
//            for (int kd = 0; kd < 6; kd ++)
//            {
//                //cout<<i<<"    ";
//                m0 = ConeTriangle(p, p1 , p2, p3 , n, kd);     //line to cone intersection
//                if (length(p - m0) < minD)
//                {
//                    m = m0;
//                    minD = length(p - m0);
//                    //cout<<fixed<<setprecision(6)<<minD<<"  3333  "<<m<<endl;
//                    flag = true;
//                }
//            }
//            for (size_t id = 0; id < 3; id++)
//            {
//                if (isq[id]>0 && length(p - q[id]) < minD)
//                {
//                    minD = length(p - q[id]);
//                    m = q[id];
//                    flag = true;
//                }
//            }
//        }
//
//    }
//    return flag;
//}

//Eigen::Vector3d GeneratingSupport::ConeTriangle(Eigen::Vector3d cone,
//                                                 Eigen::Vector3d p1,
//                                                 Eigen::Vector3d p2,
//                                                 Eigen::Vector3d p3,
//                                                 int n,
//                                                 int k)
//{
//    double Cangle = cone_apex_angle_; //Cangle is equal to cone's apex angle
//
//    Eigen::Vector3d V = cone;
//    Eigen::Vector3d X;
//
//    double lamda = k * 1.0 / n;
//    Eigen::Vector3d p1p0 = p2 - p1 + lamda * (p3 - p2);
//    Eigen::Vector3d p0   = p1 + p1p0;
//    Eigen::Vector3d p    = 0.5f * (p1 + p0);
//    Eigen::Vector3d U    = (p0 - p1) / length(p0 - p1);
//    double e = 0.5f * length(p0 - p1);
//
//    Eigen::Vector3d delta = p - V;
//    Eigen::Vector3d M1(-cos(Cangle) * cos(Cangle) * delta(0),
//                       -cos(Cangle) * cos(Cangle) * delta(1),
//                        sin(Cangle) * sin(Cangle) * delta(2));
//    Eigen::Vector3d M2(-cos(Cangle) * cos(Cangle) * U(0),
//                       -cos(Cangle) * cos(Cangle) * U(1),
//                        sin(Cangle) * sin(Cangle) * U(2));
//
//    double c0 = dot(delta , M1);
//    double c1 = dot(U     , M1);
//    double c2 = dot(U     , M2);
//
//    double test = c1 * c1 - c0 * c2;
//    double t1, t2;
//
//    if (test>0)
//    {
//        t1 = (-c1 + std::sqrt(test)) / c2;
//        t2 = (-c1 - std::sqrt(test)) / c2;
//    }
//    if (!((t1 - e) * (t1 + e) > 0))
//    {
//        //cout<<"t1 :  "<<t1<<endl;
//        X = p + t1 * U;
//    }
//    if (!((t2 - e) * (t2 + e) > 0))
//    {
//        //cout<<"t2 :  "<<t2<<endl;
//        X = p + t2 * U;
//    }
//
//    return X;
//
//}
//
//Eigen::Vector3d GeneratingSupport::PointTriangle(Eigen::Vector3d p,
//                                                Eigen::Vector3d p1,
//                                                Eigen::Vector3d p2,
//                                                Eigen::Vector3d p3)
//{
//
//    Eigen::Vector3d p2p3_p1, p1p2_p, p1p3_p, p3p2_p;
//
//    p1p2_p  = cross((p1 - p), (p2 - p));
//    p1p3_p  = cross((p1 - p), (p3 - p));
//    p3p2_p  = cross((p3 - p), (p2 - p));
//    p2p3_p1 = cross(p2 - p1, p3 - p1);
//
//    double   d = dot(p2p3_p1, p1 - p) / length(p2p3_p1);
//    double   t = -d * d / dot(p2p3_p1, p1 - p);
//    Eigen::Vector3d q = p + p2p3_p1 * t;
//
//    double isin = length(p1p2_p) + length(p1p3_p) + length(p3p2_p) - length(p2p3_p1);
//
//    if (std::abs(isin)< settings_.ZERO_EPS)
//    {
//        return q;           //q is in the triangle
//    }
//    else                    //q is out of the triangle
//    {
//        Eigen::Vector3d nearest_p, A, B;
//        if (length(q - p1) < length(q - p2))
//        {
//            if (length(q - p1) < length(q-p3))
//            {
//                nearest_p = p1;
//                A = p2;
//                B = p3;
//            }
//            else
//            {
//                nearest_p = p3;
//                A = p1;
//                B = p2;
//            }
//        }
//        else
//        {
//            if (length(q - p2) < length(q-p3))
//            {
//                nearest_p = p2;
//                A = p1;
//                B = p3;
//            }
//            else
//            {
//                nearest_p = p3;
//                A = p1;
//                B = p2;
//            }
//        }
//
//        if ( !( dot(A - nearest_p, q - nearest_p)>0 ) && !( dot(B - nearest_p, q - nearest_p)<0 ) )
//        {
//            return nearest_p;
//        }
//        else
//        {
//            double dA = dot(A - nearest_p, q - nearest_p) / length(A - nearest_p);
//            double dB = dot(B - nearest_p, q - nearest_p) / length(B - nearest_p);
//            if (dA < dB)
//            {
//                double s = (squared_length(q - nearest_p) - dA * dA) / length(A - nearest_p);
//                return s * A + nearest_p * (1 - s);
//            }
//            else
//            {
//                double s = (squared_length(q - nearest_p) - dB * dB) / length(B - nearest_p);
//                return s * B + nearest_p * (1 - s);
//            }
//        }
//    }
//}


#endif