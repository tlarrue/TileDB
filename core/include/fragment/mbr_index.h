//
//  mbrIndex.h
//  rtree
//
//  Created by Tara Larrue on 5/1/17.
//  Copyright © 2017 Tara Larrue. All rights reserved.
//
// This file defines the class MbrIndex

#ifndef __MBR_INDEX_H__
#define __MBR_INDEX_H__

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <boost/geometry/geometries/adapted/boost_fusion.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <cmath>
#include <vector>
#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

/**@{*/
/** Return code. */
#define TILEDB_MI_OK          0
#define TILEDB_MI_ERR        -1
/**@}*/

/** Default error message. */
#define TILEDB_MI_ERRMSG std::string("[TileDB::MbrIndex] Error: ")

struct point2d {
    double x,y;
};

struct point3d{
    double x,y,z;
};

BOOST_FUSION_ADAPT_STRUCT(point2d, (double, x) (double, y))
BOOST_FUSION_ADAPT_STRUCT(point3d, (double, x) (double, y) (double, z))
BOOST_GEOMETRY_REGISTER_BOOST_FUSION_CS(cs::cartesian)


class MbrIndex {
    
    // 2-D types
    typedef bg::model::box<point2d> box2d;
    typedef std::pair<box2d, unsigned long> value2d;
    typedef bgi::rtree<value2d, bgi::quadratic<16>> pack_rtree2d;
    
    // 3-D types
    typedef bg::model::box<point3d> box3d;
    typedef std::pair<box3d, unsigned long> value3d;
    typedef bgi::rtree<value3d, bgi::quadratic<16>> pack_rtree3d;
    
public:
    /** Constructor. */
    MbrIndex();
    
    /** Public Functions */
    int buildTree(std::vector<void*> &mbrs, int dim);
    int intersectQuery(const double* subarray, std::vector<int64_t> &result);
    
private:
    /** Private Attributes */
    //void* rtree_;
    pack_rtree2d rtree2d_;
    pack_rtree3d rtree3d_;
    int dim_num_;
    
    /** Private Functions */
    void build2d(std::vector<void*> &mbrs);
    
    void build3d(std::vector<void*> &mbrs);
    
    void intersect2d(const double* subarray, std::vector<int64_t> &result);
    
    void intersect3d(const double* subarray, std::vector<int64_t> &result);
};


#endif 
