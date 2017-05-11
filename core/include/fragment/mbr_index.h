//
//  mbr_index.h
//  rtree_test
//
//  Created by Tara Larrue on 5/7/17.
//  Copyright © 2017 Tara Larrue. All rights reserved.
//

#ifndef __MBR_INDEX_H__
#define __MBR_INDEX_H__

#include <stdio.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cmath>
#include <vector>
#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

/* ********************************* */
/*             CONSTANTS             */
/* ********************************* */

/**@{*/
/** Return code. */
#define TILEDB_MI_OK          0
#define TILEDB_MI_ERR        -1
/**@}*/

/** Default error message. */
#define TILEDB_MI_ERRMSG std::string("[TileDB::MbrIndex] Error: ")


class MbrIndex {
    
    // 2-D types
    typedef bg::model::point<double, 2, bg::cs::cartesian> point2d;
    typedef bg::model::box<point2d> box2d;
    typedef std::pair<box2d, unsigned long> value2d;
    typedef bgi::rtree<value2d, bgi::rstar<16,4>> pack_rtree2d;
    
    // 3-D types
    typedef bg::model::point<double, 3, bg::cs::cartesian> point3d;
    typedef bg::model::box<point3d> box3d;
    typedef std::pair<box3d, unsigned long> value3d;
    typedef bgi::rtree<value3d, bgi::rstar<16,4>> pack_rtree3d;
    
public:
    /** Constructor. */
    MbrIndex();
    
    /** Public Functions */
    int buildTree(std::vector<void*> &mbrs, int dim);
    int intersectQuery(const double* subarray, std::vector<int64_t> &result);
    
private:
    /** Private Attributes */
    pack_rtree2d rtree2d_;
    pack_rtree3d rtree3d_;
    int dim_num_;
    
    /** Private Functions */
    void build2d(std::vector<void*> &mbrs);
    void build3d(std::vector<void*> &mbrs);
    void intersect2d(const double* subarray, std::vector<int64_t> &result);
    void intersect3d(const double* subarray, std::vector<int64_t> &result);
    bool testTree2d(std::vector<void*> &mbrs);
};

#endif
