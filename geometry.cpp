//
// Created by lurker on 5/12/17.
//

#include "geometry.h"

geometry::geometry() {
    /*
     * initialize all information on meta.
     *
     *  primitive types.
     */
    _meta.numberoftriangles = 0;
    _meta.numberofpointattributes = 0;
    _meta.numberofsegments = 0;
    _meta.numberoftriangleattributes = 0;
    _meta.numberofcorners = 0;
    _meta.numberofregions = 0;
    _meta.numberofedges = 0;
    _meta.numberofholes = 0;
    _meta.numberofpoints = 0;

    /*
     * initialize all pointers.
     */
    _meta.pointlist = nullptr;
    _meta.pointattributelist = nullptr;
    _meta.pointmarkerlist = nullptr;

    _meta.trianglelist = nullptr;
    _meta.triangleattributelist = nullptr;
    _meta.trianglearealist = nullptr;

    _meta.segmentlist = nullptr;
    _meta.segmentmarkerlist = nullptr;

    _meta.edgelist = nullptr;
    _meta.edgemarkerlist = nullptr;

    _meta.regionlist = nullptr;
    _meta.holelist = nullptr;
    _meta.normlist = nullptr;
    _meta.neighborlist = nullptr;
}

geometry::~geometry() {
    /*
     * release all information of the pointers.
     *
     * todo: use unique_ptr to describe the mesh.
     */
    trifree((void *) _meta.pointlist);
    trifree((void *) _meta.pointattributelist);
    trifree((void *) _meta.pointmarkerlist);

    trifree((void *) _meta.trianglelist);
    trifree((void *) _meta.triangleattributelist);

    trifree((void *) _meta.segmentlist);
    trifree((void *) _meta.segmentmarkerlist);;

    trifree((void *) _meta.edgelist);
    trifree((void *) _meta.edgemarkerlist);

    trifree((void *) _meta.regionlist);
    trifree((void *) _meta.holelist);
    trifree((void *) _meta.normlist);
    trifree((void *) _meta.neighborlist);
}

void geometry::set_points(Array<double> &_points) {
    auto len = _points.get_size();
    assert(len % 2 == 0);
    _meta.pointlist = (double *) malloc(len * sizeof(double));
    assert(_meta.pointlist != nullptr);
    _meta.numberofpoints = (int) (len / 2);
    memcpy(_meta.pointlist, _points.get_data(), len * sizeof(double));

}

void geometry::set_facets(Array<int> &_facets) {
    auto len = _facets.get_size();
    assert(len % 2 == 0);
    _meta.segmentlist = (int *) malloc(len * sizeof(int));
    assert(_meta.segmentlist != nullptr);
    _meta.numberofsegments = (int) (len / 2);
    memcpy(_meta.segmentlist, _facets.get_data(), len * sizeof(int));
}

void geometry::build(std::string switches, geometry &out) {
    geometry vor;
    triangulate((char *) switches.c_str(), &_meta, &out._meta, &vor._meta);
    out.setProperty();
}


void geometry::refine(std::string switches, geometry &out) {
    geometry vor;
    triangulate((char *) switches.c_str(), &_meta, &out._meta, &vor._meta);
    out.setProperty();
}

scalar_t geometry::getX(int index) {
    assert(index >= 0 && index < numberofpoints);
    return this->_meta.pointlist[2 * index];
}

scalar_t geometry::getY(int index) {
    assert(index >= 0 && index < numberofpoints);
    return this->_meta.pointlist[2 * index + 1];
}

scalar_t geometry::getArea(int index) {
    int u = this->_meta.trianglelist[3 * index];
    int v = this->_meta.trianglelist[3 * index + 1];
    int w = this->_meta.trianglelist[3 * index + 2];

    scalar_t det = fabs((this->_meta.pointlist[2 * v] - this->_meta.pointlist[2 * u]) *
                        (this->_meta.pointlist[2 * w + 1] - this->_meta.pointlist[2 * u + 1]) -
                        (this->_meta.pointlist[2 * v + 1] - this->_meta.pointlist[2 * u + 1]) *
                        (this->_meta.pointlist[2 * w] - this->_meta.pointlist[2 * u]));

    return det / 2.0;
}

int geometry::getTriangle(int index) {
    return this->_meta.trianglelist[index];
}


scalar_t geometry::getDistance(int fromIndex, int toIndex) {
    return sqrt(SQR(getX(fromIndex) - getX(toIndex) + SQR(getY(fromIndex) - getY(toIndex))));
}

void geometry::getNearField(vector<vector<int>> &near) {
    near.clear();
    near.resize((unsigned long) numberoftriangles);
    assert(this->_meta.neighborlist != nullptr);

    unordered_set<int> Hashset;
    queue<int> Queue;
    unordered_set<int> cur_triangle;

    for (index_t id = 0; id < numberoftriangles; ++id) {
        Hashset.clear();
        cur_triangle.clear();

        Queue.push(id);

        for (int pid = 0; pid < 3; ++pid) {
            cur_triangle.insert(getTriangle(3 * id + pid));
        }

        while (!Queue.empty()) {
            int top = Queue.front();
            Queue.pop();

            if (Hashset.find(top) == Hashset.end()) {
                Hashset.insert(top);

                for (int nid = 0; nid < 3; ++nid) {
                    int neighbor = _meta.neighborlist[3 * top + nid];
                    if (neighbor != -1) {
                        // check if it is adjacent.
                        bool adj = false;
                        for (int pid = 0; pid < 3; ++pid) {
                            if (cur_triangle.find(getTriangle(3 * neighbor + pid)) != cur_triangle.end()) {
                                adj = true;
                                break;
                            }
                        }

                        if (adj) {
                            Queue.push(neighbor);
                        }
                    }
                }

            }
        }


        for (int neigh : Hashset) {
            near[id].push_back(neigh);
        }
    }

}

void geometry::getAllArea() {
    arealist.resize((unsigned long) numberoftriangles);
    for (index_t id = 0; id < numberoftriangles; ++id) {
        arealist[id] = getArea(id);
    }

}
