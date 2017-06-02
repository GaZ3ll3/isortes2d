/*
 * generate triangle mesh for unit square.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu. All Rights Reserved.
 *
 */
#include "utils.h"
#include "geometry_triangle.h"


extern void reorg(geometry_triangle *in);

extern void domain_information(geometry_triangle *in);

geometry_triangle *build_unit_square_mesh(int M) {

    Array<scalar_t> boundary_points((size_t) (2 * SQR(M + 1)));
    Array<index_t> edge_connection((size_t) (6 * SQR(M) + 4 * M));
    scalar_t dx = 1.0 / (scalar_t) M;

    int bid = 0, eid = 0;

    for (int _row = 0; _row < M + 1; ++_row) {
        for (int _col = 0; _col < M + 1; ++_col) {
            boundary_points[bid++] = _row * dx;
            boundary_points[bid++] = _col * dx;
        }
    }

    for (int _row = 0; _row < M; ++_row) {
        for (int _col = 0; _col < M + 1; ++_col) {
            edge_connection[eid++] = _row * (M + 1) + _col;
            edge_connection[eid++] = (_row + 1) * (M + 1) + _col;
        }
    }

    for (int _col = 0; _col < M; ++_col) {
        for (int _row = 0; _row < M + 1; ++_row) {
            edge_connection[eid++] = _row * (M + 1) + _col;
            edge_connection[eid++] = _row * (M + 1) + _col + 1;
        }
    }

    for (int _row = 0; _row < M; ++_row) {
        for (int _col = 0; _col < M; ++_col) {
            edge_connection[eid++] = _row * (M + 1) + (_col + 1);
            edge_connection[eid++] = (_row + 1) * (M + 1) + _col;
        }
    }

    geometry_triangle *setting = new geometry_triangle();
    setting->set_points(boundary_points);
    setting->set_facets(edge_connection);
    geometry_triangle *domain = new geometry_triangle();
    setting->build("pznjQ", *domain);
    free(setting);

    reorg(domain);


    domain->getAllArea();
    domain_information(domain);

    return domain;
}

extern void domain_information(geometry_triangle *in) {
#ifdef DISP
    std::cout << std::setw(15) << "vertices" << std::setw(8) << in->_meta.numberofpoints << std::endl;
    std::cout << std::setw(15) << "triangles" << std::setw(8) << in->_meta.numberoftriangles << std::endl;
    std::cout << std::setw(15) << "edges" << std::setw(8) << in->_meta.numberofedges << std::endl;
    std::cout << std::setw(15) << "segments" << std::setw(8) << in->_meta.numberofsegments << std::endl;
    std::cout << std::endl;
#endif
}

extern void reorg(geometry_triangle *in) {
    /*
     * traverse all triangles and reorder.
     */
    scalar_t x1, y1, x2, y2, x3, y3;
    scalar_t angle1, angle2, angle3; // cosine value
    scalar_t s1, s2, s3; // side


    for (int tid = 0; tid < in->numberoftriangles; ++tid) {
        int u = in->getTriangle(3 * tid);
        int v = in->getTriangle(3 * tid + 1);
        int w = in->getTriangle(3 * tid + 2);

        x1 = in->getX(u);
        y1 = in->getY(u);
        x2 = in->getX(v);
        y2 = in->getY(v);
        x3 = in->getX(w);
        y3 = in->getY(w);

        s1 = sqrt(SQR(x2 - x3) + SQR(y2 - y3));
        s2 = sqrt(SQR(x1 - x3) + SQR(y1 - y3));
        s3 = sqrt(SQR(x1 - x2) + SQR(y1 - y2));


        angle1 = (SQR(s2) + SQR(s3) - SQR(s1)) / (2 * s2 * s3);
        angle2 = (SQR(s1) + SQR(s3) - SQR(s2)) / (2 * s1 * s3);
        angle3 = (SQR(s1) + SQR(s2) - SQR(s3)) / (2 * s1 * s2);

        if (fabs(angle1) < EPS) {
            in->_meta.trianglelist[3 * tid] = v;
            in->_meta.trianglelist[3 * tid + 1] = w;
            in->_meta.trianglelist[3 * tid + 2] = u;
        } else if (fabs(angle2) < EPS) {
            in->_meta.trianglelist[3 * tid] = w;
            in->_meta.trianglelist[3 * tid + 1] = u;
            in->_meta.trianglelist[3 * tid + 2] = v;
        }
    }
}