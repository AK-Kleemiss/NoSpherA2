#include "cube.h"
#include "isosurface.h"

// --------------------------------------------------------------------------
// 1) Minimal Edge Table (partial)
//    For each of the 256 possible bit configurations, edgeTable tells you 
//    which edges are intersected by the isosurface. Each bit in edgeTable[cubeIndex] 
//    corresponds to one of the 12 edges of the cube. If the bit is set, that edge is intersected.
//    Normally, you'll see all 256 entries. Here are just a few for illustration.
//
static const int edgeTable[256] = {
    /*  0 */ 0x0000, /*  1 */ 0x0109, /*  2 */ 0x0232, /*  3 */ 0x033B,
    /*  4 */ 0x0464, /*  5 */ 0x056D, /*  6 */ 0x0656, /*  7 */ 0x075F,
    // ...
    // (fill in the other 248 entries from standard references)
};

// --------------------------------------------------------------------------
// 2) Minimal Triangle Table (partial)
//    triTable[cubeIndex] is an array of up to 16 integers, grouped in triples. 
//    Each triple indicates the edge indices that form one triangle. A value of -1 
//    indicates the end of the list for that cubeIndex.
//
//    Example below just shows a tiny portion for illustration.
//
static const int triTable[256][16] = {
    /*  0 */ { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    /*  1 */ {  0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    /*  2 */ {  0,  1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    /*  3 */ {  1,  8,  3,  9,  8,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    // ...
    // (fill in the other 252 entries from standard references)
};

// --------------------------------------------------------------------------
// Edge to corner pairing: these pairs define which 2 corners (by index) form each edge.
// The standard corner numbering convention for a cube is:
//      (x, y, z)   -> corner 0
//      (x+1,y, z)  -> corner 1
//      (x+1,y+1,z) -> corner 2
//      (x, y+1, z) -> corner 3
//      (x, y, z+1) -> corner 4
//      (x+1,y,z+1) -> corner 5
//      (x+1,y+1,z+1)-> corner 6
//      (x, y+1,z+1)-> corner 7
//
// The typical mapping for edges is shown below:

static const int cornerIndexA[12] = {
    0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3
};
static const int cornerIndexB[12] = {
    1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7
};

// --------------------------------------------------------------------------
// A tiny helper to linearly interpolate a point along an edge between two corners
// based on the isosurface value. 
//
std::array<double, 3> interpolateIso(const std::array<double, 3>& p1, const std::array<double, 3>& p2, double valP1, double valP2, double isoVal)
{
    // If valP1 == valP2, avoid divide-by-zero; just return midpoint
    if (std::fabs(valP2 - valP1) < 1e-12) {
        return { 0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1]), 0.5 * (p1[2] + p2[2])};
    }
    double t = (isoVal - valP1) / (valP2 - valP1);
    return {
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2])
    };
}

// --------------------------------------------------------------------------
// Marching Cubes on a 3D array of double values to find isosurface = 0.5
//
std::vector<Triangle> marchingCubes(const cube& volumeData)
{
    std::vector<Triangle> triangles;

    int nx = static_cast<int>(volumeData.get_size(0));
    if (nx < 2) return triangles; // not enough data
    int ny = static_cast<int>(volumeData.get_size(1));
    if (ny < 2) return triangles; // not enough data
    int nz = static_cast<int>(volumeData.get_size(2));
    if (nz < 2) return triangles; // not enough data

    double isoVal = 0.5; // fixed isosurface level

    // We will iterate through each "voxel" (cube) in the volume
    for (int x = 0; x < nx - 1; x++) {
        for (int y = 0; y < ny - 1; y++) {
            for (int z = 0; z < nz - 1; z++) {

                // 1) Collect corner positions and values
                //    corners: (x,   y,   z)
                //             (x+1, y,   z)
                //             (x+1, y+1, z)
                //             (x,   y+1, z)
                //             (x,   y,   z+1)
                //             (x+1, y,   z+1)
                //             (x+1, y+1, z+1)
                //             (x,   y+1, z+1)
                std::array<double, 3> cornerPos[8] = {
                    volumeData.get_pos(x,y,z),
                    volumeData.get_pos(x + 1,y,z),
                    volumeData.get_pos(x + 1,y + 1,z),
                    volumeData.get_pos(x,y + 1,z),
                    volumeData.get_pos(x,y,z + 1),
                    volumeData.get_pos(x + 1,y,z + 1),
                    volumeData.get_pos(x + 1,y + 1,z + 1),
                    volumeData.get_pos(x,y + 1,z + 1)
                };

                double cornerVal[8] = {
                    volumeData.get_value(x,y,z),
                    volumeData.get_value(x + 1,y,z),
                    volumeData.get_value(x + 1,y + 1,z),
                    volumeData.get_value(x,y + 1,z),
                    volumeData.get_value(x,y,z + 1),
                    volumeData.get_value(x + 1,y,z + 1),
                    volumeData.get_value(x + 1,y + 1,z + 1),
                    volumeData.get_value(x,y + 1,z + 1)
                };

                // 2) Build the "cubeIndex" from which corners are >= isoVal
                int cubeIndex = 0;
                for (int i = 0; i < 8; i++) {
                    if (cornerVal[i] >= isoVal) {
                        cubeIndex |= (1 << i);
                    }
                }

                // If this cube is entirely below or entirely above isoVal => no intersection
                if (edgeTable[cubeIndex] == 0) {
                    continue;
                }

                // 3) Find the points where the isosurface intersects the edges of this cube
                std::array<double, 3> vertList[12];
                // For each of the 12 edges, check if it is intersected
                int edges = edgeTable[cubeIndex]; // bitmask of edges
                for (int i = 0; i < 12; i++) {
                    if (edges & (1 << i)) {
                        int cA = cornerIndexA[i];
                        int cB = cornerIndexB[i];
                        vertList[i] = interpolateIso(
                            cornerPos[cA], cornerPos[cB],
                            cornerVal[cA], cornerVal[cB],
                            isoVal
                        );
                    }
                }

                // 4) Create triangles from the vertList using the triTable
                //    Each set of three edge indices forms one triangle
                for (int i = 0; i < 16; i += 3) {
                    int e0 = triTable[cubeIndex][i + 0];
                    int e1 = triTable[cubeIndex][i + 1];
                    int e2 = triTable[cubeIndex][i + 2];

                    if (e0 == -1 || e1 == -1 || e2 == -1) {
                        break; // no more triangles for this cubeIndex
                    }

                    Triangle tri;
                    tri.v1 = vertList[e0];
                    tri.v2 = vertList[e1];
                    tri.v3 = vertList[e2];
                    triangles.push_back(tri);
                }
            }
        }
    }

    return triangles;
}