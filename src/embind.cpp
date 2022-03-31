#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"

#include "interpolate_vectors.h"

using namespace emscripten;
using namespace geometrycentral;
using namespace geometrycentral::surface;

struct GeoMesh {
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;
};

std::vector<Vector3> toList(ManifoldSurfaceMesh& mesh,
                            const VertexData<Vector3>& data) {
    std::vector<Vector3> vectorList;
    for (Vertex v : mesh.vertices()) {
        vectorList.push_back(data[v]);
    }
    return vectorList;
}

VertexData<Vector3> toData(ManifoldSurfaceMesh& mesh,
                           const std::vector<Vector3>& list) {
    VertexData<Vector3> data(mesh);
    for (size_t iV = 0; iV < mesh.nVertices(); iV++) {
        data[iV] = list[iV];
    }
    return data;
}

std::vector<Vector3> generateSmoothBoundaryField(GeoMesh& geo) {
    VertexData<Vector3> bField =
        generateSmoothBoundaryVectorField(*geo.mesh, *geo.geom);
    return toList(*geo.mesh, bField);
}

std::vector<Vector3> generateWavyBoundaryField(GeoMesh& geo, size_t frequency) {
    VertexData<Vector3> bField =
        generateWavyBoundaryVectorField(*geo.mesh, *geo.geom, frequency);
    return toList(*geo.mesh, bField);
}

std::vector<Vector3>
interpolateHarmonicFunction(GeoMesh& geo,
                            const std::vector<Vector3>& boundaryData) {
    return toList(*geo.mesh,
                  interpolateByHarmonicFunction(
                      *geo.mesh, *geo.geom, toData(*geo.mesh, boundaryData)));
}

std::vector<Vector3>
interpolateConnectionLaplacian(GeoMesh& geo,
                               const std::vector<Vector3>& boundaryData,
                               bool estimateNormalDirection) {
    return toList(*geo.mesh,
                  interpolateByConnectionLaplacian(
                      *geo.mesh, *geo.geom, toData(*geo.mesh, boundaryData),
                      estimateNormalDirection));
}

std::vector<Vector3>
interpolateStereographicProjection(GeoMesh& geo,
                                   const std::vector<Vector3>& boundaryData) {
    return toList(*geo.mesh,
                  interpolateByStereographicProjection(
                      *geo.mesh, *geo.geom, toData(*geo.mesh, boundaryData)));
}

std::vector<Vector3>
interpolateHarmonicMapToSphere(GeoMesh& geo,
                               const std::vector<Vector3>& boundaryData) {
    return toList(*geo.mesh,
                  interpolateByHarmonicMapToSphere(
                      *geo.mesh, *geo.geom, toData(*geo.mesh, boundaryData)));
}

// Stolen from Ricky Reusser https://observablehq.com/d/d0df0c04ce5c94FCC
template <typename T>
void copyToVector(const val& typedArray, std::vector<T>& vec) {
    unsigned int length = typedArray["length"].as<unsigned int>();
    val memory          = val::module_property("buffer");
    vec.reserve(length);
    val memoryView = typedArray["constructor"].new_(
        memory, reinterpret_cast<uintptr_t>(vec.data()), length);
    memoryView.call<void>("set", typedArray);
}

// Mostly stolen from Ricky Reusser https://observablehq.com/d/d0df0c04ce5c94fc
EMSCRIPTEN_BINDINGS(my_module) {
    value_array<Vector3>("Vector3")
        .element(&Vector3::x)
        .element(&Vector3::y)
        .element(&Vector3::z);
    value_array<Vector2>("Vector2").element(&Vector2::x).element(&Vector2::y);

    register_vector<Vector3>("VectorVector3");
    register_vector<size_t>("VectorSizeT");
    register_vector<std::vector<size_t>>("VectorVectorSizeT");

    class_<GeoMesh>("GCMesh")
        .function("polygons", optional_override([](GeoMesh& self) {
                      return self.mesh->getFaceVertexList();
                  }))
        .function("vertexCoordinates",
                  optional_override([](const GeoMesh& self) {
                      std::vector<Vector3> vCoords;
                      for (Vertex v : self.mesh->vertices())
                          vCoords.push_back(self.geom->inputVertexPositions[v]);
                      return vCoords;
                  }));

    function("readMesh",
             optional_override([](std::string str, std::string type = "") {
                 std::stringstream in;
                 in << str;

                 GeoMesh gMesh;
                 std::tie(gMesh.mesh, gMesh.geom) =
                     readManifoldSurfaceMesh(in, type);
                 return gMesh;
             }));

    function("generateSmoothBoundaryField", &generateSmoothBoundaryField);
    function("generateWavyBoundaryField", &generateWavyBoundaryField);

    function("interpolateHarmonicFunction", &interpolateHarmonicFunction);
    function("interpolateConnectionLaplacian", &interpolateConnectionLaplacian);
    function("interpolateStereographicProjection",
             &interpolateStereographicProjection);
    function("interpolateHarmonicMapToSphere", &interpolateHarmonicMapToSphere);
}
