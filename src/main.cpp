#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "interpolate_vectors.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;
VertexData<Vector3> boundaryData;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Interpolate by harmonic function")) {
        VertexData<Vector3> interpolatedNormals =
            interpolateByHarmonicFunction(*mesh, *geom, boundaryData);
        psMesh->addVertexVectorQuantity("harmonic interpolation",
                                        interpolatedNormals);
    }
    if (ImGui::Button("Interpolate by connection laplacian")) {
        VertexData<Vector3> interpolatedNormals =
            interpolateByConnectionLaplacian(*mesh, *geom, boundaryData, false);
        psMesh->addVertexVectorQuantity(
            "pure connection laplacian interpolation", interpolatedNormals);
    }
    if (ImGui::Button(
            "Interpolate by connection laplacian + scalar laplacian")) {
        VertexData<Vector3> interpolatedNormals =
            interpolateByConnectionLaplacian(*mesh, *geom, boundaryData, true);
        psMesh->addVertexVectorQuantity(
            "connection laplacian + scalar laplacian interpolation",
            interpolatedNormals);
    }
    if (ImGui::Button("Interpolate by stereographic projection")) {
        VertexData<Vector3> interpolatedNormals =
            interpolateByStereographicProjection(*mesh, *geom, boundaryData);
        psMesh->addVertexVectorQuantity("stereographic interpolation",
                                        interpolatedNormals);
    }
    if (ImGui::Button("Interpolate by harmonic map to sphere")) {
        VertexData<Vector3> interpolatedNormals =
            interpolateByHarmonicMapToSphere(*mesh, *geom, boundaryData);
        psMesh->addVertexVectorQuantity("harmonic sphere map",
                                        interpolatedNormals);
    }
    ImGui::Separator();
    if (ImGui::Button("Set smooth boundary conditions")) {
        boundaryData = VertexData<Vector3>(*mesh, Vector3::zero());
        geom->requireVertexNormals();
        for (BoundaryLoop b : mesh->boundaryLoops()) {
            for (Vertex i : b.adjacentVertices()) {
                Vector3 Ni = geom->vertexNormals[i];
                boundaryData[i] =
                    (geom->vertexPositions[i].normalize() + 2 * Ni).normalize();
            }
        }
        geom->unrequireVertexNormals();
        psMesh->addVertexVectorQuantity("boundary data", boundaryData);
    }
    if (ImGui::Button("Set wavy boundary conditions")) {
        boundaryData = VertexData<Vector3>(*mesh, Vector3::zero());
        geom->requireVertexNormals();
        for (BoundaryLoop b : mesh->boundaryLoops()) {
            size_t iHe = 0;
            for (Halfedge ij : b.adjacentHalfedges()) {
                Vertex i   = ij.tailVertex();
                Vertex j   = ij.tipVertex();
                Vector3 Ni = geom->vertexNormals[i];
                Vector3 Ti =
                    (geom->vertexPositions[j] - geom->vertexPositions[i])
                        .normalize();
                Vector3 Bi = cross(Ti, Ni);

                // boundaryData[i] = (geom->vertexPositions[i].normalize() + 2 *
                // Ni).normalize();
                double t        = iHe / ((double)b.degree()) * 2 * M_PI;
                boundaryData[i] = cos(t) * Ni + sin(t) * Bi;

                iHe++;
            }
        }
        geom->unrequireVertexNormals();
        psMesh->addVertexVectorQuantity("boundary data", boundaryData);
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Geometry program");
    args::Positional<std::string> inputFilename(parser, "mesh",
                                                "Mesh to be processed.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = "../../meshes/bunny_small.obj";
    // Make sure a mesh name was given
    if (inputFilename) {
        filename = args::get(inputFilename);
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);

    // Center mesh
    Vector3 avgPos = Vector3::zero();
    for (Vertex v : mesh->vertices()) {
        avgPos += geom->vertexPositions[v];
    }
    avgPos /= mesh->nVertices();
    for (Vertex v : mesh->vertices()) {
        geom->vertexPositions[v] -= avgPos;
    }

    // Pick some weird boundary data
    boundaryData = VertexData<Vector3>(*mesh, Vector3::zero());
    geom->requireVertexNormals();
    for (BoundaryLoop b : mesh->boundaryLoops()) {
        for (Vertex i : b.adjacentVertices()) {
            Vector3 Ni = geom->vertexNormals[i];
            boundaryData[i] =
                (geom->vertexPositions[i].normalize() + 2 * Ni).normalize();
        }
    }
    geom->unrequireVertexNormals();

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh("mesh", geom->vertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));
    psMesh->addVertexVectorQuantity("boundary data", boundaryData);

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
