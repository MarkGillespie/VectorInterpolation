#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

VertexData<Vector3>
interpolateByHarmonicFunction(ManifoldSurfaceMesh& mesh,
                              VertexPositionGeometry& geom,
                              const VertexData<Vector3>& boundaryData);

VertexData<Vector3>
interpolateByConnectionLaplacian(ManifoldSurfaceMesh& mesh,
                                 VertexPositionGeometry& geom,
                                 const VertexData<Vector3>& boundaryData,
                                 bool estimateNormalDirection = false);

VertexData<Vector3>
interpolateByStereographicProjection(ManifoldSurfaceMesh& mesh,
                                     VertexPositionGeometry& geom,
                                     const VertexData<Vector3>& boundaryData);

VertexData<Vector3>
interpolateByHarmonicMapToSphere(ManifoldSurfaceMesh& mesh,
                                 VertexPositionGeometry& geom,
                                 const VertexData<Vector3>& boundaryData);
VertexData<Vector3>
interpolateByHarmonicBundleSection(ManifoldSurfaceMesh& mesh,
                                   VertexPositionGeometry& geom,
                                   const VertexData<Vector3>& boundaryData);

//== Misc helpers
Vector2 projectStereographic(Vector3 v);
Vector3 unprojectStereographic(Vector2 v);

double computeSphericalDirichletEnergy(ManifoldSurfaceMesh& mesh,
                                       VertexPositionGeometry& geom,
                                       const VertexData<Vector3>& f);
double computeSphericalDirichletEnergy(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const VertexData<Vector3>& f,
    const HalfedgeData<Eigen::Matrix3d>& connection);

VertexData<Vector3>
computeSphericalDirichletGradient(ManifoldSurfaceMesh& mesh,
                                  VertexPositionGeometry& geom,
                                  const VertexData<Vector3>& f);
VertexData<Vector3> computeSphericalDirichletGradient(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const VertexData<Vector3>& f,
    const HalfedgeData<Eigen::Matrix3d>& connection);

VertexData<Vector3> preconditionSphericalDirichletGradient(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const VertexData<Vector3>& f, const VertexData<Vector3>& grad,
    PositiveDefiniteSolver<double>& interiorLaplacianSolver,
    BlockDecompositionResult<double>& LaplacianDecomp);

VertexData<Vector3> preconditionSphericalDirichletGradientWithConnection(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const VertexData<Vector3>& f, const VertexData<Vector3>& grad,
    PositiveDefiniteSolver<double>& connectionLaplacianSolver,
    BlockDecompositionResult<double>& LaplacianDecomp);

VertexData<Vector3> takeSphericalStep(ManifoldSurfaceMesh& mesh,
                                      VertexPositionGeometry& geom,
                                      const VertexData<Vector3>& f,
                                      const VertexData<Vector3>& stepDir,
                                      double stepSize);

double norm(const VertexData<Vector3>& f);
double dot(const VertexData<Vector3>& a, const VertexData<Vector3>& b);

VertexData<Vector3>
generateSmoothBoundaryVectorField(ManifoldSurfaceMesh& mesh,
                                  VertexPositionGeometry& geom);
VertexData<Vector3>
generateWavyBoundaryVectorField(ManifoldSurfaceMesh& mesh,
                                VertexPositionGeometry& geom,
                                size_t frequency = 1);

// Make a matrix whose columns are the input vectors
Eigen::Matrix3d toFrame(Vector3 v1, Vector3 v2, Vector3 v3);
Eigen::Vector3d toEigen(Vector3 v);
Vector3 fromEigen(Eigen::Vector3d v);
double angle(Eigen::Vector3d a, Eigen::Vector3d b);

SparseMatrix<double>
buildConnectionLaplacian(ManifoldSurfaceMesh& mesh,
                         VertexPositionGeometry& geom,
                         const HalfedgeData<Eigen::Matrix3d>& connection);

void checkSphericalDirichletGradient(ManifoldSurfaceMesh& mesh,
                                     VertexPositionGeometry& geom,
                                     VertexData<Vector3> boundaryData);
