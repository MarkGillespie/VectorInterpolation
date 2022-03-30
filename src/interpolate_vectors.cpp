#include "interpolate_vectors.h"

VertexData<Vector3>
interpolateByHarmonicFunction(ManifoldSurfaceMesh& mesh,
                              VertexPositionGeometry& geom,
                              const VertexData<Vector3>& boundaryData) {
    geom.requireCotanLaplacian();
    geom.requireVertexIndices();

    VertexData<size_t>& vIdx = geom.vertexIndices;

    Vector<bool> isInterior = Vector<bool>::Constant(mesh.nVertices(), true);
    size_t nB               = 0;
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isInterior(vIdx[v]) = false;
            nB++;
        }
    }

    BlockDecompositionResult<double> decomp =
        blockDecomposeSquare(geom.cotanLaplacian, isInterior);

    Vector<double> boundaryX(nB), boundaryY(nB), boundaryZ(nB);
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            Vector3 boundaryVec = boundaryData[v].normalize();

            boundaryX(decomp.newInds(vIdx[v])) = boundaryVec.x;
            boundaryY(decomp.newInds(vIdx[v])) = boundaryVec.y;
            boundaryZ(decomp.newInds(vIdx[v])) = boundaryVec.z;
        }
    }

    PositiveDefiniteSolver<double> solver(decomp.AA);

    Vector<double> interiorX = solver.solve(-decomp.AB * boundaryX);
    Vector<double> interiorY = solver.solve(-decomp.AB * boundaryY);
    Vector<double> interiorZ = solver.solve(-decomp.AB * boundaryZ);

    VertexData<Vector3> solution = boundaryData;
    for (Vertex v : mesh.vertices()) {
        if (isInterior(vIdx[v])) {
            solution[v] = Vector3{interiorX(decomp.newInds(vIdx[v])),
                                  interiorY(decomp.newInds(vIdx[v])),
                                  interiorZ(decomp.newInds(vIdx[v]))}
                              .normalize();
        }
    }

    geom.unrequireCotanLaplacian();
    geom.unrequireVertexIndices();

    return solution;
}

VertexData<Vector3> interpolateByConnectionLaplacian(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const VertexData<Vector3>& boundaryData, bool estimateNormalDirection) {
    geom.requireVertexConnectionLaplacian();
    geom.requireVertexTangentBasis();
    geom.requireVertexIndices();

    VertexData<size_t>& vIdx = geom.vertexIndices;

    Vector<bool> isInterior = Vector<bool>::Constant(mesh.nVertices(), true);
    size_t nB               = 0;
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isInterior(vIdx[v]) = false;
            nB++;
        }
    }

    BlockDecompositionResult<std::complex<double>> decomp =
        blockDecomposeSquare(geom.vertexConnectionLaplacian, isInterior);

    Vector<std::complex<double>> intrinsicBoundary(nB);
    Vector<double> boundaryNormals;
    if (estimateNormalDirection) boundaryNormals = Vector<double>(nB);
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            Vector3 extVec = boundaryData[v];
            Vector2 intrinsicVec{dot(extVec, geom.vertexTangentBasis[v][0]),
                                 dot(extVec, geom.vertexTangentBasis[v][1])};

            intrinsicBoundary(decomp.newInds(vIdx[v])) = intrinsicVec;
            if (estimateNormalDirection) {
                Vector3 N = cross(geom.vertexTangentBasis[v][0],
                                  geom.vertexTangentBasis[v][1]);
                boundaryNormals(decomp.newInds(vIdx[v])) = dot(extVec, N);
            }
        }
    }

    Vector<std::complex<double>> rhs = -decomp.AB * intrinsicBoundary;
    Vector<std::complex<double>> solution =
        solvePositiveDefinite(decomp.AA, rhs);

    Vector<double> solutionNormals;
    if (estimateNormalDirection) {
        geom.requireCotanLaplacian();
        BlockDecompositionResult<double> scalarDecomp =
            blockDecomposeSquare(geom.cotanLaplacian, decomp.isA);
        Vector<double> scalarRHS = -scalarDecomp.AB * boundaryNormals;
        solutionNormals = solvePositiveDefinite(scalarDecomp.AA, scalarRHS);
        geom.unrequireCotanLaplacian();
    }

    VertexData<Vector3> solutionData = boundaryData;
    for (Vertex v : mesh.vertices()) {
        if (isInterior(vIdx[v])) {
            Vector2 intrinsicVec =
                Vector2::fromComplex(solution(decomp.newInds(vIdx[v])));
            double z;

            if (estimateNormalDirection) {
                z              = solutionNormals(decomp.newInds(vIdx[v]));
                double vecNorm = sqrt(fmax(0, intrinsicVec.norm2() + z * z));
                intrinsicVec /= vecNorm;
                z /= vecNorm;
            } else {
                if (intrinsicVec.norm() >= 1) {
                    intrinsicVec = intrinsicVec.normalize();
                    z            = 0;
                } else {
                    double z2 = 1 - intrinsicVec.norm();
                    z         = sqrt(fmax(z2, 0));
                }
            }

            Vector3 T0 = geom.vertexTangentBasis[v][0];
            Vector3 T1 = geom.vertexTangentBasis[v][1];
            Vector3 N  = cross(T0, T1);
            solutionData[v] =
                intrinsicVec[0] * T0 + intrinsicVec[1] * T1 + z * N;
        }
    }

    geom.requireVertexConnectionLaplacian();
    geom.requireVertexTangentBasis();
    geom.unrequireVertexIndices();

    return solutionData;
}

VertexData<Vector3>
interpolateByStereographicProjection(ManifoldSurfaceMesh& mesh,
                                     VertexPositionGeometry& geom,
                                     const VertexData<Vector3>& boundaryData) {
    geom.requireCotanLaplacian();
    geom.requireVertexIndices();

    VertexData<size_t>& vIdx = geom.vertexIndices;

    Vector<bool> isInterior = Vector<bool>::Constant(mesh.nVertices(), true);
    size_t nB               = 0;
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isInterior(vIdx[v]) = false;
            nB++;
        }
    }

    BlockDecompositionResult<double> decomp =
        blockDecomposeSquare(geom.cotanLaplacian, isInterior);

    Vector<double> boundaryX(nB), boundaryY(nB);
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            Vector2 projectedVec = projectStereographic(boundaryData[v]);

            boundaryX(decomp.newInds(vIdx[v])) = projectedVec.x;
            boundaryY(decomp.newInds(vIdx[v])) = projectedVec.y;
        }
    }

    PositiveDefiniteSolver<double> solver(decomp.AA);

    Vector<double> interiorX = solver.solve(-decomp.AB * boundaryX);
    Vector<double> interiorY = solver.solve(-decomp.AB * boundaryY);

    VertexData<Vector3> solution = boundaryData;
    for (Vertex v : mesh.vertices()) {
        if (isInterior(vIdx[v])) {
            Vector2 vPos{interiorX(decomp.newInds(vIdx[v])),
                         interiorY(decomp.newInds(vIdx[v]))};
            solution[v] = unprojectStereographic(vPos);
        }
    }

    geom.unrequireCotanLaplacian();
    geom.unrequireVertexIndices();

    return solution;
}

VertexData<Vector3>
interpolateByHarmonicMapToSphere(ManifoldSurfaceMesh& mesh,
                                 VertexPositionGeometry& geom,
                                 const VertexData<Vector3>& boundaryData) {
    // just get some initial guess
    VertexData<Vector3> f =
        interpolateByHarmonicFunction(mesh, geom, boundaryData);

    geom.requireCotanLaplacian();
    VertexData<size_t>& vIdx = geom.vertexIndices;
    Vector<bool> isInterior  = Vector<bool>::Constant(mesh.nVertices(), true);
    size_t nB                = 0;
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isInterior(vIdx[v]) = false;
            nB++;
        }
    }

    BlockDecompositionResult<double> decomp =
        blockDecomposeSquare(geom.cotanLaplacian, isInterior);
    PositiveDefiniteSolver<double> Lsolve(decomp.AA);

    VertexData<Vector3> grad = computeSphericalDirichletGradient(mesh, geom, f);
    size_t steps             = 0;
    while (norm(grad) > 1e-8 && steps < 25) {
        VertexData<Vector3> stepDir = preconditionSphericalDirichletGradient(
            mesh, geom, f, grad, Lsolve, decomp);

        double stepSize = 5;

        VertexData<Vector3> newF =
            takeSphericalStep(mesh, geom, f, stepDir, stepSize);

        double oldEnergy = computeSphericalDirichletEnergy(mesh, geom, f);
        double newEnergy = computeSphericalDirichletEnergy(mesh, geom, newF);

        while (stepSize > 1e-8 &&
               (!std::isfinite(newEnergy) || oldEnergy < newEnergy)) {
            stepSize /= 2;

            newF      = takeSphericalStep(mesh, geom, f, stepDir, stepSize);
            newEnergy = computeSphericalDirichletEnergy(mesh, geom, newF);
        }

        f    = newF;
        grad = computeSphericalDirichletGradient(mesh, geom, f);

        steps++;
    }

    geom.unrequireCotanLaplacian();

    return f;
}

Vector2 projectStereographic(Vector3 v) {
    v = v.normalize();
    return Vector2{v.x, v.y} / (1 - v.z);
}

Vector3 unprojectStereographic(Vector2 v) {
    double denom = 1 + v.norm2();
    return Vector3{2 * v.x, 2 * v.y, v.norm2() - 1} / denom;
}

double computeSphericalDirichletEnergy(ManifoldSurfaceMesh& mesh,
                                       VertexPositionGeometry& geom,
                                       const VertexData<Vector3>& f) {
    geom.requireEdgeCotanWeights();

    double energy = 0;
    for (Edge e : mesh.edges()) {
        Vector3 fi = f[e.halfedge().tailVertex()];
        Vector3 fj = f[e.halfedge().tipVertex()];
        energy += geom.edgeCotanWeights[e] * pow(angle(fi, fj), 2) / 2.;
    }

    geom.unrequireEdgeCotanWeights();

    return energy;
}

VertexData<Vector3>
computeSphericalDirichletGradient(ManifoldSurfaceMesh& mesh,
                                  VertexPositionGeometry& geom,
                                  const VertexData<Vector3>& f) {
    geom.requireEdgeCotanWeights();

    VertexData<Vector3> grad(mesh, Vector3::zero());
    for (Vertex i : mesh.vertices()) {
        Vector3 fi = f[i];
        if (!i.isBoundary()) {
            for (Halfedge ij : i.outgoingHalfedges()) {
                Vector3 fj = f[ij.tipVertex()];

                // Normal vector of triangle containing the origin, fj, and fi
                Vector3 N = cross(fj, fi);

                // Tangent vector pointing from fi to fj
                Vector3 T = (-cross(N, fi));
                T         = T.normalize();
                if (T != T) {
                    T = Vector3::zero();
                }

                grad[i] -= geom.edgeCotanWeights[ij.edge()] * angle(fi, fj) * T;
            }
        }
    }

    geom.unrequireEdgeCotanWeights();

    return grad;
}

VertexData<Vector3> preconditionSphericalDirichletGradient(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const VertexData<Vector3>& f, const VertexData<Vector3>& grad,
    PositiveDefiniteSolver<double>& interiorLaplacianSolver,
    BlockDecompositionResult<double>& LaplacianDecomp) {

    size_t nInteriorVertices = LaplacianDecomp.origIndsA.rows();

    geom.requireVertexIndices();
    VertexData<size_t>& vIdx = geom.vertexIndices;

    //== Solve L * x = grad;
    Vector<double> gradX(nInteriorVertices), gradY(nInteriorVertices),
        gradZ(nInteriorVertices);

    for (Vertex v : mesh.vertices()) {
        if (!v.isBoundary()) {
            gradX(LaplacianDecomp.newInds(vIdx[v])) = grad[v].x;
            gradY(LaplacianDecomp.newInds(vIdx[v])) = grad[v].y;
            gradZ(LaplacianDecomp.newInds(vIdx[v])) = grad[v].z;
        }
    }

    Vector<double> stepX = interiorLaplacianSolver.solve(-gradX);
    Vector<double> stepY = interiorLaplacianSolver.solve(-gradY);
    Vector<double> stepZ = interiorLaplacianSolver.solve(-gradZ);

    //== Collect results and project onto tangent space of sphere

    VertexData<Vector3> step(mesh, Vector3::zero());
    for (Vertex v : mesh.vertices()) {
        if (!v.isBoundary()) {
            step[v] = Vector3{stepX(LaplacianDecomp.newInds(vIdx[v])),
                              stepY(LaplacianDecomp.newInds(vIdx[v])),
                              stepZ(LaplacianDecomp.newInds(vIdx[v]))};
            step[v] -= dot(step[v], f[v]) * f[v];
        }
    }

    geom.unrequireVertexIndices();

    return step;
}

VertexData<Vector3> takeSphericalStep(ManifoldSurfaceMesh& mesh,
                                      VertexPositionGeometry& geom,
                                      const VertexData<Vector3>& f,
                                      const VertexData<Vector3>& stepDir,
                                      double stepSize) {

    VertexData<Vector3> result = f;

    for (Vertex i : mesh.vertices()) {
        if (!i.isBoundary()) {
            Vector3 p    = f[i];
            Vector3 u    = stepDir[i] * stepSize;
            Vector3 axis = cross(p, u).normalize();
            double theta = u.norm();

            // rotate p around axis by theta
            result[i] = p * cos(theta) + cross(axis, p) * sin(theta) +
                        axis * dot(axis, p) * (1 - cos(theta));
        }
    }

    return result;
}

double norm(const VertexData<Vector3>& f) {
    double result = 0;
    for (size_t i = 0; i < f.size(); i++) {
        result += f[i].norm();
    }
    return result;
}

double dot(const VertexData<Vector3>& a, const VertexData<Vector3>& b) {
    double result = 0;
    for (size_t i = 0; i < a.size(); i++) {
        result += dot(a[i], b[i]);
    }
    return result;
}

VertexData<Vector3>
generateSmoothBoundaryVectorField(ManifoldSurfaceMesh& mesh,
                                  VertexPositionGeometry& geom) {
    VertexData<Vector3> boundaryData(mesh, Vector3::zero());
    geom.requireVertexNormals();
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex i : b.adjacentVertices()) {
            Vector3 Ni = geom.vertexNormals[i];
            boundaryData[i] =
                (geom.vertexPositions[i].normalize() + 2 * Ni).normalize();
        }
    }
    geom.unrequireVertexNormals();
    return boundaryData;
}
VertexData<Vector3> generateWavyBoundaryVectorField(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, size_t frequency) {
    VertexData<Vector3> boundaryData(mesh, Vector3::zero());
    geom.requireVertexNormals();
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        size_t iHe = 0;
        for (Halfedge ij : b.adjacentHalfedges()) {
            Vertex i   = ij.tailVertex();
            Vertex j   = ij.tipVertex();
            Vector3 Ni = geom.vertexNormals[i];
            Vector3 Ti =
                (geom.vertexPositions[j] - geom.vertexPositions[i]).normalize();
            Vector3 Bi = cross(Ti, Ni);

            double t        = iHe / ((double)b.degree()) * 2 * M_PI * frequency;
            boundaryData[i] = cos(t) * Ni + sin(t) * Bi;

            iHe++;
        }
    }
    geom.unrequireVertexNormals();
    return boundaryData;
}
