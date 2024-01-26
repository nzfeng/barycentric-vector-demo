#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// == Polyscope data
polyscope::SurfaceMesh* psMesh;

std::tuple<Vector3, Vector3> boundingBox(VertexPositionGeometry& geometry) {

    SurfaceMesh& mesh = geometry.mesh;
    const double inf = std::numeric_limits<double>::infinity();
    Vector3 bboxMin = {inf, inf, inf};
    Vector3 bboxMax = {-inf, -inf, -inf};
    for (Vertex v : mesh.vertices()) {
        Vector3& pos = geometry.vertexPositions[v];
        for (int i = 0; i < 3; i++) {
            if (pos[i] <= bboxMin[i]) bboxMin[i] = pos[i];
            if (pos[i] >= bboxMax[i]) bboxMax[i] = pos[i];
        }
    }
    return std::make_tuple(bboxMin, bboxMax);
}

// void displayEdgeMidpointVectors(const EdgeData<BarycentricVector>& X, const std::string& name) {

//     Vector3 bboxMin, bboxMax;
//     std::tie(bboxMin, bboxMax) = boundingBox(*geometry);
//     double radius = (bboxMax - bboxMin).norm();
//     double charLength = radius / 50.;

//     TraceGeodesicResult tracedGeodesic;
//     TraceOptions traceOptions;
//     traceOptions.includePath = true;
//     std::vector<std::array<size_t, 2>> edgeInds;
//     std::vector<Vector3> positions;

//     double maxLength = 0.;
//     for (Edge e : mesh->edges()) {
//         maxLength = std::max(maxLength, X[e].norm(*geometry));
//     }

//     for (Edge e : mesh->edges()) {
//         Vector2 vec = Vector2::fromComplex(X[i]) / maxLength * charLength;
//         tracedGeodesic = traceGeodesic(*geometry, SurfacePoint(e, 0.5), vec, traceOptions);
//         std::vector<SurfacePoint>& pathPoints = tracedGeodesic.pathPoints;
//         int offset = positions.size();
//         for (SurfacePoint& pt : pathPoints) {
//             positions.push_back(pt.interpolate(geometry->vertexPositions));
//         }
//         for (size_t j = 0; j < pathPoints.size() - 1; j++) {
//             edgeInds.push_back({offset + j, offset + j + 1});
//         }
//     }
//     polyscope::registerCurveNetwork(name, positions, edgeInds);
// }

void testFunction() {

    std::cerr << "Testing..." << std::endl;

    // Put random vectors on edges.
    EdgeData<BarycentricVector> initVectors(*mesh);
    for (Edge e : mesh->edges()) {
        double t = unitRand();
        Vector2 edgeCoords = {t, -t};
        BarycentricVector w(e, edgeCoords);
        initVectors[e] = w;
    }

    // Rotate vectors.
    // Perform some sanity checks:
    //  - the length of rotated vectors has been preserved
    //  - the dot product between the original and rotated vectors should be 0
    //  - the dot product between the original and twice-rotated vectors should be -1
    //  - etc.
    EdgeData<BarycentricVector> vectorsRotated90(*mesh);
    EdgeData<BarycentricVector> vectorsRotated180(*mesh);
    EdgeData<BarycentricVector> vectorsRotated270(*mesh);
    EdgeData<BarycentricVector> vectorsRotated360(*mesh);
    double epsilon = 1e-5;
    for (Edge e : mesh->edges()) {
        // std::cerr << "Edge " << e << "\t" << mesh->nEdges() << std::endl;
        BarycentricVector w0 = initVectors[e];
        BarycentricVector w90 = w0.rotated90(*geometry);
        BarycentricVector w180 = w90.rotated90(*geometry);
        BarycentricVector w270 = w180.rotated90(*geometry);
        BarycentricVector w360 = w270.rotated90(*geometry);
        vectorsRotated90[e] = w90;
        vectorsRotated180[e] = w180;
        vectorsRotated270[e] = w270;
        vectorsRotated360[e] = w360;

        double origLength = w0.norm(*geometry);

        std::cerr << dot(*geometry, w0, w90) << std::endl;
        assert(std::abs(dot(*geometry, w0, w90)) < epsilon);

        std::cerr << dot(*geometry, w0, w180) << std::endl;
        assert(std::abs(dot(*geometry, w0, w180) + 1.) < epsilon);

        std::cerr << dot(*geometry, w0, w270) << std::endl;
        assert(std::abs(dot(*geometry, w0, w270)) < epsilon);

        std::cerr << dot(*geometry, w0, w360) << std::endl;
        assert(std::abs(dot(*geometry, w0, w360) - 1.) < epsilon);

        assert(std::abs(origLength - w90.norm(*geometry)) < epsilon);
        assert(std::abs(origLength - w180.norm(*geometry)) < epsilon);
        assert(std::abs(origLength - w270.norm(*geometry)) < epsilon);
        assert(std::abs(origLength - w360.norm(*geometry)) < epsilon);
    }

    // // Visualize the vectors, before and after -- just draw the barycentric vectors as single vectors emanating from
    // // edge midpoints, even though that's not really the case.
    // displayEdgeMidpointVectors(initVectors, "initial vectors");
    // displayEdgeMidpointVectors(vectorsRotated90, "rotated 90");
    // displayEdgeMidpointVectors(vectorsRotated180, "rotated 180");
    // displayEdgeMidpointVectors(vectorsRotated270, "rotated 270");
    // displayEdgeMidpointVectors(vectorsRotated360, "rotated 360");

    std::cerr << "Done testing." << std::endl;
}

void myCallback() {}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser(
        "A program for debugging additions to the BarycentricVector class in geometry-central.");
    args::Positional<std::string> meshFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    polyscope::init();

    polyscope::state::userCallback = myCallback;

    // Load mesh
    if (meshFilename) {
        std::string MESH_FILEPATH = args::get(meshFilename);
        std::string MESHNAME = polyscope::guessNiceNameFromPath(MESH_FILEPATH);
        std::tie(mesh, geometry) = readManifoldSurfaceMesh(MESH_FILEPATH);
        psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->vertexPositions, mesh->getFaceVertexList());
        psMesh->setAllPermutations(polyscopePermutations(*mesh));
    }

    testFunction();

    polyscope::show();

    return EXIT_SUCCESS;
}