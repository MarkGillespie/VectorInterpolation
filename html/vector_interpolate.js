import {
  Vector3,
  Matrix4,
} from "https://unpkg.com/three@0.125.1/build/three.module.js";

// import { Geoptic } from "./geoptic.js/src/geoptic.js";
import { Geoptic } from "./geoptic.js/build/geoptic.module.min.js";
// import { standardizeVector3Array } from "./geoptic.js/src/standardize_data_array.js";

let geo = undefined;
let boundaryField = undefined;
let psBaseMesh = undefined;

// create geoptic manager
let geoptic = new Geoptic({
  // parent: document.getElementById("geoptic-panel"),
  path: "geoptic.js",
  picks: false,
});
geoptic.setGroundPlaneEnabled(false);

geoptic.message("waiting for webassembly to load");

if (Module.runtimeInitialized) {
  // Once the wasm has loaded, we can start our app
  geoptic.message("webassembly loaded");

  // Load the meshes and set up our state
  init(hand);

  // Start animating with geoptic
  // This will call geoptic.userCallback() every frame
  geoptic.animate();
} else {
  Module.onRuntimeInitialized = (_) => {
    // Once the wasm has loaded, we can start our app
    geoptic.message("webassembly loaded");

    // Load the meshes and set up our state
    init(hand);

    // Start animating with geoptic
    // This will call geoptic.userCallback() every frame
    geoptic.animate();
  };
}

function vec3ToTHREE(v) {
  return new Vector3(v[0], v[1], v[2]);
}

// Set up UI panel
let io = geoptic.commandGui.addFolder("IO");
geoptic.commandGuiFields["Load New Base Mesh"] = function () {
  geoptic.loadMesh((text) => {
    geo.delete();
    geo = Module.readMesh(text, "obj");

    psBaseMesh = geoptic.registerSurfaceMesh(
      "Base Mesh",
      geo.vertexCoordinates(),
      geo.polygons()
    );

    boundaryField = Module.generateSmoothBoundaryField(geo);
    psBaseMesh.addVertexVectorQuantity("Boundary field", boundaryField);
  });
};
io.add(geoptic.commandGuiFields, "Load New Base Mesh");
io.close();

let boundaryIO = geoptic.commandGui.addFolder("Boundary");
geoptic.commandGuiFields["Generate smooth boundary field"] = function () {
  boundaryField = Module.generateSmoothBoundaryField(geo);
  psBaseMesh.addVertexVectorQuantity("Boundary field", boundaryField);
};
boundaryIO.add(geoptic.commandGuiFields, "Generate smooth boundary field");
geoptic.commandGuiFields["Generate wavy boundary field"] = function () {
  boundaryField = Module.generateWavyBoundaryField(geo, 1);
  psBaseMesh.addVertexVectorQuantity("Boundary field", boundaryField);
};
boundaryIO.add(geoptic.commandGuiFields, "Generate wavy boundary field");
boundaryIO.open();

let interpolationIO = geoptic.commandGui.addFolder("Interpolation");
geoptic.commandGuiFields["interpolate by harmonic function"] = function () {
  let interpolant = Module.interpolateHarmonicFunction(geo, boundaryField);
  psBaseMesh.addVertexVectorQuantity("interpolant (harmonic)", interpolant);
};
interpolationIO.add(
  geoptic.commandGuiFields,
  "interpolate by harmonic function"
);
geoptic.commandGuiFields["interpolate by connection Laplacian"] = function () {
  let interpolant = Module.interpolateConnectionLaplacian(
    geo,
    boundaryField,
    false
  );
  psBaseMesh.addVertexVectorQuantity(
    "interpolant (connection Laplacian)",
    interpolant
  );
};
interpolationIO.add(
  geoptic.commandGuiFields,
  "interpolate by connection Laplacian"
);
geoptic.commandGuiFields[
  "interpolate by connection Laplacian + scalar Laplacian"
] = function () {
  let interpolant = Module.interpolateConnectionLaplacian(
    geo,
    boundaryField,
    true
  );
  psBaseMesh.addVertexVectorQuantity(
    "interpolant (connection Laplacian++)",
    interpolant
  );
};
interpolationIO.add(
  geoptic.commandGuiFields,
  "interpolate by connection Laplacian + scalar Laplacian"
);
geoptic.commandGuiFields["interpolate by stereographic projection"] =
  function () {
    let interpolant = Module.interpolateStereographicProjection(
      geo,
      boundaryField
    );
    psBaseMesh.addVertexVectorQuantity(
      "interpolant (stereographic projection)",
      interpolant
    );
  };
interpolationIO.add(
  geoptic.commandGuiFields,
  "interpolate by stereographic projection"
);
geoptic.commandGuiFields["interpolate by harmonic map to sphere"] =
  function () {
    let interpolant = Module.interpolateHarmonicMapToSphere(geo, boundaryField);
    psBaseMesh.addVertexVectorQuantity(
      "interpolant (harmonic map to sphere)",
      interpolant
    );
  };
interpolationIO.add(
  geoptic.commandGuiFields,
  "interpolate by harmonic map to sphere"
);
interpolationIO.open();

function init(text) {
  if (geo) geo.delete();

  geoptic.message("reading mesh ...");
  // give browser time to print the message
  setTimeout(() => {
    geo = Module.readMesh(text, "obj");

    geoptic.message("registering meshes with geoptic ...");
    setTimeout(() => {
      psBaseMesh = geoptic.registerSurfaceMesh(
        "Base Mesh",
        geo.vertexCoordinates(),
        geo.polygons()
      );

      boundaryField = Module.generateSmoothBoundaryField(geo);
      psBaseMesh.addVertexVectorQuantity("Boundary field", boundaryField);

      // update metadata
      geoptic.message("Done");

      // turn off spinner
      document.getElementById("spinner").style.display = "none";
    }, 0);
  }, 0);
}
