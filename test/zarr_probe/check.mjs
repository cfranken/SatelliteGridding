// Ordering check via zarrita (spec §4.3) — authoritative for the browser dashboard.
// Reads the store over HTTP through the consolidated metadata (.zmetadata), exactly
// as the front end does. zarrita 0.7 ships a default Blosc codec, so no manual
// codec registration is needed (numcodecs is pulled in transitively).
import * as zarr from "zarrita";

const url = process.argv[2] || "http://localhost:8123/probe.zarr";
const store = await zarr.withConsolidatedMetadata(new zarr.FetchStore(url));
const arr = await zarr.open(zarr.root(store).resolve("probe"), { kind: "array" });

const { data, shape } = await zarr.get(arr, [0, null, null]);   // Float32Array, shape [180, 360]
if (shape[0] !== 180 || shape[1] !== 360) {
  throw new Error(`JS shape FAIL: [${shape}] != [180,360]`);
}
const NX = 360;
for (const [iy, ix] of [[0, 0], [44, 59], [89, 180], [179, 359]]) {
  const got = data[iy * NX + ix];
  const exp = iy * 1000 + ix;
  if (got !== exp) throw new Error(`JS order FAIL @(${iy},${ix}): ${got} != ${exp}`);
}
console.log("JS ORDERING OK");
