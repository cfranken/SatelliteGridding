// Exercises the exact read path the browser dashboard uses against the REAL
// store: withConsolidatedMetadata -> FetchStore -> get([ti,null,null]).
// Asserts physically sane values for a 2024 day. Usage:
//   node check_dashboard.mjs http://127.0.0.1:8137/oco2_daily.zarr
import * as zarr from "zarrita";

const url = process.argv[2] || "http://127.0.0.1:8137/oco2_daily.zarr";
const NX = 360, NY = 180;
const ti = 1277;                       // 2024-07-01  (days since 2021-01-01)
const dateStr = new Date(Date.UTC(2021,0,1) + ti*86400000).toISOString().slice(0,10);

const store = await zarr.withConsolidatedMetadata(new zarr.FetchStore(url));

function stats(data){
  let n=0, min=Infinity, max=-Infinity, sum=0;
  for (const v of data){ if (Number.isNaN(v)) continue; n++; sum+=v; if(v<min)min=v; if(v>max)max=v; }
  return { n, min, max, mean: n? sum/n : NaN };
}
// lon/lat -> flat index (row 0 = +89.5N, col 0 = -179.5E)
const idx = (lat,lon) => Math.round(89.5-lat)*NX + Math.round(lon+179.5);

console.log(`store: ${url}`);
console.log(`time index ${ti} = ${dateStr}\n`);

for (const name of ["xco2","sif757","sifdaily","xco2_unc"]) {
  const a = await zarr.open(zarr.root(store).resolve(name), { kind: "array" });
  const { data, shape } = await zarr.get(a, [ti, null, null]);
  if (shape[0] !== NY || shape[1] !== NX) throw new Error(`${name} shape ${shape}`);
  const s = stats(data);
  const pct = (100*s.n/(NX*NY)).toFixed(1);
  console.log(`${name.padEnd(9)} populated ${String(s.n).padStart(6)}/${NX*NY} (${pct}%)  `
    + `min ${s.min.toFixed(3)}  mean ${s.mean.toFixed(3)}  max ${s.max.toFixed(3)}`);
}

// spot-check a couple of land cells in the consolidated array shape
const xco2 = await zarr.open(zarr.root(store).resolve("xco2"), { kind: "array" });
const { data: xd } = await zarr.get(xco2, [ti, null, null]);
const samples = [["Amazon",-3,-60],["Sahara",23,10],["Siberia",62,100],["US-Midwest",41,-93]];
console.log("\nxco2 spot cells (ppm):");
for (const [nm,la,lo] of samples){ const v = xd[idx(la,lo)]; console.log(`  ${nm.padEnd(11)} (${la},${lo})  ${Number.isNaN(v)?"— (no data)":v.toFixed(2)}`); }

// sanity gates
const sx = stats(xd);
if (!(sx.n > 1000)) throw new Error("FAIL: xco2 has too few populated cells");
if (!(sx.mean > 400 && sx.mean < 440)) throw new Error(`FAIL: xco2 mean ${sx.mean} out of [400,440] ppm`);
console.log("\nDASHBOARD DATA PATH OK");
