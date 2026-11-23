// Run with: deno test --allow-read tests/review_page/core_test.js
import { assertEquals, assert } from "jsr:@std/assert@1";
const src = await Deno.readTextFile(new URL("../../workflow/scripts/review_page/core.js", import.meta.url));
(0, eval)(src);
const Core = globalThis.Core;

Deno.test("canonLabel merges spellings", () => {
  assertEquals(Core.canonLabel("non functional", []), "Non-functional");
  assertEquals(Core.canonLabel(" a/d ", ["A/D", "B"]), "A/D");
  assertEquals(Core.canonLabel("yes", ["Yes"]), "Yes");
  assertEquals(Core.canonLabel("  ", ["A"]), "");
  assertEquals(Core.canonLabel("P1", ["A"]), "P1");
});

Deno.test("tsv round trip, tabs and newlines in text are flattened", () => {
  const rows = [["a", "b"], ["x\ty", "line1\nline2"], [1.5, null]];
  assertEquals(Core.tsvParse(Core.tsvSerialize(rows)), [["a", "b"], ["x y", "line1 line2"], ["1.5", ""]]);
});

Deno.test("crc32 known value", () => {
  assertEquals(Core.crc32(new TextEncoder().encode("123456789")), 0xcbf43926);
});

Deno.test("xlsx written then read, with special characters", async () => {
  const rows = [["id", "Subtype", "n", "Note"], ["S1-1", "A", 12, "a < b & \"c\""], ["S1-2", "", 0.25, ""], ["S1-3", "Non-functional", null, "?"]];
  const bytes = Core.xlsxBuild([{ name: "suggestions", rows, header: true, freeze: { rows: 1, cols: 1 }, fillCols: [1], widths: [12, 12, 8, 20] },
                                { name: "legend", rows: [["column", "meaning"], ["x", "y"]] }]);
  const back = await Core.xlsxRead(bytes);
  assertEquals(back[0], ["id", "Subtype", "n", "Note"]);
  assertEquals(back[1], ["S1-1", "A", "12", "a < b & \"c\""]);
  assertEquals(back[2][1], "");
  assertEquals(back[3][3], "?");
  // R openxlsx (used by Stage 4) cannot parse <t xml:space="preserve">, so the writer must not emit it
  const sheetXml = new TextDecoder().decode((await Core.unzip(bytes))["xl/worksheets/sheet1.xml"]);
  assert(!sheetXml.includes("xml:space"));
  assert(sheetXml.includes("<is><t>S1-1</t></is>"));
  // valid zip: central directory entries and end record are present
  const dv = new DataView(bytes.buffer);
  assertEquals(dv.getUint32(bytes.length - 22, true), 0x06054b50);
});

Deno.test("booleans are written as boolean cells", async () => {
  const bytes = Core.xlsxBuild([{ name: "s", rows: [["a", "b"], [true, false]] }]);
  const xml = new TextDecoder().decode((await Core.unzip(bytes))["xl/worksheets/sheet1.xml"]);
  assert(xml.includes('t="b"><v>1</v>') && xml.includes('t="b"><v>0</v>'));
});

Deno.test("reads an xlsx written by openxlsx (deflate, shared strings)", async () => {
  const p = Deno.env.get("KORV_TEST_XLSX");
  if (!p) return;
  const rows = await Core.xlsxRead(await Deno.readFile(p));
  assert(rows.length > 1);
  assert(rows[0].includes("new_rep_97_ID"));
  assert(rows[0].includes("Subtype"));
});

Deno.test("treeLayout", () => {
  // ((A:1,B:2):1,C:4);  root has 2 children
  const nodes = [
    { parent: -1, bl: 0, tip: null }, { parent: 0, bl: 1, tip: null },
    { parent: 1, bl: 1, tip: "A" }, { parent: 1, bl: 2, tip: "B" }, { parent: 0, bl: 4, tip: "C" },
  ];
  const lc = Core.treeLayout(nodes, "clado");
  assertEquals(lc.order.map((v) => nodes[v].tip), ["A", "B", "C"]);
  assertEquals([lc.x[2], lc.x[3], lc.x[4]], [lc.maxX, lc.maxX, lc.maxX]);   // tips aligned
  assertEquals(lc.x[lc.root], 0);
  assertEquals(lc.y[1], 0.5);
  const ld = Core.treeLayout(nodes, "dist");
  assertEquals([ld.x[2], ld.x[3], ld.x[4]], [2, 3, 4]);
  assertEquals(lc.tipsUnder(1).sort(), [2, 3]);
});

Deno.test("subtypeRank puts intermediates next to their nearer parent, facing the other parent", () => {
  const cfg = { order: ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N"],
                intermediates: { QL: { parents: ["L", "A"], nearer: "L" }, QD: { parents: ["D", "L"], nearer: "D" } } };
  const labels = ["Non-functional", "QL", "L", "possible new subtype (nearest A)", "A/D", "D", "QD", "A", "", "Not_KoRV", "B", "L/D", "D/F"];
  const sorted = labels.slice().sort((a, b) => Core.subtypeCompare(a, b, cfg));
  // A/D after A (nearer = first parent, facing D); QD after D (facing L); QL and L/D before L (facing A or D); ties by label
  assertEquals(sorted, ["A", "A/D", "B", "D", "D/F", "QD", "L/D", "QL", "L", "possible new subtype (nearest A)", "Not_KoRV", "Non-functional", ""]);
});

Deno.test("subtypeRank without a configuration still orders sensibly", () => {
  assertEquals(["Non-functional", "X", "possible new subtype"].sort((a, b) => Core.subtypeCompare(a, b, null)), ["possible new subtype", "X", "Non-functional"]);
  assertEquals(Core.subtypeRank("A/D", { order: ["A", "B"] }), 3);   // D is not in the order: an unknown label, after "possible new subtype"
});
