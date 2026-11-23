// Pure logic for the review page: no DOM access, so it can be tested outside a browser.
// Exposes a single global `Core`.
const Core = (function () {
  "use strict";
  const enc = new TextEncoder();
  const dec = new TextDecoder();

  // ---- labels -------------------------------------------------------------------------------
  // Key used to decide that two spellings are the same label.
  const labelKey = (s) => String(s).toLowerCase().replace(/[\s_-]+/g, "");
  // Return the spelling to use for a typed label: an existing spelling that differs only in case,
  // spaces, '-' or '_', "Non-functional" for any spelling of that, otherwise the trimmed text.
  function canonLabel(text, existing) {
    const t = String(text).trim();
    if (!t) return "";
    const k = labelKey(t);
    if (k === "nonfunctional") return "Non-functional";
    for (const e of existing) if (labelKey(e) === k) return e;
    return t;
  }

  // ---- TSV -----------------------------------------------------------------------------------
  const tsvCell = (v) => (v === null || v === undefined ? "" : String(v)).replace(/[\t\r\n]+/g, " ");
  const tsvSerialize = (rows) => rows.map((r) => r.map(tsvCell).join("\t")).join("\n") + "\n";
  function tsvParse(text) {
    const lines = text.replace(/^﻿/, "").split(/\r?\n/).filter((l) => l.length);
    return lines.map((l) => l.split("\t"));
  }

  // ---- CRC32 and a store-only ZIP writer -----------------------------------------------------
  const crcTable = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
      t[n] = c >>> 0;
    }
    return t;
  })();
  function crc32(bytes) {
    let c = 0xffffffff;
    for (let i = 0; i < bytes.length; i++) c = crcTable[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
    return (c ^ 0xffffffff) >>> 0;
  }
  // files: [{name, data: Uint8Array}] -> Uint8Array of a valid ZIP archive (no compression)
  function zipStore(files) {
    const parts = [], central = [];
    let offset = 0;
    const u16 = (v) => [v & 255, (v >>> 8) & 255];
    const u32 = (v) => [v & 255, (v >>> 8) & 255, (v >>> 16) & 255, (v >>> 24) & 255];
    for (const f of files) {
      const name = enc.encode(f.name), crc = crc32(f.data), size = f.data.length;
      const local = new Uint8Array([0x50, 0x4b, 3, 4, ...u16(20), ...u16(0x0800), ...u16(0), ...u16(0), ...u16(0x21),
        ...u32(crc), ...u32(size), ...u32(size), ...u16(name.length), ...u16(0)]);
      parts.push(local, name, f.data);
      central.push(new Uint8Array([0x50, 0x4b, 1, 2, ...u16(20), ...u16(20), ...u16(0x0800), ...u16(0), ...u16(0), ...u16(0x21),
        ...u32(crc), ...u32(size), ...u32(size), ...u16(name.length), ...u16(0), ...u16(0), ...u16(0), ...u16(0), ...u32(0), ...u32(offset)]), name);
      offset += local.length + name.length + size;
    }
    const cdSize = central.reduce((s, p) => s + p.length, 0);
    const end = new Uint8Array([0x50, 0x4b, 5, 6, ...u16(0), ...u16(0), ...u16(files.length), ...u16(files.length), ...u32(cdSize), ...u32(offset), ...u16(0)]);
    const all = [...parts, ...central, end];
    const out = new Uint8Array(all.reduce((s, p) => s + p.length, 0));
    let o = 0;
    for (const p of all) { out.set(p, o); o += p.length; }
    return out;
  }

  // ---- XLSX writer ---------------------------------------------------------------------------
  const xmlEsc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]))
    .replace(/[\u0000-\u0008\u000b\u000c\u000e-\u001f]/g, "");
  const colLetter = (n) => { let s = ""; n++; while (n > 0) { const m = (n - 1) % 26; s = String.fromCharCode(65 + m) + s; n = Math.floor((n - 1) / 26); } return s; };
  // sheets: [{name, rows: [[cell...]], header: bool, freeze: {rows, cols}, fillCols: [indices], widths: [chars]}]
  // cell: number -> numeric cell, boolean -> boolean cell, null/undefined/"" -> empty, other -> inline string
  function xlsxBuild(sheets) {
    const files = [];
    const add = (name, text) => files.push({ name, data: enc.encode(text) });
    add("[Content_Types].xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>' +
      sheets.map((s, i) => `<Override PartName="/xl/worksheets/sheet${i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>`).join("") + "</Types>");
    add("_rels/.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>');
    add("xl/workbook.xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' +
      sheets.map((s, i) => `<sheet name="${xmlEsc(s.name)}" sheetId="${i + 1}" r:id="rId${i + 1}"/>`).join("") + "</sheets></workbook>");
    add("xl/_rels/workbook.xml.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' +
      sheets.map((s, i) => `<Relationship Id="rId${i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet${i + 1}.xml"/>`).join("") +
      `<Relationship Id="rId${sheets.length + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>`);
    // styles: 0 plain, 1 bold header with bottom border, 2 yellow fill (columns filled in by the expert)
    add("xl/styles.xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FFFFF2CC"/><bgColor indexed="64"/></patternFill></fill></fills><borders count="2"><border><left/><right/><top/><bottom/><diagonal/></border><border><left/><right/><top/><bottom style="thin"><color auto="1"/></bottom><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="1" xfId="0" applyFont="1" applyBorder="1"/><xf numFmtId="0" fontId="0" fillId="2" borderId="0" xfId="0" applyFill="1"/></cellXfs></styleSheet>');
    sheets.forEach((s, i) => {
      const ncol = Math.max(...s.rows.map((r) => r.length), 1);
      let x = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">';
      if (s.freeze) {
        const fr = s.freeze.rows || 0, fc = s.freeze.cols || 0;
        x += `<sheetViews><sheetView workbookViewId="0"><pane ${fc ? `xSplit="${fc}" ` : ""}${fr ? `ySplit="${fr}" ` : ""}topLeftCell="${colLetter(fc)}${fr + 1}" activePane="${fr && fc ? "bottomRight" : fr ? "bottomLeft" : "topRight"}" state="frozen"/></sheetView></sheetViews>`;
      }
      if (s.widths) x += "<cols>" + s.widths.map((w, j) => `<col min="${j + 1}" max="${j + 1}" width="${w}" customWidth="1"/>`).join("") + "</cols>";
      x += "<sheetData>";
      s.rows.forEach((row, ri) => {
        x += `<row r="${ri + 1}">`;
        row.forEach((v, ci) => {
          const ref = colLetter(ci) + (ri + 1);
          const style = s.header && ri === 0 ? 1 : (s.fillCols && s.fillCols.includes(ci) && ri > 0 ? 2 : 0);
          const st = style ? ` s="${style}"` : "";
          if (v === null || v === undefined || v === "") { if (style) x += `<c r="${ref}"${st}/>`; }
          else if (typeof v === "number" && isFinite(v)) x += `<c r="${ref}"${st}><v>${v}</v></c>`;
          else if (typeof v === "boolean") x += `<c r="${ref}"${st} t="b"><v>${v ? 1 : 0}</v></c>`;
          else x += `<c r="${ref}"${st} t="inlineStr"><is><t>${xmlEsc(String(v).trim())}</t></is></c>`;
        });
        x += "</row>";
      });
      x += "</sheetData></worksheet>";
      add(`xl/worksheets/sheet${i + 1}.xml`, x);
    });
    return zipStore(files);
  }

  // ---- XLSX reader (first sheet, values as strings) --------------------------------------------
  const u16le = (b, o) => b[o] | (b[o + 1] << 8);
  const u32le = (b, o) => (b[o] | (b[o + 1] << 8) | (b[o + 2] << 16) | (b[o + 3] << 24)) >>> 0;
  async function inflateRaw(bytes) {
    const ds = new DecompressionStream("deflate-raw");
    const stream = new Blob([bytes]).stream().pipeThrough(ds);
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }
  // returns {name: Uint8Array}
  async function unzip(bytes) {
    let e = bytes.length - 22;
    while (e >= 0 && u32le(bytes, e) !== 0x06054b50) e--;
    if (e < 0) throw new Error("not a zip file");
    const n = u16le(bytes, e + 10); let p = u32le(bytes, e + 16);
    const out = {};
    for (let i = 0; i < n; i++) {
      if (u32le(bytes, p) !== 0x02014b50) throw new Error("bad zip directory");
      const method = u16le(bytes, p + 10), csize = u32le(bytes, p + 20), nlen = u16le(bytes, p + 28), xlen = u16le(bytes, p + 30), clen = u16le(bytes, p + 32), loff = u32le(bytes, p + 42);
      const name = dec.decode(bytes.subarray(p + 46, p + 46 + nlen));
      const lnlen = u16le(bytes, loff + 26), lxlen = u16le(bytes, loff + 28), start = loff + 30 + lnlen + lxlen;
      const raw = bytes.subarray(start, start + csize);
      out[name] = method === 0 ? raw : method === 8 ? await inflateRaw(raw) : null;
      p += 46 + nlen + xlen + clen;
    }
    return out;
  }
  const xmlUnesc = (s) => s.replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&quot;/g, '"').replace(/&apos;/g, "'")
    .replace(/&#x([0-9a-fA-F]+);/g, (_, h) => String.fromCodePoint(parseInt(h, 16))).replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(+d)).replace(/&amp;/g, "&");
  const colIndex = (letters) => { let n = 0; for (const ch of letters) n = n * 26 + ch.charCodeAt(0) - 64; return n - 1; };
  async function xlsxRead(bytes) {
    const z = await unzip(bytes);
    const text = (k) => (z[k] ? dec.decode(z[k]) : null);
    let sheetPath = "xl/worksheets/sheet1.xml";
    const wb = text("xl/workbook.xml"), rels = text("xl/_rels/workbook.xml.rels");
    if (wb && rels) {
      const rid = (wb.match(/<sheet\b[^>]*r:id="([^"]+)"/) || [])[1];
      const tgt = rid && (rels.match(new RegExp('<Relationship\\b[^>]*Id="' + rid + '"[^>]*>')) || [""])[0].match(/Target="([^"]+)"/);
      if (tgt) sheetPath = "xl/" + tgt[1].replace(/^\/?(xl\/)?/, "");
    }
    const sheet = text(sheetPath);
    if (!sheet) throw new Error("no worksheet found");
    const shared = [];
    const sst = text("xl/sharedStrings.xml");
    if (sst) for (const m of sst.matchAll(/<si\b[^>]*>([\s\S]*?)<\/si>/g)) {
      let s = ""; for (const t of m[1].matchAll(/<t\b[^>]*>([\s\S]*?)<\/t>/g)) s += t[1];
      shared.push(xmlUnesc(s));
    }
    const rows = [];
    for (const rm of sheet.matchAll(/<row\b[^>]*?(\/>|>([\s\S]*?)<\/row>)/g)) {
      if (rm[1] === "/>" || !rm[2]) continue;
      const row = [];
      for (const cm of rm[2].matchAll(/<c\b([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g)) {
        const attrs = cm[1], body = cm[2] || "";
        const ref = (attrs.match(/\br="([A-Z]+)\d+"/) || [])[1], type = (attrs.match(/\bt="([^"]+)"/) || [])[1];
        let v = "";
        if (type === "inlineStr") { for (const t of body.matchAll(/<t\b[^>]*>([\s\S]*?)<\/t>/g)) v += t[1]; v = xmlUnesc(v); }
        else { const vm = body.match(/<v>([\s\S]*?)<\/v>/); if (vm) v = type === "s" ? shared[+vm[1]] : xmlUnesc(vm[1]); }
        row[ref ? colIndex(ref) : row.length] = v;
      }
      for (let i = 0; i < row.length; i++) if (row[i] === undefined) row[i] = "";
      rows.push(row);
    }
    return rows;
  }

  // ---- tree layout ---------------------------------------------------------------------------
  // nodes: [{parent: -1 for the root, bl: branch length, tip: label or null}], index = node id.
  // mode "clado": equal-length branches with the tips aligned; "dist": cumulative branch length.
  // Returns {children, order (tip node ids top to bottom), x, y (tip index units), maxX, root, tipsUnder(node)}
  function treeLayout(nodes, mode) {
    const n = nodes.length, children = Array.from({ length: n }, () => []);
    let root = -1;
    nodes.forEach((nd, i) => { if (nd.parent < 0) root = i; else children[nd.parent].push(i); });
    const order = [], post = [], stack = [[root, 0]];
    while (stack.length) {            // iterative depth-first walk
      const top = stack[stack.length - 1], [v, k] = top;
      if (k === 0 && children[v].length === 0) order.push(v);
      if (k < children[v].length) { top[1]++; stack.push([children[v][k], 0]); }
      else { post.push(v); stack.pop(); }
    }
    const y = new Float64Array(n), h = new Float64Array(n), x = new Float64Array(n);
    order.forEach((v, i) => { y[v] = i; });
    for (const v of post) {
      if (children[v].length) {
        y[v] = (y[children[v][0]] + y[children[v][children[v].length - 1]]) / 2;
        h[v] = Math.max(...children[v].map((c) => h[c])) + 1;
      }
    }
    if (mode === "dist") {
      const pre = [root];
      x[root] = 0;
      while (pre.length) { const v = pre.pop(); for (const c of children[v]) { x[c] = x[v] + (nodes[c].bl || 0); pre.push(c); } }
    } else {
      const H = h[root];
      for (let v = 0; v < n; v++) x[v] = H - h[v];
    }
    let maxX = 0; for (let v = 0; v < n; v++) maxX = Math.max(maxX, x[v]);
    function tipsUnder(v) {
      const out = [], st = [v];
      while (st.length) { const u = st.pop(); if (children[u].length) st.push(...children[u]); else out.push(u); }
      return out;
    }
    return { children, order, x, y, maxX, root, tipsUnder };
  }

  // Position of a subtype label in the row order used when the page is sorted by subtype.
  // cfg = {order: [subtypes], intermediates: {QL: {parents: [L, A], nearer: "L"}}}. An intermediate (or a label
  // written "X/Y") goes next to its nearer parent, on the side that faces the other parent, so it sits between the two
  // as closely as a single row order allows. Then come "possible new subtype", other labels, Non-functional, no label.
  function subtypeRank(label, cfg) {
    const order = (cfg && cfg.order) || [], inter = (cfg && cfg.intermediates) || {}, n = order.length;
    if (!label) return n + 3;
    if (label === "Non-functional") return n + 2;
    if (/^possible new/.test(label)) return n + 0.5;
    const at = order.indexOf(label);
    if (at >= 0) return at;
    let spec = inter[label];
    if (!spec && label.indexOf("/") > 0) { const p = label.split("/").map((t) => t.trim()); if (p.length === 2) spec = { parents: p, nearer: p[0] }; }
    if (spec && spec.parents && spec.parents.length === 2) {
      const near = spec.nearer || spec.parents[0], far = spec.parents[0] === near ? spec.parents[1] : spec.parents[0];
      const i = order.indexOf(near), j = order.indexOf(far);
      if (i >= 0 && j >= 0 && i !== j) return i + (j > i ? 0.5 : -0.5);
    }
    return n + 1;
  }
  const subtypeCompare = (a, b, cfg) => subtypeRank(a, cfg) - subtypeRank(b, cfg) || String(a).localeCompare(String(b));

  return { labelKey, canonLabel, subtypeRank, subtypeCompare, tsvSerialize, tsvParse, crc32, zipStore, xlsxBuild, xlsxRead, unzip, colLetter, treeLayout };
})();
if (typeof globalThis !== "undefined") globalThis.Core = Core;
