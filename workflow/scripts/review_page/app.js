// Interactive review page: tree and alignment side by side with linked selection and subtype assignment.
(function () {
  "use strict";
  const D = JSON.parse(document.getElementById("korv-data").textContent);
  const $ = (s, r) => (r || document).querySelector(s);
  const esc = (s) => String(s === null || s === undefined ? "" : s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  // ---- data ---------------------------------------------------------------------------------------
  const HEAD = 52, LABW = 364;
  const ids = D.alignment.order, ncol = D.alignment.ncol;
  const refInfo = {}; D.references.forEach((r) => { refInfo[r.id] = r; });
  const isRef = (id) => !!refInfo[id];
  const seqById = {}; D.sequences.forEach((s) => { seqById[s.id] = s; });
  const reps = D.sequences.map((s) => s.id);
  const nodes = D.tree.nodes;
  const colIdx = {}; D.export_columns.forEach((c, i) => { colIdx[c] = i; });
  const FILL_COLS = ["Subtype", "Endogenous", "Unsure", "Note"];
  const CX = { b: "#80a0f0", r: "#f01505", m: "#c048c0", g: "#15c015", o: "#f09048", c: "#15a4a4", p: "#f08080", y: "#c0c000", ".": "#ffffff" };
  const vrSet = new Set(D.rules.variable_region_cols), nqSet = new Set(D.rules.nqqhnn_cols);
  const x = (id, col) => { const s = seqById[id]; return s ? s.export[col] : null; };
  const reads = (id) => (x(id, "Total_number_of_reads") || 0);
  const suggested = (id) => x(id, "suggested_subtype") || "";
  const confidence = (id) => x(id, "confidence") || "";

  // ---- state --------------------------------------------------------------------------------------
  const S = { order: "suggested", treeMode: "clado", colourBy: "suggested", alnColour: "clustalx", rowH: 14, cellW: 11,
    sel: new Set(), anchor: null, focus: null, assign: {}, undo: [], redo: [],
    hide: true, f: { conf: "", flagged: false, unassigned: false, differs: false, text: "", minReads: 0 } };
  const blank = () => ({ Subtype: "", Endogenous: "", Unsure: "", Note: "" });
  const asg = (id) => S.assign[id] || null;
  const assigned = (id) => (S.assign[id] && S.assign[id].Subtype) || "";

  // ---- tree layout and row orders ----------------------------------------------------------------------
  const layouts = {};
  const layout = (m) => layouts[m || S.treeMode] || (layouts[m || S.treeMode] = Core.treeLayout(nodes, m || S.treeMode));
  const treeIds = layout("clado").order.map((v) => nodes[v].tip);
  const treeIndex = {}; treeIds.forEach((id, i) => { treeIndex[id] = i; });
  const tipNode = {}; nodes.forEach((n, i) => { if (n.tip) tipNode[n.tip] = i; });
  const lo = new Int32Array(nodes.length), hi = new Int32Array(nodes.length);
  (function ranges() {
    const ch = layout("clado").children;
    const go = (v) => {
      if (!ch[v].length) { lo[v] = hi[v] = treeIndex[nodes[v].tip]; return; }
      lo[v] = 1e9; hi[v] = -1;
      for (const c of ch[v]) { go(c); lo[v] = Math.min(lo[v], lo[c]); hi[v] = Math.max(hi[v], hi[c]); }
    };
    go(layout("clado").root);
  })();
  let rowOrder = ids.slice(), rowIndex = {}, rowGroup = null;
  // Subtype order: the configured order, or the alphabetical order of the reference subtypes
  const SUB_CFG = D.subtype_order || { order: [...new Set(D.references.map((r) => r.subtype))].sort(), intermediates: {} };
  const subtypeGroup = (id, mode) => { const l = isRef(id) ? refInfo[id].subtype : mode === "assigned" ? (assigned(id) || suggested(id)) : suggested(id); return /^possible new/.test(l) ? "possible new subtype" : l; };
  function setOrder(mode) {
    S.order = mode; rowGroup = null;
    const repsOnly = reps.slice(), refs = ids.filter(isRef);
    if (mode === "suggested" || mode === "assigned") {
      // rows grouped by subtype: references first, then the representatives in tree order
      rowGroup = {}; ids.forEach((id) => { rowGroup[id] = subtypeGroup(id, mode); });
      rowOrder = ids.slice().sort((a, b) => Core.subtypeCompare(rowGroup[a], rowGroup[b], SUB_CFG) || (isRef(b) - isRef(a)) || treeIndex[a] - treeIndex[b]);
    }
    else if (mode === "tree") rowOrder = treeIds.slice();
    else if (mode === "file") rowOrder = ids.slice();
    else if (mode === "reads") rowOrder = repsOnly.sort((a, b) => reads(b) - reads(a)).concat(refs);
    else if (mode === "workbook") rowOrder = repsOnly.concat(refs);
    else if (mode === "status") rowOrder = repsOnly.sort((a, b) => (assigned(a) ? 1 : 0) - (assigned(b) ? 1 : 0)).concat(refs);
    if (S.hide && filterOn()) rowOrder = rowOrder.filter(matches);   // filtered rows are hidden, references stay unless searched for
    rowIndex = {}; rowOrder.forEach((id, i) => { rowIndex[id] = i; });
  }
  // the tree and the alignment share rows (and one scroll position) only when nothing is hidden
  const treeAligned = () => S.order === "tree" && !(S.hide && filterOn());

  // ---- colours ------------------------------------------------------------------------------------------
  const FALLBACK = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d", "#1f78b4", "#00bcd4", "#33a02c", "#6a3d9a", "#e6ab02", "#8c564b"];
  function colourOf(label) {
    if (!label) return "#bbbbbb";
    const g = /^possible new/.test(label) ? "possible new subtype" : label;
    if (D.palette[g]) return D.palette[g];
    let h = 0; for (const ch of g) h = (h * 31 + ch.charCodeAt(0)) >>> 0;
    return FALLBACK[h % FALLBACK.length];
  }
  const tipColour = (id) => isRef(id) ? colourOf(refInfo[id].subtype) : S.colourBy === "assigned" ? (assigned(id) ? colourOf(assigned(id)) : "#ffffff") : colourOf(suggested(id));

  // ---- filters ------------------------------------------------------------------------------------------
  function matches(id) {
    const f = S.f;
    if (isRef(id)) return !f.text || (id + refInfo[id].name).toLowerCase().includes(f.text);
    if (f.text && !id.toLowerCase().includes(f.text) && !suggested(id).toLowerCase().includes(f.text) && !assigned(id).toLowerCase().includes(f.text)) return false;
    if (f.conf && confidence(id) !== f.conf) return false;
    if (f.flagged && !seqById[id].flags.length) return false;
    if (f.unassigned && assigned(id)) return false;
    if (f.differs && (!assigned(id) || assigned(id) === suggested(id))) return false;
    if (reads(id) < f.minReads) return false;
    return true;
  }
  const filterOn = () => { const f = S.f; return !!(f.conf || f.flagged || f.unassigned || f.differs || f.text || f.minReads); };

  // ---- assignments, undo, persistence --------------------------------------------------------------------
  const vocabulary = () => {
    const seen = new Map();
    const add = (l) => { if (l && !seen.has(l)) seen.set(l, 0); };
    D.vocabulary_seed.forEach(add); reps.forEach((id) => add(suggested(id)));
    reps.forEach((id) => { const a = assigned(id); if (a) { add(a); seen.set(a, (seen.get(a) || 0) + 1); } });
    return [...seen.entries()].filter(([l]) => !/^possible new/.test(l)).sort((a, b) => b[1] - a[1]).map(([l, n]) => ({ label: l, n }));
  };
  // Quick labels have a fixed order, so a key always means the same label: 1-9 are the first nine subtypes of the
  // references (alphabetical), 0 is Non-functional. Further labels (the rest of the references and any label you
  // have used) get unnumbered buttons in alphabetical order.
  const quickLabels = () => {
    const seed = D.vocabulary_seed.filter((l) => l !== "Non-functional"), numbered = seed.slice(0, 9);
    const extra = [...new Set(seed.slice(9).concat(reps.map(assigned).filter(Boolean)))].filter((l) => l !== "Non-functional" && !numbered.includes(l)).sort();
    return { numbered, extra };
  };
  const existingLabels = () => vocabulary().map((v) => v.label);
  function commit(idList, patchFn, note) {
    const before = idList.map((id) => JSON.stringify(S.assign[id] || null)), after = [];
    idList.forEach((id) => {
      const cur = Object.assign(blank(), S.assign[id] || {}); const nxt = patchFn(cur, id);
      const empty = !nxt.Subtype && !nxt.Endogenous && !nxt.Unsure && !nxt.Note;
      if (empty) delete S.assign[id]; else S.assign[id] = nxt;
      after.push(JSON.stringify(S.assign[id] || null));
    });
    if (before.some((b, i) => b !== after[i])) { S.undo.push({ ids: idList.slice(), before, after, note }); if (S.undo.length > 200) S.undo.shift(); S.redo = []; }
    afterChange();
  }
  function applyState(entry, which) {
    entry.ids.forEach((id, i) => { const v = JSON.parse(entry[which][i]); if (v) S.assign[id] = v; else delete S.assign[id]; });
    afterChange();
  }
  const undo = () => { const e = S.undo.pop(); if (e) { S.redo.push(e); applyState(e, "before"); toast("Undid: " + (e.note || "change")); } };
  const redo = () => { const e = S.redo.pop(); if (e) { S.undo.push(e); applyState(e, "after"); toast("Redid: " + (e.note || "change")); } };
  const storeKey = "korv-review:" + D.analysis_name + ":" + D.fingerprint;
  let saveTimer = null;
  function scheduleSave() {
    clearTimeout(saveTimer);
    saveTimer = setTimeout(() => { try { localStorage.setItem(storeKey, JSON.stringify({ assign: S.assign, t: Date.now() })); } catch (e) { /* storage may be blocked */ } }, 300);
  }
  function restore() {
    try { const v = JSON.parse(localStorage.getItem(storeKey) || "null"); if (v && v.assign && Object.keys(v.assign).length) { S.assign = v.assign; return Object.keys(v.assign).length; } } catch (e) { /* ignore */ }
    return 0;
  }
  const nAssigned = () => reps.filter((id) => assigned(id)).length;
  function afterChange() { scheduleSave(); renderTop(); refreshRows(); renderSide(); }
  function refreshRows() { setOrder(S.order); sizeCanvas(); renderTree(); requestDraw(); }

  // ---- selection -------------------------------------------------------------------------------------------
  function select(idList, mode, focusId) {
    if (mode === "set") S.sel = new Set(idList);
    else if (mode === "add") idList.forEach((i) => S.sel.add(i));
    else if (mode === "toggle") idList.forEach((i) => (S.sel.has(i) ? S.sel.delete(i) : S.sel.add(i)));
    S.focus = focusId || idList[idList.length - 1] || (S.sel.size ? [...S.sel][0] : null);
    updateTreeHighlight(); renderSide(); requestDraw();
  }
  function clickRow(id, ev, source) {
    if (ev.shiftKey && S.anchor && rowIndex[S.anchor] !== undefined) {
      const a = rowIndex[S.anchor], b = rowIndex[id], [s, e] = a < b ? [a, b] : [b, a];
      select(rowOrder.slice(s, e + 1), ev.ctrlKey || ev.metaKey ? "add" : "set", id);
    } else if (ev.ctrlKey || ev.metaKey) { select([id], "toggle", id); S.anchor = id; }
    else { select([id], "set", id); S.anchor = id; }
    reveal(source);
  }
  function scrollToRow(el, idx) {
    const top = HEAD + idx * S.rowH;
    if (top < el.scrollTop + HEAD || top + S.rowH > el.scrollTop + el.clientHeight) el.scrollTop = Math.max(0, top - el.clientHeight / 2);
  }
  function reveal(source) {
    if (!S.focus) return;
    if (source !== "aln" && rowIndex[S.focus] !== undefined) scrollToRow($("#alnScroll"), rowIndex[S.focus]);
    if (source !== "tree" && !treeAligned()) scrollToRow($("#treeScroll"), treeIndex[S.focus]);
  }

  // ---- tree view ---------------------------------------------------------------------------------------------
  let nodeGroups = [], bandsG = null;
  function renderTree() {
    const svg = $("#treeSvg"), L = layout(), W = Math.max(160, $("#treeScroll").clientWidth - 2), n = treeIds.length;
    const showText = !treeAligned(), M = { l: 8, r: showText ? 130 : 12 };
    const sx = (W - M.l - M.r) / (L.maxX || 1), H = HEAD + n * S.rowH + 24;
    const X = (v) => M.l + L.x[v] * sx, Y = (v) => HEAD + (L.y[v] + 0.5) * S.rowH;
    svg.setAttribute("width", W); svg.setAttribute("height", H);
    let out = '<g id="bands"></g>';
    nodes.forEach((nd, v) => {
      const kids = L.children[v]; let d = "";
      if (nd.parent >= 0) d += `M${X(nd.parent).toFixed(1)},${Y(v).toFixed(1)}L${X(v).toFixed(1)},${Y(v).toFixed(1)}`;
      if (kids.length) d += `M${X(v).toFixed(1)},${Y(kids[0]).toFixed(1)}L${X(v).toFixed(1)},${Y(kids[kids.length - 1]).toFixed(1)}`;
      if (d) out += `<g class="nd" data-n="${v}"><path class="br" d="${d}"/><path class="hit" d="${d}"/></g>`;
    });
    L.order.forEach((v) => {
      const id = nodes[v].tip, cx = X(v), cy = Y(v), dim = matches(id) ? 1 : 0.22, col = tipColour(id);
      const r = Math.max(2.5, Math.min(4.5, S.rowH / 3));
      const flagged = !isRef(id) && seqById[id].flags.length, novel = !isRef(id) && confidence(id) === "novel";
      const stroke = novel ? "#e41a1c" : isRef(id) ? "#222" : "#666", sw = novel ? 2 : isRef(id) ? 1.5 : 1;
      out += flagged
        ? `<polygon class="tipdot" data-id="${esc(id)}" opacity="${dim}" points="${cx},${cy - r - 1} ${cx + r + 1},${cy + r} ${cx - r - 1},${cy + r}" fill="${col}" stroke="${stroke}" stroke-width="${sw}"/>`
        : `<circle class="tipdot" data-id="${esc(id)}" opacity="${dim}" cx="${cx}" cy="${cy}" r="${r}" fill="${col}" stroke="${stroke}" stroke-width="${sw}"/>`;
      if (showText) out += `<text x="${cx + 8}" y="${cy + 3}" opacity="${dim}" ${isRef(id) ? 'font-weight="bold"' : ""}>${esc(isRef(id) ? "REF " + refInfo[id].name : id)}</text>`;
    });
    if (S.treeMode === "dist") {
      const raw = (W - M.l - M.r) * 0.2 / sx, p = Math.pow(10, Math.floor(Math.log10(raw))), len = [1, 2, 5, 10].map((m) => m * p).filter((v) => v <= raw).pop() || p;
      out += `<g><path d="M${M.l},${H - 14}L${M.l + len * sx},${H - 14}" stroke="#333"/><text x="${M.l}" y="${H - 18}">${len} substitutions per site</text></g>`;
    }
    svg.innerHTML = out; bandsG = $("#bands", svg);
    nodeGroups = [...svg.querySelectorAll("g.nd")];
    updateTreeHighlight();
  }
  function updateTreeHighlight() {
    if (!bandsG) return;
    const pre = new Int32Array(treeIds.length + 1);
    treeIds.forEach((id, i) => { pre[i + 1] = pre[i] + (S.sel.has(id) ? 1 : 0); });
    nodeGroups.forEach((g) => { const v = +g.dataset.n, cnt = pre[hi[v] + 1] - pre[lo[v]]; g.classList.toggle("sel", cnt > 0 && cnt === hi[v] - lo[v] + 1); });
    const W = +$("#treeSvg").getAttribute("width"); let bands = "";
    S.sel.forEach((id) => { if (treeIndex[id] !== undefined) bands += `<rect class="band" x="0" y="${HEAD + treeIndex[id] * S.rowH}" width="${W}" height="${S.rowH}"/>`; });
    bandsG.innerHTML = bands;
    document.querySelectorAll("#treeSvg .tipdot").forEach((d) => d.setAttribute("stroke-width", S.sel.has(d.dataset.id) ? 3 : (d.getAttribute("stroke") === "#e41a1c" ? 2 : 1)));
  }
  $("#treeSvg").addEventListener("click", (ev) => {
    const dot = ev.target.closest(".tipdot"), g = ev.target.closest("g.nd");
    if (dot) return clickRow(dot.dataset.id, ev, "tree");
    if (g) {
      const v = +g.dataset.n, tips = layout("clado").tipsUnder(v).map((t) => nodes[t].tip).sort((a, b) => treeIndex[a] - treeIndex[b]);
      const m = ev.ctrlKey || ev.metaKey || ev.shiftKey ? "add" : "set";
      select(tips, m, tips[0]); S.anchor = tips[0]; reveal("tree");
      toast(tips.length + " sequences in clade");
    }
  });

  // ---- alignment canvas --------------------------------------------------------------------------------------
  const cv = $("#alnCanvas"), ctx = cv.getContext("2d"); let raf = 0;
  const requestDraw = () => {   // animation frames are paused in background tabs, so a timer backs them up
    if (raf) return; raf = 1;
    const run = () => { if (!raf) return; raf = 0; draw(); };
    requestAnimationFrame(run); setTimeout(run, 50);
  };
  function sizeCanvas() {
    const sc = $("#alnScroll"), dpr = window.devicePixelRatio || 1, w = Math.max(100, sc.clientWidth), h = Math.max(100, sc.clientHeight);
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); cv.style.width = w + "px"; cv.style.height = h + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const inner = $("#alnInner");
    inner.style.width = Math.max(w, LABW + ncol * S.cellW) + "px"; inner.style.height = Math.max(h, HEAD + rowOrder.length * S.rowH) + "px";
  }
  function cellColour(id, c, row) {
    if (S.alnColour === "none") return "#ffffff";
    const ch = row[c];
    if (S.alnColour === "differs") {
      if (ch === "-") return "#ffffff";
      if (isRef(id)) return "#cfe3f5";
      const ref = D.alignment.seq[seqById[id].nearest_ref_row];
      return ref && ref[c] === ch ? "#f2f2f2" : "#f6b26b";
    }
    return CX[D.alignment.fill[id][c]] || "#ffffff";
  }
  function draw() {
    const w = cv.clientWidth, h = cv.clientHeight, sc = $("#alnScroll"), st = sc.scrollTop, sl = sc.scrollLeft, rh = S.rowH, cw = S.cellW;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, h);
    const r0 = Math.max(0, Math.floor((st) / rh)), r1 = Math.min(rowOrder.length - 1, Math.ceil((st + h - HEAD) / rh));
    const c0 = Math.max(0, Math.floor(sl / cw)), c1 = Math.min(ncol - 1, Math.ceil((sl + w - LABW) / cw));
    // alignment body (clipped to the area right of the labels)
    ctx.save(); ctx.beginPath(); ctx.rect(LABW, HEAD, w - LABW, h - HEAD); ctx.clip();
    const fontPx = Math.max(6, Math.min(rh - 3, cw - 1)); ctx.font = fontPx + "px ui-monospace, Menlo, monospace"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
    for (let r = r0; r <= r1; r++) {
      const id = rowOrder[r], row = D.alignment.seq[id], y = HEAD + r * rh - st, dim = matches(id) ? 1 : 0.28;
      ctx.globalAlpha = dim;
      for (let c = c0; c <= c1; c++) {
        const px = LABW + c * cw - sl, col = cellColour(id, c, row);
        if (col !== "#ffffff") { ctx.fillStyle = col; ctx.fillRect(px, y, cw, rh); }
        if (cw >= 7 && row[c] !== "-") { ctx.fillStyle = "#111"; ctx.fillText(row[c], px + cw / 2, y + rh / 2 + 0.5); }
      }
      ctx.globalAlpha = 1;
      if (S.sel.has(id)) { ctx.fillStyle = "rgba(255,200,0,0.32)"; ctx.fillRect(LABW, y, w - LABW, rh); }
      if (rowGroup && r > 0 && rowGroup[id] !== rowGroup[rowOrder[r - 1]]) { ctx.fillStyle = "#6b7684"; ctx.fillRect(LABW, y - 1, w - LABW, 2); }   // subtype group boundary
    }
    // variable-region box and NQQHNN shading
    const vc = D.rules.variable_region_cols;
    if (vc.length) {
      const bx = LABW + (vc[0] - 1) * cw - sl, bw = vc.length * cw;
      ctx.strokeStyle = "#000"; ctx.lineWidth = 1.5; ctx.strokeRect(bx, HEAD - st - 1, bw, rowOrder.length * rh + 2);
    }
    ctx.restore();
    // header: ruler, conservation, region markers
    ctx.fillStyle = "#fff"; ctx.fillRect(LABW, 0, w - LABW, HEAD);
    ctx.save(); ctx.beginPath(); ctx.rect(LABW, 0, w - LABW, HEAD); ctx.clip();
    ctx.font = "10px sans-serif"; ctx.textAlign = "center"; ctx.textBaseline = "alphabetic";
    for (let c = c0; c <= c1; c++) {
      const px = LABW + c * cw - sl, cons = D.alignment.conservation[c], g = Math.round(255 - 215 * cons);
      ctx.fillStyle = `rgb(${g},${g},${g})`; ctx.fillRect(px, HEAD - 14, cw - (cw > 4 ? 1 : 0), 12);
      if (nqSet.has(c + 1)) { ctx.fillStyle = "#9ec5e8"; ctx.fillRect(px, HEAD - 18, cw, 3); }
      if (vrSet.has(c + 1)) { ctx.fillStyle = "#000"; ctx.fillRect(px, HEAD - 18, cw, 3); }
      if ((c + 1) % 10 === 0 && cw >= 4) { ctx.fillStyle = "#555"; ctx.fillText(String(c + 1), px + cw / 2, HEAD - 22); }
    }
    ctx.restore();
    // labels (fixed left column)
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, LABW, h); ctx.textBaseline = "middle";
    ctx.font = "bold 10px sans-serif"; ctx.textAlign = "left"; ctx.fillStyle = "#6b7684";
    ctx.fillText("sequence", 16, HEAD - 8); ctx.fillText("assigned", 152, HEAD - 8); ctx.fillText("suggested", 222, HEAD - 8); ctx.textAlign = "right"; ctx.fillText("reads", LABW - 8, HEAD - 8);
    ctx.save(); ctx.beginPath(); ctx.rect(0, HEAD, LABW, h - HEAD); ctx.clip();
    const lf = Math.max(8, Math.min(12, rh - 3));
    for (let r = r0; r <= r1; r++) {
      const id = rowOrder[r], y = HEAD + r * rh - st, mid = y + rh / 2, dim = matches(id) ? 1 : 0.35, ref = isRef(id);
      ctx.globalAlpha = dim;
      if (S.sel.has(id)) { ctx.fillStyle = "#fff3b0"; ctx.fillRect(0, y, LABW, rh); }
      ctx.fillStyle = tipColour(id); ctx.fillRect(2, y + 2, 5, rh - 4);
      if (!ref && seqById[id].flags.length) { ctx.fillStyle = "#b7791f"; ctx.font = lf + "px sans-serif"; ctx.textAlign = "left"; ctx.fillText("▲", 9, mid); }
      ctx.font = (ref ? "bold " : "") + lf + "px sans-serif"; ctx.textAlign = "left"; ctx.fillStyle = "#1f2933";
      ctx.fillText(ref ? "REF " + refInfo[id].name : id, 22, mid);
      if (!ref) {
        const a = assigned(id);
        if (a) { ctx.font = "bold " + (lf - 1) + "px sans-serif"; const at = a.length > 9 ? a.slice(0, 8) + "…" : a, tw = ctx.measureText(at).width + 8; ctx.fillStyle = colourOf(a); ctx.fillRect(150, y + 1, tw, rh - 2); ctx.fillStyle = "#fff"; ctx.fillText(at, 154, mid); }
        const sg = suggested(id), conf = confidence(id);
        ctx.font = (lf - 1) + "px sans-serif"; ctx.fillStyle = conf === "novel" ? "#c53030" : conf === "high" ? "#2f855a" : "#6b7684";
        const t = sg.replace("possible new subtype", "new?"); ctx.fillText(t.length > 17 ? t.slice(0, 16) + "…" : t, 222, mid);
        ctx.textAlign = "right"; ctx.fillStyle = "#6b7684"; ctx.fillText(String(reads(id)), LABW - 8, mid);
      }
      ctx.globalAlpha = 1;
      ctx.strokeStyle = "#eef0f3"; ctx.beginPath(); ctx.moveTo(0, y + rh - 0.5); ctx.lineTo(LABW, y + rh - 0.5); ctx.stroke();
      if (rowGroup && r > 0 && rowGroup[id] !== rowGroup[rowOrder[r - 1]]) { ctx.fillStyle = "#6b7684"; ctx.fillRect(0, y - 1, LABW, 2); }
    }
    ctx.restore();
    ctx.strokeStyle = "#d5d9e0"; ctx.beginPath(); ctx.moveTo(LABW - 0.5, 0); ctx.lineTo(LABW - 0.5, h); ctx.stroke();
  }

  // pointer handling on the canvas
  let dragging = false;
  function rowAt(ev) { const r = cv.getBoundingClientRect(), y = ev.clientY - r.top, sc = $("#alnScroll"); if (y < HEAD) return null; const i = Math.floor((y - HEAD + sc.scrollTop) / S.rowH); return i >= 0 && i < rowOrder.length ? rowOrder[i] : null; }
  cv.addEventListener("mousedown", (ev) => { const id = rowAt(ev); if (!id) return; dragging = ev.clientX - cv.getBoundingClientRect().left < LABW && !ev.shiftKey && !ev.ctrlKey && !ev.metaKey; clickRow(id, ev, "aln"); });
  window.addEventListener("mouseup", () => { dragging = false; });
  cv.addEventListener("mousemove", (ev) => {
    const id = rowAt(ev); if (!id) { $("#status").textContent = ""; return; }
    if (dragging && S.anchor) { const a = rowIndex[S.anchor], b = rowIndex[id], [s, e] = a < b ? [a, b] : [b, a]; select(rowOrder.slice(s, e + 1), "set", id); }
    const r = cv.getBoundingClientRect(), sc = $("#alnScroll"), x0 = ev.clientX - r.left;
    let t = (isRef(id) ? "REF " + refInfo[id].name : id);
    if (x0 >= LABW) { const c = Math.floor((x0 - LABW + sc.scrollLeft) / S.cellW); if (c >= 0 && c < ncol) t += "  |  column " + (c + 1) + "  |  residue " + D.alignment.seq[id][c] + (vrSet.has(c + 1) ? "  |  variable region" : nqSet.has(c + 1) ? "  |  NQQHNN zone" : ""); }
    else if (!isRef(id) && x0 >= ASG_X0 && x0 <= ASG_X1) t += "  |  double-click to edit the assignment, Delete clears the selection";
    else if (!isRef(id)) t += "  |  suggested " + suggested(id) + " (" + confidence(id) + ")  |  " + esc(x(id, "basis") || "");
    $("#status").textContent = t;
  });
  const ASG_X0 = 148, ASG_X1 = 226;       // the Assigned cell of a label row, in canvas pixels
  cv.addEventListener("dblclick", (ev) => {
    const id = rowAt(ev); if (!id || isRef(id)) return;
    const x0 = ev.clientX - cv.getBoundingClientRect().left;
    if (x0 >= ASG_X0 && x0 <= ASG_X1) return openEditor(id, false);
    const a = suggested(id); if (!/^possible new/.test(a) && !assigned(id)) commit([id], (c) => Object.assign(c, { Subtype: canon(a) }), "accept suggestion");
  });
  // inline editor: a text field laid over the Assigned cell. Enter saves (an empty value clears), Tab saves and moves
  // to the next unassigned row, Escape cancels. With toSelection the value is applied to every selected sequence.
  let editor = null;
  function closeEditor() { if (editor) { const e = editor; editor = null; e.remove(); } }
  function openEditor(id, toSelection) {
    closeEditor();
    if (rowIndex[id] === undefined || isRef(id)) return;
    const sc = $("#alnScroll"); scrollToRow(sc, rowIndex[id]);
    const pane = $("#alnPane"), pr = pane.getBoundingClientRect(), cr = cv.getBoundingClientRect(), rh = S.rowH;
    const y = cr.top - pr.top + HEAD + rowIndex[id] * rh - sc.scrollTop, h = Math.max(20, rh + 2);
    if (y < HEAD - 2 || y > cv.clientHeight) return;
    let dl = $("#dlInline"); if (!dl) { dl = document.createElement("datalist"); dl.id = "dlInline"; document.body.appendChild(dl); }
    dl.innerHTML = vocabulary().map((v) => `<option value="${esc(v.label)}">`).join("");
    const inp = document.createElement("input"); inp.type = "text"; inp.className = "inline-edit"; inp.setAttribute("list", "dlInline");
    inp.value = assigned(id); inp.placeholder = toSelection ? "label for " + [...S.sel].filter((i) => !isRef(i)).length + " selected" : "Subtype";
    Object.assign(inp.style, { left: cr.left - pr.left + ASG_X0 + "px", top: y - (h - rh) / 2 + "px", width: "190px", height: h + "px" });
    pane.appendChild(inp); editor = inp; inp.focus(); inp.select();
    const targets = () => (toSelection ? [...S.sel].filter((i) => !isRef(i)) : [id]);
    const save = () => {
      const val = canon(inp.value), ids = targets(); closeEditor();
      commit(ids, (c) => { c.Subtype = val; return c; }, val ? "assign " + val : "clear");
    };
    inp.addEventListener("keydown", (ev) => {
      ev.stopPropagation();
      if (ev.key === "Enter") { ev.preventDefault(); save(); }
      else if (ev.key === "Escape") { ev.preventDefault(); closeEditor(); }
      else if (ev.key === "Tab") { ev.preventDefault(); save(); nextUnassigned(ev.shiftKey ? -1 : 1); if (S.focus && !assigned(S.focus)) openEditor(S.focus, false); }
    });
    inp.addEventListener("blur", () => setTimeout(closeEditor, 120));
  }

  // scroll handling and sync between the two panes
  let syncing = false;
  $("#alnScroll").addEventListener("scroll", () => { closeEditor(); requestDraw(); if (treeAligned() && !syncing) { syncing = true; $("#treeScroll").scrollTop = $("#alnScroll").scrollTop; syncing = false; } });
  $("#treeScroll").addEventListener("scroll", () => { if (treeAligned() && !syncing) { syncing = true; $("#alnScroll").scrollTop = $("#treeScroll").scrollTop; syncing = false; } });
  new ResizeObserver(() => { sizeCanvas(); requestDraw(); }).observe($("#alnScroll"));
  new ResizeObserver(() => { renderTree(); }).observe($("#treeScroll"));

  // ---- helpers for labels -------------------------------------------------------------------------------------
  const canon = (t) => Core.canonLabel(t, existingLabels());

  // ---- top bar -------------------------------------------------------------------------------------------------
  function renderTop() {
    const bar = $("#topbar"), n = reps.length, a = nAssigned(), missing = n - a;
    if (!bar.dataset.built) {
      bar.dataset.built = "1";
      bar.innerHTML = `<h1>${esc(D.analysis_name)}</h1>
      <div class="grp"><span id="progress"></span><button id="bNext" title="N">Next unassigned</button></div>
      <div class="grp"><button id="bUndo" title="Ctrl+Z">Undo</button><button id="bRedo">Redo</button></div>
      <div class="grp"><label>Order <select id="fOrder"><option value="suggested">suggested subtype</option><option value="assigned">assigned subtype</option><option value="tree">tree</option><option value="file">alignment file</option><option value="reads">reads</option><option value="workbook">workbook</option><option value="status">unassigned first</option></select></label>
        <label>Tree <select id="fTree"><option value="clado">cladogram</option><option value="dist">distance</option></select></label>
        <label>Tips <select id="fColourBy"><option value="suggested">suggested</option><option value="assigned">assigned</option></select></label>
        <label>Alignment <select id="fAln"><option value="clustalx">ClustalX</option><option value="differs">differs from nearest ref</option><option value="none">plain</option></select></label></div>
      <div class="grp"><label title="Confidence tier of the pipeline's suggestion">Show <select id="fConf"><option value="">all</option><option value="high">high</option><option value="check">check</option><option value="novel">novel</option></select></label>
        <label title="Only sequences with one or more flags (marked with a triangle): internal stop, deletion in a conserved region, missing ends, NQQHNN deletion, very short variable region, possible chimera"><input type="checkbox" id="fFlag">has flags</label><label title="Only sequences that have no Subtype yet"><input type="checkbox" id="fUn">unassigned</label><label title="Only sequences where the Subtype you assigned is different from the pipeline's suggestion, to review your overrides"><input type="checkbox" id="fDiff">differs from suggestion</label><label title="Hide the rows that do not match the filters. Untick to dim them instead and keep the tree and alignment aligned"><input type="checkbox" id="fHide" checked>hide others</label>
        <input type="text" id="fText" placeholder="search (/)" ><label title="Only sequences with at least this many reads">reads &ge; <input type="number" id="fMin" value="0" min="0" style="width:64px"></label></div>
      <div class="grp"><button id="bRowDn">Rows -</button><button id="bRowUp">Rows +</button><button id="bColDn">Cols -</button><button id="bColUp">Cols +</button></div>
      <div class="grp"><button id="bImport">Import</button><button class="primary" id="bXlsx">Export xlsx</button><button id="bTsv">TSV</button><button id="bJson">Session</button><button id="bHelp">Help</button><button id="bPanel" title="Show or hide the side panel">Panel</button></div>`;
      const on = (id, ev, fn) => $("#" + id).addEventListener(ev, fn);
      on("bNext", "click", nextUnassigned); on("bUndo", "click", undo); on("bRedo", "click", redo);
      on("fOrder", "change", (e) => { setOrder(e.target.value); relayout(); treeAligned() && ($("#treeScroll").scrollTop = $("#alnScroll").scrollTop); });
      on("fTree", "change", (e) => { S.treeMode = e.target.value; renderTree(); });
      on("fColourBy", "change", (e) => { S.colourBy = e.target.value; renderTree(); requestDraw(); });
      on("fAln", "change", (e) => { S.alnColour = e.target.value; requestDraw(); });
      const refilter = () => { refreshRows(); };
      on("fConf", "change", (e) => { S.f.conf = e.target.value; refilter(); });
      on("fFlag", "change", (e) => { S.f.flagged = e.target.checked; refilter(); });
      on("fUn", "change", (e) => { S.f.unassigned = e.target.checked; refilter(); });
      on("fDiff", "change", (e) => { S.f.differs = e.target.checked; refilter(); });
      on("fHide", "change", (e) => { S.hide = e.target.checked; refilter(); });
      on("fText", "input", (e) => { S.f.text = e.target.value.trim().toLowerCase(); refilter(); });
      on("fMin", "input", (e) => { S.f.minReads = +e.target.value || 0; refilter(); });
      on("bRowDn", "click", () => zoom(-2, 0)); on("bRowUp", "click", () => zoom(2, 0)); on("bColDn", "click", () => zoom(0, -2)); on("bColUp", "click", () => zoom(0, 2));
      on("bImport", "click", importDialog); on("bXlsx", "click", () => exportFile("xlsx")); on("bTsv", "click", () => exportFile("tsv")); on("bJson", "click", exportSession); on("bHelp", "click", helpDialog); on("bPanel", "click", () => { const s = $("#side"); s.style.display = s.style.display === "none" ? "" : "none"; });
    }
    const p = $("#progress"); p.textContent = `assigned ${a} of ${n}` + (missing ? `, ${missing} missing` : ", complete"); p.classList.toggle("done", missing === 0);
    $("#bUndo").disabled = !S.undo.length; $("#bRedo").disabled = !S.redo.length;
  }
  function zoom(dr, dc) { S.rowH = Math.max(6, Math.min(26, S.rowH + dr)); S.cellW = Math.max(3, Math.min(24, S.cellW + dc)); relayout(); }
  function relayout() { sizeCanvas(); renderTree(); requestDraw(); }
  function nextUnassigned(dir) {
    const start = S.focus && rowIndex[S.focus] !== undefined ? rowIndex[S.focus] : -1, n = rowOrder.length, step = dir === -1 ? -1 : 1;
    for (let k = 1; k <= n; k++) { const id = rowOrder[((start + step * k) % n + n) % n]; if (!isRef(id) && !assigned(id) && matches(id)) { select([id], "set", id); S.anchor = id; reveal(); return; } }
    toast("Every visible sequence has a Subtype");
  }
  function stepRow(d) {
    const start = S.focus && rowIndex[S.focus] !== undefined ? rowIndex[S.focus] : -1, n = rowOrder.length;
    for (let k = 1; k <= n; k++) { const id = rowOrder[((start + d * k) % n + n) % n]; if (matches(id)) { select([id], "set", id); S.anchor = id; reveal(); return; } }
  }

  // ---- side panel ----------------------------------------------------------------------------------------------
  const touched = {};
  function renderSide() {
    const side = $("#side"), selIds = [...S.sel].filter((i) => !isRef(i)), nsel = selIds.length, f = S.focus, vocab = vocabulary();
    const common = (k) => { const vals = new Set(selIds.map((i) => (asg(i) || blank())[k])); return vals.size === 1 ? [...vals][0] : ""; };
    const mixed = (k) => selIds.length > 1 && new Set(selIds.map((i) => (asg(i) || blank())[k])).size > 1;
    let h = "";
    if (!side.dataset.seen) h += `<div id="banner">Assignments are kept in this browser while you work, but <b>Export xlsx</b> is the real save. Send or keep the exported file.</div>`;
    h += `<h2>Assign</h2><div class="muted" style="margin-bottom:4px">${nsel ? nsel + " selected sequence" + (nsel === 1 ? "" : "s") : "Select sequences in the tree or the alignment"}</div>`;
    h += `<div class="row"><label>Subtype</label><input type="text" id="aSub" list="dlSub" value="${esc(common("Subtype"))}" placeholder="${mixed("Subtype") ? "(mixed)" : ""}" ${nsel ? "" : "disabled"}></div>`;
    h += `<div class="row"><label>Endogenous</label><input type="text" id="aEnd" list="dlEnd" value="${esc(common("Endogenous"))}" placeholder="${mixed("Endogenous") ? "(mixed)" : ""}" ${nsel ? "" : "disabled"}></div>`;
    h += `<div class="row"><label>Unsure</label><input type="text" id="aUns" value="${esc(common("Unsure"))}" placeholder="${mixed("Unsure") ? "(mixed)" : ""}" ${nsel ? "" : "disabled"}><button id="aQ" ${nsel ? "" : "disabled"}>?</button></div>`;
    h += `<div class="row"><label>Note</label><input type="text" id="aNote" value="${esc(common("Note"))}" placeholder="${mixed("Note") ? "(mixed)" : ""}" ${nsel ? "" : "disabled"}></div>`;
    h += `<datalist id="dlSub">${vocab.map((v) => `<option value="${esc(v.label)}">`).join("")}</datalist><datalist id="dlEnd">${[...new Set(reps.map((i) => (asg(i) || {}).Endogenous).filter(Boolean))].map((v) => `<option value="${esc(v)}">`).join("")}</datalist>`;
    h += `<div class="btns"><button class="primary" id="aApply" ${nsel ? "" : "disabled"}>Apply to ${nsel || ""}</button><button id="aAccept" ${nsel ? "" : "disabled"} title="A">Accept suggestion</button><button id="aClear" ${nsel ? "" : "disabled"}>Clear</button></div>`;
    const ql = quickLabels(), dis = nsel ? "" : "disabled";
    h += `<div class="btns">${ql.numbered.map((l, i) => `<button class="q" data-l="${esc(l)}" ${dis} title="key ${i + 1}">${i + 1} ${esc(l)}</button>`).join("")}<button class="q" data-l="Non-functional" ${dis} title="key 0">0 Non-functional</button></div>`;
    if (ql.extra.length) h += `<div class="btns">${ql.extra.map((l) => `<button class="q" data-l="${esc(l)}" ${dis}>${esc(l)}</button>`).join("")}</div>`;
    h += `<h2>Select</h2><div class="btns"><button id="sUn">Unassigned</button><button id="sFlag">Flagged</button><button id="sNovel">Novel</button><button id="sVis">All visible</button><button id="sNone">None</button><button id="sLabels">Labels...</button></div>`;
    if (f && (seqById[f] || refInfo[f])) h += detailHtml(f);
    side.innerHTML = h; side.dataset.seen = "1";
    Object.keys(touched).forEach((k) => delete touched[k]);
    const bind = (id, k) => { const e = $("#" + id); if (e) e.addEventListener("input", () => { touched[k] = true; }); };
    bind("aSub", "Subtype"); bind("aEnd", "Endogenous"); bind("aUns", "Unsure"); bind("aNote", "Note");
    const on = (id, fn) => { const e = $("#" + id); if (e) e.addEventListener("click", fn); };
    on("aQ", () => { $("#aUns").value = $("#aUns").value === "?" ? "" : "?"; touched.Unsure = true; });
    on("aApply", applyForm); on("aAccept", acceptSuggestion);
    on("aClear", clearSelection);
    side.querySelectorAll("button.q").forEach((b) => b.addEventListener("click", () => assignLabel(b.dataset.l)));
    ["aSub", "aEnd", "aUns", "aNote"].forEach((id) => { const e = $("#" + id); if (e) e.addEventListener("keydown", (ev) => { if (ev.key === "Enter") applyForm(); }); });
    on("sUn", () => select(reps.filter((i) => !assigned(i) && matches(i)), "set")); on("sFlag", () => select(reps.filter((i) => seqById[i].flags.length), "set"));
    on("sNovel", () => select(reps.filter((i) => confidence(i) === "novel"), "set")); on("sVis", () => select(rowOrder.filter(matches), "set")); on("sNone", () => select([], "set"));
    on("sLabels", labelsDialog);
  }
  function detailHtml(id) {
    if (isRef(id)) return `<h2>Reference</h2><table><tr><td>name</td><td>${esc(refInfo[id].name)}</td></tr><tr><td>subtype</td><td><span class="chip" style="background:${colourOf(refInfo[id].subtype)}">${esc(refInfo[id].subtype)}</span></td></tr></table>`;
    const s = seqById[id], e = s.export, fx = (v, d) => (typeof v === "number" ? v.toFixed(d === undefined ? 3 : d) : esc(v || ""));
    const top = (e.library_top_hits || "").split("; ").filter(Boolean).map(esc).join("<br>");
    const ps = Object.entries(s.per_sample).sort((a, b) => b[1] - a[1]).map(([k, v]) => `${esc(k)}: ${v}`).join("<br>");
    const rows = [
      ["suggested", `<span class="chip" style="background:${colourOf(e.suggested_subtype)}">${esc(e.suggested_subtype)}</span> ${esc(e.confidence)}`],
      ["basis", esc(e.basis)],
      ["flags", s.flags.length ? `<span class="flag">${s.flags.map(esc).join(", ")}</span>` : "none"],
      ["reads", `${esc(e.Total_number_of_reads)} in ${esc(e.N_samples_present)} sample(s)`],
      ["nearest reference", `${esc(e.nearest_reference)} (${fx(e.nearest_reference_identity)})`],
      ["other subtype", `${esc(e.other_subtype_reference)} (${fx(e.other_subtype_identity)}), margin ${fx(e.margin_to_other_subtype)}`],
      ["conserved identity", fx(e.conserved_identity)],
      ["variable region", `${fx(e.variable_region_identity)} to ${esc(e.variable_region_nearest_reference)}; ${esc(e.variable_region_residues)} residues, ${esc(e.variable_region_missing_bp)} bp missing`],
      ["NQQHNN zone", `${fx(e.nqqhnn_identity)}; ${esc(e.missing_nqqhnn_bp)} bp missing`],
      ["stops / ends", `${esc(e.n_internal_stops)} internal stops; missing 5' ${esc(e.terminal_missing_5_bp)} bp, 3' ${esc(e.terminal_missing_3_bp)} bp; conserved deletion ${esc(e.deletion_conserved_bp)} bp; inserted ${esc(e.inserted_residues)}`],
      ["library", `${e.library_exact_match ? "exact match; " : ""}best ${fx(e.library_best_identity)} ${esc(e.library_labels_best)}${e.library_endogenous ? " (endogenous " + esc(e.library_endogenous) + ")" : ""}`],
      ["library hits", top || "none"], ["per sample", ps || "none"],
    ];
    return `<h2>${esc(id)}</h2><table>${rows.map((r) => `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join("")}</table><h2>Protein</h2><div class="seq">${esc(e.Sequence)}</div>`;
  }
  function applyForm() {
    const selIds = [...S.sel].filter((i) => !isRef(i)); if (!selIds.length) return;
    const val = { Subtype: $("#aSub").value, Endogenous: $("#aEnd").value, Unsure: $("#aUns").value, Note: $("#aNote").value };
    const keys = FILL_COLS.filter((k) => touched[k]); if (!keys.length) return toast("Nothing to apply: type a value first");
    if (touched.Subtype) val.Subtype = canon(val.Subtype), val.Subtype !== $("#aSub").value.trim() && val.Subtype && toast("Using existing spelling: " + val.Subtype);
    if (touched.Endogenous) val.Endogenous = canonEndo(val.Endogenous);
    commit(selIds, (c) => { keys.forEach((k) => { c[k] = String(val[k]).trim(); }); return c; }, "assign " + selIds.length);
  }
  function canonEndo(t) { const have = [...new Set(reps.map((i) => (asg(i) || {}).Endogenous).filter(Boolean))]; return Core.labelKey(t) ? Core.canonLabel(t, have) : ""; }
  function clearSelection() {
    const selIds = [...S.sel].filter((i) => !isRef(i) && asg(i));
    if (!selIds.length) return toast("The selected sequences have no assignment to clear");
    commit(selIds, () => blank(), "clear " + selIds.length); toast("Cleared " + selIds.length + " assignment(s)");
  }
  function assignLabel(l) { const selIds = [...S.sel].filter((i) => !isRef(i)); if (selIds.length) commit(selIds, (c) => Object.assign(c, { Subtype: l }), "assign " + l); }
  function acceptSuggestion() {
    const selIds = [...S.sel].filter((i) => !isRef(i)), ok = selIds.filter((i) => !/^possible new/.test(suggested(i))), skipped = selIds.length - ok.length;
    if (!ok.length) return toast(skipped ? "Possible new subtypes need an explicit label" : "Select sequences first");
    commit(ok, (c, id) => Object.assign(c, { Subtype: canon(suggested(id)) }), "accept " + ok.length);
    if (skipped) toast(skipped + " possible new subtype(s) skipped: assign them yourself");
  }

  // ---- dialogs -----------------------------------------------------------------------------------------------------
  function toast(msg) { document.querySelectorAll(".toast").forEach((e) => e.remove()); const t = document.createElement("div"); t.className = "toast"; t.textContent = msg; document.body.appendChild(t); setTimeout(() => t.remove(), 3200); }
  function dialog(title, bodyHtml, buttons) {
    const d = $("#dialog"); d.hidden = false;
    d.innerHTML = `<div class="box"><h3>${esc(title)}</h3>${bodyHtml}<div class="btns">${buttons.map((b, i) => `<button data-i="${i}" class="${b.primary ? "primary" : ""}">${esc(b.label)}</button>`).join("")}</div></div>`;
    d.querySelectorAll(".btns button").forEach((b) => b.addEventListener("click", () => { const fn = buttons[+b.dataset.i].fn; d.hidden = true; d.innerHTML = ""; fn && fn(); }));
    return d;
  }
  const close = () => { const d = $("#dialog"); d.hidden = true; d.innerHTML = ""; };
  function helpDialog() {
    dialog("How to use this page", `<p>Click a tip, a branch (selects the clade), a row label or any alignment cell to select sequences. Selections show in both views. Shift+click selects a range, Ctrl/Cmd+click toggles, dragging over the labels selects a range.</p>
      <p>Type a Subtype in the side panel and press Enter or <b>Apply</b>, or use <b>Accept suggestion</b> (possible new subtypes always need an explicit label). You can also edit in the table: <b>double-click the Assigned cell</b> of a row (or press <kbd>Enter</kbd> to edit the selected rows), type a label and press Enter (<kbd>Tab</kbd> saves and moves to the next unassigned row, <kbd>Esc</kbd> cancels, an empty label clears it). Double-click elsewhere in an unassigned row to accept its suggestion. <kbd>Delete</kbd> or <kbd>Backspace</kbd> clears the assignment of the selected rows.</p><p><b>Filters</b> (they combine, and by default hide the rows that do not match):<br><b>Show</b>: the confidence tier of the suggestion (high, check or novel).<br><b>has flags</b>: sequences with one or more flags, marked with a triangle (internal stop, deletion in a conserved region, missing ends, NQQHNN deletion, very short variable region, possible chimera). The flags are listed in the side panel.<br><b>unassigned</b>: sequences with no Subtype yet.<br><b>differs from suggestion</b>: sequences where the Subtype you assigned is not the one the pipeline suggested, to review your overrides.<br><b>hide others</b>: hide the rows that do not match; untick to dim them instead.<br><b>search</b> matches IDs, suggested and assigned labels; <b>reads</b> sets a minimum read count.</p>
      <p><b>Export xlsx</b> writes the suggestions workbook with Subtype, Endogenous, Unsure and Note filled in. Give that file to the pipeline as <code>subtype_assignments</code>. <b>Import</b> reads such a workbook, a TSV or a session file, so you can continue or merge work.</p>
      <p>Keys: <kbd>J</kbd>/<kbd>K</kbd> next/previous row, <kbd>N</kbd> next unassigned, <kbd>A</kbd> accept suggestion, <kbd>1</kbd>-<kbd>9</kbd> quick labels (always the first nine subtypes, A to I), <kbd>0</kbd> Non-functional, <kbd>Enter</kbd> edit assignment in the table, <kbd>Delete</kbd> clear assignment, <kbd>E</kbd> <kbd>U</kbd> focus Endogenous/Unsure, <kbd>/</kbd> search, <kbd>C</kbd> toggle cladogram/distance, <kbd>Ctrl+Z</kbd> undo, <kbd>Esc</kbd> clear selection.</p>
      <p class="muted">Triangles mark sequences with flags, a red ring marks novel suggestions. Black box: variable region. Blue bar: NQQHNN zone.</p>`, [{ label: "Close", primary: true }]);
  }
  function labelsDialog() {
    const v = vocabulary().filter((x) => x.n > 0);
    if (!v.length) return toast("No labels assigned yet");
    let inputs = [];
    const d = dialog("Labels in use", `<p class="muted">Rename a label everywhere, or give two labels the same name to merge them.</p><table><tr><th>Label</th><th>Sequences</th><th>Rename to</th></tr>${v.map((x) => `<tr><td>${esc(x.label)}</td><td>${x.n}</td><td><input type="text" data-from="${esc(x.label)}" value="${esc(x.label)}"></td></tr>`).join("")}</table>`,
      [{ label: "Cancel" }, { label: "Apply", primary: true, fn: () => {
        const map = {}; inputs.forEach((i) => { const to = i.value.trim(); if (to && to !== i.dataset.from) map[i.dataset.from] = to; });
        const idl = reps.filter((id) => map[assigned(id)]);
        if (idl.length) commit(idl, (c) => Object.assign(c, { Subtype: map[c.Subtype] }), "rename labels");
      } }]);
    inputs = [...d.querySelectorAll("input[data-from]")];     // detached inputs keep their values after the dialog closes
  }

  // ---- export and import -----------------------------------------------------------------------------------------------
  function exportRows() {
    const cols = D.export_columns;
    const rows = [cols.slice()];
    D.sequences.forEach((s) => {
      const a = Object.assign(blank(), S.assign[s.id] || {});
      rows.push(cols.map((c) => (FILL_COLS.includes(c) ? a[c] : s.export[c])));
    });
    return rows;
  }
  function download(name, bytes, type) {
    const url = URL.createObjectURL(new Blob([bytes], { type })), a = document.createElement("a");
    a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function checkBeforeExport(proceed) {
    const missing = reps.filter((id) => !assigned(id)), warn = [];
    const labs = [...new Set(reps.map(assigned).filter(Boolean))], keys = {};
    labs.forEach((l) => { (keys[Core.labelKey(l)] = keys[Core.labelKey(l)] || []).push(l); });
    Object.values(keys).filter((v) => v.length > 1).forEach((v) => warn.push("Labels differing only in case or spacing: " + v.join(", ")));
    const newish = reps.filter((id) => /^possible new/i.test(assigned(id))); if (newish.length) warn.push(newish.length + " Subtype value(s) still start with \"possible new subtype\"");
    if (!missing.length && !warn.length) return proceed(false);
    const body = (missing.length ? `<p><b>${missing.length} sequence(s) have no Subtype</b>, so the pipeline will refuse this file until they are filled in:</p><p class="muted">${missing.slice(0, 12).map(esc).join(", ")}${missing.length > 12 ? ", ..." : ""}</p>` : "") + warn.map((w) => `<p>${esc(w)}</p>`).join("");
    dialog("Check before export", body, [{ label: "Cancel" }, { label: missing.length ? "Export as incomplete" : "Export anyway", primary: true, fn: () => proceed(missing.length > 0) }]);
  }
  function exportFile(kind) {
    checkBeforeExport((incomplete) => {
      const rows = exportRows(), name = D.analysis_name + "_subtype_assignments" + (incomplete ? "_INCOMPLETE" : "") + "." + kind;
      if (kind === "tsv") return download(name, Core.tsvSerialize(rows), "text/tab-separated-values");
      const widths = D.export_columns.map((c, i) => Math.min(48, Math.max(c.length + 2, ...rows.slice(1, 40).map((r) => String(r[i] === null || r[i] === undefined ? "" : r[i]).length + 1))));
      const bytes = Core.xlsxBuild([
        { name: "suggestions", rows, header: true, freeze: { rows: 1, cols: 3 }, fillCols: FILL_COLS.map((c) => colIdx[c]).filter((i) => i !== undefined), widths },
        { name: "legend", rows: [["column", "meaning"]].concat(D.legend.map((l) => [l.column, l.meaning])), header: true, widths: [45, 120] }]);
      download(name, bytes, "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");
      toast("Saved " + name);
    });
  }
  function exportSession() { download(D.analysis_name + "_review_session.json", JSON.stringify({ fingerprint: D.fingerprint, analysis_name: D.analysis_name, assign: S.assign }, null, 1), "application/json"); toast("Session saved"); }
  function importDialog() {
    const d = dialog("Import assignments", `<p>Choose an exported workbook (.xlsx), a TSV or CSV with a <code>new_rep_97_ID</code> (or <code>seq_md5</code>) column and a Subtype column, or a session file (.json).</p><input type="file" id="impFile" accept=".xlsx,.tsv,.txt,.csv,.json">`, [{ label: "Cancel" }]);
    $("#impFile").addEventListener("change", async (ev) => {
      const f = ev.target.files[0]; if (!f) return;
      try { close(); await importFile(f); } catch (e) { dialog("Import failed", `<p>${esc(e.message || e)}</p>`, [{ label: "Close", primary: true }]); }
    });
  }
  async function importFile(f) {
    let incoming = {}, other = 0, note = "";
    if (/\.json$/i.test(f.name)) {
      const j = JSON.parse(await f.text());
      if (j.fingerprint && j.fingerprint !== D.fingerprint) note = "This session file comes from a different run (the sequences differ). ";
      Object.entries(j.assign || {}).forEach(([id, v]) => { if (seqById[id]) incoming[id] = Object.assign(blank(), v); else other++; });
    } else {
      let rows;
      if (/\.xlsx$/i.test(f.name)) rows = await Core.xlsxRead(new Uint8Array(await f.arrayBuffer()));
      else { const t = await f.text(); rows = /\.csv$/i.test(f.name) ? t.split(/\r?\n/).filter(Boolean).map((l) => l.split(",")) : Core.tsvParse(t); }
      const head = rows[0].map((h) => String(h).trim()), ci = (n) => head.indexOf(n);
      const kId = ci("new_rep_97_ID") >= 0 ? ci("new_rep_97_ID") : ci("seq_md5"), byMd5 = ci("new_rep_97_ID") < 0;
      if (kId < 0 || ci("Subtype") < 0) throw new Error("The file needs a new_rep_97_ID (or seq_md5) column and a Subtype column.");
      const md5 = {}; D.sequences.forEach((s) => { md5[s.export.seq_md5] = s.id; });
      rows.slice(1).forEach((r) => {
        const id = byMd5 ? md5[r[kId]] : r[kId];
        if (!id || !seqById[id]) { other++; return; }
        const v = blank(); FILL_COLS.forEach((k) => { if (ci(k) >= 0) v[k] = String(r[ci(k)] || "").trim(); });
        if (v.Subtype || v.Endogenous || v.Unsure || v.Note) incoming[id] = v;
      });
    }
    const idl = Object.keys(incoming), conflicts = idl.filter((id) => assigned(id) && incoming[id].Subtype && assigned(id) !== incoming[id].Subtype);
    const msg = `${note}<p>${idl.length} sequence(s) with assignments found${other ? `, ${other} row(s) not in this run were ignored` : ""}.</p>` + (conflicts.length ? `<p><b>${conflicts.length}</b> differ from what you have now:</p><p class="muted">${conflicts.slice(0, 10).map((id) => `${esc(id)}: you ${esc(assigned(id))}, file ${esc(incoming[id].Subtype)}`).join("<br>")}${conflicts.length > 10 ? "<br>..." : ""}</p>` : "");
    const run = (takeTheirs) => commit(idl.filter((id) => takeTheirs || !conflicts.includes(id)), (c, id) => { const v = incoming[id]; FILL_COLS.forEach((k) => { if (v[k] || takeTheirs) c[k] = v[k] || c[k]; }); if (v.Subtype) c.Subtype = canon(v.Subtype) || v.Subtype; return c; }, "import");
    if (!idl.length) return dialog("Import", msg, [{ label: "Close", primary: true }]);
    dialog("Import assignments", msg, conflicts.length
      ? [{ label: "Cancel" }, { label: "Keep mine where they differ", fn: () => run(false) }, { label: "Take the file's", primary: true, fn: () => run(true) }]
      : [{ label: "Cancel" }, { label: "Import", primary: true, fn: () => run(true) }]);
  }

  // ---- keyboard ----------------------------------------------------------------------------------------------------------
  document.addEventListener("keydown", (ev) => {
    const tag = (ev.target.tagName || "").toLowerCase(), typing = tag === "input" || tag === "textarea" || tag === "select";
    if (!$("#dialog").hidden) { if (ev.key === "Escape") close(); return; }
    if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === "z" && !typing) { ev.preventDefault(); ev.shiftKey ? redo() : undo(); return; }
    if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === "y" && !typing) { ev.preventDefault(); redo(); return; }
    if (ev.key === "Escape") { if (typing) ev.target.blur(); else select([], "set"); return; }
    if (typing || ev.ctrlKey || ev.metaKey || ev.altKey) return;
    const k = ev.key.toLowerCase();
    if (ev.key === "Delete" || ev.key === "Backspace") { ev.preventDefault(); clearSelection(); return; }
    if (ev.key === "Enter") { const id = S.focus && !isRef(S.focus) ? S.focus : null; if (id) { ev.preventDefault(); openEditor(id, S.sel.size > 1); } return; }
    if (k === "j") stepRow(1); else if (k === "k") stepRow(-1); else if (k === "n") nextUnassigned(ev.shiftKey ? -1 : 1);
    else if (k === "a") acceptSuggestion(); else if (k === "e") { const e = $("#aEnd"); if (e) { ev.preventDefault(); e.focus(); } }
    else if (k === "u") { const e = $("#aUns"); if (e) { ev.preventDefault(); e.focus(); } }
    else if (k === "/") { ev.preventDefault(); $("#fText").focus(); }
    else if (k === "c") { S.treeMode = S.treeMode === "clado" ? "dist" : "clado"; $("#fTree").value = S.treeMode; renderTree(); }
    else if (k === "?") helpDialog();
    else if (/^[0-9]$/.test(k)) { const ql = quickLabels(), l = k === "0" ? "Non-functional" : ql.numbered[+k - 1]; if (l) assignLabel(l); }
  });

  // ---- start -----------------------------------------------------------------------------------------------------------------
  setOrder("suggested");
  const restored = restore();
  renderTop(); sizeCanvas(); renderTree(); renderSide(); requestDraw();
  if (restored) toast("Restored " + restored + " assignment(s) saved in this browser");
  window.__korv = { rows: () => rowOrder, S, D, select, commit, exportRows, importFile, layout, matches, draw, renderTree, setOrder, relayout };   // handle for tests
})();
