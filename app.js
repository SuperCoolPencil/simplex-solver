/**
 * UI controller — wires up the DOM to the solver.
 */
(function () {
  const F = (v) => Frac.from(v);

  // ── State ──
  let method = "simplex";

  // ── DOM refs ──
  const $method = document.querySelectorAll(".method-card");
  const $numVars = document.getElementById("num-vars");
  const $numCons = document.getElementById("num-constraints");
  const $optType = document.getElementById("opt-type");
  const $generate = document.getElementById("btn-generate");
  const $form = document.getElementById("dynamic-form");
  const $solve = document.getElementById("btn-solve");
  const $output = document.getElementById("output-container");
  const $clear = document.getElementById("btn-clear");
  const $example = document.getElementById("btn-load-example");
  const $outputSection = document.getElementById("output-section");

  // ── Method selection ──
  $method.forEach((btn) => {
    btn.addEventListener("click", () => {
      $method.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      method = btn.dataset.method;
    });
  });

  // ── Generate dynamic form ──
  function generateForm() {
    const n = +$numVars.value;
    const m = +$numCons.value;
    let html = "";

    // Objective
    html += `<div class="form-label">Objective coefficients</div><div class="obj-row">`;
    for (let j = 0; j < n; j++) {
      html += `<input type="text" id="c${j}" placeholder="c${j + 1}" value="0">`;
      html += `<span class="var-label">x${j + 1}${j < n - 1 ? " +" : ""}</span>`;
    }
    html += `</div>`;

    // Constraints
    html += `<div class="form-label">Constraints</div>`;
    for (let i = 0; i < m; i++) {
      html += `<div class="constraint-row">`;
      for (let j = 0; j < n; j++) {
        html += `<input type="text" id="a${i}_${j}" placeholder="a${i + 1}${j + 1}" value="0">`;
        html += `<span class="var-label">x${j + 1}${j < n - 1 ? " +" : ""}</span>`;
      }
      html += `<select id="ct${i}"><option value="le">≤</option><option value="ge">≥</option><option value="eq">=</option></select>`;
      html += `<input type="text" id="b${i}" placeholder="b${i + 1}" value="0">`;
      html += `</div>`;
    }
    $form.innerHTML = html;
  }

  $generate.addEventListener("click", generateForm);
  generateForm(); // initial

  // ── Load example ──
  $example.addEventListener("click", () => {
    $numVars.value = 3;
    $numCons.value = 3;
    $optType.value = "max";
    // Select ILPP method
    $method.forEach((b) => b.classList.remove("active"));
    document.querySelector('[data-method="ilpp"]').classList.add("active");
    method = "ilpp";
    generateForm();
    // Max Z = 3x1 + x2 + 3x3
    document.getElementById("c0").value = "3";
    document.getElementById("c1").value = "1";
    document.getElementById("c2").value = "3";
    // -x1 + 2x2 + x3 <= 4
    document.getElementById("a0_0").value = "-1";
    document.getElementById("a0_1").value = "2";
    document.getElementById("a0_2").value = "1";
    document.getElementById("ct0").value = "le";
    document.getElementById("b0").value = "4";
    // 4x2 - 3x3 <= 2
    document.getElementById("a1_0").value = "0";
    document.getElementById("a1_1").value = "4";
    document.getElementById("a1_2").value = "-3";
    document.getElementById("ct1").value = "le";
    document.getElementById("b1").value = "2";
    // x1 - 3x2 + 2x3 <= 3
    document.getElementById("a2_0").value = "1";
    document.getElementById("a2_1").value = "-3";
    document.getElementById("a2_2").value = "2";
    document.getElementById("ct2").value = "le";
    document.getElementById("b2").value = "3";
  });

  // ── Clear ──
  $clear.addEventListener("click", () => { $output.innerHTML = ""; });

  // ── Solve ──
  $solve.addEventListener("click", () => {
    $output.innerHTML = "";
    try {
      const { cj, varNames, cb, basis, xb, table, maximize, n } = buildTableau();
      let allSteps = [];

      if (method === "simplex") {
        const res = solveSimplex([...cj], [...varNames], [...cb], [...basis], [...xb], table.map(r => [...r]), maximize);
        allSteps = res.steps;
        if (res.optimal) appendSolution(res, maximize, n);
      } else if (method === "dual") {
        const res = solveDualSimplex([...cj], [...varNames], [...cb], [...basis], [...xb], table.map(r => [...r]));
        allSteps = res.steps;
        if (res.optimal) appendSolution(res, maximize, n);
      } else if (method === "ilpp") {
        // First solve continuous
        const cjCopy = [...cj], vnCopy = [...varNames], cbCopy = [...cb], bCopy = [...basis], xbCopy = [...xb], tCopy = table.map(r => [...r]);
        const simpRes = solveSimplex(cjCopy, vnCopy, cbCopy, bCopy, xbCopy, tCopy, maximize);
        allSteps.push(...simpRes.steps);
        if (simpRes.optimal) {
          appendSolution(simpRes, maximize, n, "Continuous Optimum");
          // Now ILPP
          const ilpRes = solveILPP(simpRes.cj, simpRes.varNames, simpRes.cb, simpRes.basis, simpRes.xb, simpRes.table);
          allSteps.push(...ilpRes.steps);
          if (ilpRes.steps.some(s => s.success)) appendSolution(ilpRes, maximize, n, "Integer Optimum");
        }
      }

      renderSteps(allSteps);
    } catch (e) {
      $output.innerHTML = `<div class="info-card">Error: ${e.message}</div>`;
    }
    $outputSection.scrollIntoView({ behavior: "smooth" });
  });

  // ── Build initial tableau from form ──
  function buildTableau() {
    const n = +$numVars.value;
    const m = +$numCons.value;
    const maximize = $optType.value === "max";

    const obj = [];
    for (let j = 0; j < n; j++) obj.push(F(document.getElementById(`c${j}`).value));

    const constraints = [], ctypes = [], rhs = [];
    for (let i = 0; i < m; i++) {
      const row = [];
      for (let j = 0; j < n; j++) row.push(F(document.getElementById(`a${i}_${j}`).value));
      constraints.push(row);
      ctypes.push(document.getElementById(`ct${i}`).value);
      rhs.push(F(document.getElementById(`b${i}`).value));
    }

    // Add slack/surplus
    const varNames = [];
    for (let j = 0; j < n; j++) varNames.push(`x${j + 1}`);
    const cj = [...obj];

    let sCount = 0;
    for (let i = 0; i < m; i++) {
      if (ctypes[i] === "le") {
        sCount++;
        varNames.push(`S${sCount}`);
        cj.push(Frac.ZERO);
        for (let k = 0; k < m; k++) constraints[k].push(k === i ? Frac.ONE : Frac.ZERO);
      } else if (ctypes[i] === "ge") {
        sCount++;
        varNames.push(`S${sCount}`);
        cj.push(Frac.ZERO);
        for (let k = 0; k < m; k++) constraints[k].push(k === i ? Frac.ONE.neg() : Frac.ZERO);
      }
    }

    const basis = [], cb = [], xb = [...rhs];
    let sIdx = n;
    for (let i = 0; i < m; i++) {
      if (ctypes[i] === "le" || ctypes[i] === "ge") {
        basis.push(varNames[sIdx]);
        cb.push(Frac.ZERO);
        sIdx++;
      } else {
        const aName = `A${i + 1}`;
        varNames.push(aName);
        const bigM = maximize ? new Frac(-1000) : new Frac(1000);
        cj.push(bigM);
        for (let k = 0; k < m; k++) constraints[k].push(k === i ? Frac.ONE : Frac.ZERO);
        basis.push(aName);
        cb.push(bigM);
      }
    }

    return { cj, varNames, cb, basis, xb, table: constraints, maximize, n };
  }

  // ── Render steps ──
  function renderSteps(steps) {
    steps.forEach((s, idx) => {
      const el = document.createElement("div");
      if (s.type === "tableau") {
        el.className = "tableau-card";
        el.style.animationDelay = `${idx * 0.05}s`;
        el.innerHTML = buildTableauHTML(s);
      } else if (s.type === "pivot") {
        el.className = "pivot-info";
        el.innerHTML = `<span><span class="label">Leaving:</span> <span class="val">${s.leaving}</span></span><span><span class="label">Entering:</span> <span class="val">${s.entering}</span></span><span><span class="label">Pivot:</span> <span class="val">${s.element}</span></span>`;
      } else if (s.type === "cut") {
        el.className = "tableau-card";
        el.innerHTML = `<div class="tableau-title"><span class="tag tag-cut">Cut #${s.num}</span> Gomory Cut</div><div style="font-family:var(--mono);font-size:0.82rem;color:var(--text-dim);line-height:1.8">Source: <span style="color:var(--accent)">${s.source} = ${s.value}</span><br>f<sub>i</sub> = ${s.fi}<br>${s.equation}</div>`;
      } else if (s.type === "info") {
        el.className = s.success ? "solution-card" : "info-card";
        el.innerHTML = s.success ? `<h3>✓ ${s.msg}</h3>` : s.msg;
      }
      $output.appendChild(el);
    });
  }

  // ── Build tableau HTML table ──
  function buildTableauHTML(s) {
    const { title, tag, cj, varNames, cb, basis, xb, table, zjcj, ratios, pivotRow, pivotCol } = s;
    const m = basis.length;
    const nv = varNames.length;
    const tagClass = tag === "optimal" ? "tag-optimal" : "tag-iter";
    const tagLabel = tag === "optimal" ? "Optimal" : "Iteration";

    let h = `<div class="tableau-title"><span class="tag ${tagClass}">${tagLabel}</span>${title}</div>`;
    h += `<table class="tableau-table"><thead>`;
    // Cj row
    h += `<tr class="row-cj"><td></td><td></td><td></td>`;
    for (let j = 0; j < nv; j++) h += `<td>${cj[j].toHTML()}</td>`;
    if (ratios) h += `<td></td>`;
    h += `</tr>`;
    // Header
    h += `<tr><th>C<sub>B</sub></th><th style="text-align:left">Basis</th><th>X<sub>B</sub></th>`;
    for (let j = 0; j < nv; j++) h += `<th>${varNames[j]}</th>`;
    if (ratios) h += `<th>Ratio</th>`;
    h += `</tr></thead><tbody>`;
    // Body
    for (let i = 0; i < m; i++) {
      h += `<tr>`;
      h += `<td class="col-cb">${cb[i].toHTML()}</td>`;
      h += `<td class="col-basis">${basis[i]}</td>`;
      h += `<td class="col-xb">${xb[i].toHTML()}</td>`;
      for (let j = 0; j < nv; j++) {
        const cls = (i === pivotRow && j === pivotCol) ? "pivot-cell" : "";
        h += `<td class="${cls}">${table[i][j].toHTML()}</td>`;
      }
      if (ratios) {
        h += `<td class="col-ratio">${ratios[i] !== null ? ratios[i].toHTML() : "—"}</td>`;
      }
      h += `</tr>`;
    }
    // Zj-Cj row
    if (zjcj) {
      h += `<tr class="row-zjcj"><td></td><td class="col-basis" style="font-size:0.7rem">Z<sub>j</sub>−C<sub>j</sub></td><td></td>`;
      for (let j = 0; j < nv; j++) {
        const cls = zjcj[j].isNeg() ? "negative-zjcj" : (zjcj[j].gt(Frac.ZERO) ? "positive-zjcj" : "");
        h += `<td class="${cls}">${zjcj[j].toHTML()}</td>`;
      }
      if (ratios) h += `<td></td>`;
      h += `</tr>`;
    }
    h += `</tbody></table>`;
    return h;
  }

  // ── Append solution summary ──
  function appendSolution(res, maximize, nDecision, label = "Optimal Solution") {
    const { cb, basis, xb } = res;
    let z = Frac.ZERO;
    for (let i = 0; i < cb.length; i++) z = z.add(cb[i].mul(xb[i]));

    let html = `<h3>✓ ${label}</h3>`;
    html += `<div class="solution-z">Z = ${z.toHTML()}</div>`;
    for (let j = 0; j < nDecision; j++) {
      const name = `x${j + 1}`;
      const idx = basis.indexOf(name);
      const val = idx >= 0 ? xb[idx] : Frac.ZERO;
      html += `<div class="solution-line">${name} = ${val.toHTML()}</div>`;
    }

    const el = document.createElement("div");
    el.className = "solution-card";
    el.innerHTML = html;
    $output.appendChild(el);
  }
})();
