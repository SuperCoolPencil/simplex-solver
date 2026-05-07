/**
 * Simplex, Dual Simplex, and ILPP (Gomory) solver.
 * All arithmetic uses the Frac class for exact results.
 * Each function collects iteration snapshots for rendering.
 */

function computeZjCj(cj, cb, table, nVars, m) {
  const zjcj = [];
  for (let j = 0; j < nVars; j++) {
    let zj = Frac.ZERO;
    for (let i = 0; i < m; i++) zj = zj.add(cb[i].mul(table[i][j]));
    zjcj.push(zj.sub(cj[j]));
  }
  return zjcj;
}

/**
 * Run the simplex method.
 * Returns { steps: [...], optimal: bool, unbounded: bool, cj, varNames, cb, basis, xb, table }
 */
function solveSimplex(cj, varNames, cb, basis, xb, table, maximize) {
  const m = basis.length;
  const nVars = varNames.length;
  const steps = [];
  let iter = 0;

  while (true) {
    const zjcj = computeZjCj(cj, cb, table, nVars, m);

    // Check optimality
    let entering = -1;
    if (maximize) {
      let minVal = Frac.ZERO;
      for (let j = 0; j < nVars; j++) {
        if (zjcj[j].lt(minVal)) { minVal = zjcj[j]; entering = j; }
      }
      if (entering === -1) {
        steps.push({ type: "tableau", title: `Iteration ${iter}`, tag: "optimal", cj: [...cj], varNames: [...varNames], cb: [...cb], basis: [...basis], xb: [...xb], table: table.map(r => [...r]), zjcj: [...zjcj], ratios: null, pivotRow: -1, pivotCol: -1 });
        return { steps, optimal: true, unbounded: false, cj, varNames, cb, basis, xb, table };
      }
    } else {
      let maxVal = Frac.ZERO;
      for (let j = 0; j < nVars; j++) {
        if (zjcj[j].gt(maxVal)) { maxVal = zjcj[j]; entering = j; }
      }
      if (entering === -1) {
        steps.push({ type: "tableau", title: `Iteration ${iter}`, tag: "optimal", cj: [...cj], varNames: [...varNames], cb: [...cb], basis: [...basis], xb: [...xb], table: table.map(r => [...r]), zjcj: [...zjcj], ratios: null, pivotRow: -1, pivotCol: -1 });
        return { steps, optimal: true, unbounded: false, cj, varNames, cb, basis, xb, table };
      }
    }

    // Ratios
    const ratios = [];
    for (let i = 0; i < m; i++) {
      if (table[i][entering].gt(Frac.ZERO)) {
        ratios.push(xb[i].div(table[i][entering]));
      } else {
        ratios.push(null);
      }
    }

    let leaving = -1, minRatio = null;
    for (let i = 0; i < m; i++) {
      if (ratios[i] !== null && (minRatio === null || ratios[i].lt(minRatio))) {
        minRatio = ratios[i]; leaving = i;
      }
    }

    steps.push({ type: "tableau", title: `Iteration ${iter}`, tag: "iter", cj: [...cj], varNames: [...varNames], cb: [...cb], basis: [...basis], xb: [...xb], table: table.map(r => [...r]), zjcj: [...zjcj], ratios: ratios, pivotRow: leaving, pivotCol: entering });

    if (leaving === -1) {
      steps.push({ type: "info", msg: "Problem is UNBOUNDED — no valid ratio found." });
      return { steps, optimal: false, unbounded: true, cj, varNames, cb, basis, xb, table };
    }

    steps.push({ type: "pivot", leaving: basis[leaving], entering: varNames[entering], element: table[leaving][entering].toString() });

    // Pivot
    const pivot = table[leaving][entering];
    xb[leaving] = xb[leaving].div(pivot);
    for (let j = 0; j < nVars; j++) table[leaving][j] = table[leaving][j].div(pivot);
    for (let i = 0; i < m; i++) {
      if (i === leaving) continue;
      const f = table[i][entering];
      xb[i] = xb[i].sub(f.mul(xb[leaving]));
      for (let j = 0; j < nVars; j++) table[i][j] = table[i][j].sub(f.mul(table[leaving][j]));
    }
    basis[leaving] = varNames[entering];
    cb[leaving] = cj[entering];
    iter++;
  }
}

/**
 * Dual simplex method.
 */
function solveDualSimplex(cj, varNames, cb, basis, xb, table) {
  const nVars = varNames.length;
  let m = basis.length;
  const steps = [];
  let iter = 0;

  while (true) {
    const zjcj = computeZjCj(cj, cb, table, nVars, m);

    if (xb.every(v => v.gte(Frac.ZERO))) {
      steps.push({ type: "tableau", title: `Dual Simplex Iter ${iter}`, tag: "optimal", cj: [...cj], varNames: [...varNames], cb: [...cb], basis: [...basis], xb: [...xb], table: table.map(r => [...r]), zjcj: [...zjcj], ratios: null, pivotRow: -1, pivotCol: -1 });
      return { steps, optimal: true, infeasible: false, cj, varNames, cb, basis, xb, table };
    }

    // Leaving: most negative XB
    let leaving = -1, minXb = Frac.ZERO;
    for (let i = 0; i < m; i++) {
      if (xb[i].lt(minXb)) { minXb = xb[i]; leaving = i; }
    }

    steps.push({ type: "tableau", title: `Dual Simplex Iter ${iter}`, tag: "iter", cj: [...cj], varNames: [...varNames], cb: [...cb], basis: [...basis], xb: [...xb], table: table.map(r => [...r]), zjcj: [...zjcj], ratios: null, pivotRow: leaving, pivotCol: -1 });

    // Entering: max {zjcj[j] / a[leaving][j]} for a[leaving][j] < 0
    let entering = -1, maxRatio = null;
    for (let j = 0; j < nVars; j++) {
      if (table[leaving][j].lt(Frac.ZERO)) {
        const r = zjcj[j].div(table[leaving][j]);
        if (maxRatio === null || r.gt(maxRatio)) { maxRatio = r; entering = j; }
      }
    }

    if (entering === -1) {
      steps.push({ type: "info", msg: "Problem is INFEASIBLE — no negative element in leaving row." });
      return { steps, optimal: false, infeasible: true, cj, varNames, cb, basis, xb, table };
    }

    steps.push({ type: "pivot", leaving: basis[leaving], entering: varNames[entering], element: table[leaving][entering].toString() });

    const pivot = table[leaving][entering];
    xb[leaving] = xb[leaving].div(pivot);
    for (let j = 0; j < nVars; j++) table[leaving][j] = table[leaving][j].div(pivot);
    for (let i = 0; i < m; i++) {
      if (i === leaving) continue;
      const f = table[i][entering];
      xb[i] = xb[i].sub(f.mul(xb[leaving]));
      for (let j = 0; j < nVars; j++) table[i][j] = table[i][j].sub(f.mul(table[leaving][j]));
    }
    basis[leaving] = varNames[entering];
    cb[leaving] = cj[entering];
    iter++;
  }
}

/**
 * ILPP — Gomory's Cutting Plane.
 * Takes an already-optimal continuous solution and applies cuts.
 */
function solveILPP(cj, varNames, cb, basis, xb, table) {
  let m = basis.length;
  let nVars = varNames.length;
  const steps = [];
  let cutNum = 0;

  while (true) {
    // Check integrality
    let cutRow = -1, maxFrac = Frac.ZERO;
    for (let i = 0; i < m; i++) {
      const f = xb[i].fracPart();
      if (!f.isZero() && f.gt(maxFrac)) { maxFrac = f; cutRow = i; }
    }
    if (cutRow === -1) {
      steps.push({ type: "info", msg: "✓ All basic variables are integers — ILPP solved!", success: true });
      return { steps, cj, varNames, cb, basis, xb, table };
    }

    cutNum++;
    const gName = `G${cutNum}`;

    // Build cut info
    const cutTerms = [];
    for (let j = 0; j < nVars; j++) {
      const fij = table[cutRow][j].fracPart();
      if (!fij.isZero()) cutTerms.push(`${fij.neg()} ${varNames[j]}`);
    }
    steps.push({ type: "cut", num: cutNum, source: basis[cutRow], value: xb[cutRow].toString(), fi: maxFrac.toString(), equation: cutTerms.join(" + ") + ` + ${gName} = ${maxFrac.neg()}` });

    // Add G column
    cj.push(Frac.ZERO);
    varNames.push(gName);
    for (let i = 0; i < m; i++) table[i].push(Frac.ZERO);

    // Build cut row
    const newRow = [];
    for (let j = 0; j < varNames.length - 1; j++) {
      newRow.push(table[cutRow][j].fracPart().neg());
    }
    newRow.push(Frac.ONE); // G coefficient

    table.push(newRow);
    xb.push(maxFrac.neg());
    basis.push(gName);
    cb.push(Frac.ZERO);
    m++;
    nVars = varNames.length;

    // Dual simplex
    const dsResult = solveDualSimplex(cj, varNames, cb, basis, xb, table);
    steps.push(...dsResult.steps);

    if (!dsResult.optimal) {
      steps.push({ type: "info", msg: "ILPP is infeasible after adding Gomory cut." });
      return { steps, cj, varNames, cb, basis, xb, table };
    }
    // Update references
    cj = dsResult.cj; varNames = dsResult.varNames;
    cb = dsResult.cb; basis = dsResult.basis;
    xb = dsResult.xb; table = dsResult.table;
    m = basis.length; nVars = varNames.length;
  }
}
