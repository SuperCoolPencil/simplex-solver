#!/usr/bin/env python3
"""
Simplex, Dual Simplex, and ILPP (Gomory's Cutting Plane) Solver
Uses exact fractions for all arithmetic.
Prints formatted tableaux at each iteration.
"""

from fractions import Fraction
import copy
import math
import sys


# ─────────────────────────────────────────────
#  Pretty Printing
# ─────────────────────────────────────────────

def fmt(val):
    """Format a Fraction for display."""
    if val.denominator == 1:
        return str(val.numerator)
    return f"{val.numerator}/{val.denominator}"


def print_tableau(cj, var_names, cb, basis, xb, table, zj_cj, ratios=None, title=""):
    """
    Print a full simplex tableau.

    cj        : list[Fraction]  – objective coefficients for all variables
    var_names : list[str]       – variable names (x1, x2, ..., S1, S2, ..., G1, ...)
    cb        : list[Fraction]  – objective coefficients of current basis
    basis     : list[str]       – current basic variable names
    xb        : list[Fraction]  – current RHS values
    table     : list[list[Fraction]] – constraint coefficient matrix
    zj_cj     : list[Fraction]  – Zj - Cj row
    ratios    : list[str|None]  – ratio column (optional)
    title     : str             – title for this tableau
    """
    if title:
        print(f"\n{'═' * 80}")
        print(f"  {title}")
        print(f"{'═' * 80}")

    n_vars = len(var_names)
    m = len(basis)

    # Determine column widths
    col_w = {}
    col_w['CB'] = max(4, max((len(fmt(c)) for c in cb), default=4))
    col_w['B'] = max(4, max((len(b) for b in basis), default=4))
    col_w['XB'] = max(4, max((len(fmt(v)) for v in xb), default=4))
    for j, name in enumerate(var_names):
        vals = [fmt(table[i][j]) for i in range(m)]
        vals.append(fmt(cj[j]))
        vals.append(name)
        if zj_cj:
            vals.append(fmt(zj_cj[j]))
        col_w[name] = max(len(v) for v in vals)
    if ratios:
        col_w['Ratio'] = max(5, max((len(r) for r in ratios), default=5))

    def pad(s, w):
        return str(s).rjust(w)

    # ── Cj row ──
    cj_label = "Cj"
    header_left = f"  {pad('', col_w['CB'])}  {pad('', col_w['B'])}  {pad('', col_w['XB'])}"
    cj_vals = "  ".join(pad(fmt(cj[j]), col_w[var_names[j]]) for j in range(n_vars))
    print(f"\n  {pad(cj_label, col_w['CB'])}  {pad('', col_w['B'])}  {pad('', col_w['XB'])}  {cj_vals}")

    # ── Header row ──
    hdr = f"  {pad('CB', col_w['CB'])}  {pad('B', col_w['B'])}  {pad('XB', col_w['XB'])}"
    for name in var_names:
        hdr += f"  {pad(name, col_w[name])}"
    if ratios:
        hdr += f"  {pad('Ratio', col_w['Ratio'])}"
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))

    # ── Body rows ──
    for i in range(m):
        row = f"  {pad(fmt(cb[i]), col_w['CB'])}  {pad(basis[i], col_w['B'])}  {pad(fmt(xb[i]), col_w['XB'])}"
        for j in range(n_vars):
            row += f"  {pad(fmt(table[i][j]), col_w[var_names[j]])}"
        if ratios:
            row += f"  {pad(ratios[i], col_w['Ratio'])}"
        print(row)

    # ── Zj - Cj row ──
    if zj_cj:
        sep = "  " + "─" * (len(hdr) - 2)
        print(sep)
        zrow = f"  {pad('', col_w['CB'])}  {pad('Zj-Cj', col_w['B'])}  {pad('', col_w['XB'])}"
        for j in range(n_vars):
            zrow += f"  {pad(fmt(zj_cj[j]), col_w[var_names[j]])}"
        print(zrow)

    print()


# ─────────────────────────────────────────────
#  Compute Zj - Cj
# ─────────────────────────────────────────────

def compute_zj_cj(cj, cb, table, n_vars, m):
    """Compute Zj - Cj for each variable column."""
    zj_cj = []
    for j in range(n_vars):
        zj = sum(cb[i] * table[i][j] for i in range(m))
        zj_cj.append(zj - cj[j])
    return zj_cj


# ─────────────────────────────────────────────
#  Simplex Method
# ─────────────────────────────────────────────

def simplex(cj, var_names, cb, basis, xb, table, maximize=True):
    """
    Perform the simplex method.

    Returns (cj, var_names, cb, basis, xb, table) of the optimal solution,
    or None if unbounded.
    """
    m = len(basis)
    n_vars = len(var_names)
    iteration = 0

    while True:
        zj_cj = compute_zj_cj(cj, cb, table, n_vars, m)

        # Determine entering variable
        if maximize:
            min_val = min(zj_cj)
            if min_val >= 0:
                # Optimal
                print_tableau(cj, var_names, cb, basis, xb, table, zj_cj,
                              title=f"Simplex Iteration {iteration} (Optimal)")
                return cj, var_names, cb, basis, xb, table
            entering = zj_cj.index(min_val)
        else:
            max_val = max(zj_cj)
            if max_val <= 0:
                print_tableau(cj, var_names, cb, basis, xb, table, zj_cj,
                              title=f"Simplex Iteration {iteration} (Optimal)")
                return cj, var_names, cb, basis, xb, table
            entering = zj_cj.index(max_val)

        # Compute ratios
        ratios_display = []
        ratios_val = []
        for i in range(m):
            if table[i][entering] > 0:
                r = xb[i] / table[i][entering]
                ratios_display.append(fmt(r))
                ratios_val.append(r)
            else:
                ratios_display.append("--")
                ratios_val.append(None)

        print_tableau(cj, var_names, cb, basis, xb, table, zj_cj,
                      ratios=ratios_display,
                      title=f"Simplex Iteration {iteration}")

        # Find leaving variable (minimum ratio)
        min_ratio = None
        leaving = -1
        for i in range(m):
            if ratios_val[i] is not None:
                if min_ratio is None or ratios_val[i] < min_ratio:
                    min_ratio = ratios_val[i]
                    leaving = i

        if leaving == -1:
            print("  ⚠ Problem is UNBOUNDED!")
            return None

        pivot = table[leaving][entering]
        print(f"  → Pivot: {basis[leaving]} leaves, {var_names[entering]} enters")
        print(f"  → Pivot element: {fmt(pivot)} (Row {leaving + 1}, Col {var_names[entering]})")

        # Pivot operation
        # 1) Divide pivot row by pivot element
        xb[leaving] = xb[leaving] / pivot
        for j in range(n_vars):
            table[leaving][j] = table[leaving][j] / pivot

        # 2) Eliminate other rows
        for i in range(m):
            if i == leaving:
                continue
            factor = table[i][entering]
            xb[i] = xb[i] - factor * xb[leaving]
            for j in range(n_vars):
                table[i][j] = table[i][j] - factor * table[leaving][j]

        # Update basis
        basis[leaving] = var_names[entering]
        cb[leaving] = cj[entering]

        iteration += 1


# ─────────────────────────────────────────────
#  Dual Simplex Method
# ─────────────────────────────────────────────

def dual_simplex(cj, var_names, cb, basis, xb, table):
    """
    Perform the dual simplex method.
    Used when solution is dual feasible (Zj-Cj >= 0) but primal infeasible (some XB < 0).

    Returns (cj, var_names, cb, basis, xb, table) of the optimal solution,
    or None if infeasible.
    """
    m = len(basis)
    n_vars = len(var_names)
    iteration = 0

    while True:
        zj_cj = compute_zj_cj(cj, cb, table, n_vars, m)

        # Check if all XB >= 0 (primal feasible)
        all_non_negative = all(xb[i] >= 0 for i in range(m))
        if all_non_negative:
            print_tableau(cj, var_names, cb, basis, xb, table, zj_cj,
                          title=f"Dual Simplex Iteration {iteration} (Optimal)")
            return cj, var_names, cb, basis, xb, table

        # Find leaving variable: row with most negative XB
        leaving = -1
        min_xb = Fraction(0)
        for i in range(m):
            if xb[i] < min_xb:
                min_xb = xb[i]
                leaving = i

        if leaving == -1:
            print("  ⚠ No negative XB found — already feasible.")
            return cj, var_names, cb, basis, xb, table

        print_tableau(cj, var_names, cb, basis, xb, table, zj_cj,
                      title=f"Dual Simplex Iteration {iteration}")
        print(f"  → Leaving variable: {basis[leaving]} (XB = {fmt(xb[leaving])})")

        # Find entering variable: max {(Zj-Cj) / a_kj} for a_kj < 0
        entering = -1
        max_ratio = None
        for j in range(n_vars):
            if table[leaving][j] < 0:
                ratio = zj_cj[j] / table[leaving][j]
                if max_ratio is None or ratio > max_ratio:
                    max_ratio = ratio
                    entering = j

        if entering == -1:
            print("  ⚠ Problem is INFEASIBLE (no negative element in leaving row)!")
            return None

        pivot = table[leaving][entering]
        print(f"  → Entering variable: {var_names[entering]}")
        print(f"  → Pivot element: {fmt(pivot)} (Row {leaving + 1}, Col {var_names[entering]})")

        # Pivot operation
        xb[leaving] = xb[leaving] / pivot
        for j in range(n_vars):
            table[leaving][j] = table[leaving][j] / pivot

        for i in range(m):
            if i == leaving:
                continue
            factor = table[i][entering]
            xb[i] = xb[i] - factor * xb[leaving]
            for j in range(n_vars):
                table[i][j] = table[i][j] - factor * table[leaving][j]

        basis[leaving] = var_names[entering]
        cb[leaving] = cj[entering]

        iteration += 1


# ─────────────────────────────────────────────
#  ILPP — Gomory's Cutting Plane Method
# ─────────────────────────────────────────────

def fractional_part(val):
    """Return the fractional part of a Fraction (always in [0, 1))."""
    return val - Fraction(math.floor(val))


def gomory_cut(cj, var_names, cb, basis, xb, table):
    """
    Apply Gomory's cutting plane method for ILPP.

    Steps:
    1. Solve continuous LPP using simplex.
    2. If all basic variables are integer, done.
    3. Otherwise, pick the basic variable with the largest fractional part.
    4. Generate the Gomory cut and add it to the tableau.
    5. Use dual simplex to re-optimize.
    6. Repeat until all basic variables are integer.
    """
    m = len(basis)
    n_vars = len(var_names)
    cut_count = 0

    while True:
        # Check if all basic variables are integers
        all_integer = True
        max_frac = Fraction(0)
        cut_row = -1

        for i in range(m):
            f = fractional_part(xb[i])
            if f != 0:
                all_integer = False
                if f > max_frac:
                    max_frac = f
                    cut_row = i

        if all_integer:
            print("\n" + "═" * 80)
            print("  ✓ All basic variables are integers. ILPP SOLVED!")
            print("═" * 80)
            return cj, var_names, cb, basis, xb, table

        print(f"\n{'─' * 80}")
        print(f"  Generating Gomory Cut #{cut_count + 1}")
        print(f"  Source row: {basis[cut_row]} = {fmt(xb[cut_row])}")
        print(f"  Fractional part (f_i): {fmt(max_frac)}")
        print(f"{'─' * 80}")

        # Build the Gomory cut
        # Cut: -f_ij * x_j + G = -f_i   for all non-basic variables
        # where f_ij = fractional part of a_ij
        cut_count += 1
        g_name = f"G{cut_count}"

        # New variable column for G (add 0 to cj)
        cj.append(Fraction(0))
        var_names.append(g_name)

        # Add a 0 column for G in existing rows
        for i in range(m):
            table[i].append(Fraction(0))

        # Build the new cut row
        new_row = []
        for j in range(len(var_names) - 1):  # exclude the new G column
            f_ij = fractional_part(table[cut_row][j])
            new_row.append(-f_ij)
        new_row.append(Fraction(1))  # coefficient of G in the cut row

        new_xb = -max_frac

        # Print the cut equation
        cut_terms = []
        for j in range(len(var_names) - 1):
            if new_row[j] != 0:
                cut_terms.append(f"{fmt(new_row[j])} {var_names[j]}")
        cut_eq = " + ".join(cut_terms) + f" + {g_name} = {fmt(new_xb)}"
        print(f"  Cut: {cut_eq}")

        # Add the cut row to the tableau
        table.append(new_row)
        xb.append(new_xb)
        basis.append(g_name)
        cb.append(Fraction(0))
        m += 1
        n_vars = len(var_names)

        # Now apply dual simplex to restore feasibility
        print(f"\n  Applying Dual Simplex to restore feasibility...")
        result = dual_simplex(cj, var_names, cb, basis, xb, table)
        if result is None:
            print("  ⚠ ILPP is infeasible!")
            return None
        cj, var_names, cb, basis, xb, table = result
        m = len(basis)
        n_vars = len(var_names)


# ─────────────────────────────────────────────
#  Problem Input Helpers
# ─────────────────────────────────────────────

def parse_fraction(s):
    """Parse a string like '3', '-1/2', '0' into a Fraction."""
    s = s.strip()
    if '/' in s:
        parts = s.split('/')
        return Fraction(int(parts[0]), int(parts[1]))
    return Fraction(int(s))


def input_problem():
    """Interactively read an LP problem from the user."""
    print("\n╔══════════════════════════════════════════════════╗")
    print("║   Simplex / Dual Simplex / ILPP Solver           ║")
    print("║   (All computations use exact fractions)          ║")
    print("╚══════════════════════════════════════════════════╝\n")

    print("Select method:")
    print("  1. Simplex Method")
    print("  2. Dual Simplex Method")
    print("  3. ILPP (Gomory's Cutting Plane)")
    print("  4. Run example from notes")
    choice = input("\nChoice [1/2/3/4]: ").strip()

    if choice == '4':
        run_example()
        return

    n = int(input("Number of decision variables: "))
    m = int(input("Number of constraints: "))

    print(f"\nEnter objective function coefficients (c1 c2 ... c{n}):")
    obj = list(map(parse_fraction, input("  > ").split()))

    print(f"\nEnter constraint type for each constraint:")
    print("  <= : 1")
    print("  >= : 2")
    print("  =  : 3")

    constraints = []
    constraint_types = []
    rhs = []

    for i in range(m):
        print(f"\nConstraint {i + 1}:")
        print(f"  Coefficients (a1 a2 ... a{n}):")
        coeffs = list(map(parse_fraction, input("    > ").split()))
        ct = input("  Type (1: <=, 2: >=, 3: =): ").strip()
        constraint_types.append(int(ct))
        b = parse_fraction(input("  RHS value: "))
        constraints.append(coeffs)
        rhs.append(b)

    opt = input("\nMaximize or Minimize? (max/min): ").strip().lower()
    maximize = (opt == 'max')

    # Build the initial tableau
    # Add slack/surplus/artificial variables as needed
    var_names = [f"x{j + 1}" for j in range(n)]
    cj = list(obj)

    slack_count = 0
    for i in range(m):
        if constraint_types[i] == 1:  # <=
            slack_count += 1
            s_name = f"S{slack_count}"
            var_names.append(s_name)
            cj.append(Fraction(0))
            for k in range(m):
                constraints[k].append(Fraction(1) if k == i else Fraction(0))
        elif constraint_types[i] == 2:  # >=
            slack_count += 1
            s_name = f"S{slack_count}"
            var_names.append(s_name)
            cj.append(Fraction(0))
            for k in range(m):
                constraints[k].append(Fraction(-1) if k == i else Fraction(0))
        # For = constraints, no slack variable

    n_vars = len(var_names)

    # Determine initial basis (slack variables)
    basis = []
    cb = []
    xb = list(rhs)

    slack_idx = n
    for i in range(m):
        if constraint_types[i] == 1:
            basis.append(var_names[slack_idx])
            cb.append(Fraction(0))
            slack_idx += 1
        elif constraint_types[i] == 2:
            # For >= constraints with dual simplex, multiply by -1 to make <= first
            # Or use Big-M / Two-phase. For simplicity, handle via dual simplex.
            basis.append(var_names[slack_idx])
            cb.append(Fraction(0))
            slack_idx += 1
        else:
            # Need artificial variable for = constraints
            a_name = f"A{i + 1}"
            var_names.append(a_name)
            big_m = Fraction(-1000) if maximize else Fraction(1000)
            cj.append(big_m)
            for k in range(m):
                constraints[k].append(Fraction(1) if k == i else Fraction(0))
            basis.append(a_name)
            cb.append(big_m)

    n_vars = len(var_names)
    table = [row[:] for row in constraints]

    if choice == '1':
        print("\n" + "=" * 80)
        print("  SIMPLEX METHOD")
        print("=" * 80)
        result = simplex(cj, var_names, cb, basis, xb, table, maximize)
        if result:
            print_solution(result, maximize)

    elif choice == '2':
        print("\n" + "=" * 80)
        print("  DUAL SIMPLEX METHOD")
        print("=" * 80)
        # For dual simplex, we need dual feasibility first
        # Typically used when RHS has negatives after reformulation
        result = dual_simplex(cj, var_names, cb, basis, xb, table)
        if result:
            print_solution(result, maximize)

    elif choice == '3':
        print("\n" + "=" * 80)
        print("  ILPP — GOMORY'S CUTTING PLANE METHOD")
        print("=" * 80)

        # First solve continuous relaxation
        print("\n  Step 1: Solving continuous relaxation using Simplex...")
        result = simplex(cj, var_names, cb, basis, xb, table, maximize)
        if result is None:
            print("  ⚠ Problem is unbounded!")
            return

        cj, var_names, cb, basis, xb, table = result
        print_solution(result, maximize)

        # Then apply Gomory cuts
        print("\n  Step 2: Applying Gomory Cuts for integer solution...")
        result = gomory_cut(cj, var_names, cb, basis, xb, table)
        if result:
            print_solution(result, maximize, integer=True)


def print_solution(result, maximize=True, integer=False):
    """Print the final solution."""
    cj, var_names, cb, basis, xb, table = result
    prefix = "Integer " if integer else ""
    print(f"\n┌──────────────────────────────────────┐")
    print(f"│  {prefix}Optimal Solution                  │")
    print(f"└──────────────────────────────────────┘")

    z = sum(cb[i] * xb[i] for i in range(len(basis)))
    print(f"  Z = {fmt(z)}")

    # Show values of decision variables
    for name in var_names:
        if name.startswith('x'):
            if name in basis:
                idx = basis.index(name)
                print(f"  {name} = {fmt(xb[idx])}")
            else:
                print(f"  {name} = 0")


# ─────────────────────────────────────────────
#  Example from PDF Notes
# ─────────────────────────────────────────────

def run_example():
    """
    Run the exact example from the notes:

    Max Z = 3x1 + x2 + 3x3
    Subject to:
        -x1 + 2x2 +  x3 <= 4
              4x2 - 3x3 <= 2
         x1 - 3x2 + 2x3 <= 3
        x1, x2, x3 >= 0 and integer
    """
    print("\n" + "=" * 80)
    print("  EXAMPLE FROM NOTES")
    print("  Max Z = 3x1 + x2 + 3x3")
    print("  Subject to:")
    print("    -x1 + 2x2 +  x3 <= 4")
    print("          4x2 - 3x3 <= 2")
    print("     x1 - 3x2 + 2x3 <= 3")
    print("    x1, x2, x3 >= 0, integer")
    print("=" * 80)

    F = Fraction

    # Objective: Max Z = 3x1 + x2 + 3x3 + 0*S1 + 0*S2 + 0*S3
    cj = [F(3), F(1), F(3), F(0), F(0), F(0)]
    var_names = ['x1', 'x2', 'x3', 'S1', 'S2', 'S3']

    # Initial basis: S1, S2, S3
    cb = [F(0), F(0), F(0)]
    basis = ['S1', 'S2', 'S3']
    xb = [F(4), F(2), F(3)]

    # Constraint matrix (with slack columns)
    table = [
        [F(-1), F(2),  F(1),  F(1), F(0), F(0)],   # -x1 + 2x2 + x3 + S1 = 4
        [F(0),  F(4),  F(-3), F(0), F(1), F(0)],   #       4x2 - 3x3 + S2 = 2
        [F(1),  F(-3), F(2),  F(0), F(0), F(1)],   #  x1 - 3x2 + 2x3 + S3 = 3
    ]

    # Step 1: Solve continuous relaxation using Simplex
    print("\n  ── Step 1: Solving continuous LPP using Simplex ──")
    result = simplex(cj, var_names, cb, basis, xb, table, maximize=True)

    if result is None:
        print("  ⚠ Unbounded!")
        return

    cj, var_names, cb, basis, xb, table = result
    print_solution(result, maximize=True)

    # Step 2: Apply Gomory cuts for integer solution
    print("\n  ── Step 2: Applying Gomory Cuts (ILPP) ──")
    result = gomory_cut(cj, var_names, cb, basis, xb, table)

    if result:
        print_solution(result, maximize=True, integer=True)


# ─────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────

if __name__ == '__main__':
    input_problem()
