#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
latex_solver.py

Parse a differential equation given in LaTeX, classify it, and solve it.
Supported (automatic) cases:
 - ODE: try symbolic dsolve; otherwise numeric RK4 (system).
 - 1D Heat equation (u_t = alpha u_xx): Crank-Nicolson numeric.
 - 1D Wave equation (u_tt = c^2 u_xx): Leapfrog numeric.
 - Radial TDSE (i*hbar u_t = - (hbar^2/2m) u_xx + V u for u=r*psi): Crank-Nicolson complex.

Usage:
  python latex_solver.py --latex "<latex string>"
Examples:
  python latex_solver.py --latex "\\frac{d^2 y}{dt^2} + y = 0"
  python latex_solver.py --latex "u_t = u_{xx}"
  python latex_solver.py --latex "i*hbar \\frac{\\partial \\psi}{\\partial t} = -\\frac{\\hbar^2}{2m} \\frac{\\partial^2 \\psi}{\\partial r^2} + V(r) \\psi"
"""

import argparse
import sys
import math
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sympy import symbols, Function, Eq, Derivative, Symbol, simplify, pprint
from sympy.parsing.latex import parse_latex
from sympy import dsolve, solve, lambdify, S

# =========================================================
# Utilidad: Algoritmo de Thomas (sistemas tridiagonales)
# =========================================================
def thomas(a, b, c, d):
    """
    Resuelve un sistema tridiagonal Ax=d con:
      a: subdiagonal (len n-1)
      b: diagonal principal (len n)
      c: superdiagonal (len n-1)
      d: vector RHS (len n)
    Devuelve x (len n).
    """
    a = np.array(a, dtype=complex if np.iscomplexobj(b) or np.iscomplexobj(d) else float)
    b = np.array(b, dtype=a.dtype)
    c = np.array(c, dtype=a.dtype)
    d = np.array(d, dtype=a.dtype)

    n = len(b)
    if len(a) != n-1 or len(c) != n-1 or len(d) != n:
        raise ValueError("Dimensiones incompatibles para sistema tridiagonal.")

    # forward elimination
    cp = np.zeros(n-1, dtype=a.dtype)
    dp = np.zeros(n, dtype=a.dtype)
    bp = b.astype(a.dtype).copy()

    cp[0] = c[0] / bp[0]
    dp[0] = d[0] / bp[0]
    for i in range(1, n-1):
        denom = bp[i] - a[i-1]*cp[i-1]
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i-1]*dp[i-1]) / denom
    dp[n-1] = (d[n-1] - a[n-2]*dp[n-2]) / (bp[n-1] - a[n-2]*cp[n-2])

    # back substitution
    x = np.zeros(n, dtype=a.dtype)
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i]*x[i+1]
    return x

# ---------- Utilities ----------
def safe_parse_latex(latex_str):
    try:
        expr = parse_latex(latex_str)
        return expr
    except Exception as e:
        print("LaTeX parsing failed:", e)
        print("Try simpler LaTeX (use \\partial, \\frac{d}{dt}, u(x,t) style).")
        raise

def expr_to_eq(expr):
    # If parsed is an Eq/Relational, return; else treat as expr == 0
    if hasattr(expr, 'is_Relational') and expr.is_Relational:
        return expr
    return Eq(expr, 0)

def find_derivative_atoms(expr):
    return list(expr.atoms(Derivative))

def unique_independent_vars(derivs):
    vars_set = set()
    for d in derivs:
        vars_set.update(d.variables)
    return vars_set

# ---------- Classification ----------
def classify(parsed, latex_str: str | None = None):
    """
    Devuelve metadatos de clasificación de la ecuación.
    Usa objetos SymPy y, si hace falta, heurísticas sobre el LaTeX original.
    """
    # Si es igualdad, pasamos a expr = lhs - rhs
    if getattr(parsed, 'is_Equality', False) or getattr(parsed, 'is_Relational', False):
        lhs = parsed.lhs
        rhs = parsed.rhs
        expr = lhs - rhs
    else:
        expr = parsed

    derivs = find_derivative_atoms(expr)
    vars_in_derivs = unique_independent_vars(derivs)
    atoms = list(expr.atoms())

    cls = {
        'expr': expr,
        'derivs': derivs,
        'vars_in_derivs': vars_in_derivs,
        'is_pde': False,
        'is_ode': False,
        'time_vars': set(),
        'space_vars': set(),
        'order_time': 0,
        'order_space': 0,
        'guess': 'unknown',
        'notes': []
    }

    # Heurística de variables
    for v in vars_in_derivs:
        name = str(v)
        if name in ('t', 'time'):
            cls['time_vars'].add(v)
        elif name in ('x', 'r', 'radius', 'rho', 'theta', 'phi'):
            cls['space_vars'].add(v)
        else:
            # si solo hay una variable en derivadas, podría ser ODE temporal
            if len(vars_in_derivs) == 1:
                cls['time_vars'].add(v)

    # Órdenes detectados desde objetos Derivative
    for d in derivs:
        tv = [v for v in d.variables if str(v) in ('t', 'time')]
        sv = [v for v in d.variables if str(v) in ('x','r','theta','phi','rho')]
        cls['order_time'] = max(cls['order_time'], len(tv))
        cls['order_space'] = max(cls['order_space'], len(sv))

    # Si no hay Derivative en el árbol, usar heurísticas de texto LaTeX
    if not derivs and latex_str:
        s = latex_str.lower()
        if re.search(r"(\\ddot|y''|\w+''\(|\\frac\s*\{\s*d\^?2\s*\}\s*\{\s*dt\^?2\s*\}|d\^?2\s*/\s*dt\^?2)", s):
            cls['order_time'] = max(cls['order_time'], 2)
            cls['time_vars'].add(Symbol('t'))
            cls['is_ode'] = True
            cls['notes'].append("Detected second-time derivative by LaTeX heuristics.")
        elif re.search(r"(\\dot|y'|\w+'\(|\\frac\s*\{\s*d\s*\}\s*\{\s*dt\s*\}|d\s*/\s*dt)", s):
            cls['order_time'] = max(cls['order_time'], 1)
            cls['time_vars'].add(Symbol('t'))
            cls['is_ode'] = True
            cls['notes'].append("Detected first-time derivative by LaTeX heuristics.")
        if re.search(r"(u_{xx}|\\partial\^?2\s*u/\s*\\partial x\^?2|d\^?2\s*/\s*dx\^?2|u_{rr}|u_{\theta\theta}|u_{\phi\phi})", s):
            cls['order_space'] = max(cls['order_space'], 2)
            cls['space_vars'].add(Symbol('x'))
            cls['is_pde'] = True
            cls['notes'].append("Detected second-space derivative by LaTeX heuristics.")

    # Decidir ODE/PDE si no está decidido
    atom_names = {str(a) for a in atoms}
    if any(name in atom_names for name in ('x','r','theta','phi','rho')):
        cls['space_vars'].update([Symbol(n) for n in ['x','r','theta','phi','rho'] if n in atom_names])
    if cls['space_vars']:
        cls['is_pde'] = True
    elif cls['time_vars']:
        cls['is_ode'] = True

    # Adivinar tipo de PDE
    srepr = (latex_str or str(parsed)).lower()
    if any(k in srepr for k in ['hbar', 'psi', 'i*', ' i ']):
        cls['guess'] = 'schrodinger'
    elif (('u_tt' in srepr) or ('ddot' in srepr) or re.search(r'd\^?2\s*/\s*dt\^?2', srepr)):
        cls['guess'] = 'wave'
    elif (('u_t' in srepr) or ('u_{xx}' in srepr) or re.search(r'd\^?2\s*/\s*dx\^?2', srepr) or ('\\partial' in srepr)):
        cls['guess'] = 'heat_or_diffusion'
    else:
        cls['guess'] = 'unknown'

    return cls

# ---------- ODE solvers ----------
def try_symbolic_ode(expr):
    try:
        print("Trying SymPy dsolve for symbolic ODE solution...")
        sol = dsolve(Eq(expr, 0))
        print("Symbolic dsolve result:")
        pprint(sol)
        return sol
    except Exception as e:
        print("Symbolic solve failed or not supported:", e)
        return None

def numeric_ode_integrate(expr, t0=0.0, tmax=10.0, dt=0.01, y0=None):
    """
    Integra numéricamente una ODE de orden n transformándola en sistema de 1er orden (RK4).
    - Detecta función dependiente y variable independiente a partir de 'expr = 0'.
    - Aísla la derivada de orden máximo y construye f(t, y, y', ..., y^(n-1)).
    - y0: condiciones iniciales [y(0), y'(0), ..., y^(n-1)(0)].
    Devuelve: times, Y (matriz con columnas = [y, y', ...]).
    """
    # Buscar derivadas:
    derivs = list(expr.atoms(Derivative))
    if not derivs:
        raise RuntimeError("No derivatives found for numeric ODE integration.")

    # Variable independiente: tomamos la primera que aparezca en derivadas
    indep_vars = []
    for d in derivs:
        for v in d.variables:
            if v not in indep_vars:
                indep_vars.append(v)
    if not indep_vars:
        raise RuntimeError("Could not detect independent variable.")
    t = indep_vars[0]

    # Función dependiente: el objeto Function que aparece en derivadas
    funcs = [d.expr for d in derivs]  # Derivative(y(t), ...) -> expr = y(t)
    # Asegurar que es una Function
    funcs = [f for f in funcs if hasattr(f, 'free_symbols') or hasattr(f, 'func')]
    if not funcs:
        raise RuntimeError("Could not detect dependent function.")
    y_func = funcs[0]       # y(t)
    y_class = y_func.func   # y

    # Orden de la ODE = máximo número de veces que aparece 't' en variables de la derivada de y
    order = 1
    for d in derivs:
        if d.expr == y_func:
            order = max(order, len(d.variables))
    if order < 1:
        order = 1

    # Identificar la derivada de orden máximo D^{order} y aislarla: expr==0 => D^{order} y = F(...)
    highest = None
    target = Derivative(y_func, *([t]*order))
    # Intentar resolver directamente para la derivada de orden máximo
    try:
        sol_for_high = solve(Eq(expr, 0), target, dict=True)
    except Exception:
        sol_for_high = []
    if sol_for_high:
        rhs_high = sol_for_high[0][target]
    else:
        # Si falla, intenta resolver respecto a "alguna" derivada de orden 'order'
        candidates = [d for d in derivs if d.expr == y_func and len(d.variables) == order]
        if not candidates:
            raise RuntimeError("Couldn't find highest derivative to isolate.")
        try:
            sol_for_high = solve(Eq(expr, 0), candidates[0], dict=True)
            rhs_high = sol_for_high[0][candidates[0]]
        except Exception as e:
            raise RuntimeError(f"Cannot isolate highest derivative: {e}")

    # Construir lambda f(t, y0, y1, ..., y_{n-1})
    # Creamos símbolos Y0..Y_{n-1} y sustituimos:
    Y = symbols(' '.join([f'Y{i}' for i in range(order)]))
    subs_map = {y_func: Y[0]}
    for i in range(1, order):
        subs_map[Derivative(y_func, *([t]*i))] = Y[i]
    rhs_sym = rhs_high.subs(subs_map)

    f = lambdify((t,) + tuple(Y), rhs_sym, 'numpy')

    def system_rhs(tval, state):
        # state = [y, y', ..., y^{n-1}]
        out = np.zeros_like(state, dtype=float)
        # y' = y1, y1' = y2, ..., y^{n-2}' = y^{n-1}
        out[:-1] = state[1:]
        # y^{n-1}' = f(t, Y0..Y_{n-1})
        out[-1] = float(f(tval, *[float(s) for s in state]))
        return out

    # Integración RK4
    times = np.arange(t0, tmax + 1e-12, dt)
    if y0 is None:
        y0 = np.zeros(order, dtype=float)
        y0[0] = 1.0  # valor por defecto
    sol = np.zeros((len(times), order), dtype=float)
    sol[0, :] = y0

    for k in range(1, len(times)):
        tt = times[k-1]; h = dt; s = sol[k-1, :]
        k1 = system_rhs(tt, s)
        k2 = system_rhs(tt + h/2, s + 0.5*h*k1)
        k3 = system_rhs(tt + h/2, s + 0.5*h*k2)
        k4 = system_rhs(tt + h, s + h*k3)
        sol[k, :] = s + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

    return times, sol

# ---------- PDE solvers (1D) ----------
def solve_heat_1d_cn(u0, L=1.0, Nx=201, dt=1e-4, Nt=400, alpha=1.0):
    dx = L/(Nx-1)
    r = alpha*dt/(2*dx*dx)
    N = Nx
    # interior size N-2
    a = -r*np.ones(N-3); b = (1+2*r)*np.ones(N-2); c = -r*np.ones(N-3)
    aB =  r*np.ones(N-3); bB = (1-2*r)*np.ones(N-2); cB =  r*np.ones(N-3)
    u_prev = u0.copy()
    results = [u_prev.copy()]
    for n in range(Nt):
        u_in = u_prev[1:-1]
        d = bB * u_in.copy()
        d[:-1] += aB * u_in[1:]
        d[1:]  += cB * u_in[:-1]
        u_in_new = thomas(a, b, c, d)
        u_new = np.zeros_like(u_prev)
        u_new[1:-1] = np.real(u_in_new)
        # BC Dirichlet 0
        u_new[0] = 0.0; u_new[-1] = 0.0
        results.append(u_new)
        u_prev = u_new
    return np.array(results)

def solve_wave_1d_leapfrog(u0, v0, L=1.0, Nx=201, dt=1e-4, Nt=400, c=1.0):
    dx = L/(Nx-1)
    lam = c*dt/dx
    if lam > 1.0:
        print("Warning: CFL > 1; scheme unstable.")
    u_prev = u0.copy()
    u_curr = u0.copy()
    uxx = np.zeros_like(u0)
    for i in range(1, Nx-1):
        uxx[i] = (u0[i+1] - 2*u0[i] + u0[i-1]) / dx**2
    u_curr[1:-1] = u0[1:-1] + dt*v0[1:-1] + 0.5*(dt**2)*c**2*uxx[1:-1]
    results = [u_prev.copy(), u_curr.copy()]
    for n in range(1, Nt):
        u_next = np.zeros_like(u0)
        for i in range(1, Nx-1):
            u_next[i] = 2*u_curr[i] - u_prev[i] + (lam**2)*(u_curr[i+1] - 2*u_curr[i] + u_curr[i-1])
        # BC Dirichlet 0
        u_next[0] = 0.0; u_next[-1] = 0.0
        results.append(u_next.copy())
        u_prev, u_curr = u_curr, u_next
    return np.array(results)

def solve_tdse_radial_cn(u0, rmax=20.0, Nr=400, dt=0.005, Nt=400, V_func=None, hbar=1.0, m=1.0):
    dr = rmax/Nr
    r = np.linspace(0, rmax, Nr+1)
    V = np.zeros_like(r)
    if V_func is None:
        V[1:] = -1.0/r[1:]
        V[0] = V[1]
    else:
        for j in range(1, Nr+1):
            V[j] = V_func(r[j])

    u = np.array(u0, dtype=complex)
    results = [u.copy()]
    im = 1j
    coef = im * dt / (2.0 * hbar)
    k = hbar*hbar/(2.0*m*dr*dr)
    N_inner = Nr-1

    a = np.full(N_inner-1, -coef*k, dtype=complex)
    c = np.full(N_inner-1, -coef*k, dtype=complex)
    b = np.zeros(N_inner, dtype=complex)
    for j in range(1, Nr):
        b[j-1] = 1.0 + 2.0*coef*k + coef*V[j]

    aB = np.full(N_inner-1,  coef*k, dtype=complex)
    cB = np.full(N_inner-1,  coef*k, dtype=complex)
    bB = np.zeros(N_inner, dtype=complex)
    for j in range(1, Nr):
        bB[j-1] = 1.0 - 2.0*coef*k - coef*V[j]

    for n in range(Nt):
        u_in = u[1:-1]
        d = bB * u_in
        d[:-1] += aB * u_in[1:]
        d[1:]  += cB * u_in[:-1]
        u_inner_next = thomas(a, b, c, d)
        u[1:-1] = u_inner_next
        u[0] = 0.0; u[-1] = 0.0
        results.append(u.copy())
    return r, np.array(results)

# ---------- Orchestrator ----------
def auto_solve_from_latex(latex_str):
    print("Parsing LaTeX:", latex_str)
    parsed = safe_parse_latex(latex_str)
    print("Parsed SymPy expression:")
    pprint(parsed)
    classification = classify(parsed, latex_str)
    print("Classification summary:", classification['guess'],
          "is_pde?", classification['is_pde'],
          "order_time", classification['order_time'],
          "order_space", classification['order_space'])

    # ODE
    if classification['is_ode'] and not classification['is_pde']:
        print("--> Detected ODE. Trying symbolic solve first.")
        sym_sol = try_symbolic_ode(classification['expr'])
        if sym_sol is not None:
            print("Symbolic solution obtained. Writing to file symbolic_ode_solution.txt")
            with open("symbolic_ode_solution.txt", "w") as f:
                f.write(str(sym_sol))
            return
        else:
            print("Falling back to numeric RK4 integration (defaults).")
            try:
                t, sol = numeric_ode_integrate(classification['expr'])
                df = pd.DataFrame(sol)
                df.insert(0, "t", t)
                df.to_csv("ode_numeric_solution.csv", index=False)
                print("Wrote ode_numeric_solution.csv")
                # Plot rápido
                plt.figure()
                plt.plot(t, sol[:, 0])
                plt.xlabel('t'); plt.ylabel('y (first component)')
                plt.title('Numeric ODE RK4 solution')
                plt.show()
                return
            except Exception as e:
                print("Numeric ODE solve failed:", e)
                return

    # PDE
    if classification['is_pde']:
        guess = classification.get('guess', 'unknown')
        if guess.startswith('heat'):
            print("--> Heat equation detected (or similar). Running Crank-Nicolson 1D with defaults.")
            L = 1.0; Nx = 201; dt = 5e-5; Nt = 400; alpha = 1.0
            x = np.linspace(0, L, Nx)
            u0 = np.exp(-200*(x-0.5)**2)
            results = solve_heat_1d_cn(u0, L, Nx, dt, Nt, alpha)
            pd.DataFrame(results).to_csv("heat1d_solution.csv", index=False)
            print("Wrote heat1d_solution.csv. Will animate.")
            animate_1d(results, x, dt)
            return
        elif guess == 'wave':
            print("--> Wave equation detected (or similar). Running leapfrog 1D with defaults.")
            L = 1.0; Nx = 201; dt = 0.0005; Nt = 800; c = 1.0
            x = np.linspace(0, L, Nx)
            u0 = np.exp(-200*(x-0.3)**2)
            v0 = np.zeros_like(x)
            results = solve_wave_1d_leapfrog(u0, v0, L, Nx, dt, Nt, c)
            pd.DataFrame(results).to_csv("wave1d_solution.csv", index=False)
            print("Wrote wave1d_solution.csv. Will animate.")
            animate_1d(results, x, dt)
            return
        elif guess == 'schrodinger':
            print("--> Schrödinger-like equation detected. Running radial TDSE CN with defaults.")
            rmax = 20.0; Nr = 400; dt = 0.005; Nt = 400
            r = np.linspace(0, rmax, Nr+1)
            r0 = 8.0
            u0 = (r * np.exp(-0.1*(r-r0)**2)).astype(complex)
            rr, data = solve_tdse_radial_cn(u0, rmax, Nr, dt, Nt, V_func=None)
            # Guardar densidad |psi|^2
            psisq = []
            for n in range(data.shape[0]):
                psi = data[n] / rr
                psi[0] = 0.0
                psisq.append(np.abs(psi)**2)
            psisq = np.array(psisq)
            pd.DataFrame(psisq).to_csv("tdse_radial_density.csv", index=False)
            print("Wrote tdse_radial_density.csv. Will animate.")
            animate_1d(psisq, rr, dt, ylabel='|psi|^2', is_complex=False)
            return
        else:
            print("PDE detected but type unknown or unsupported for automatic solving.")
            return

    print("Could not classify equation or unsupported; try giving a more specific LaTeX or use --examples.")

# ---------- Animation helper ----------
def animate_1d(result_array, x_grid, dt, ylabel='u', is_complex=False):
    fig, ax = plt.subplots()
    line, = ax.plot(x_grid, result_array[0, :], lw=2)
    ax.set_xlim(x_grid[0], x_grid[-1])
    ymin = np.min(result_array); ymax = np.max(result_array)
    if ymin == ymax:
        ymin -= 1.0; ymax += 1.0
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('x')
    ax.set_ylabel(ylabel)
    title = ax.set_title('t = 0.0')
    def update(i):
        line.set_ydata(result_array[i, :])
        title.set_text(f"t = {i*dt:.4f}")
        return line, title
    FuncAnimation(fig, update, frames=result_array.shape[0], interval=30, blit=False)
    plt.show()

# ---------- CLI ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--latex", type=str, default=None, help="LaTeX string for the equation")
    parser.add_argument("--examples", action='store_true', help="Show example LaTeX inputs")
    args = parser.parse_args()
    if args.examples:
        print("Examples you can paste (wrap in quotes):")
        print(r'  "\\frac{d^2}{dt^2} y(t) + \\omega^2 y(t) = 0"  # simple harmonic oscillator')
        print(r'  "u_t = D u_{xx}"  # heat/diffusion')
        print(r'  "u_{tt} = c^2 u_{xx}"  # wave')
        print(r'  "i*hbar \\frac{\\partial \\psi}{\\partial t} = -\\frac{\\hbar^2}{2m} \\frac{\\partial^2 \\psi}{\\partial r^2} + V(r) \\psi"  # radial TDSE')
        sys.exit(0)
    if args.latex is None:
        print("No LaTeX provided. Use --examples to see sample inputs, or run the demo by providing --latex.")
        sys.exit(0)
    auto_solve_from_latex(args.latex)
