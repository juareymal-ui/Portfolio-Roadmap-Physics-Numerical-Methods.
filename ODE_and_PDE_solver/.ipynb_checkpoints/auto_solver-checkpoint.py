#!/usr/bin/env python3
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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sympy import symbols, Function, Eq, Derivative, Symbol, simplify, pprint
from sympy.parsing.latex import parse_latex
from sympy import dsolve
from sympy import S

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
    # If the parsed expression is an Eq (Relational), return it.
    # Otherwise treat as expr == 0
    from sympy import Eq as SymEq, Relational
    if isinstance(expr, RelEq := getattr(sys.modules['sympy'], 'Relational', None)):
        # handle generically
        return expr
    # Better: check .is_Relational
    if hasattr(expr, 'is_Relational') and expr.is_Relational:
        return expr
    # treat as expr = 0
    return Eq(expr, 0)

def find_derivative_atoms(expr):
    return list(expr.atoms(Derivative))

def unique_independent_vars(derivs):
    vars_set = set()
    for d in derivs:
        vars_set.update(d.variables)
    return vars_set

# ---------- Classification ----------
def classify(parsed):
    """
    Return a dictionary describing what we think the equation is.
    """
    # If equal relation: lhs, rhs
    if getattr(parsed, 'is_Equality', False) or getattr(parsed, 'is_Relational', False):
        lhs = parsed.lhs
        rhs = parsed.rhs
        expr = lhs - rhs
    else:
        expr = parsed

    derivs = find_derivative_atoms(expr)
    # collect derivative variable symbols (e.g. t, x)
    vars_in_derivs = unique_independent_vars(derivs)

    # Look for partial derivatives vs total derivatives by name
    # We'll search textual for "partial" / "Derivative" etc.
    text = str(parsed)

    classification = {
        'expr': expr,
        'derivs': derivs,
        'vars_in_derivs': vars_in_derivs,
        'is_pde': False,
        'is_ode': False,
        'time_vars': set(),
        'space_vars': set(),
        'order_time': 0,
        'order_space': 0,
        'notes': []
    }

    # Heuristics: treat symbols named t as time, x/r as space, theta/phi as angular
    for v in vars_in_derivs:
        name = str(v)
        if name in ('t', 'time'):
            classification['time_vars'].add(v)
        elif name in ('x', 'r', 'radius', 'rho'):
            classification['space_vars'].add(v)
        elif name in ('theta', 'phi'):
            classification['space_vars'].add(v)
        else:
            # if only single var present, guess it's time
            classification['time_vars'].add(v)

    # count orders
    max_time_order = 0
    max_space_order = 0
    for d in derivs:
        # d.variables is tuple of variables in order they were derived: e.g. (t, t) is d2/dt2
        tv = [v for v in d.variables if str(v) in ('t', 'time')]
        sv = [v for v in d.variables if str(v) in ('x','r','theta','phi','rho')]
        if len(tv) > max_time_order:
            max_time_order = len(tv)
        if len(sv) > max_space_order:
            max_space_order = len(sv)
    classification['order_time'] = max_time_order
    classification['order_space'] = max_space_order

    # decide ODE vs PDE
    if len(classification['space_vars'])>0:
        classification['is_pde'] = True
    else:
        if len(classification['time_vars'])>0:
            classification['is_ode'] = True

    # guess PDE type
    s = str(parsed)
    if 'I' in s or 'I*' in s or 'I' in s or 'i' in s or 'hbar' in s or 'psi' in s:
        classification['guess'] = 'schrodinger'
    elif 'u_tt' in s or 'd2/dt2' in s or 'ddot' in s:
        classification['guess'] = 'wave'
    elif 'u_t' in s or 'partial' in s or 'u_xx' in s or 'd2/dx2' in s:
        classification['guess'] = 'heat_or_diffusion'
    else:
        classification['guess'] = 'unknown'

    return classification

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

def numeric_ode_integrate(expr, dependent_fn=None, t0=0.0, tmax=10.0, dt=0.01, y0=None):
    """
    Convert higher-order ODE into first-order system and integrate by RK4.
    dependent_fn: symbolic Function object (like y(t)) if available
    y0: initial values list for y, y', ...
    """
    from sympy import lambdify, Function, Derivative
    # detect dependent function and independent var
    derivs = list(expr.atoms(Derivative))
    if len(derivs)==0:
        raise RuntimeError("No derivatives found for numeric ODE integration.")
    # find the function symbol (e.g. y(t))
    funcs = [a.expr for a in derivs if hasattr(a, 'expr')]
    # fallback: search for Function in expr free symbols
    # pick independent variable as the first variable found in derivatives
    indep_vars = set()
    for d in derivs:
        indep_vars.update(d.variables)
    indep = list(indep_vars)[0]
    t = indep
    # find dependent function symbol name by scanning expr free symbols that are functions
    funcs = [f for f in expr.atoms(Function)]
    if len(funcs)==0:
        # try to find symbol that depends on t
        funcs = []
        for a in expr.atoms(Symbol):
            pass
        raise RuntimeError("Could not detect dependent function for numeric ODE. Provide simple y(t) style.")
    y_func = funcs[0]  # e.g. y(t)
    # determine order
    order = 0
    for d in derivs:
        if d.expr==y_func.func and len(d.variables)>order:
            order = len(d.variables)
    # Build first-order system: define state vector [y, y', y'', ...]
    # Solve for highest derivative: write expr = 0 => y^(n) = f(t, y, y', ...)
    # rearrange to isolate highest derivative symbolically
    # For simplicity, use sympy.solve for the derivative symbol
    highest_deriv = None
    for d in derivs:
        if d.expr==y_func.func and len(d.variables)==order:
            highest_deriv = d
            break
    if highest_deriv is None:
        raise RuntimeError("Couldn't find highest derivative form.")

    # Symbol for the highest derivative function (Derivative(y(t), (t,order)))
    from sympy import solve
    sol_for_high = solve(Eq(expr,0), highest_deriv)
    if len(sol_for_high)==0:
        raise RuntimeError("Cannot algebraically isolate highest derivative.")
    rhs_high = sol_for_high[0]
    # Create lambdified RHS function that takes (t, y0, y1, ... y_{n-1})
    Ys = symbols(' '.join([f'y{i}' for i in range(order)]))
    # replace y, y', ... with Ys in rhs_high
    repl = {}
    for i in range(order):
        repl[Derivative(y_func, *([t]*i))] = Ys[i] if i>0 else y_func.func(Ys[0]) # not ideal
    # Simpler: build lambda by lambdify of rhs_high with y(t), y'(t), ... passed as symbols
    lam = lambdify((t,)+tuple(Derivative(y_func, *([t]*i)) for i in range(order)), rhs_high, 'numpy')
    # We will create an evaluation wrapper that given (t, state) returns derivatives
    def rhs_time(tval, state):
        # state has length = order: state[0]=y, state[1]=y', etc
        args = [tval] + list(state)
        return lam(*args)

    # Now build full first-order system: dy0/dt = y1; dy1/dt = y2; ...; dy_{n-1}/dt = rhs_time
    def system_rhs(tval, state):
        out = np.zeros_like(state, dtype=float)
        n = len(state)
        for i in range(n-1):
            out[i] = state[i+1]
        out[-1] = rhs_time(tval, *state)
        return out

    # integrate by RK4
    times = np.arange(t0, tmax+1e-12, dt)
    if y0 is None:
        y0 = np.zeros(order)
        y0[0] = 1.0  # arbitrary
    sol = np.zeros((len(times), order))
    sol[0,:] = y0
    for k in range(1,len(times)):
        tt = times[k-1]
        h = dt
        s = sol[k-1,:]
        k1 = system_rhs(tt, s)
        k2 = system_rhs(tt + h/2, s + 0.5*h*k1)
        k3 = system_rhs(tt + h/2, s + 0.5*h*k2)
        k4 = system_rhs(tt + h, s + h*k3)
        sol[k,:] = s + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return times, sol

# ---------- PDE solvers (1D) ----------
def solve_heat_1d_cn(u0, L=1.0, Nx=201, dt=1e-4, Nt=400, alpha=1.0):
    dx = L/(Nx-1)
    r = alpha*dt/(2*dx*dx)
    N = Nx
    # interior size N-2
    a = -r*np.ones(N-3); b = (1+2*r)*np.ones(N-2); c = -r*np.ones(N-3)
    aB = r*np.ones(N-3); bB = (1-2*r)*np.ones(N-2); cB = r*np.ones(N-3)
    # use complex in case
    u_prev = u0.copy()
    results = [u_prev.copy()]
    for n in range(Nt):
        u_in = u_prev[1:-1]
        d = bB * u_in.copy()
        d[:-1] += aB * u_in[1:]
        d[1:] += cB * u_in[:-1]
        # solve tri
        u_in_new = thomas(a,b,c,d)
        u_new = np.zeros_like(u_prev)
        u_new[1:-1] = np.real(u_in_new)
        results.append(u_new)
        u_prev = u_new
    return np.array(results)

def solve_wave_1d_leapfrog(u0, v0, L=1.0, Nx=201, dt=1e-4, Nt=400, c=1.0):
    dx = L/(Nx-1)
    lam = c*dt/dx
    if lam>1.0:
        print("Warning: CFL > 1; scheme unstable.")
    u_prev = u0.copy()
    # u_curr by Taylor
    u_curr = u0.copy()
    uxx = np.zeros_like(u0)
    for i in range(1, Nx-1):
        uxx[i] = (u0[i+1] - 2*u0[i] + u0[i-1]) / dx**2
    u_curr[1:-1] = u0[1:-1] + dt*v0[1:-1] + 0.5*(dt**2)*c**2*uxx[1:-1]
    results = [u_prev.copy(), u_curr.copy()]
    for n in range(1,Nt):
        u_next = np.zeros_like(u0)
        for i in range(1,Nx-1):
            u_next[i] = 2*u_curr[i] - u_prev[i] + (lam**2)*(u_curr[i+1] - 2*u_curr[i] + u_curr[i-1])
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
        for j in range(1,Nr+1):
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
    for j in range(1,Nr):
        b[j-1] = 1.0 + 2.0*coef*k + coef*V[j]
    aB = np.full(N_inner-1, coef*k, dtype=complex)
    cB = np.full(N_inner-1, coef*k, dtype=complex)
    bB = np.zeros(N_inner, dtype=complex)
    for j in range(1,Nr):
        bB[j-1] = 1.0 - 2.0*coef*k - coef*V[j]
    for n in range(Nt):
        u_in = u[1:-1]
        d = bB * u_in
        d[:-1] += aB * u_in[1:]
        d[1:] += cB * u_in[:-1]
        u_inner_next = thomas(a,b,c,d)
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
    classification = classify(parsed)
    print("Classification summary:", classification['guess'],
          "is_pde?", classification['is_pde'],
          "order_time", classification['order_time'],
          "order_space", classification['order_space'])
    # If ODE:
    if classification['is_ode'] and not classification['is_pde']:
        print("--> Detected ODE. Trying symbolic solve first.")
        sym_sol = try_symbolic_ode(classification['expr'])
        if sym_sol is not None:
            print("Symbolic solution obtained. Writing to file symbolic_solution.txt")
            with open("symbolic_ode_solution.txt","w") as f:
                f.write(str(sym_sol))
            return
        else:
            print("Falling back to numeric RK4 integration (defaults).")
            try:
                t, sol = numeric_ode_integrate(classification['expr'])
                # save first component (y) if higher order
                df = pd.DataFrame(sol)
                df.insert(0, "t", t)
                df.to_csv("ode_numeric_solution.csv", index=False)
                print("Numeric ODE solution saved to ode_numeric_solution.csv")
                # quick plot
                import matplotlib.pyplot as plt
                plt.plot(t, sol[:,0])
                plt.xlabel('t'); plt.ylabel('y (first component)')
                plt.title('Numeric ODE RK4 solution')
                plt.show()
                return
            except Exception as e:
                print("Numeric ODE solve failed:", e)
                return
    # If PDE:
    if classification['is_pde']:
        guess = classification.get('guess','unknown')
        if guess.startswith('heat'):
            print("--> Heat equation detected (or similar). Running Crank-Nicolson 1D with defaults.")
            L = 1.0; Nx = 201; dt = 5e-5; Nt = 400; alpha = 1.0
            x = np.linspace(0,L,Nx)
            u0 = np.exp(-200*(x-0.5)**2)
            results = solve_heat_1d_cn(u0, L, Nx, dt, Nt, alpha)
            pd.DataFrame(results).to_csv("heat1d_solution.csv", index=False)
            print("Wrote heat1d_solution.csv. Will animate.")
            animate_1d(results, x, dt)
            return
        elif guess=='wave':
            print("--> Wave equation detected (or similar). Running leapfrog 1D with defaults.")
            L=1.0; Nx=201; dt=0.0005; Nt=800; c=1.0
            x = np.linspace(0,L,Nx)
            u0 = np.exp(-200*(x-0.3)**2)
            v0 = np.zeros_like(x)
            results = solve_wave_1d_leapfrog(u0, v0, L, Nx, dt, Nt, c)
            pd.DataFrame(results).to_csv("wave1d_solution.csv", index=False)
            print("Wrote wave1d_solution.csv. Will animate.")
            animate_1d(results, x, dt)
            return
        elif guess=='schrodinger':
            print("--> Schrödinger-like equation detected. Running radial TDSE CN with defaults.")
            rmax=20.0; Nr=400; dt=0.005; Nt=400
            r = np.linspace(0, rmax, Nr+1)
            r0 = 8.0
            u0 = (r * np.exp(-0.1*(r-r0)**2)).astype(complex)
            rr, data = solve_tdse_radial_cn(u0, rmax, Nr, dt, Nt, V_func=None)
            # save density |psi|^2
            psisq = []
            for n in range(data.shape[0]):
                psi = data[n] / rr
                psi[0] = 0.0
                psisq.append(np.abs(psi)**2)
            pd.DataFrame(np.array(psisq)).to_csv("tdse_radial_density.csv", index=False)
            print("Wrote tdse_radial_density.csv. Will animate.")
            animate_1d(np.array(psisq), rr, dt, ylabel='|psi|^2', is_complex=False)
            return
        else:
            print("PDE detected but type unknown or unsupported for automatic solving.")
            return

    print("Could not classify equation or unsupported; try giving a more specific LaTeX or use --examples.")

# ---------- Animation helper ----------
def animate_1d(result_array, x_grid, dt, ylabel='u', is_complex=False):
    # result_array: shape (Nt+1, Nx)
    fig, ax = plt.subplots()
    line, = ax.plot(x_grid, result_array[0,:], lw=2)
    ax.set_xlim(x_grid[0], x_grid[-1])
    ymin = np.min(result_array); ymax = np.max(result_array)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('x')
    ax.set_ylabel(ylabel)
    title = ax.set_title('t=0.0')
    def update(i):
        line.set_ydata(result_array[i,:])
        title.set_text(f"t = {i*dt:.4f}")
        return line, title
    ani = FuncAnimation(fig, update, frames=result_array.shape[0], interval=30, blit=False)
    plt.show()

# ---------- CLI ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--latex", type=str, default=None, help="LaTeX string for the equation")
    parser.add_argument("--examples", action='store_true', help="Show example LaTeX inputs")
    args = parser.parse_args()
    if args.examples:
        print("Examples you can paste (wrap in quotes):")
        print(r'  "\\frac{d^2 y}{dt^2} + \\omega^2 y = 0"  # simple harmonic oscillator')
        print(r'  "u_t = D u_{xx}"  # heat/diffusion')
        print(r'  "u_{tt} = c^2 u_{xx}"  # wave')
        print(r'  "i*hbar \\frac{\\partial \\psi}{\\partial t} = -\\frac{\\hbar^2}{2m} \\frac{\\partial^2 \\psi}{\\partial r^2} + V(r) \\psi"  # radial TDSE')
        sys.exit(0)
    if args.latex is None:
        print("No LaTeX provided. Use --examples to see sample inputs, or run the demo by providing --latex.")
        sys.exit(0)
    auto_solve_from_latex(args.latex)
