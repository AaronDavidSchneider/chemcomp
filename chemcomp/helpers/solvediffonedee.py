from .tridag import *


def solvediffonedee_components(x, y, v, d, g, s, bcl, bcr, dt=None, int=False, upwind=False,
                               extracond=None, retall=False):
    """
    Vectorised version of solvediffonedee (works with sigma_g_components_mol)

    Solve the 1-D linear advection-diffusion equation with boundary conditions,
    either in stationary state, or as an implicit time step.

      dy(x,t)    d /                          d / y(x,t) \ \
      ------- + -- | y(x,t)*v(x) - d(x)*g(x)*-- | ------ | | = s(x)
         dt     dx \                         dx \  g(x)  / /

    If dt is not given, then it is assumed to be infinite, meaning that the
    time-derivative is zero. Geometric terms (e.g. r^2 for use in spherical
    coordinates) can be included by treating y as yorig*r^2 where yorig is
    the actual variable you want to advect/diffuse. You must then also
    take g = gorig*r^2 and s = sorig*r^2. For cylindrical coordinates it
    would be r instead of r^2. Cell interfaces will be exactly between the
    cell center positions x. For the boundaries they are at the positions
    x[0] and x[-1] respectively.

    NOTE: The cells and interfaces are indexed as follows:

      cells:       |  0  |  1  |   2  |   3  |   4  |   5  |       for n=5
      interfaces:        0     1      2      3      4

    The boundary conditions are of the type:

        dy
      p -- + q*y = r
        dx

    ARGUMENTS:

      x          Array of grid cell center points in coordinate x (n-elements array)
      y          Current values of y at grid points (important if dt given)
      v          Velocity values at cell centers, except if int=True: then at
                 cell interfaces, in which case v has only n-1 elements.
      d          Diffusion constant at cell centers, except if int=True (see v)
      g          The g-function values at cell-centers (see equation above).
      s          The s-function values at cell-centers (see equation above).
      bcl        The left boundary condition: tuple of 3 values (p,q,r,type)
                 If type=0: apply it to y, if type=1: apply bc to y/g
      bcr        The right boundary condition: tuple of 3 values (p,q,r,type,ir)
                 modified to be able to specify at which point the bc sits
      dt         If not None, then this is the time step
      int        Determines if v and d are located at cell interface (if True) or center.
      upwind     If true, then treat advection term as upwind
      extracond  Sometimes it may be useful to provide additional internal conditions.
                 They are as bcl but with an extra ir: ic=(p,q,r,type,ir). The extracond
                 (if not None) is a list of tuples: [ic1,ic2,ic3].
    """
    #
    # Get the array size
    #
    n_components = y.shape
    #
    # Set up the a,b,c bands of the tridag matrix and the rhs r
    #
    a = np.zeros(n_components)
    b = np.zeros(n_components)
    c = np.zeros(n_components)
    r = np.zeros(n_components)
    #
    # Grid stuff
    #
    xi = 0.5 * (x[1:, np.newaxis] + x[:-1, np.newaxis])
    dxc = xi[1:] - xi[:-1]
    dxc = np.array([xi[0] - x[0, np.newaxis], *dxc, x[-1, np.newaxis] - xi[-1]])
    dxi = x[1:, np.newaxis] - x[:-1, np.newaxis]
    #
    # If v and d are at cell centers, we have to compute the interface values
    #
    if int:
        vi = v
        di = d
    else:
        vi = 0.5 * (v[1:] + v[:-1])
        di = 0.5 * (d[1:] + d[:-1])
    #
    # Find g values at interfaces
    #
    gi = 0.5 * (g[1:] + g[:-1])
    #
    # Add velocity term to matrix
    #
    if upwind:
        a[1:-1] -= (vi[:-1] >= 0) * vi[:-1] / dxc[1:-1]
        b[1:-1] -= (vi[:-1] < 0) * vi[:-1] / dxc[1:-1]
        b[1:-1] += (vi[1:] >= 0) * vi[1:] / dxc[1:-1]
        c[1:-1] += (vi[1:] < 0) * vi[1:] / dxc[1:-1]
    else:
        a[1:-1] -= 0.5 * vi[:-1] / dxc[1:-1]
        b[1:-1] -= 0.5 * vi[:-1] / dxc[1:-1]
        b[1:-1] += 0.5 * vi[1:] / dxc[1:-1]
        c[1:-1] += 0.5 * vi[1:] / dxc[1:-1]
    #
    # Add diffusion term to matrix
    #
    a[1:-1] -= di[:-1] * gi[:-1] / (g[:-2] * dxi[:-1] * dxc[1:-1])
    b[1:-1] += di[:-1] * gi[:-1] / (g[1:-1] * dxi[:-1] * dxc[1:-1])
    b[1:-1] += di[1:] * gi[1:] / (g[1:-1] * dxi[1:] * dxc[1:-1])
    c[1:-1] -= di[1:] * gi[1:] / (g[2:] * dxi[1:] * dxc[1:-1])
    #
    # Add chemcomp term
    #
    r[1:-1] += s[1:-1]
    #
    # If dt is set, add terms for time derivative
    #
    if dt is not None:
        b[1:-1] += 1.0 / dt
        r[1:-1] += y[1:-1] / dt
    #
    # Left boundary condition
    #
    if bcl[3] == 0:
        #
        # Apply it to y
        #
        a[0] = 0.
        b[0] = bcl[1] - bcl[0] / dxi[0]
        c[0] = bcl[0] / dxi[0]
        r[0] = bcl[2]
    else:
        #
        # Apply it to y/g
        #
        a[0] = 0.
        b[0] = (bcl[1] - bcl[0] / dxi[0]) / g[0]
        c[0] = (bcl[0] / dxi[0]) / g[1]
        r[0] = bcl[2]
    #
    # Right boundary condition
    #
    if bcr[3] == 0:
        #
        # Apply it to y
        #
        a[-1] = -bcr[0] / dxi[-1]
        b[-1] = bcr[1] + bcr[0] / dxi[-1]
        c[-1] = 0.
        r[-1] = bcr[2]
    else:
        #
        # Apply it to y/g
        #
        a[-1] = (-bcr[0] / dxi[-1]) / g[-2]
        b[-1] = (bcr[1] + bcr[0] / dxi[-1]) / g[-1]
        c[-1] = 0.
        r[-1] = bcr[2]
    #
    # Internal conditions, if present
    #
    if extracond is not None:
        for ic in extracond:
            ir = ic[4]
            ddx = x[ir[0] + 1] - x[ir[0] - 1]
            if ic[3] == 0:
                #
                # Apply it to y
                #
                a[ir] = -ic[0] / ddx
                b[ir] = ic[1]
                c[ir] = ic[0] / ddx
                r[ir] = ic[2]
            else:
                #
                # Apply it to y/g
                #
                a[ir] = (-ic[0] / ddx) / g[1]
                b[ir] = (ic[1]) / g[0]
                c[ir] = (ic[0] / ddx) / g[1]
                r[ir] = ic[2]
    #
    # Normalize
    #
    fact = 1 / np.abs(b)
    a *= fact
    b *= fact
    c *= fact
    r *= fact
    #
    # Now solve
    #
    sol = tridag_components(a, b, c, r)
    #
    # Return
    #
    if retall:
        return sol, a, b, c, r
    else:
        return sol


def solvediffonedee(x, y, v, d, g, s, bcl, bcr, dt=None, int=False, upwind=False, extracond=None, retall=False):
    """
    Solve the 1-D linear advection-diffusion equation with boundary conditions,
    either in stationary state, or as an implicit time step.

      dy(x,t)    d /                          d / y(x,t) \ \
      ------- + -- | y(x,t)*v(x) - d(x)*g(x)*-- | ------ | | = s(x)
         dt     dx \                         dx \  g(x)  / /

    If dt is not given, then it is assumed to be infinite, meaning that the
    time-derivative is zero. Geometric terms (e.g. r^2 for use in spherical
    coordinates) can be included by treating y as yorig*r^2 where yorig is
    the actual variable you want to advect/diffuse. You must then also
    take g = gorig*r^2 and s = sorig*r^2. For cylindrical coordinates it
    would be r instead of r^2. Cell interfaces will be exactly between the
    cell center positions x. For the boundaries they are at the positions
    x[0] and x[-1] respectively.

    NOTE: The cells and interfaces are indexed as follows:

      cells:       |  0  |  1  |   2  |   3  |   4  |   5  |       for n=5
      interfaces:        0     1      2      3      4

    The boundary conditions are of the type:

        dy
      p -- + q*y = r
        dx

    ARGUMENTS:

      x          Array of grid cell center points in coordinate x (n-elements array)
      y          Current values of y at grid points (important if dt given)
      v          Velocity values at cell centers, except if int=True: then at
                 cell interfaces, in which case v has only n-1 elements.
      d          Diffusion constant at cell centers, except if int=True (see v)
      g          The g-function values at cell-centers (see equation above).
      s          The s-function values at cell-centers (see equation above).
      bcl        The left boundary condition: tuple of 3 values (p,q,r,type)
                 If type=0: apply it to y, if type=1: apply bc to y/g
      bcr        The right boundary condition: tuple of 3 values (p,q,r,type)
      dt         If not None, then this is the time step
      int        Determines if v and d are located at cell interface (if True) or center.
      upwind     If true, then treat advection term as upwind
      extracond  Sometimes it may be useful to provide additional internal conditions.
                 They are as bcl but with an extra ir: ic=(p,q,r,type,ir). The extracond
                 (if not None) is a list of tuples: [ic1,ic2,ic3].
    """
    #
    # Get the array size
    #
    n = x.size
    #
    # Set up the a,b,c bands of the tridag matrix and the rhs r
    #
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    r = np.zeros(n)
    #
    # Grid stuff
    #
    xi = 0.5 * (x[1:] + x[:-1])
    dxc = xi[1:] - xi[:-1]
    dxc = np.hstack((xi[0] - x[0], dxc, x[-1] - xi[-1]))
    dxi = x[1:] - x[:-1]
    #
    # If v and d are at cell centers, we have to compute the interface values
    #
    if int:
        vi = v
        di = d
    else:
        vi = 0.5 * (v[1:] + v[:-1])
        di = 0.5 * (d[1:] + d[:-1])
    #
    # Find g values at interfaces
    #
    gi = 0.5 * (g[1:] + g[:-1])
    #
    # Add velocity term to matrix
    #
    if upwind:
        a[1:-1] -= (vi[:-1] >= 0) * vi[:-1] / dxc[1:-1]
        b[1:-1] -= (vi[:-1] < 0) * vi[:-1] / dxc[1:-1]
        b[1:-1] += (vi[1:] >= 0) * vi[1:] / dxc[1:-1]
        c[1:-1] += (vi[1:] < 0) * vi[1:] / dxc[1:-1]
    else:
        a[1:-1] -= 0.5 * vi[:-1] / dxc[1:-1]
        b[1:-1] -= 0.5 * vi[:-1] / dxc[1:-1]
        b[1:-1] += 0.5 * vi[1:] / dxc[1:-1]
        c[1:-1] += 0.5 * vi[1:] / dxc[1:-1]
    #
    # Add diffusion term to matrix
    #
    a[1:-1] -= di[:-1] * gi[:-1] / (g[:-2] * dxi[:-1] * dxc[1:-1])
    b[1:-1] += di[:-1] * gi[:-1] / (g[1:-1] * dxi[:-1] * dxc[1:-1])
    b[1:-1] += di[1:] * gi[1:] / (g[1:-1] * dxi[1:] * dxc[1:-1])
    c[1:-1] -= di[1:] * gi[1:] / (g[2:] * dxi[1:] * dxc[1:-1])
    #
    # Add chemcomp term
    #
    r[1:-1] += s[1:-1]
    #
    # If dt is set, add terms for time derivative
    #
    if dt is not None:
        b[1:-1] += 1.0 / dt
        r[1:-1] += y[1:-1] / dt
    #
    # Left boundary condition
    #
    if bcl[3] == 0:
        #
        # Apply it to y
        #
        a[0] = 0.
        b[0] = bcl[1] - bcl[0] / dxi[0]
        c[0] = bcl[0] / dxi[0]
        r[0] = bcl[2]
    else:
        #
        # Apply it to y/g
        #
        a[0] = 0.
        b[0] = (bcl[1] - bcl[0] / dxi[0]) / g[0]
        c[0] = (bcl[0] / dxi[0]) / g[1]
        r[0] = bcl[2]
    #
    # Right boundary condition
    #
    if bcr[3] == 0:
        #
        # Apply it to y
        #
        a[-1] = -bcr[0] / dxi[-1]
        b[-1] = bcr[1] + bcr[0] / dxi[-1]
        c[-1] = 0.
        r[-1] = bcr[2]
    else:
        #
        # Apply it to y/g
        #
        a[-1] = (-bcr[0] / dxi[-1]) / g[-2]
        b[-1] = (bcr[1] + bcr[0] / dxi[-1]) / g[-1]
        c[-1] = 0.
        r[-1] = bcr[2]
    #
    # Internal conditions, if present
    #
    if extracond is not None:
        for ic in extracond:
            ir = ic[4]
            ddx = x[ir + 1] - x[ir - 1]
            if ic[3] == 0:
                #
                # Apply it to y
                #
                a[ir] = -ic[0] / ddx
                b[ir] = ic[1]
                c[ir] = ic[0] / ddx
                r[ir] = ic[2]
            else:
                #
                # Apply it to y/g
                #
                a[ir] = (-ic[0] / ddx) / g[1]
                b[ir] = (ic[1]) / g[0]
                c[ir] = (ic[0] / ddx) / g[1]
                r[ir] = ic[2]
    #
    # Normalize
    #
    fact = 1.0 / np.abs(b)
    a *= fact
    b *= fact
    c *= fact
    r *= fact
    #
    # Now solve
    #
    sol = tridag(a, b, c, r)
    #
    # Return
    #
    if retall:
        return sol, a, b, c, r
    else:
        return sol


def getfluxonedee(x, y, v, d, g, int=False, upwind=False, oned=False):
    """
    Compute the flux at the cell interfaces, which are located at
    xi  = 0.5 * (x[1:] + x[:-1]). The fluxes are computed in the same
    way as solvediffonedee(), so these routines are mutually consistent.
    works with sigma_g_components
    """
    #
    # Get the array size
    #
    n = x.size
    #
    # Set up the b and c bands of the tridag matrix. We do not need a, because
    # we only look at the right cell wall
    #
    b = np.zeros(n - 1)
    c = np.zeros(n - 1)
    #
    # Grid stuff
    #
    dxi = x[1:] - x[:-1]
    #
    # If v and d are at cell centers, we have to compute the interface values
    #
    if int:
        vi = v
        di = d
    else:
        vi = 0.5 * (v[1:] + v[:-1])
        di = 0.5 * (d[1:] + d[:-1])
    #
    # Find g values at interfaces
    #
    gi = 0.5 * (g[1:] + g[:-1])
    #
    # Add velocity term to matrix
    #
    if upwind:
        b += (vi >= 0) * vi
        c += (vi < 0) * vi
    else:
        b += 0.5 * vi
        c += 0.5 * vi
    #
    # Add diffusion term to matrix
    #
    b += di * gi / (g[0:-1] * dxi)
    c -= di * gi / (g[1:] * dxi)
    #
    # Now compute the flux
    #
    if not oned:
        flux = b[:, np.newaxis, np.newaxis] * y[:-1] + c[:, np.newaxis, np.newaxis] * y[1:]
    else:
        flux = b * y[:-1] + c * y[1:]
    #
    # Return
    #
    return flux


def solvehyperviscosity(x, y, d, dt, int=False, retall=False):
    """
    Solve the 1-D linear hyperviscosity equation with keeping the two boundary
    cells at each side fixed. Integrate as an implicit time step.

      dy(x,t)    d /      d^3 y(x,t) \ \
      ------- + -- | d(x)*---------- | | = 0
         dt     dx \         dx^3    / /

    NOTE: The cells and interfaces are indexed as follows:

      cells:       |  0  |  1  |   2  |   3  |   4  |   5  |       for n=5
      interfaces:        0     1      2      3      4

    ARGUMENTS:

      x          Array of grid cell center points in coordinate x (n-elements array)
      y          Current values of y at grid points (important if dt given)
      d          Diffusion constant at cell centers, except if int=True (see v)
      dt         If not None, then this is the time step
      int        Determines if v and d are located at cell interface (if True) or center.
    """
    #
    # Get the array size
    #
    n = x.size
    #
    # Set up the aa,a,b,c,cc bands of the pentdag matrix and the rhs r
    #
    aa = np.zeros(n)
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    cc = np.zeros(n)
    r = np.zeros(n)
    #
    # Grid stuff
    #
    xi = 0.5 * (x[1:] + x[:-1])
    dxc = xi[1:] - xi[:-1]
    dxc = np.hstack((xi[0] - x[0], dxc, x[-1] - xi[-1]))
    dxi = x[1:] - x[:-1]
    #
    # If d is at cell centers, we have to compute the interface values
    #
    if int:
        di = d
    else:
        di = 0.5 * (d[1:] + d[:-1])
    #
    # Add hyperdiffusion term to matrix
    #
    aa[2:-2] += di[1:-2] / (dxi[0:-3] * dxc[1:-3] * dxi[1:-2] * dxc[2:-2])
    a[2:-2] -= di[1:-2] / (dxi[0:-3] * dxc[1:-3] * dxi[1:-2] * dxc[2:-2])
    a[2:-2] -= di[1:-2] / (dxi[1:-2] * dxc[1:-3] * dxi[1:-2] * dxc[2:-2])
    a[2:-2] -= di[1:-2] / (dxi[1:-2] * dxc[2:-2] * dxi[1:-2] * dxc[2:-2])
    b[2:-2] += di[1:-2] / (dxi[1:-2] * dxc[1:-3] * dxi[1:-2] * dxc[2:-2])
    b[2:-2] += di[1:-2] / (dxi[2:-1] * dxc[2:-2] * dxi[1:-2] * dxc[2:-2])
    b[2:-2] += di[1:-2] / (dxi[1:-2] * dxc[2:-2] * dxi[1:-2] * dxc[2:-2])
    c[2:-2] -= di[1:-2] / (dxi[2:-1] * dxc[2:-2] * dxi[1:-2] * dxc[2:-2])
    a[2:-2] -= di[2:-1] / (dxi[1:-2] * dxc[2:-2] * dxi[2:-1] * dxc[2:-2])
    b[2:-2] += di[2:-1] / (dxi[2:-1] * dxc[2:-2] * dxi[2:-1] * dxc[2:-2])
    b[2:-2] += di[2:-1] / (dxi[1:-2] * dxc[2:-2] * dxi[2:-1] * dxc[2:-2])
    b[2:-2] += di[2:-1] / (dxi[2:-1] * dxc[3:-1] * dxi[2:-1] * dxc[2:-2])
    c[2:-2] -= di[2:-1] / (dxi[2:-1] * dxc[2:-2] * dxi[2:-1] * dxc[2:-2])
    c[2:-2] -= di[2:-1] / (dxi[2:-1] * dxc[3:-1] * dxi[2:-1] * dxc[2:-2])
    c[2:-2] -= di[2:-1] / (dxi[3:] * dxc[3:-1] * dxi[2:-1] * dxc[2:-2])
    cc[2:-2] += di[2:-1] / (dxi[3:] * dxc[3:-1] * dxi[2:-1] * dxc[2:-2])
    #
    # If dt is set, add terms for time derivative
    #
    b += 1.0 / dt
    r += y / dt
    #
    # Normalize
    #
    fact = 1.0 / np.abs(b)
    aa *= fact
    a *= fact
    b *= fact
    c *= fact
    cc *= fact
    r *= fact
    #
    # Now solve
    #
    sol = pentdag(aa, a, b, c, cc, r)
    #
    # Return
    #
    if retall:
        return sol, aa, a, b, c, cc, r
    else:
        return sol
